# oligoscan — project handoff

**An exhaustive, uncapped alignment scanner for short oligos (siRNA / ASO)
against cDNA and whole-genome references.**

Name: **oligoscan** (decided 2026-11-12). License: **source-available,
free for academic / non-profit research, commercial use requires a licence**
(**PolyForm Noncommercial 1.0.0**, chosen 2026-11-16, replacing MIT; see 10). Target venue:
**Bioinformatics, Application Note**. Version 0.9.0.

Status as of 2026-11-12: **Stage 1 (mismatches), Stage 2 (indels), per-locus
selection and the packed reference (`.osr`) are built and validated; the
mismatch path is identical to bowtie1 `-v 3` on the full Ensembl human cDNA**
(§4), and full-reference timing is measured (§6), including `.osr` loading in
Docker (~1 s). The public release tree is packaged (§10), and bowtie2's
recall shortfall is measured: **93.8% at ≤4 edits, 58.0% at ≤6, and
non-monotonic in the budget**, with oligoscan at **100.0%** bulge recall
against bowtie2's 70.2% (§8a). **Repeated on a published library** (2,361
Huesken guides): identical to bowtie1 at **517,754 = 517,754 sites, zero
discordance**, and bowtie2 at **75.0%** (§8b). Measurement C was re-run
twice on 2026-11-12 with a fixed harness: all five library gates are
**byte-identical across independent passes**, and one earlier claim is
**corrected — at <=3 MM on a full library bowtie1 is 9-14x FASTER than
oligoscan** (8.0-12.7 s vs 112-113 s); oligoscan's case is flatness in k
and in the indel budget, not throughput at k=3 (§8b). The manuscript is
drafted with every number measured. Next: §3.4 and the abstract need that
correction, then the open items in §9a (affiliation, repo URL/DOI, two
figures, a 2+-bulge benchmark stratum), then the GPU path (§9).

This document consolidates everything discussed and built so far, so a new
project can start from it without re-deriving anything.

---

## 1. Where this came from

The original proposal: build "bowtie3" / "bowtie_oligo" — raise bowtie1's
mismatch cap from 3 to 6 over the whole oligo, integrate bowtie1 + bowtie2 so
indels up to 2 edits are covered, publish to GitHub for other researchers, then
adopt it in the siRNA/ASO off-target platform.

Clarified later, and this is the actual scope: **not a source merge of bowtie1
and bowtie2 — a complete rewrite purpose-built for short-oligo (<50 nt)
alignment against both cDNA and whole-genome references.**

### Premise check (done, with sources)

| claim | verdict |
|---|---|
| bowtie1 unmaintained | **Substantially right.** v1.3.1, last push **May 2022**; no commits or PRs in the past year; 39 open issues. ~4.5 years dormant, not 10 — and it is the one the exhaustive `-v` path depends on. |
| bowtie2 unmaintained | **Not correct.** **2.5.5 shipped 16 Feb 2026**, after 2.5.4 (May 2024) and 2.5.3 (Jan 2024). Actively maintained. Its limit is that multiseed is a heuristic **by design** — measured plateauing at 84.9% recall with `-L 7`, `-D 500`, `-R 20` all failing to move it. Upstream maintenance would not change that. |

### Why not fork bowtie

- bowtie1's 3-MM cap is **not arbitrary**. `-v` is exhaustive FM-index
  backtracking and its cost grows combinatorially in k; raising it to 6 makes
  bowtie1 unusably slow rather than merely slower.
- bowtie2's recall gap is architectural, so gluing the two together inherits
  **both** failure modes and leaves a C++ fork to maintain and re-validate
  against upstream forever.

---

## 2. The algorithmic argument

A 21-mer carries 42 bits; locating one position in 560 Mb needs ≈29.7 bits. At
4–6 mismatches only 15–17 bases remain informative, so **every filtering scheme
loses specificity exactly where this problem lives**:

- **q-gram lemma** — a match with ≤k errors shares ≥ L−q+1−kq common q-grams.
  For L=21 this is *negative* for any useful q (q=5, k=4 → −3). No filter
  exists.
- **Pigeonhole partition** — k+1 blocks, one of which must match exactly.
  Selectivity against 560 Mb needs blocks ≥ ~13 nt, i.e. reads of
  L ≥ 13(k+1) = 65 nt at k=4. A 21-mer cannot supply it.

With nothing to filter *with*, the honest method is to examine every position —
and that turns out to be cheap: one byte per base, inner loop a byte compare
against a scalar (one AVX2/AVX-512 instruction per 32–64 positions), tiled so
the L passes stay in L2 and the reference streams once rather than L times.

What that buys **structurally**, not by tuning:

- **no mismatch cap** — k is a comparison at the end, not a search budget
- **exhaustive by construction** — recall is 100% by definition
- **runtime independent of k** and of reference repetitiveness
- no `-k` truncation, no `-L` seed length, no `--gbar`, no `-D`/`-R` effort
  limits — an entire class of bug disappears

---

## 3. What is built

Source on the device at `C:\WebUI\oligoscan\` (moved 2026-11-11 from
`C:\WebUI\siRNA_OffTarget_Bowtie2_092026\oligoscan\`, where the originals
remain). Layout: `src/` (scanner; `src/archive/` frozen builds), `tests/`,
`scripts/`, `data/` (host mount for reference + indexes, git- and
docker-ignored), `Dockerfile` + `docker-compose.yml` (ubuntu 24.04 with
bowtie1, bowtie2 and oligoscan; the build runs `--selftest` and the selection
suite and fails on any failure — builds and passes on Jeff's machine), `scripts/build_indexes.sh`
(bowtie1 + bowtie2 indexes with provenance).

| file | purpose |
|---|---|
| `src/oligoscan.c` | the scanner: mismatches, indels, per-locus selection, `pack` |
| `tests/test_oligoscan_select.py` | self-contained suite, brute-force oracles, no bowtie or external data |
| `tests/validate_oligoscan.py` | bowtie1 `-v 3` equivalence harness |
| `tests/example_oligos.tsv` | 40 real TRIM63 21-mers |
| `scripts/build_indexes.sh` | bowtie1 + bowtie2 indexes, and `oligoscan pack`, with provenance |
| `scripts/bench_gates.sh` | runtime across gates (§6) |
| `scripts/bench_paper.py` | head-to-head vs bowtie1/bowtie2 for the manuscript (§8a) |
| `examples/` | 54 kb demo reference with sites of known kind + `make_demo.py` |
| `README.md`, `BUILD.md`, `LICENSE`, `CITATION.cff` | release files |
| `src/archive/` | frozen pre-selection build (`86a9ea4d`) and the superseded Stage 2 suite — **local only** |

```
oligoscan <ref.fa|ref.osr> <oligos.tsv> <maxSeedMM> <maxNonseedMM> <maxTotalEdits> \
          [anti|sense|both] [maxIndelBases=0] [maxIndelEvents=1] [--all-descriptions]
oligoscan pack ref.fa ref.osr   # one-time; ~20x faster loading, identical output
oligoscan --selftest            # numerics + selection check, no inputs needed
oligoscan --version
```

### Stage 1 — mismatches

- **Transcript-boundary safety is structural**, not a check: transcripts are
  concatenated with `PAD = 64` sentinel bases that match nothing, so any window
  straddling a junction accrues mismatches on every sentinel column and cannot
  pass the gate (the indel path also checks the junction explicitly).
- **Structured gate**: mismatches counted in three **disjoint** regions of the
  antisense oligo — position 1 (5′ terminal base), positions 2–8 (seed),
  positions 9..L (body) — gated on seed ≤ S, non-seed ≤ N, total ≤ T. Keeping
  them disjoint is the point: a seed-perfect site with 5 body mismatches and a
  seed-mangled site with 5 body matches share a total and share nothing of the
  risk.
- **Orientation handled once**: the region map is built in *query-column*
  space, so the scan never reasons about strand again — historically where every
  off-by-one in this problem has lived.
- Output TSV: oligo, strand, seq_id, 0-based position within the sequence,
  seed_mm, nonseed_mm, total_mm, gap_events, oligo_bulge, target_bulge, class,
  edit_total.

### Oligo length — verified 2026-11-12

**Supported range is 10–64 nt, not 21.** Each oligo is scanned at its own
length and lengths may be mixed freely in one file: verified exact against a
numpy brute force at 10, 12, 15, 19, 21, 25, 30, 40, 50, 60 and 64 nt
simultaneously (447 = 447 sites, 0 missed, 0 extra), and each per-length run
is identical to that oligo's rows in the mixed run. 64 is Myers' one-word
limit; 10 is information-theoretic (20 bits against a transcriptome's ~30).

**Three silent input failures found and fixed the same day.** Each returned
results for a *different oligo than the user supplied*, with no indication:

| input | old behaviour | now |
|---|---|---|
| >64 nt | silently truncated to the first 64 bases — a 70-mer and its 64-nt prefix gave byte-identical output | skipped, with a reason on stderr |
| <10 nt | silently dropped | skipped, with a reason |
| degenerate base (N, R, Y…) | the character was **deleted**, shortening the oligo and shifting every base after it | skipped, with a reason |

A summary line now counts the skips, the accepted length range is reported
(`11 queries, 10-64 nt`), and an all-oligos-skipped run exits non-zero
instead of reporting zero sites. This matters for a published tool: a
truncated or shifted query produces a plausible-looking site list that is
simply wrong, and nothing downstream can detect it.

*Caveat for ASOs*: the seed region (positions 2–8) is an siRNA/RISC concept.
For a chemistry where it does not apply, leave `maxSeedMM` unconstrained and
gate on body and total instead. Noted in the README.

### Stage 2 — indels

**Myers' bit-vector edit distance** (JACM 46(3), 1999) gives the edit distance
of the best alignment of the whole pattern ending at every text position, in
O(1) 64-bit words per position for patterns ≤64 nt. Exact and exhaustive, so it
is a **complete candidate filter with no recall loss** — unlike a seeded
extender, which is precisely where bowtie2 loses its 15%.

Candidates are rare, so each then gets a small **banded DP with traceback** to
recover the actual alignment: which oligo positions mismatch (so seed and body
are attributed separately), and where the gaps are and which way round — pattern
base unpaired = **oligo bulge**, text base unpaired = **target bulge**.

---

## 4. Validation

### Mismatch path — byte-identical to bowtie1 `-v 3`

bowtie1's `-v` mode is exhaustive backtracking within its cap, so it is a true
oracle there — and the only one available, since nothing else is exhaustive
above 3 MM.

80.5 Mb reference with repeat families, paralogs and low-complexity runs;
40 real TRIM63 21-mers:

```
bowtie1   : 21 sites
oligoscan : 21 sites
only bowtie1 : 0 · only oligoscan : 0 · shared : 21
mismatch-count disagreements: 0
PASS — identical to bowtie1 -v 3
```

A set comparison, not a count. Independent cross-check: scaling site counts to
560 Mb gives 3.5 / 70 / 715 / 6,050 per oligo at ≤3/4/5/6 MM against the
platform's own documented estimates of 5 / 67 / 696 / 5,733.

### Full human cDNA — identical to bowtie1 `-v 3` (2026-11-12)

Ensembl `Homo_sapiens.GRCh38.cdna.all.fa` (328,868 transcripts, 743 Mb of
sequence; md5 `64ca420f…`), 40 TRIM63 21-mers, seed≤3 / body≤3 / total≤3,
antisense, 22 threads, in Docker (`tests/validate_oligoscan.py`):

```
bowtie1   : 1,917 sites   7.6 s (incl. index load)
oligoscan : 1,917 sites  21.4 s = 19.7 s FASTA parse + 1.70 s scan
shared 1,917 · only bowtie1 0 · only oligoscan 0 · count differs 0
PASS
```

The scan is ~4.5× faster than bowtie1's whole run even at ≤3 MM (bowtie1's
best case); >90% of oligoscan's wall time is parsing the 800 MB FASTA, which
makes the `pack` step (§9) the next real speed win. Linear scaling gives ~75 s
mismatch-only for 1,750 oligos on 22 threads. No disagreement involved a
non-ACGT reference base (the reference has ~100k), so the N-handling
difference (oligoscan counts N as a mismatch; bowtie1 should not report alignments over
ambiguous reference bases) did not arise in this set. Index builds:
bowtie-build 517 s / 4.5 GB peak, bowtie2-build 548 s / 4.8 GB.

### Indel path — planted truth

bowtie1 cannot do gapped alignment at all, so there is no external oracle;
truth is planted. 48 sites = 12 guide positions (9…L−1) × {oligo, target}
bulge × {0, 1} extra mismatches, plus 8 perfect and 8 random controls, in 3 Mb.

| check | result |
|---|---|
| every planted bulge found | **48/48** |
| each planted site is **one** gap event | yes |
| seed clean where mismatches were planted outside it | yes |
| 1-nt sites still found at `maxIndel=2` | yes |
| bulge direction correct | 47/48 — the exception is a terminal representational ambiguity (see §5) |

### Two bugs the self-test caught, and would otherwise have hidden

A transcription error in the Myers recurrence is invisible in a hit count, so
`--selftest` checks it against a brute-force O(mn) DP **before** anything else
runs. It earned its keep immediately: the first build disagreed with the DP on
**3,531 of 4,000** random trials.

**Bug 1 — the horizontal shift-in.** The code had `(HP << 1) | 1`, which is the
**global** alignment boundary condition. Semi-global pattern-in-text search has
the first DP row all zeros, so the delta shifted into row 0 is **0**, not +1.
Rather than keep re-deriving it, all 32 recurrence variants were enumerated
against the DP:

```
every variant shifting in 0  -> EXACT MATCH (0 disagreements / 930,098 columns)
every variant shifting in 1  -> fails
```

The `Eq`-vs-`X` and `D0`-vs-`X` choices turn out not to matter; the shift-in is
the whole bug. The canonical Hyyrö form is used, with the reason in a comment
so it cannot be "simplified" back.

**Bug 2 — the traceback window was one band too long.** The window must *end*
at the alignment end, because traceback requires the optimum to finish in the
last column. `wn = L + 2*band` left `band` trailing text bases that had to be
consumed as insertions, inflating every cost by exactly that much. Corrected to
`wn = L + band`, window spanning `[endpos+1-L-band, endpos+1)`.

After the fix:

```
Myers bit-vector vs brute-force O(mn) DP        ok  (0/4000 disagreements)
traceback cost vs DP, region split sums total   ok  (0/2000 disagreements)
```

---

## 5. The one real finding — no per-locus selection layer

Three apparent test failures turned out to have **one** root cause. The scanner
emits a row per candidate *end position*, carrying whatever alignment the DP
prefers there, with **nothing deciding between competing descriptions of the
same locus**. Symptoms:

1. A target bulge at guide position 19 reported as an oligo bulge — equal cost,
   both valid, a terminal representational ambiguity.
2. Perfect sites *also* reported as inferior gapped alignments — within the
   gate, but worse descriptions of a locus already explained better.
3. Sites appearing to go "missing" between the mismatch and indel paths.

Symptom 3 is worth recording precisely, because it looked like a recall bug and
is not. Every "lost" site is present **one base away**, described more cheaply:

```
siRNA_0015 BG0037 59453   mismatch path: 4 MM
    indel path 59452  target_bulge  mm=2 gaps=1 edit=3
siRNA_0021 BG0101 291709  mismatch path: 4 MM
    indel path 291710 oligo_bulge   mm=1 gaps=1 edit=2
```

The DP correctly found a **better** explanation of the same locus. Nothing is
lost; the start position and the class change, so a naive key comparison reads
it as a loss.

**The fix is a layer the production platform already has**, learned the hard
way: best-explanation-per-locus with a **total** tie-break, plus suppression of
a gapped alignment when a strictly better ungapped one exists at the same
locus. In the platform this is `verify_bulge_candidates`' `best` map,
`_better_ungapped`, and the total-rank fix of 2026-10-27 that made the choice
deterministic across batch sizes. The design is converging on the same answer
from a different direction, which is reassuring.

---

## 5a. Per-locus selection — built 2026-11-11

- At every Myers candidate end, **both** the gapped traceback and the plain
  ungapped placement are offered, so a traceback tie can never hide an
  ungapped description.
- All gate-passing descriptions are collected; a global post-pass keeps the
  best per locus by greedy non-maximum suppression in a **total** rank order:
  edit total → gap bases → gap events → seed MM → pos-1 MM → body MM → start →
  end. Equal cost goes to the ungapped description.
- Same locus = same oligo, same sequence, start **or** end within
  R = max(1, maxIndel).
- Selection runs **after** the gate (a locus is reported if any description
  passes — errs toward sensitivity). *Open decision: confirm this matches the
  platform's `verify_bulge_candidates` semantics.*
- Global post-pass ⇒ tile boundaries cannot split a locus, and output is
  sorted and byte-identical across thread counts. Also removed a latent
  overflow (indel path flushed its 1 MB buffer only once per tile).
- Mismatch-only mode (`maxIndel = 0`) is untouched and still reports every
  offset, so bowtie1 `-v` equivalence is unaffected.
- `--all-descriptions` disables selection for debugging.
- **Tie policy — decided 2026-11-12 (Jeff): prefer the ungapped description
  on equal cost.** This is the opposite of the deployed platform, whose
  `_better_ungapped` suppresses a bulge only when an ungapped alignment is
  *strictly* better (it keeps the bulge on ties). Edit totals are identical
  either way — only the class label differs — so summary edit counts are
  unaffected. Consequence for any comparison with the platform: a platform
  bulge at a tie (typically guide positions 18–20) appears in oligoscan as
  `mismatch_only` at the same `edit_total`; count it as covered.

**Validation** (`test_oligoscan_select.py`, 2.1 Mb synthetic reference, all pass):

| check | result |
|---|---|
| mismatch path vs numpy brute-force Hamming | 166 = 166, counts identical |
| planted bulges: rows per locus | exactly 1 for all 48 + 8 perfect controls |
| reported edit vs gap-constrained brute-force DP optimum | 0 differ |
| 300 random sites (0–4 MM, 0–2 indel bases): within-budget loci reported | 256/256, all at the optimum cost |
| two winners within R | none |
| mismatch-path sites covered by an indel winner at ≤ cost | 171/171 (the §5 "missing site" symptom) |
| site split across a junction | not reported |
| 1 vs 2 threads | byte-identical |

The guide-position-19 bulges now come out as `mismatch_only` at equal cost
— a 1-nt bulge that near the 3′ terminus is equally a single terminal
mismatch — which is the §5 symptom 1 ambiguity, now resolved deterministically.

Cost: none measurable — 20 Mb × 40 oligos, maxIndel 2: 1.66 s vs 1.69 s
before; memory unchanged. Descriptions collapse ~3.3× into loci.

---

## 6. Performance

### Benchmark machine (confirmed 2026-11-12)

**Intel Core Ultra 9 285** (Arrow Lake-S, 3 nm): 24 cores — 8 performance +
16 efficient — 24 threads, no SMT, 2.5 GHz base / 5.6 GHz turbo, 36 MB L3.
128 GB RAM. Windows + Docker Desktop, which exposed 22 of the 24 threads to
the container; all measurements used those 22.

**This CPU has no AVX-512** (`/proc/cpuinfo`: `avx avx2 avx_vnni bmi2 popcnt
sse4_2`; Arrow Lake dropped AVX-512). The mismatch scan therefore vectorised
to **32 bytes per instruction, not 64**, and the measured figures are a
floor rather than a ceiling: a part with AVX-512 (Xeon, EPYC, Zen 5) should
go faster on the same code. Worth stating in the methods section — it
forecloses "you benchmarked on unusually wide hardware".

### Measured on the full Ensembl human cDNA (2026-11-12)

743 Mb of sequence, 40 TRIM63 21-mers, antisense, 22 threads, Docker on
Windows (`scripts/bench_gates.sh`):

| gate (seed / body / total, maxIndel / events) | sites or loci | scan |
|---|---|---|
| 3 / 3 / 3 | 1,917 | 1.73 s |
| 1 / 3 / 4 | 7,607 | 1.71 s |
| 0 / 6 / 6 (seed-perfect) | 47,218 | 1.68 s |
| 1 / 5 / 6 | 207,262 | 1.68 s |
| unconstrained / 6 | 681,388 | 1.81 s |
| 1 / 3 / 4, indel 1 / 1 | 28,996 | 5.24 s |
| 1 / 3 / 4, indel 2 / 1 | 33,923 | 5.16 s |
| 1 / 3 / 4, indel 2 / 2 | 63,289 | 5.19 s |

Flat across the gate and across the indel budget, confirmed at full scale.
Projected for 1,750 oligos on 22 threads: ~75 s mismatch-only, ~230 s with
indels. Site counts here are **per transcript**: Ensembl cDNA holds many
isoforms of each gene, so one genomic site is counted once per isoform that
contains it — a gene-level collapse belongs in reporting, not in the scanner.

**Load time was the real cost**: ~20 s per run to parse the 800 MB FASTA
inside Docker, against 2.3 s raw read / 4.4 s parse in a plain Linux VM. Cause:
stdio's default 4 KB buffer turns the file into ~200k small reads, each slow
through Docker Desktop's Windows bind mount (bowtie1 reads its index in large
blocks, so it is unaffected). Fixed two ways — a 16 MB read buffer, and the
packed `.osr` reference (207.5 MB vs 838.8 MB). **In Docker, loading from
`.osr` takes 0.9–1.1 s instead of ~20 s**, with identical site counts at all
8 gates. End to end at ≤3 MM on **40 oligos**, oligoscan is now ~2.8 s
against bowtie1's 7.6 s. **This inverts at library scale and must not be
quoted as a general result**: bowtie1's time here is almost all fixed index
load, so on 2,361 guides bowtie1 is 8.0-12.7 s against oligoscan's
112-113 s (§8b). Scan times rose slightly at the
loosest gate (2.08 s vs 1.81 s at 681k sites) because mismatch output is now
sorted; negligible at realistic gates.

### Measured, 80.5 Mb × 40 queries, 2 cores

| gate | sites | per oligo | scan time |
|---|---|---|---|
| seed ≤0, nonseed ≤3, total ≤3 | 2 | 0.1 | **2.17 s** |
| seed ≤1, nonseed ≤3, total ≤4 | 175 | 4.4 | 2.2 s |
| seed ≤0, nonseed ≤6, total ≤6 | 2,216 | 55.4 | 2.2 s |
| seed ≤1, nonseed ≤5, total ≤6 | 9,819 | 245.5 | 2.26 s |
| unconstrained seed, total ≤6 | 34,583 | 864.6 | **2.40 s** |

**Runtime is flat across the gate** — 2.17 s at 3 MM to 2.40 s at 6 MM, and the
growth is output writing, not search. This is the property bowtie structurally
cannot have.

| mode | 80.5 Mb × 40 queries, 2 cores |
|---|---|
| mismatches only | **2.31 s** |
| maxIndel = 1 | 9.04 s |
| maxIndel = 2 | **9.05 s** |

Myers costs ~4× the byte-compare scan and is **flat in the indel budget**:
1 nt and 2 nt cost the same.

### Extrapolated to production (1,750 oligos, 20 cores)

| reference | mismatch-only | with indels |
|---|---|---|
| human cDNA 560 Mb | ≈66 s | ≈275 s |
| GRCh38 3.1 Gb | ≈366 s | — |

Runtime is independent of reference repeat content — the property that matters
for GRCh38, where ~50% is repeat and every seeded method degrades. Memory is
one byte per base (3.1 GB for GRCh38); 2-bit packing would quarter it.

---

## 7. The gate is what makes 6 MM usable

At ≤6 MM a 21-mer keeps 15 matched bases — about the information needed to
locate *one* position in 560 Mb — so an unconstrained 6-MM gate sits at the
noise floor by construction: ~865 sites/oligo on 80 Mb, ~6,000 on 560 Mb,
**10.5 million across 1,750 oligos**. No aligner makes those real; a better one
just finds the noise faster.

The seed constraint is what recovers signal, measured at total ≤6:

| seed gate | sites/oligo | vs unconstrained |
|---|---|---|
| unconstrained | 864.6 | 1.0× |
| ≤3 | 804.4 | 1.1× |
| ≤2 | 602.8 | 1.4× |
| ≤1 | 245.5 | 3.5× |
| **≤0 (seed-perfect)** | **55.4** | **15.6× smaller** |

So the deliverable is not "the cap is now 6" — it is **"the cap is gone, and
the gate is what you choose"**.

---

## 8. Fresh motivating evidence (2026-11-11)

A production indel run put hard numbers on why the current architecture needs
replacing, not patching. Settings: MAX MM ALONGSIDE AN INDEL = 2,
MAX INDEL LENGTH = 2.

```
var_oligo_bulge.fa    1900 queries, 17–18 nt  →      17 MB of hits
var_target_bulge.fa   2000 queries, 12–13 nt  →  14,219 MB of hits
```

Same query count, **800× the output**. Cause: bowtie1 charges each `N` in a
target-bulge variant as a mismatch, so

```
need = MAX MM ALONGSIDE AN INDEL + MAX INDEL LENGTH = 2 + 2 = 4  >  -v cap of 3
```

and the pass falls back to 12-nt fragment queries. A 12-mer carries 24 bits
where 560 Mb needs ~29.7, so it cannot constrain the search: ~140 million
candidate alignments for one batch of 100 oligos, which Python then rejects one
at a time. 30 minutes in, that batch had produced zero output rows.

**Three separate workarounds in the deployed module exist for the same two
limitations**, none of them biological:

| workaround | exists because |
|---|---|
| fragment-then-verify for mismatches | the gate needs `-v 4+`, bowtie1 caps at 3 |
| `N`-substituted bulge variants | bowtie1 cannot do gaps; a gap becomes a forced mismatch against the same capped budget |
| a second aligner for the sub-sequence screen | historical; bowtie2's *gapped* path was measured and retired (below) |

And the measured reason bowtie2 was abandoned for bulges — from
`indel_offtarget.py`, 720 planted bulge sites, 1-nt bulge in a 21-mer, recall by
guide position:

```
                    p9    p12   p16   p19
--gbar 4 (default)  100%   99%   86%   41%
--gbar 1            100%   99%   87%   54%
```

`--gbar` defaults to 4, categorically forbidding a gap within 4 nt of either
read end, which silently excluded the 3′ end. Forcing `--gbar 1` moved overall
recall 81.4% → 84.9%, then it plateaus. oligoscan's Myers pass is exhaustive by
construction: 100%.

**oligoscan collapses all three workarounds into one pass with no cap and no
cliff.**

---

## 8a. Head-to-head against bowtie2 — measured 2026-11-12

Full Ensembl human cDNA, 10 TRIM63 21-mers, antisense, 22 threads on a Core
Ultra 9 285 (§6), Docker (`scripts/bench_paper.py`, run
`paper_20261112_164554`). bowtie2 in its best
honest configuration: flat edit-distance scoring, `--gbar 1`, `-a` not `-k`,
`-N 1`, `-i C,1,0`, `-D 200 -R 10`, and the provably-largest seed length.

| edit budget | oligoscan sites | bowtie2 alignments | **bowtie2 recall** | oligoscan | bowtie2 |
|---|---|---|---|---|---|
| ≤4 | 5,600 | 5,967 | **93.8%** | 1.5 s | 1,444.5 s |
| ≤6 | 252,831 | 213,691 | **58.0%** | 2.5 s | 1,709.6 s |

So bowtie2 misses **6.2% of sites at ≤4 edits and 42.0% at ≤6**, while taking
**963× and 684× longer end to end** (wall clock, oligoscan including its
reference load).

Three things make this stronger than the bare percentages:

1. **It is a conservative bound.** oligoscan ran mismatch-only here while
   bowtie2 was allowed gaps — 12.1% of its budget-4 alignments and 32.9% of
   its budget-6 alignments carry an indel. bowtie2 therefore searched a
   *strictly larger* space than oligoscan and still missed this much; its
   true mismatch-only recall can only be the same or worse.
2. **bowtie2 misses sites inside bowtie1's exhaustive range.** 4 of the
   budget-4 misses and 18 of the budget-6 misses are sites of only 3 edits —
   sites `bowtie -v 3` finds by construction, and which §4 confirms
   oligoscan finds exactly.
3. **Its recall is non-monotonic in the search budget.** Of the 5,600 sites
   within ≤4 edits, bowtie2 covered 5,254 at budget 4 but only 4,671 at
   budget 6: **widening the budget lost 735 sites it had already found** and
   gained 152, a net loss of 583. Only 1 of the edit-3 misses is common to
   both budgets. So bowtie2's miss set is not a fixed, characterisable
   deficit that a wider search repairs — it *moves with the parameters*.
   This is the sharpest available statement of the architectural point, and
   far more damning than a flat recall figure: you cannot tune your way out,
   and you cannot even say in advance what you will lose.

Re-confirmation at ≤3 edits (measurement C): bowtie1 752 sites, oligoscan
752 sites, 7.7 s vs 1.4 s — a second independent equivalence check on a
different oligo subset.

### B. 1-nt bulge recall — final, measured 2026-11-12 (run `paper_20261112_182117`)

702 verified planted sites (10 guides × guide positions 9–20 × {oligo,
target} bulge × {0,1,2} further substitutions) in 19,851 real cDNA
transcripts / 41 Mb, 18 further candidates excluded at construction as
ambiguous. Gate: ≤2 substitutions alongside a 1-nt bulge. Loose-gate sanity
pass: 702/702 reachable.

| recall by guide position | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | **all** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| oligoscan | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | **100.0%** |
| bowtie2 `--gbar 1` | 67 | 67 | 67 | 68 | 67 | 67 | 67 | 70 | 75 | 70 | 74 | 91 | **70.2%** |
| bowtie2 `--gbar 4` (default) | 67 | 67 | 67 | 68 | 67 | 67 | 67 | 70 | 70 | 60 | 72 | 91 | **68.8%** |

**The dominant variable is error load, not guide position** — and this is
the finding to build the figure on:

| recall by further substitutions alongside the bulge | 0 | 1 | 2 |
|---|---|---|---|
| oligoscan | 100% | 100% | 100% |
| bowtie2 `--gbar 1` | 100% | 100% | **11%** |
| bowtie2 `--gbar 4` | 100% | 96% | **11%** |

bowtie2 is *complete* on a clean bulge and collapses to **11%** once two
substitutions accompany it. That is 3 edits in a 21-mer — exactly the point
§2 predicts from first principles: 18 of 21 bases informative leaves no
window long enough to seed on, so no `-L`, `-D` or `-R` setting recovers it.
Across guide position bowtie2 is nearly flat at 67–75% (rising to 91% at
position 20), so `--gbar` is a *secondary* effect worth only 1.4 points
overall. This reframes §8's historical p19 numbers: the gap is not mainly
the gap barrier, it is total error load.

oligoscan is **100.0% (702/702)**, at every guide position and every error
load.

Calls in the unplanted real-cDNA sequences: 45, every one at `edit_total` 3,
i.e. sitting exactly at the gate boundary. These are **not** false
positives — a 21-mer within this gate genuinely occurs in 41 Mb of real
transcript sequence, and the count is the expected order of magnitude from
§6's site densities.

### Measurement B had to be redesigned — twice

The first bulge-recall design scored **both** tools at 100% and so
discriminated nothing: it planted clean bulges with no other error, one
guide, 24 sites, a small reference — bowtie2's effort limits never bound.
Worse, 4 of 24 "misses" by oligoscan turned out to be a **planting** fault,
not a tool fault: in a homopolymer run a substitution can make the sequence
re-describable more cheaply, so "one 1-nt bulge + two substitutions"
(3 edits) actually read as "two separate 1-nt bulges" (2 edits). oligoscan
correctly reported the cheaper description; the planted "truth" was not the
truth of the sequence.

B now follows the platform's own 720-site design — N guides × guide
positions 9…L−1 × {oligo, target} bulge × {0,1,2} further substitutions,
in a real-cDNA background, gate ≤2 substitutions alongside a 1-nt bulge —
and **verifies each planted site against an independent O(mn) DP**, dropping
any candidate whose intended description is not the optimum (logged to
`B_excluded.tsv` with its reason). A loose-gate sanity pass must find every
verified site before any stratum is reported. bowtie2 is run at both
`--gbar 1` (tuned) and `--gbar 4` (its default, which categorically forbids
a gap within 4 nt of either read end), so the figure separates "bowtie2
misconfigured" from "bowtie2 cannot".

**Second fault, found in the real-reference run (2026-11-12).** The one
site oligoscan "missed", `PL_siRNA_0010_oligo_bulge_13_2`, is reported at a
loose gate as *one oligo bulge + one target bulge, 0 mismatches, 2 edits* —
cheaper than the intended 3. The verification had checked each planted core
**in isolation**, but the core is embedded in 300 nt of real cDNA on each
side and the alignment can **borrow a base from the flank** to reach that
cheaper description. A core that is unambiguous on its own can therefore
acquire a cheaper one once embedded. Fixed: flanks are now chosen *before*
verification and the optimum is checked over core + 8 nt of real flank on
each side. On synthetic data the fixed check gives 360/360 with no
exclusions.

Both faults ran the same way: oligoscan found the genuinely cheapest
description and the benchmark's "truth" was wrong. Worth stating in the
manuscript's methods — a planted-truth benchmark that does not verify its
own truth in context will score a correct exhaustive tool as incomplete.

**Re-run done** (`paper_20261112_182117`): with flank-aware verification,
702 verified sites, 18 excluded, oligoscan **100.0%**. The figures above are
final.

---

## 8b. Published library (Huesken) — measured 2026-11-12

2,361 unique 19-nt antisense guides, full Ensembl human cDNA, 22 threads
(runs `paper_20261112_184456` and `_192516`).

**Exactness at library scale — the strongest validation result so far.**
Against `bowtie -v 3` over all 2,361 guides and all 328,868 transcripts:

```
bowtie1   : 517,754 sites
oligoscan : 517,754 sites
only bowtie1 0 · only oligoscan 0 · shared 517,754
```

Half a million sites, two independent implementations, zero discordance in
either direction. oligoscan: 189.16 s scan + 1.1 s load. (bowtie1 reported
2,360 of 2,361 reads aligning; one guide has no ≤3-edit site.)

**bowtie2 recall is worse on this library than on the 21-mers.** 25-guide
seeded random subsample, budget 4, `-L 7` (the bound for a 19-mer):

| | |
|---|---|
| oligoscan | 52,303 sites, 1.24 s scan |
| bowtie2 | 45,765 alignments, ~46 min |
| **bowtie2 recall** | **75.05%** |

Of 13,047 missed sites: 12,085 at 4 edits, 942 at 3, and **20 at only 2
edits**. On the 21-mers the shallowest miss was at 3 edits; on a shorter
published library it drops sites at 2.

**Throughput, whole library against whole transcriptome — re-measured
twice 2026-11-12** (runs `paper_20261112_201840` and `_204208`, fixed
harness, 22 threads). Both passes are reported because the same job had
previously measured 136 s and 189 s an hour apart; the tight re-measurement
shows those earlier figures were contaminated by machine load.

| gate | sites | sites/oligo | scan, pass 1 / pass 2 |
|---|---|---|---|
| seed ≤3, body ≤3, total ≤3 | 517,754 | 219 | 110.8 / 112.3 s |
| seed ≤1, body ≤3, total ≤4 | 2,248,742 | 952 | 131.0 / 127.1 s |
| seed-perfect, total ≤6 | 10,280,009 | 4,354 | 146.5 / 143.6 s |
| ≤3 edits + 1 nt indel | 10,646,838 | 4,510 | 398.8 / 392.8 s |
| ≤3 edits + 2 nt indel | 22,475,487 | 9,520 | 429.2 / 415.0 s |

**Determinism at library scale.** All five outputs are **byte-identical
between the two passes**, the 22.5-million-row indel file included. This is
the strongest form of the §5a reproducibility claim: it holds across
processes, not just across thread counts.

**bowtie1, the number the crash had cost:** **12.7 s (pass 1) and 8.0 s
(pass 2)** for all 2,361 guides and 517,754 sites — the spread is index
load against a warm or cold page cache. Reconstructed beforehand from file
mtimes as 12.5 s / 7.6 s, which the timers confirmed to within 0.4 s; the
mtime method is sound for recovering a lost run.

**This corrects a claim.** At ≤3 MM on a full library bowtie1 is
**9-14× faster** than oligoscan (8.0-12.7 s vs 112-113 s). The earlier
"faster even at bowtie1's best case" reading came from 40 oligos, where
bowtie1's fixed ~7 s index load dominates a tiny query set; scale the
library 59× and bowtie1 barely moves while oligoscan scales linearly in
query count. oligoscan's case was never throughput at k=3, where bowtie1
is exhaustive, mature and faster. It is:

1. **Flat in k.** Sites rise 19.9× from the ≤3 gate to seed-perfect ≤6
   (517,754 → 10,280,009) while scan time rises 1.28× (112 → 144 s).
   Fitting the two ends gives a fixed search cost of ~110 s and ~3.2 µs per
   site written: **the search is flat and the output is not** — now
   quantified rather than asserted. bowtie1 cannot run any of these gates.
2. **Flat in the indel budget.** 1 nt → 2 nt doubles the sites (10.6M →
   22.5M) for **6% more time** (393 → 415 s). The indel path costs 3.5× the
   mismatch path (111 → 393 s), matching the ~4× measured at 80 Mb.
3. **Able to run at all.** At budget 4 bowtie1 is absent by construction and
   bowtie2 took ~46 min for 25 guides at 75% recall — ~72 h extrapolated to
   2,361 against oligoscan's 128 s at 100%.

**Gates do not transfer between libraries.** Seed-perfect ≤6 returns 4,354
sites/oligo on these 19-mers and ≤3 edits + 2 nt indel returns 9,520 — at
or near the §7 noise floor. A 19-mer carries 38 bits against a 21-mer's 42,
so a gate that is selective at 21 nt is loose at 19. Every gate table in the
manuscript must name its library, and none may be compared across them.

### What the crash cost, and the lesson

Both runs died when Docker Desktop's engine wedged (`docker ps` → HTTP 500).
Cause: **my harness ran the unconstrained-seed ≤6 gate on the full library,
writing a 14.2 GB TSV per run — 28 GB concurrently through a Windows bind
mount.** §7 of this document already said that gate is at the noise floor
and unusable; running it at library scale was avoidable.

Fixed in `bench_paper.py`: the full-library loop now uses gates someone
would actually set (seed ≤1/body ≤3/total ≤4, and seed-perfect ≤6), with a
warning if any single output exceeds 2 GB.

**The lost number has been recovered** — see the re-measurement above.
The root cause was a harness design fault worth naming: `runtime` was
accumulated in memory and `C_runtime.tsv` written **once, at the end of the
run**, so a crash three-quarters of the way through lost *every* timing,
including a 46-minute bowtie2 measurement. Fixed: a `rt()` helper now
rewrites `C_runtime.tsv` after every single measurement. Any harness that
measures expensive things should persist each result as it is produced.

---

## 8c. Open bug — 2-nt bulges split into two gap events (found 2026-11-12)

**A quarter of genuine single 2-nt bulges are silently dropped at
`maxIndelEvents=1`**, the setting the production pipeline uses.

40 independent 2-nt oligo bulges (adjacent deletions) planted in random
sequence, each verified by an independent DP to have a minimum-cost
alignment with exactly **one** gap run:

| build | found at maxI=2 maxE=1 | reported as 2 events |
|---|---|---|
| 0.9.0 as shipped | **30 / 40** | 10 |
| patched | **40 / 40** | 0 |

**Cause** — `traceback()` in `src/oligoscan.c`. When several moves are
optimal it tries the diagonal first, so a free match can interrupt an open
gap and split one 2-base run into two 1-base events. The edit total is
identical, so nothing looks wrong; the site is simply rejected at the gate
(`a->gap_events > c->maxE`). No 1-event alternative is offered, because the
traceback yields one alignment per candidate end.

**Fix** (written and validated in sandbox, NOT yet applied): on a tie prefer
*continuing* an open run before taking the diagonal — the description with
fewest gap events, which is also what per-locus selection already ranks
highest. Costs are untouched: only moves already equal to `D[j][i]` are
taken. `--selftest` 0/4000 and 0/2000; the full selection suite passes
identically on both builds.

**The suite does not catch this** — check C plants 0-2 indel bases but runs
at `maxEvents=2`, where a split run still passes. Needs a case at
`maxEvents=1` with a genuine 2-base run. Same hole as the missing 2+-bulge
benchmark stratum (§9a item 6), now with a real defect sitting in it.

Verified working alongside it, so the scope of the bug is exactly the event
count: one 2-nt target bulge, one 2-nt oligo bulge, two separate 1-nt
bulges, two 1-nt bulges + 3 mismatches, and one 2-nt bulge + 3 mismatches
are all found at the DP optimum. Note `edit_total` counts gap bases, so
3 mismatches + a 2-nt bulge needs `maxTotalEdits 5`, not 3.

---

## 8d. Measurement B's 11% was a scoring artifact (found 2026-11-13)

**§3.3's central claim is void and must be re-measured.** bowtie2 charges a
gap of length N as `open + N*extend`, so with `--rdg 1,1 --rfg 1,1` a 1-nt
bulge costs **2**, while oligoscan counts it as **1**:

| planted site | oligoscan edits | bowtie2 cost | within B's budget 3? |
|---|---|---|---|
| bulge only | 1 | 2 | yes |
| bulge + 1 substitution | 2 | 3 | yes, at the limit |
| bulge + 2 substitutions | 3 | **4** | **no — forbidden** |

So the stratum where bowtie2 "collapses to 11%" is the stratum its
`--score-min C,-3,0` excluded by construction. Confirmed in isolation by
`bench_three.py --selfcheck-only`, which plants three sites of known cost in
a toy reference:

```
--rdg 0,1 budget 3: mm3=y, bulge0=y, bulge2=y   <- usable
--rdg 1,1 budget 3: mm3=y, bulge0=y, bulge2=NO
--rdg 1,1 budget 4: mm3=y, bulge0=y, bulge2=y   <- usable
```

**Contaminated**: both 11% cells; the 70.2% / 68.8% totals (means over three
strata, one impossible); the by-guide-position table (pools all strata); and
the conclusion that "the dominant variable is error load, not guide
position", which rests entirely on the excluded stratum. The planned central
figure (bulge recall against error load) cannot be drawn until B is re-run.
The abstract does not carry this number, so it is unaffected.

**Fixed**: `bt2_cmd` in `bench_paper.py` now uses `--rdg 0,1 --rfg 0,1`, so a
length-N gap costs exactly N and the script's "flat edit-distance scoring"
claim is true rather than approximate. bowtie2 accepts a zero open penalty.
Measurement A's existing numbers were taken under the old scoring; cheaper
gaps can only give bowtie2 more ways to cover a site, so its recall there
will move up slightly on a re-run.

### The pattern worth naming

This is the second configuration error found on 2026-11-13, after the `-L`
bound (§8e). Both ran the same way: a bowtie2 parameter was set from a rule
reasoned out on paper, never tested against what the tool actually did, and
the error made bowtie2 look worse. The fix in both cases raised bowtie2's
score. **Any claim that bowtie2 was given "its best honest configuration"
has to rest on measured behaviour, not derivation** -- hence the startup
self-check in `bench_three.py`, which verifies the scoring against planted
sites of known cost before any comparison runs, and exits rather than
produce a comparison on a broken footing.

---

## 8e. The bowtie2 seed-length bound is the wrong choice (found 2026-11-13)

100 Huesken guides, budget 4, identical oligoscan truth set of 204,213 sites;
only `-L` changed:

| bowtie2 at budget 4 | `-L 7` (the provable bound) | `-L 8` (forced) |
|---|---|---|
| alignments | 179,881 | 231,045 |
| **recall** | **75.46%** | **97.3%** |
| missed at 2 edits | 52 | **0** |
| missed at 3 edits | 2,501 | 72 |
| missed at 4 edits | 47,560 | 5,511 |
| wall clock | 3,087 s | **2,391 s** |

One extra base of seed is worth **22 points of recall and runs 23% faster**.
`seed_len_bound` returns the longest `-L` that can *provably* nucleate a
k-error alignment (7 for a 19-mer at k=4); `-L 8` is theoretically unsound
yet finds nearly everything the sound setting misses.

**Mechanism**: with `-i C,1,0` seeding at every offset, a 7-nt seed (14 bits)
against 743 Mb (~30 bits needed) yields enormous numbers of seed hits, and
`-D 200 -R 10` exhausts itself extending noise before reaching the true
sites. The longer seed is less exhaustive in principle, far more selective in
practice. The speedup corroborates it: less wasted work, not more work.

**Consequence for the paper.** Every bowtie2 figure was taken at the bound,
and §3.2 claims bowtie2 was given "its best honest configuration". That claim
is false as it stands. Either sweep `-L` and report bowtie2's *best*, or drop
the claim. The finding survives either way and gets sharper: **bowtie2's
recall swings 75% -> 97% on one parameter, the value its own theory
prescribes is the worse one, and nothing in the tool tells a user which to
pick.** An aligner that must be tuned blindly, per read length and per
budget, to a value contradicting its documented bound, is not one to trust
for off-target enumeration.

**Gapped fractions recounted** on the 100-guide library (the manuscript's
12.1% / 32.9% came from the TRIM63 tiling set): budget 4 `-L 7` **15.2%**,
budget 4 `-L 8` **14.8%**, budget 6 `-L 7` **34.9%**. The conservativeness
argument survives intact.

---

## 8f. `example_oligos.tsv` is a tiling walk, not a library (found 2026-11-13)

All 40 "TRIM63 21-mers" are consecutive 21-nt windows sliding 1 nt at a time
across ~60 nt of one transcript: `0002 = C + 0001[:-1]`, `0016[2:] ==
0014[:-2]`, and so on throughout. **They are one sequence sampled at 40
offsets, not 40 queries.**

Surfaced by a count mismatch: `paper_20261112_164554` reported 252,831 sites
at total <=6 where `paper_20261113_124832` reported 172,920 for what looked
like the same run. Cause was benign -- the harness changed from "first N" to
a seeded random subsample, so the two drew different oligos -- but the old
draw was `siRNA_0001..0010`, ten *consecutive* windows, effectively **n = 1**,
which is why its site count ran 47% above the file average.

**Consequence**: every 21-mer result -- §3.1 exactness, §3.2's 93.8% / 58.0%,
§3.4's runtime table -- comes from 40 overlapping windows of one 60-nt region
of one gene. **Huesken must become the primary library.** TRIM63 stays only
where correlation is harmless: §3.1's bowtie1 set-equality check.

---

## 8g. GPU backend — built and validated 2026-11-13

**oligoscan-gpu 0.2.0-mismatch: byte-identical to the CPU build and 11.6-12.0x
faster on the scan.** 1,000 Huesken guides, 19 nt, full Ensembl human cDNA,
20 CPU threads / RTX 4000 Ada (run `three_20261113_195946`):

| gate | guides | sites | CPU scan | GPU kernel | kernel | end to end | bowtie1 |
|---|---|---|---|---|---|---|---|
| total <=1 | 1,000 | 15,932 | 50.56 s | **4.29 s** | 11.8x | 8.7x | 7.8 s |
| total <=2 | 1,000 | 33,058 | 50.64 s | **4.35 s** | 11.6x | 8.5x | 8.0 s |
| total <=3 | 1,000 | 233,898 | 51.38 s | **4.29 s** | 12.0x | 8.2x | 8.1 s |
| total <=4 | 2,361 | 4,795,910 | 129.16 s | **10.10 s** | 12.8x | 6.0x | cannot express |
| total <=5 | 2,361 | 37,288,299 | 194.31 s | **27.17 s** | 7.2x | 1.9x | cannot express |

Wide gates are the full library (run `three_20261113_202715`), narrow ones
1,000 guides alongside bowtie1 (`three_20261113_195946`). **Every row is
byte-identical to the CPU build**, the <=5 row only after the ambiguous-base
fix below.

**Exactness now has three independent witnesses** — an FM-index backtracker
(bowtie1), a CPU byte-compare scan, and a CUDA popcount kernel -- agreeing on
233,898 sites with no shared code between them.

**The GPU kernel is LESS flat than the CPU scan**, which is worth stating
plainly rather than claiming the property transfers. Over the narrow gates it
holds (4.29 / 4.35 / 4.29 s while sites rise 15x), but from <=4 to <=5 the
site count rises 7.8x, the CPU scan 1.5x and the GPU kernel **2.7x**. Cause:
every hit is appended through a single global `atomicAdd`, so emission cost
sits inside the kernel, where on the CPU counting and writing are separate
passes. Per-block counters plus a prefix sum would fix it; until then the
kernel's cost grows with the number of hits.

**The bottleneck has left the GPU altogether at wide gates.** At <=5 the
kernel is 27 s of a 102.6 s run; the other ~75 s is host-side -- copying 37M
hits back, sorting, sentinel filtering, and `printf` of 2.26 GB. End to end
that is 1.9x while the kernel is 7.2x. The next real optimisation is
therefore the OUTPUT PATH (wider buffers, parallel formatting, or a binary
format plus a converter), not a better kernel. Same conclusion §6 reached
about the CPU from the other direction: once the search is cheap, what
remains is writing the answer down.

**bowtie1's last advantage is nearly gone, but state it precisely.** End to
end at 1,000 guides the GPU wins (5.9-6.4 s vs 7.8-8.1 s), but both are
fixed-cost dominated at opposite ends: bowtie1 ~= 7.8 s + ~0/oligo (index
load), GPU ~= 1.6 s + 4.3 ms/oligo (reference load + host->device copy, then
scan). **They cross near 1,470 oligos**; above that bowtie1 leads again at
<=3, and on the full 2,361-guide library it would be ~8 s against the GPU's
~12 s. The claim for the paper is therefore *parity with bowtie1 in the one
regime bowtie1 can express, with no cap and flat cost in k* — not a speed win.

### The indel path on GPU — full library, exact, 14.2x (2026-11-16)

**The Myers kernel is done and it is the bigger win.** Full 2,361-guide
Huesken library, 19 nt antisense, whole Ensembl human cDNA, gate total <=3
edits / 3 indel bases / 3 gap events, `oligoscan-cuda --gpu` against the CPU
build:

| | sites | selection | scan | throughput |
|---|---|---|---|---|
| CPU | 2,487,939 | 0.84 s | 357.69 s | 5.04 Gbase-compare/s |
| GPU | 2,487,939 | 0.78 s | **25.17 s** | **71.68 Gbase-compare/s** |

**Output files byte-identical** (`cmp`, whole file, not a count check -- the
discipline §8g's region-accounting bug earned). Selection is host-side in
both and is unchanged, as expected.

**14.2x, against 11.6-12.8x on the mismatch path, and the asymmetry is
structural rather than lucky.** Myers' bit-vector does ~20 ALU ops per
64-bit word of pattern-profile against one text column; the mismatch kernel
does three masked popcounts per window. The indel kernel therefore carries
far more arithmetic per byte loaded, sits further from the memory wall, and
runs closer to the `__popc`/ALU ceiling. Two consequences worth keeping:

- The larger speedup landed on the path that costs ~3.5x more on CPU. The
  indel path is where the runtime lives, so this is the useful end.
- The mismatch kernel's 33-36% of peak is a *memory* ceiling, not an
  instruction one. Shared-memory tiling and multiple positions per thread
  should still help it; the indel kernel has less headroom left.

**This is the scan phase, not end to end.** The run reported no total wall
time, so no end-to-end figure is claimed here. At this gate the CPU's
output path was ~106 s of 483 s (§8g), and that cost is untouched by the
kernel -- the output path remains the next real optimisation.

Retires the "Myers kernel for the indel path" item below.

### Ambiguous bases: a 4-row deficit in 37 million (found and fixed 2026-11-13)

`.osr` records sentinel runs for **two** different things and the CPU treats
them differently: inter-transcript **padding** (a window overlapping it must
not be reported -- the structural junction safety) and **non-ACGT reference
bases** (~100k; the CPU charges each as a mismatch and still reports the
site). The first GPU filter discarded both, losing **4 sites of 37,288,299**
at total<=5 -- 0.00001%, invisible to any count-based check, and genuinely
wrong: real off-target sites silently absent from a user's report.

Padding lies between transcripts and ambiguous bases inside one, so the
transcript table separates them. There is also a second half: a packed
ambiguous base reads as 0 (= 'A'), so the kernel can *undercount* mismatches
there and pass a site the CPU rejects. Both directions are handled by
recomputing any hit whose window touches an ambiguous run -- charging each
ambiguous position as a mismatch in its own region -- and re-applying the
gate on the host. Such hits are rare, so the cost is nil.

Found only because the full library at a wide gate was diffed byte for byte
and a 4-row discrepancy was chased rather than shrugged at. `<=4` passed;
only `<=5` was loose enough to admit a site with an `N` in the window.

### Design notes that cost real time to find

- **Sentinels are stored as 0 in `.osr`, indistinguishable from `A`.**
  `pack_ref` writes `r.b[i] == SENT ? 0 : r.b[i]` and records the runs
  separately; the CPU re-inserts them when expanding to 1 byte/base, which is
  why its junction safety is structural. A kernel reading packed data
  directly cannot tell padding from adenines, so a window spanning a junction
  can score as a near-perfect match. Padding-as-A can only make a window look
  *better*, so it yields false positives and never false negatives: the host
  discards any hit whose window overlaps a sentinel run, which is exactly
  equivalent to masking on the device and costs nothing because candidates
  are rare. The CPU's separate check that a window not extend past its
  transcript's end is also replicated — padding makes it redundant except for
  the final transcript.
- **The three regions are disjoint and `N` bounds the BODY alone.** The CPU
  gate is `if (cs > S || cn > N || cs + cn + c1 > T)`, so the 5′ terminal base
  (region 0) is constrained only through the total, and `nonseed_mm` in the
  output is the body alone. The first kernel folded position 1 into
  "non-seed", which changed both the gate and that column. The site sets were
  already identical (35,133 = 35,133, same positions) and only one column
  differed — **a count-based check would have passed it.** Caught by
  `tests/test_gpu_equivalence.py`, which diffs whole files.
- Packing is little-endian within the byte (base `i` at bits `(i&3)*2`), so a
  64-bit load yields 32 bases with base `k` at bits `2k..2k+1`, which is what
  the XOR/popcount primitive wants. Windows at arbitrary offsets are a shift
  of the containing word combined with the next.

### Measured budget and remaining headroom

The card reports 48 SMs at 2.17 GHz. At ~16 `__popc` per SM per clock the
ceiling is 1.67e12 popcounts/s, or **5.57e11 windows/s** at three masked
popcounts each, against the CPU's measured 1.61e10 — a **~35x ceiling**.
At 11.6-12.8x the kernel is near 33-36% of peak, which is what an unoptimised
first version should give: one thread per position, uncoalesced overlapping
window loads, no shared-memory staging, no ILP. Tiling the reference through
shared memory and giving each thread several positions should recover a good
part of the gap. **Note the gate dependence**: at the widest library-scale
gates ~106 s of the 483 s is output writing, which no kernel helps, so the
realised speedup falls with gate width (≈24x at <=2 edits, ≈4x at <=4). Any
GPU claim in the paper needs that qualification or a reviewer reproducing it
at a wide gate will contradict it.

### Scope and layout

`src/oligoscan_gpu.cu` (mismatch scan only, oligos 10-32 nt),
`Dockerfile.gpu` (nvidia/cuda:12.6.2-devel, also carries the CPU build and
bowtie so the three-way benchmark runs in one image), an `oligoscan-gpu`
service in `docker-compose.yml`, and `tests/test_gpu_equivalence.py`. Same
tree deliberately: the CPU build is the oracle, and a separate repository
would guarantee drift. `bench_three.py --gpu oligoscan-gpu` runs the GPU on
every mismatch-only gate and **diffs its output against the CPU's**, so a GPU
time is never reported without confirmation that the run was correct.
The Myers kernel for the indel path landed 2026-11-16 (see above) and is
validated on the full library. Still to do: per-block counters plus a prefix
sum to replace the global `atomicAdd`, the host-side output path (the real
bottleneck at wide gates), direct 2-bit scanning to quarter scan
memory, and a multi-stage image (`nvcc -cudart static` into
`nvidia/cuda:*-base`) to cut ~7 GB to ~300 MB for release.

---

## 8h. RIsearch2 — measured 2026-11-16 (run `risearch_20261116_133212`)

100 Huesken guides, full Ensembl human cDNA, 22 threads, RIsearch2 **v2.1**
(`scripts/bench_risearch.py`).

**Provenance matters here.** The binary is the PUBLISHED C tool (Alkan et al.,
NAR 2017, 45(8):e60), GPL3, pulled from the authors' own conda channel
(`rthtools`, md5 `82b021a622ef6dc9ea187242c0154e9c`) and pinned in the
Dockerfile. **`github.com/RTH-tools/risearch` is no longer that tool**: it now
hosts a Rust rewrite, 3.0.0-alpha.5, different first authors, BUSL-1.1, whose
citation is a paper *in preparation*. Benchmarking that and citing Alkan 2017
would be wrong twice. `VERSIONS` in the image records the version and md5.

### Why this needed its own harness

bowtie1/bowtie2 and oligoscan all answer "where does this oligo align within
k edits". RIsearch2 answers "where could this oligo form a favourable duplex".
Three differences, each handled explicitly rather than averaged over:

- **G:U wobble.** RIsearch2 pairs G:U by default, so a "perfect 7-mer seed"
  can be almost all wobbles — the first test run returned `PWWWWPP...` noise.
  oligoscan counts G:U as a mismatch. `--noGUseed` is mandatory, and even then
  wobbles still score during extension, so the tools stay non-identical.
- **Local vs full-length.** oligoscan aligns the WHOLE oligo and reports an
  edit count; RIsearch2 reports the energy-optimal LOCAL duplex. Verified on a
  planted site: a 19/19 alignment with 2 mismatches came back as a 16-nt
  duplex with 1 unpaired base at shifted coordinates. Matching is therefore at
  LOCUS level (transcript + >= L/2 overlap), which is **generous to
  RIsearch2** — a stricter rule lowers its recall, never raises it.
- **Energy, not edits.** `-e` governs what is reported and there is no
  principled value for enumeration, so it is swept, not chosen.

Coordinates verified against planted sites: RIsearch2 target positions are
1-based inclusive, oligoscan's `pos` is 0-based. Only `+` strand hits count
('-' means the guide binds the transcript's reverse complement, not an mRNA
target).

### Runtime — RIsearch2 is fast, and the paper must not pretend otherwise

| | oligoscan | RIsearch2 | ratio |
|---|---|---|---|
| seed-perfect, total <=3, no indels | **5.23 s** | ~45 s | 8.6x |
| seed-perfect, total <=4, 2 nt indel | **14.58 s** | ~45 s | 3.2x |

RIsearch2 is **flat in the energy threshold** (44.56 / 45.98 / 46.32 / 46.20 /
47.00 s from -30 to -10) while its output grows 600x. Flatness in a search
parameter is therefore NOT unique to oligoscan and must not be claimed as
such; oligoscan's flatness is in the EDIT BUDGET, which is a different axis.
RIsearch2's 45 s excludes a one-time 291.6 s index build producing an **11.9
GB** suffix array (16.00 bytes/base, measured); oligoscan needs no index.

**This kills the speed argument against RIsearch2.** 3-9x is not 500x. The
case against it is completeness, and only completeness.

### Recall of oligoscan's site set, by edit count

Gate seed-perfect / total <=3 / no indels — 6,686 sites
(per edit_total: 0:1287  1:80  2:708  3:4611):

| -e | 0 | 1 | 2 | 3 | **all** | RIsearch2-only |
|---|---|---|---|---|---|---|
| -30 | 60.8% | 55.0% | 43.2% | 15.1% | **27.3%** | 6,669 |
| -25 | 99.4% | 88.8% | 65.4% | 43.6% | **57.2%** | 95,720 |
| -20 | 100% | 95.0% | 89.5% | 72.9% | **80.2%** | 700,555 |
| -15 | 100% | 100% | 99.7% | 90.7% | **93.6%** | 2,599,571 |
| -10 | 100% | 100% | 100% | 99.1% | **99.4%** | 5,271,710 |

Gate seed-perfect / total <=4 / 2 nt indel, 1 event — 196,433 sites
(per edit_total: 0:1287  1:113  2:1846  3:21807  4:171380):

| -e | 0 | 1 | 2 | 3 | 4 | **all** | RIsearch2-only |
|---|---|---|---|---|---|---|---|
| -30 | 60.8% | 38.9% | 27.0% | 6.3% | 1.2% | **2.5%** | 3,681 |
| -25 | 99.4% | 62.8% | 49.5% | 21.1% | 9.2% | **11.5%** | 76,836 |
| -20 | 100% | 73.5% | 64.4% | 37.9% | 24.9% | **27.2%** | 652,043 |
| -15 | 100% | 86.7% | 70.7% | 53.5% | 40.2% | **42.4%** | 2,521,758 |
| -10 | 100% | 86.7% | 71.1% | 57.3% | 48.2% | **49.8%** | 5,178,169 |

The edit_total=1 column rests on 80 and 113 sites respectively; do not
over-read it.

### Three findings, in order of how much they matter

**1. An energy cutoff discriminates by GC content, not by complementarity.**
At `-e -30`, RIsearch2 reports **60.8% of PERFECTLY complementary sites** --
505 of 1,287 missed. Nothing is wrong with its search: those are AU-rich
19-mers whose duplexes never reach -30 kcal/mol. An off-target screen built on
an energy cutoff therefore carries a **systematic compositional bias**, and no
amount of tuning removes it -- lowering the threshold to admit AU-rich true
sites admits everything else too. This is the cleanest criticism available and
does not depend on any recall percentage.

**2. You cannot loosen your way to completeness.** On the gapped gate,
-20 -> -10 recovers ~44,300 sites and adds **4,526,126** unmatched
interactions: about **102 extra interactions per additional true site**. The
recall curve flattens (+15.2 then +7.4 points) exactly as the noise explodes.
Architecturally the same conclusion as bowtie2's non-monotonicity (8a), from
the opposite direction.

**3. At its own default (`-e -20`) on a realistic gapped gate, RIsearch2
reports 27.2% of the sites in the budget.** Even fully relaxed it reaches
49.8%.

### Open, before any of this is published

- ~~Re-run at the authors' recommended siRNA pipeline settings~~ —
  **decided 2026-11-16 (Jeff): not doing this, and the reasoning holds up.**
  The parameter that would have been the real exposure is `-e`, and it was
  not chosen but SWEPT from -30 to -10 with the default row reported, so no
  single threshold carries the result; at the most generous setting RIsearch2
  still reaches 49.8%. The headline finding — 39.2% of PERFECTLY
  complementary sites missing at -30 — is a property of energy cutoffs and
  does not depend on tuning at all. And `-s 2:8/7` is not an arbitrary guess
  the way `-L 7` (8e) and `--rdg 1,1` (8d) were: it mirrors the seed-perfect
  gate the compared sites were drawn from, so both tools are required to
  match the seed exactly. That is a principled choice, not a reasoned-out one.
  **What was done instead:** one sentence of disclosure in the Results saying
  the threshold was swept rather than chosen, that the seed mirrors the gate,
  and that the authors distribute an siRNA pipeline whose settings are not the
  ones used here. Costs nothing and converts a possible ambush into a stated
  choice. Reopen only if a referee asks.
- Locus-level matching is weaker than the exact-coordinate matching used
  against bowtie2 (3.6). Say so in the methods rather than letting a reviewer
  find it.
- One library, one length (19 nt), 100 guides, `-l 20` fixed.

---

## 8i. The 71 bowtie2-only sites were an oligoscan bug — found and fixed 2026-11-16 (0.9.3)

§3.6's 100-guide run left 6,382 bowtie2-only alignments at `total<=4 + 4 nt
indel`. Most are a per-locus-merge artifact: 7,008 of 7,085 (independent
recount) sit within 8 nt of an oligoscan site, so the two tools found the same
locus and drew its boundary differently. **71 were not.** Those 71 had no
oligoscan site anywhere on the transcript — not a different description of the
same site, a missing site. A paper whose central claim is exhaustiveness could
not carry them.

They were real, and oligoscan was wrong.

### What they were

All 71 share one shape:

| property | value |
|---|---|
| oligoscan sites on that (oligo, transcript) | **0** |
| bowtie2 `NM` | 4 (62), 3 (9) |
| true semi-global edit distance (independent DP) | **3 or 4 — all within gate** |
| strand | reverse, 71/71 |
| alignment start, transcript coordinate | **1 (65), 2 (5), 3 (1)** |
| ambiguous reference bases (`XN`) | 0 in all 71 |
| distinct oligos / transcripts | 29 / 63 |

Every one is an alignment flush against a transcript's 5′ end whose **first
oligo base is unpaired**.

### The mechanism

Two guards, each correct alone, with a gap between them.

`on_candidate()` ran the traceback over a window extending `L + band` bases
back from the alignment's end position. For an alignment starting at
transcript position 1 that window reaches ~5 bases into the sentinel padding
ahead of the transcript. `traceback()` returns **one** description — the
lexicographic minimum of (cost, gap bases, gap events). So at a single end
position two descriptions could tie:

| description | start | cost | gap bases |
|---|---|---|---|
| ungapped, spans the junction | sentinel padding | 3 | **0** ← wins the tie |
| real, inside the transcript, 1 nt oligo bulge | transcript position 0 | 3 | 1 |

The ungapped one wins on gap bases. `push_aln()` then rejects it because
`tx_at(start) < 0` — its start is in padding, which is exactly the junction
safety the sentinels exist to provide. The real description was never
generated. **The site is lost, not merely described differently.**

Instrumented, on a two-record reproducer (`Hu_0998` × `ENST00000429864.1`):

```
[DROP] start=367 endpos=385 in SENTINEL: mm=3 gapev=0 ob=0 tb=0 total_edits=3
```

Three edits, well inside a gate of four, discarded.

### Why the mismatch-only path was never affected

The `maxI == 0` fast path screens on the alignment *start* (`tx_at(abspos) <
0` → skip) and never consults a traceback, so no tie-break can mislead it. An
ungapped 19-mer window that spans a junction genuinely is not a site, so
skipping it is right. **§3.1's exact agreement with `bowtie -v 3` across
517,754 sites is therefore untouched by this bug and by its fix** — that
comparison runs at `maxI = 0`.

### The fix

Clamp the traceback window to the transcript containing the end position, so
every description the traceback can return is one `push_aln()` will accept:

```c
long ti_end = tx_at(c->r, endpos);
if (ti_end < 0) return;                    /* endpos is itself padding */
long tx0 = (long)c->r->tx[ti_end].start;
if (wstart < tx0) wstart = tx0;
int wn = (int)((long)endpos + 1 - wstart);
```

This also subsumes the old `if (wstart < 0) return;`, which silently dropped
candidates near the start of the reference array. The GPU backend routes
through the same `on_candidate()` (§8g), so one change covers both and the
byte-identical property is preserved.

### Validation

- `--selftest`: 0/4000 and 0/2000 disagreements, unchanged.
- Reproducer: site now reported, `pos 0, 2 mismatches, 1 nt oligo bulge,
  edit_total 3` — matching an independent Python DP exactly.
- **All 71 recovered, 0 pairs lost**, on a 64-record reference holding the 63
  implicated transcripts plus a leading decoy (so the padding is faithful):
  original binary recovered 1/71, fixed binary 71/71, and the fix added 70
  (oligo, transcript) pairs while removing none.
- `tests/test_junction_edge.py` — 30 planted sites (oligo bulge at pattern
  positions 1, 2, 3, 10, 19 × transcript offsets 0, 1, 2, 5, 40, plus ungapped
  controls). Passes on 0.9.3, fails on 0.9.2 with exactly one missing case:
  `BULGE_p1_off0` — bulge at pattern position 1, site at transcript position 0.
  Wired into the Dockerfile, so the image will not build if it regresses.

### Re-run confirmed — 0.9.3, same bowtie2 SAMs, 2026-11-16

`scripts/rerun_after_fix.cmd`. **Orphans 0 / 0 / 0**, and the result is
stronger than a closed hole:

| gate | oligoscan sites | bowtie2 | bt2 recall, exact | bt2 recall, locus | bt2-only, locus |
|---|---|---|---|---|---|
| total <=2 + 2 nt | 5,773 | 5,734 | 99.32% | 99.32% | **0** |
| total <=3 + 3 nt | 101,981 | 99,189 | 97.17% | 97.26% | **0** |
| total <=4 + 4 nt | 1,552,873 | 1,358,868 | 87.05% | 87.33% | **0** |

**Every bowtie2 alignment is now contained in oligoscan's site set at locus
level — 1,358,868 of 1,358,868 at the widest gate, precision 100.00% at all
three.** That is the exhaustiveness claim carrying direct evidence instead of
an argument from the algorithm, and it is the single most useful sentence the
benchmark has produced. Before the fix the same criterion left 9 and 71
alignments unmatched, so the claim could not have been made.

What moved:

| | 0.9.2 | 0.9.3 |
|---|---|---|
| sites, total <=3 | 101,960 | **101,981** (+21) |
| sites, total <=4 | 1,552,674 | **1,552,873** (+199) |
| bt2 recall (exact), <=3 | 97.18% | **97.17%** |
| bt2 recall (exact), <=4 | 87.06% | **87.05%** |
| oligoscan-only (exact), <=4 | 200,891 | **201,022** |
| bt2-only (exact), <=4 | 7,085 | **7,017** |
| bt2-only (locus), <=4 | 77 | **0** |

Note 71 orphan *alignments* became 199 extra *sites* and 9 became 21: a
recovered locus often carries several descriptions.

Scan times 14.05 / 15.26 / 17.23 s, so flatness is unchanged and marginally
better than published: **a 269x rise in sites costs 1.23x more time** (was
stated as 1.28x). bowtie2's times are reused unchanged at 2460.1 / 2381.2 /
2411.8 s, i.e. 140-175x oligoscan's.

### Every number that must change before submission

- Abstract: 87.1% -> **87.0%**; 1.28-fold -> **1.23-fold**.
- §3.6 table and the same table in `README.md`: all three gates.
- Figure 2 and its legend: bar heights, the 200,891 -> **201,022** shortfall,
  13% -> **12.9%**, 1.28x -> **1.23x**, and 137-178x -> **140-175x**.
- Add the locus-criterion result, which did not exist before: it is a better
  headline than the recall percentage.
- The exact-criterion bt2-only residual is now **7,017 at <=4 and 96 at <=3**,
  all of them boundary differences on a shared locus, none of them missing
  sites. Say so explicitly; it is the question the 71 made a referee ask.

### The pattern worth naming

This is the third bug of the same family (§8c, §8d, this one): the cost was
right and the *description* was wrong, in a way that then changed what passed
a filter. Myers gives a cost; everything downstream depends on the traceback
turning that cost into the right alignment, and a tie-break that looks like a
cosmetic preference becomes a correctness bug the moment a later stage filters
on the description. Both earlier instances were found by cross-checking
against another tool, and so was this one. The selftests check the recurrence
against brute force — they cannot see a site that was never generated.

---

## 8j. CPU vs GPU across oligo count and gate — measured 2026-11-16 (`gpu_scaling`)

`scripts/bench_gpu_scaling.py`, one grid in the GPU image so every cell shares
a machine, a reference and an I/O path. oligoscan 0.9.3, 22 CPU threads,
RTX 4000 Ada. **All 12 cells byte-identical, CPU against GPU, md5 of the whole
output.**

| oligos | gate | sites | CPU scan | GPU scan | scan × | CPU wall | GPU wall | wall × |
|---|---|---|---|---|---|---|---|---|
| 100 | total ≤3 | 20,068 | 4.23 s | 0.44 s | 9.6 | 5.18 s | 1.80 s | 2.9 |
| 100 | total ≤4 | 192,077 | 4.61 s | 0.43 s | 10.7 | 5.61 s | 2.04 s | 2.8 |
| 100 | ≤3 + 3 nt | 94,102 | 14.76 s | 1.54 s | 9.6 | 15.62 s | 2.47 s | 6.3 |
| 100 | ≤4 + 4 nt | 1,510,546 | 17.80 s | 6.23 s | **2.9** | 18.79 s | 7.31 s | 2.6 |
| 1,000 | total ≤3 | 235,006 | 44.91 s | 4.29 s | 10.5 | 45.85 s | 6.31 s | 7.3 |
| 1,000 | total ≤4 | 2,050,337 | 52.92 s | 4.33 s | 12.2 | 53.78 s | 10.00 s | 5.4 |
| 1,000 | ≤3 + 3 nt | 1,052,765 | 171.93 s | 12.56 s | 13.7 | 172.76 s | 13.49 s | 12.8 |
| 1,000 | ≤4 + 4 nt | 15,666,765 | 203.06 s | 62.49 s | **3.2** | 203.98 s | 63.51 s | 3.2 |
| 2,361 | total ≤3 | 517,754 | 116.38 s | 10.12 s | 11.5 | 117.25 s | 13.20 s | 8.9 |
| 2,361 | total ≤4 | 4,795,910 | 133.18 s | 10.14 s | 13.1 | 134.11 s | 21.88 s | 6.1 |
| 2,361 | ≤3 + 3 nt | 2,488,108 | 415.71 s | 28.80 s | 14.4 | 416.53 s | 29.84 s | 14.0 |
| 2,361 | ≤4 + 4 nt | 37,591,667 | 485.42 s | 152.18 s | **3.2** | 486.34 s | 153.19 s | 3.2 |

### 1. The speedup is not one number, and the paper must stop quoting it as one

Scan speedup runs **9.6× to 14.4× at three of the four gates and collapses to
2.9–3.2× at ≤4 + 4 nt**, at every oligo count. That is not noise and it is not
the output path — it is inside the kernel. Every hit is appended through a
single global `atomicAdd`, so emission cost lives in the kernel, and at
≤4 + 4 nt the site counts (1.5M, 15.7M, 37.6M) make emission the dominant
term. On the CPU, counting and writing are separate passes, so the CPU degrades
far more gently.

**Any GPU claim in the paper must name its gate.** A reviewer who reproduces
"14×" at a wide gate gets 3× and is right to say so. The per-block counters
plus prefix sum already noted as to-do in §8g are no longer a nicety: they
target precisely the cell where the speedup collapses.

### 2. Scan speedup rises with oligo count; wall-clock speedup is a different number

Per-launch fixed cost amortises, so scan speedup climbs with the library:
9.6 → 10.5 → 11.5 at total ≤3, and 9.6 → 13.7 → 14.4 at ≤3 + 3 nt. It is
converging near 13–14×, so the full library is where the GPU earns its place.

Wall clock tells a blunter story at small scale: at 100 oligos the ≤3 gate is
**9.6× on the scan but 2.9× end to end**, because ~1 s of reference load and
process start is most of a 5 s run. Quote the scan figure to describe the
kernel and the wall figure to describe a user's experience, and never swap
them.

### 3. The two GPU binaries do not define "scan" the same way — do not compare the column across them

CPU wall minus scan is a steady ~0.9 s (reference load) in every row, and so is
`oligoscan-cuda --gpu`'s. But `oligoscan-gpu`'s grows with hit count: 1.4 s,
1.6 s, 5.7 s, 11.7 s at 20k, 192k, 2.05M, 4.80M sites — about 2.4 µs per site.
Its host-side copy-back, sort and sentinel filter sit **outside** its scan
timer, where `oligoscan-cuda`'s equivalent work sits inside. Hence the odd row
at 2,361 × total ≤4: 13.1× on scan but 6.1× end to end.

**Wall clock is comparable across every row in this table. Scan time is not
comparable between the mismatch and indel backends.** This should be stated
wherever the table appears, or the mismatch rows will read as better than they
are.

### 4. Two incidental confirmations

- **2,361 × total ≤3 gives 517,754 sites**, the same number `bowtie -v 3`
  produced, now reproduced by a third implementation.
- **2,361 × ≤3 + 3 nt gives 2,488,108 sites, where 2026-11-13 measured
  2,487,939** — 169 more, the 0.9.3 junction-edge fix (§8i) showing up at
  library scale, on the GPU as well as the CPU, exactly as predicted from the
  shared `on_candidate()`. The CPU scan also reads 415.71 s against 357.69 s
  then; the 14.2× figure in §8g is superseded by **14.4×** here, measured on
  one machine in one run rather than across two.

**Supersedes** the scattered GPU timings in §8g for anything the paper quotes:
those were three runs at different oligo counts and gates, and this is one grid.

---

## 8k. Are the missed sites real? — measured 2026-11-16 (`dg_20261116`)

The question §3.2 invites and could not previously answer: are the ~197,000
sites bowtie2 does not report worth reporting, or is exhaustive search
enumerating noise. Answered as a population comparison, not an anecdote.
`scripts/dg_missed_vs_shared.py`, 100 guides at `total<=4 + 4 nt indel`,
oligoscan 0.9.3, ViennaRNA `RNA.duplexfold` (Turner 2004).

### Part 1 — structure, and it needs no thermodynamics at all

Classifying the 1,552,873 sites against the bowtie2 SAM at locus level:

| | shared with bowtie2 | **oligoscan only** |
|---|---|---|
| n | 1,356,223 | 196,650 |
| intact RISC seed (positions 2–8) | 38.6% | **53.2%** |
| carries a bulge | 86.3% | **97.2%** |
| intact seed **and** a bulge | 36.5% | **53.2%** |
| by edit total | 2:4,316 3:93,480 4:1,257,008 | 2:38 3:2,728 4:193,884 |

**104,580 of the missed sites carry a perfectly intact seed, and the missed
set is ENRICHED for intact seeds relative to the shared set — 53.2% against
38.6%.** These are not the irrelevant tail; they disproportionately carry the
one feature most associated with miRNA-like off-target repression.

The two seed rows are identical because **every one of those 104,580 sites
also carries a bulge.** That is the mechanism, observed rather than argued: a
bulge breaks the exactly matching stretch, so a seeded aligner has nothing
long enough to nucleate on. It is the Introduction's pigeonhole argument,
measured. This result costs nothing — it is in the scan output already.

### Part 2 — thermodynamics, and the first answer was the inconvenient one

Matched sampling, strata on **(guide, edit_total)** so neither guide GC
content nor error load can explain a difference; 180 matched strata,
1,845 sites per group:

| group | n | mean dG | median | 10th pct | strongest | ≤ −25 | ≤ −30 |
|---|---|---|---|---|---|---|---|
| shared | 1,845 | −20.62 | −20.40 | −28.20 | −36.30 | 22.2% | 5.9% |
| oligoscan only | 1,845 | **−18.58** | −18.30 | −26.40 | −35.00 | 13.4% | 3.8% |

**The missed sites are weaker**, by +2.04 kcal/mol, ~0.36 pooled SD,
Mann-Whitney p = 3.6e-27. Reported as measured; the convenient result would
have been no difference.

**But most of that is composition, not biology.** The missed set is far more
heavily gapped (3–4 gap bases: 45% against 18%), and a bulge costs energy.
Matched on gap load as well, the difference largely dissolves:

| gap bases | n shared | n only | shared dG | only dG | difference | p |
|---|---|---|---|---|---|---|
| 1 | 626 | 307 | −21.48 | −20.95 | +0.53 | 0.13 (n.s.) |
| 2 | 629 | 695 | −19.63 | −18.60 | +1.03 | 8e-4 |
| 3 | 287 | 659 | −19.32 | −18.00 | +1.32 | 1.4e-4 |
| 4 | 44 | 161 | −17.49 | −16.28 | +1.21 | 0.22 (n.s.) |

At one gap base there is **no significant difference at all**. So the missed
sites are weaker chiefly *because they are gapped* — which is a restatement of
why a seeded aligner misses them, not independent evidence that they do not
matter. A residual of 0.5–1.3 kcal/mol remains within a gap class, plausibly
because bowtie2's own scoring correlates with stability and it finds the
better alignments within a class first.

### Part 3 — absolute numbers, which the matched design cannot give

The matched comparison discards strata where bowtie2 found nothing, so it
cannot be read for population counts. An **unmatched** random sample of 4,000
of the 196,650 missed sites:

| | fraction | extrapolated count |
|---|---|---|
| dG ≤ −20 | 35.2% | ~69,200 |
| dG ≤ −25 | 8.2% | ~16,200 |
| dG ≤ −30 | 0.7% | ~1,300 |
| seed-intact **and** dG ≤ −25 | 8.4% of 104,580 | **~8,800** |
| seed-intact **and** dG ≤ −30 | 0.6% of 104,580 | ~630 |

Strongest missed site −36.9 kcal/mol; within the matched comparison the
strongest missed and shared sites are −35.0 and −36.3. The missed
distribution is not truncated at the top.

**Version check, and it passed.** The matched comparison was first run with a
local ViennaRNA 2.7.2 and then in the image with **2.8.0**: every figure in
the matched table is identical to the decimal place, so the two agree on this
computation. The unmatched figures moved slightly (−25: 8.5% → 8.2%) because
the in-image run draws a different 4,000-site sample, not because the energy
model changed. **The in-image run is canonical** — it is the one a reader can
reproduce — and the manuscript carries its numbers.

### What the paper can now say, and what it cannot

**Can:** on one hundred guides, roughly **9,100 sites combine a perfectly
intact RISC seed with a predicted duplex at or below −25 kcal/mol, and no
seeded aligner reports any of them.** That is the answer to "are you
enumerating noise".

**Cannot:** that the missed sites are as strong as the shared ones. They are
modestly weaker, and the paper must say so in the same breath — the defence is
that the deficit is mostly gap composition and that seed-mediated repression
does not require a strong full-length duplex, not that the deficit is absent.

**Caveats to carry into the text.** `duplexfold` scores the most stable
intermolecular structure of two sequences in isolation; it models neither RISC
geometry nor target accessibility, so dG is a hypothesis about binding, not a
prediction of silencing. G:U wobbles pair in the energy model but count as
mismatches in oligoscan's gate, so the gate can exclude a site the energy
model would like. 100 guides, one library, one length.

**This supersedes the planned single-gene case study**, which would have been
n = 1 against n = 196,650. `scripts/case_study.py` remains for a named
illustration if one is wanted.

---

## 8l. bowtie2 cannot use a whole machine at a wide gate — measured 2026-11-17

The 1,000-guide run of 2026-11-16 failed, and the way it failed matters more
than the failure.

**What happened.** bowtie2 exited 137 (SIGKILL — out of memory) after 69.7 s,
wrote a 0-byte SAM, and `bench_three.py` parsed that empty file and reported
**"0 alignments, recall 0.0%"**. A tool that was killed was scored as a tool
that found nothing, and the error flattered oligoscan. Both gates did it.

**Cause, measured rather than guessed.** bowtie2's `-a` alignment buffers are
per thread and large at this gate. Running bowtie2 alone under
`/usr/bin/time -v`, 4 threads, `-S /dev/null`:

| guides | wall | peak RSS |
|---|---|---|
| 100 | 1:06:18 | 35.3 GB |
| 300 | 2:41:09 | 36.8 GB |
| 1,000 | 8:41:49 | 38.3 GB |

**~9 GB per thread, nearly independent of read count.** So 22 threads wants
roughly 190 GB; the machine offers 110 GB (`.wslconfig memory=112GB`,
container sees 110). It died during early allocation, which is why 69.7 s.

The 100-guide run at `-p 22` survived only because each thread saw ~5 reads
and never met a pathological one; at 1,000 reads every thread does.

**This is a publishable practical limitation, not just an ops note.** At a
realistic gate on 19-mers, bowtie2 cannot use all 22 cores of a 128 GB
workstation — it is capped near 11 threads by memory alone. oligoscan at the
same gate uses every core and peaks under 2 GB.

**Fixes applied.**

- `bench_three.py` now builds the site index **after** bowtie2 exits, not
  before. It used to hold 16.2M sites (several GB) across bowtie2's whole
  run. This was not the cause — 110 GB swallows it — but holding it was
  pointless and it narrowed the margin.
- A non-zero exit or empty SAM is now reported as `**FAILED** (exit N) --
  NOT a zero-recall result`, with the stderr tail, and exit 137 prints the
  memory advice. **Never score a dead tool.**
- New `--bt2-threads`, so bowtie2 can be throttled while oligoscan keeps
  every core and its timings stay comparable with full-thread runs.

**Cost of the diagnosis, and a lesson.** The scaling test wrote to
`/dev/null` to measure survival, so the 1,000-guide run that **succeeded**
in 8:41:49 produced no SAM. A survival test that writes a real output costs
nothing extra and would have been the run itself. Write the output.

### A second thing it turned up: the two 1,000-guide sets are different

oligoscan reported 16,191,128 sites here against 15,666,765 in the GPU grid
(§8j) at nominally the same gate. Benign, and worth recording so nobody
chases it twice: **`bench_three.py` takes a seeded random subsample**
(`--subsample-seed 20261112`) while `bench_gpu_scaling.py` takes the first
*n* lines. The two 1,000-guide sets **overlap by only 423 guides**. The gate
strings are identical in both logs. Different guides, not different answers;
nothing in the manuscript is affected.

A random subsample is the better default — the head of a library can carry
ordering artefacts — so `bench_gpu_scaling.py` should adopt it if the GPU
grid is ever re-run. Until then the Methods must say that the GPU table and
the bowtie2 table do not rest on the same 1,000 guides.

---

## 8m. Exhaustiveness, against an independent oracle — 2026-11-17 (`oracle_20261117`)

The paper's central claim had three legs and one was weaker than it looked:

- `bowtie -v 3` agreement — **mismatches only, capped at 3**;
- planted truth — proves oligoscan finds what was *put there*, not everything;
- `--selftest` — Myers against a DP **in the same source file, same author**.
  A shared misconception passes both.

Nothing independent covered the **indel path at the full gate**. Now something
does.

### What is checked

`scripts/exhaustiveness_oracle.py`, pure Python/numpy, a different language
and a different algorithm (full O(*mn*) DP, not a bit-vector), sharing no code
with oligoscan. For every (oligo, transcript, end position) it decides by
brute force whether a full-length alignment ending there exists within the
gate, and at what cost, then requires:

> oligoscan reports a description ending at a position **if and only if** one
> exists within the gate, **and** the cheapest description it reports there
> costs exactly the true minimum.

Completeness, soundness, and optimality of the description. The third clause
matters: it is the one that catches a traceback finding the right END but
describing it suboptimally, which is exactly how 71 real sites were lost in
0.9.2 (§8i) — the description's cost exceeded the gate and the site was
dropped though a cheaper description existed.

**Why a plain DP can serve as the oracle.** The gate bounds total edits, gap
bases and gap events. Each gap base costs one edit and each event spans at
least one base, so when `maxIndelBases >= maxTotal` and `maxIndelEvents >=
maxTotal` the only binding constraint is the edit total. At `maxIndelBases =
0` the oracle switches to **Hamming** distance instead, since the
mismatch-only fast path is a different code path and an edit-distance oracle
would wrongly credit it with gapped alignments.

The oracle replicates only two documented SPEC decisions — a non-ACGT
reference base mismatches everything, and alignments never span transcripts —
because those are claims about behaviour, not implementation details.

### Result: 3.0 Mb, 1,541 sequences, 100 guides, seven gates

| gate | end positions in gate | (oligo, transcript) pairs | missed | invented | cost disagreements |
|---|---|---|---|---|---|
| total ≤2, no indel | 18 | 16 | 0 | 0 | 0 |
| total ≤3, no indel | 76 | 74 | 0 | 0 | 0 |
| total ≤4, no indel | 880 | 863 | 0 | 0 | 0 |
| total ≤2 + 2 nt | 61 | 26 | 0 | 0 | 0 |
| total ≤3 + 3 nt | 585 | 421 | 0 | 0 | 0 |
| total ≤4 + 4 nt | 9,866 | 6,308 | 0 | 0 | 0 |
| **total ≤6 + 6 nt** | **982,156** | **115,033** | **0** | **0** | **0** |

**993,642 end positions checked across the seven gates, zero disagreements in
all three directions**, from roughly 2.1 billion (position, oligo) decisions.

### The negative control, which is what makes this worth anything

A test that cannot fail proves nothing. Run against the **0.9.2** binary, the
same oracle **fails 4 of 7 gates**, reporting missed end positions on
`EDGE_start_leading_bulge` — the transcript planted precisely for that failure
mode. It passes 0.9.3. The oracle detects a real historical bug and clears the
fix.

### Wired in

`tests/test_exhaustiveness.py` is self-contained (builds its own 138-sequence
reference, no downloaded data) and runs `--quick` in the Dockerfile, so the
image will not build if exhaustiveness regresses. The full seven-gate sweep
over 3 Mb is the version to quote.

### What this does NOT establish

A finite test, not a proof. The honest sentence is *verified exhaustive over
this space*, not *proven exhaustive*. It does not cover reference sizes where
a 32-bit index might overflow.

### Length sweep — the word boundaries are clean (added same day)

The gap flagged above is now closed. Myers packs the pattern into 64-bit
words, so 32/33 and 63/64 are where an off-by-one in the word handling would
live. Two passes:

**Edge cases at every length**, against a 95 kb reference carrying, per
length, an exact site, a site at position 0 with the first oligo base
unpaired, a site flush at the 3′ end, a transcript that IS the site, a
transcript shorter than the oligo, and a two-bulge site:

| length | 10 | 12 | 19 | 20 | 21 | 31 | 32 | 33 | 48 | 63 | 64 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| verdict | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS |

**Realistic hit counts**, with guides drawn FROM the 3 Mb reference (a random
long oligo never occurs in real sequence, so random guides would have tested
nothing but their own planted sites), 12 guides per length, gate total ≤4 +
4 nt:

| length | end positions | missed | invented | cost disagreements |
|---|---|---|---|---|
| 19 | 1,634 | 0 | 0 | 0 |
| 32 | 175 | 0 | 0 | 0 |
| **33** | 495 | 0 | 0 | 0 |
| 63 | 182 | 0 | 0 | 0 |
| **64** | 226 | 0 | 0 | 0 |

**Mixed lengths in one scan** (20 guides across all five lengths): 1,305 end
positions, 0 disagreements, and the output is **byte-identical to the union of
the five per-length runs** — which verifies the documented claim that lengths
may be mixed freely, rather than leaving it as an assertion.

**65 nt is refused with a reason**, not truncated: `skipping L65: longer than
64 nt`, then `no usable oligos; nothing to scan`.

So the claim now covers lengths 10 to 64 including both word boundaries, seven
gates, and mixed-length runs. What remains untested is reference scale, not
oligo length.

---

## 8n. Cas-OFFinder — measured 2026-11-17 (`tools_matrix`)

The last comparator that had a paragraph instead of numbers. 100 Huesken
guides, whole human cDNA, mismatch gates, one container, one reference,
one I/O path. oligoscan reads the packed `.osr`, Cas-OFFinder the FASTA:
each tool in its deployed form.

| gate | sites | oligoscan CPU | oligoscan GPU | Cas-OFFinder 2.4.1 | 3.0.0b3 (beta) |
|---|---|---|---|---|---|
| total ≤2 | **3,036** | 5.5 s | 1.8 s | 156.1 s | 169.4 s |
| total ≤3 | **21,995** | 6.1 s | 1.7 s | 210.0 s | 225.9 s |
| total ≤4 | **204,213** | 7.0 s | 2.0 s | 313.3 s | 333.2 s |
| peak RAM | — | 928 MB | 327 MB | 5,945 MB | 5,829 MB |
| output, gate 4 | — | 12.4 MB | 12.4 MB | 116.2 MB | 117.7 MB |

**All four tools report the SAME sites at every gate — identical sets at
exact transcript coordinates, nothing unique to anyone.**

### 1. The best exactness evidence in the paper

bowtie1 is an oracle only to three mismatches. **Cas-OFFinder is uncapped,
like oligoscan**, so this tests the same regime by an entirely different
method — and agrees on 204,213 sites at gate 4. Independent witnesses to
exactness are now: an FM-index backtracker (bowtie1, ≤3), an OpenCL
exhaustive searcher (Cas-OFFinder, uncapped), and a brute-force DP oracle
(§8m). No two share code.

Worth noting separately: **2.4.1 and the 3.0.0b3 beta agree with each other
exactly**, which the authors' own warning that v3 "may differ from
Cas-OFFinder 2" did not guarantee.

### 2. Flatness is NOT a property of exhaustiveness — it is a property of this design

The finding that is easy to miss. Across a **67-fold** rise in sites:

| | time multiplier |
|---|---|
| oligoscan CPU | **1.27x** |
| oligoscan GPU | **1.11x** |
| Cas-OFFinder | **2.01x** |

Cas-OFFinder is exhaustive and its time still **doubles**. So the paper must
stop implying that flatness follows from scanning every position; it does
not. It follows from *this* scan's structure — a fixed cost per reference
position that the error budget does not change — and Cas-OFFinder is the
proof that another exhaustive tool need not have it. That makes the flatness
claim both narrower and much better supported.

(RIsearch2 is flat in its ENERGY threshold, §8h, which is a different axis
again. Three tools, three different flatness behaviours.)

### 3. Resources

- **RAM: 6.4x less than Cas-OFFinder on CPU, 18.2x less on GPU.** The GPU arm
  uses *less* host memory than the CPU arm (327 MB against 928 MB) because
  the reference lives on the device.
- **Output 9.4x smaller** at gate 4 (12.4 MB against 116.2 MB). Partly a
  format difference — Cas-OFFinder writes the matched genomic sequence for
  every hit — so report it as such rather than as pure overhead.
- **Speed 28-45x on CPU, 87-157x on GPU**, and the ratio *grows* with the
  gate because oligoscan is flat and Cas-OFFinder is not.

### 4. State the one-time cost, or a referee will

oligoscan reads a packed `.osr`: **15.4 s one-time, 839 MB FASTA -> 207 MB**.
Cas-OFFinder needs no preparation at all. This must sit beside any timing
quoted from the table. For scale, RIsearch2's suffix array costs 291.6 s and
11.9 GB (§8h), so oligoscan's preparation is an order of magnitude cheaper
than the only other indexed comparator — but it is not zero, and the
comparison is only honest if it is named. `--also-fasta` adds an
oligoscan-on-FASTA row (≈24 s) for anyone who wants every tool parsing the
same file.

### Two harness bugs this run exposed, both silent

- **Cas-OFFinder keeps the WHOLE fasta header as the sequence name**, where
  oligoscan keeps the first token. Ensembl headers carry a description, so
  every key differed and the first run reported **3,036 sites each, 0
  shared** — found by Jeff, who spotted that identical counts with zero
  overlap is arithmetically impossible.
- **2.4.1 and 3.0.0b3 write DIFFERENT column layouts.** v3 is
  `Id, Bulge Type, crRNA, DNA, Chromosome, Location, Direction, Mismatches,
  Bulge Size`; v2 is `crRNA, chrom, location, DNA, strand, mismatches[, id]`.
  Parsing v3 with the v2 layout put "Chromosome" where the strand belongs, so
  v3 reported **0 sites beside a 1.3 MB output file**. The reader is now
  header-aware.

Both failures produced a *plausible-looking* number rather than an error.
That is the third time in this project that a silent parse has flattered
oligoscan; the matching layer, not the scanner, is where the risk lives.

### Bulges — the comparison is bounded, not completed (2026-11-17)

2.4.1 has no bulge support. 3.0.0b3 has it, and three things stand in the
way of a straight comparison. All three were found by the conventions
self-check, which is the reason it exists.

**1. It models ONE bulge of ONE kind per alignment.** The output carries a
single `Bulge Type` (X/DNA/RNA) and a single `Bulge Size`. oligoscan allows
several gap events and both directions at once, so oligoscan's multi-event
sites are outside Cas-OFFinder's model *by construction* and must never be
counted as sites it "missed". A tool cannot miss what it cannot represent.

**2. The bulge-type mapping, confirmed from planted sites:**

| Cas-OFFinder | oligoscan |
|---|---|
| `DNA` bulge | `target_bulge` |
| `RNA` bulge | `oligo_bulge` |
| `X` | ungapped |

**3. Coordinates shift in bulge mode, and differently per bulge type.** Three
sites all planted at position **120** were reported at **121** (DNA), **122**
(exact) and **123** (RNA), where oligoscan reports 120 for all three. In
mismatch-only mode the coordinates agreed exactly (§8n). So exact-coordinate
matching is meaningless here and the comparison must be at locus level with a
tolerance, as it is against RIsearch2.

**4. It emits many descriptions of one site.** A single PERFECT match
produced **eight** rows — `DNA 2` at 118, `DNA 1` at 120 and 121, `X 0` at
122, `RNA 1` at 123 and 124, `RNA 2` at 124 and 126 — because with bulges
allowed a perfect match can also be read as a shifted alignment carrying one.
Raw rows are not sites; they must be collapsed before counting.

### The runtime bound, which is itself the finding

**At 100 guides with a bulge budget of 2, Cas-OFFinder 3.0.0b3 did not
complete in over two hours.** oligoscan finished the same gate in ~18 s.
Budgets of 4 were abandoned earlier the same day after hours without
finishing — note that `DNA 4 + RNA 4 + 4 mismatches` permits **twelve** edits,
since Cas-OFFinder treats bulge size and mismatches as separate allowances,
so it is a far wider search than oligoscan's gate of four.

This mirrors bowtie2: gapped search is where every other tool collapses.
The cost of enabling gaps, measured:

| tool | mismatch-only, gate 4 | with gaps | factor |
|---|---|---|---|
| oligoscan CPU | 7.0 s | 17.8 s | **2.5x** |
| Cas-OFFinder 3.0.0b3 | 313 s | > 2 h (100 guides, bulge 2) | **> 23x** |
| bowtie2 | — | 2,412 s | dominated by `-a` |

oligoscan pays a constant factor — Myers' recurrence is more arithmetic per
position than a popcount, on the same fixed amount of work. A search that
enumerates bulge configurations pays a combinatorial price instead. **That is
the structural argument of the paper, measured against the one tool that
shares its exhaustiveness.**

`--timeout` caps a run and records "DID NOT COMPLETE", because an
unfinished search has no defined site set and must not be scored as one.
### Bulges — completed 2026-11-17, 50 guides

Three failed designs preceded this one; the arithmetic that makes it work is
Jeff's. **At edit cost <= 4 every alignment has at most 4 mismatches by
arithmetic**, so with Cas-OFFinder at `mm <= 4, bulge <= 1` and oligoscan at
`total <= 5, 1 nt indel, 1 event`, each tool can express everything the other
can and the (mismatch, bulge) split stops mattering. A 1 nt bulge also cannot
be split across two gap events, which removes the multi-event difference that
produced a 436-locus residual at bulge 2.

| | oligoscan CPU | oligoscan GPU | Cas-OFFinder 3.0.0b3 |
|---|---|---|---|
| wall | **31.5 s** | 37.5 s | **5,180.3 s** |
| peak RAM | 918 MB | 964 MB | — |
| output | 162.7 MB | 162.7 MB | 3,557.9 MB |

**164x faster on CPU, with 21.9x less output.** Site agreement at cost <= 4:

| | loci |
|---|---|
| shared | **328,431** |
| oligoscan only | 120 |
| Cas-OFFinder only | 2,952 |

**99.96% of oligoscan's set and 99.11% of Cas-OFFinder's.** The residual is
**concentrated at the cost boundary** — 2,930 of 2,952 and 117 of 120 sit at
cost 4 exactly — which is the signature of an equivalence edge rather than a
broken convention. A residual spread evenly across costs would have meant
something was still wrong.

**The GPU is SLOWER than the CPU here** (37.5 s against 31.5 s). Not an
anomaly: at 50 guides the per-launch fixed cost exceeds the kernel saving,
the same crossover §8j measured. Report it; it tells users when the GPU is
worth using.

### The two-gap-event example — NOT part of the comparison above

**Where this came from, because it is easy to misread.** This case was met
during the **bulge-2** exploration, not in the 50-guide benchmark, and it
takes **no part** in the 328,431 / 120 / 2,952 figures. In the bulge-1 run
Cas-OFFinder cannot emit a 2 nt bulge at all and oligoscan runs at one
event, so neither tool can express the site and both exclude it,
symmetrically. The Cas-OFFinder output IS filtered before comparison —
rows to `bulge <= 1` and `cost <= 5`, then loci to `cost <= 4` — so a
6-edit description never enters.

It is recorded because it shows what the single-bulge MODEL can represent,
which is a separate finding from recall.
`Hu_2070` against `ENST00000679802.1` near position 1088:

```
rc(oligo) : AGCTGCTG   TGTCC   CGGTGG
target    : AGCTGCTG C TGTCC A CGGTGG
```

**Two edits — two separate single-base target bulges.** The best ungapped
window there carries **eight** mismatches.

- A single-bulge model **cannot represent this alignment at all**.
  Cas-OFFinder describes the site as 4 mismatches + a 2 nt bulge: **6 edits**.
- oligoscan constrained to `events 1`, to match that model, also cannot, and
  correctly reports its cheapest single-event reading (5-6).
- oligoscan at `events 2` finds it as **2 edits**.

A site two edits from the guide, reported as six by the only other
exhaustive tool. This is the multi-event blind spot measured against bowtie2
(§8d, §8e), now shown against an exhaustive tool with a concrete sequence —
and it is a MODEL limitation, not a search failure, which is the honest and
more interesting claim.

**A caution this episode earned.** The diagnostic that first flagged these
loci used unconstrained edit distance, while both tools ran at one gap event.
That is not a like-for-like oracle, and it briefly looked as though oligoscan
had missed a 2-edit site. When comparing against a constrained run, constrain
the oracle the same way.


---

## 9. What remains

1. ~~Per-locus selection layer~~ — **done 2026-11-11**, see §5a.
2. ~~Validate on the full human cDNA against bowtie1 at ≤3 MM~~ — **PASS
   2026-11-12**, see §4; 4–6 MM and indel timing on the full reference done
   (§6). Still open: comparison against the deployed bowtie1 bulge path at
   1 nt (`bt1_search.py`; snapshot with md5s in `data/platform_snapshot/`;
   compare antisense, seed≤1, body≤3, ≤2 substitutions + 1-nt bulge, one
   event, no seed/terminal bulges; oligoscan gate S=1 N=3 T=3 maxI=1 maxE=1).
3. ~~Packed reference~~ — **done 2026-11-12**: `oligoscan pack ref.fa ref.osr`
   (2 bits/base + sequence table + sentinel runs for padding and non-ACGT;
   version byte checked on load). Loading expands to 1 byte/base, so the scan
   is unchanged and output byte-identical to FASTA (checked on synthetic
   edge cases, the selection suite §H, and the full cDNA in mismatch and indel
   mode). `build_indexes.sh` now packs too. Mismatch output is now **sorted**
   like indel output, so every mode is identical across runs and threads.
   Still open: scanning the 2-bit data directly (XOR/popcount) to quarter
   scan memory — shares its kernel design with the GPU path.
4. **GPU path** (decided 2026-11-11; next): target is Jeff's
   NVIDIA RTX 4000 Ada (20 GB, 48 SMs, sm_89). Mismatch kernel = 2-bit window
   XOR + popcount with seed / pos-1 / body masks; indel kernel = Myers per
   thread, candidates back to CPU for traceback + selection. Estimated
   ~2–5 s mismatch-only and ~5–15 s with indels for 560 Mb × 1,750 oligos
   (vs ~66 s / ~275 s projected on 20 CPU cores) — estimates, unmeasured.
   Second compose service on `nvidia/cuda:12.x-devel-ubuntu24.04`; CPU build
   stays the default and the oracle (GPU output must be byte-identical).
5. ~~Decide the name and license~~ — **done 2026-11-12**, see §10.
6. **Manuscript** — drafted 2026-11-12 as a Claude doc, *oligoscan —
   Bioinformatics Application Note (draft)*: abstract, introduction with the
   information-theoretic argument, implementation, four result tables,
   discussion, limitations, references, and an open-items checklist. All
   benchmark numbers in it are measured. See §9a for what remains.
7. Only then wire into Tabs 2/3 behind a toggle, so results can be A/B'd
   against the current implementation before anything is replaced.

---

## 9a. What the manuscript still needs (as of 2026-11-12)

`scripts/bench_paper.py` produces most of it in one run and writes
`results.md` ready to paste. It needs the real bowtie1 and bowtie2 indexes
and takes a few minutes (bowtie2 `-a` is the slow part):

    docker compose run --rm oligoscan python3 /opt/oligoscan/scripts/bench_paper.py \
      --ref /data/ref/Homo_sapiens.GRCh38.cdna.all.fa \
      --scan-ref /data/ref/Homo_sapiens.GRCh38.cdna.all.osr \
      --oligos /opt/oligoscan/tests/example_oligos.tsv \
      --bt1 /data/index/bowtie1/Homo_sapiens.GRCh38.cdna.all \
      --bt2 /data/index/bowtie2/Homo_sapiens.GRCh38.cdna.all

**Still missing, in order of how much a reviewer will care:**

1. ~~bowtie2 recall at 4–6 edits~~ — **measured 2026-11-12, see §8a.**
   93.8% at ≤4, 58.0% at ≤6, non-monotonic, 963×/684× slower.
   Measurement B also measured (§8a): oligoscan **100.0%** (702/702),
   bowtie2 **70.2%** tuned / **68.8%** at defaults, and **11% once two
   substitutions accompany the bulge**. All benchmark numbers are final.
2. **Published oligo set — obtained 2026-11-12.**
   `tests/huesken_guides.tsv`: **2,361 unique 19-nt antisense guides** from
   Huesken et al. (2005), the set RIsearch2 used, so the numbers land on
   common ground with the existing literature.
   *Provenance matters for the citation*: the canonical set is 2,431 siRNAs
   at 21 nt; this copy comes from the OligoFormer repository
   (`lulab/OligoFormer`, `data/Hu.csv`), which truncated the 3′ overhang at
   positions 20–21 and filtered 70 sequences. **Cite it as 2,361 Huesken
   sequences via that redistribution, never as "the Huesken 2,431 set."**
   Verified on import: all 2,361 unique, all 19 nt, all ACGT, and all
   **antisense** — the reverse complement of each matched its own mRNA
   context in 300/300 spot checks, which is the orientation oligoscan's
   `anti` mode expects. Library-like composition (GC 11–95%, median 53%;
   balanced base and 5′-base distribution), unlike the 40 TRIM63 oligos.
   Note 19-mers carry 38 bits against a 21-mer's 42, so site counts run
   higher than the 21-mer tables; each table must name its library.
3. ~~CPU model~~ — **confirmed 2026-11-12**: Core Ultra 9 285, 24C/24T, no
   AVX-512, 128 GB RAM; see §6.
4. **Positioning against Cas-OFFinder** — the real related-work risk, since
   it is *also* exhaustive over mismatches and bulges. The honest
   differences: PAM-free, transcriptome-aware with structural
   junction safety, a region-structured seed/body gate rather than one
   total, one row per locus by a total order, no OpenCL dependency, and no
   cap at all. A reviewer will ask; a sentence is the minimum, a benchmark
   is stronger.
5. **Figures**: scan time flat across the gate (§6), and bulge recall
   against error load (§8a) — the cliff at two substitutions is the
   paper's central figure. Both from TSVs the harness already writes.
6. **A 2+-bulge benchmark stratum** — the gap found 2026-11-12. The
   engine is correct here: `maxIndelBases` gates total gap bases across all
   events and `maxIndelEvents` gates the event count, so one 2-nt bulge and
   two separate 1-nt bulges are distinct configurations, and selection
   check C plants 0-2 indel bases at independent random positions and both
   directions, validated 256/256 against a gap-constrained DP oracle. But
   **every published recall number is 1-nt bulge only** — measurement B's
   702 sites, both stratum tables and the planned figure. The stated
   motivation (§8: `var_target_bulge.fa` at MAX INDEL LENGTH = 2 produced
   14,219 MB of hits and zero usable rows in 30 min) is precisely the
   2-bulge case, and it has no table. Add strata {one 1-nt bulge, one 2-nt
   bulge, two separate 1-nt bulges} × {0,1,2 further substitutions}, same
   flank-aware verification, gate maxI=2 maxE=2. bowtie1 is absent by
   construction; bowtie2 collapsed to 11% at one bulge + two substitutions
   and should be near zero here.
7. ~~Repository URL~~ — **settled 2026-11-16: `https://github.com/breezedu/OligoScan`**. Still open: a Zenodo DOI
   for the tagged release, and the commercial-licensing contact that
   `LICENSE` and the README still carry as a placeholder. Note the
   repository capitalises the S (`OligoScan`) while the tool, its binary
   and every usage string are lowercase `oligoscan`; GitHub redirects
   either way, but pick one spelling for the paper and stay with it.
8. Author affiliation; funding and competing-interest statements.
9. A primary citation for seed-region complementarity driving off-targeting
   (Birmingham 2006 / Jackson 2006).

**Not needed for this venue**: wet-lab or expression validation. It would be
needed for NAR Methods Online, and the natural test would be whether sites
oligoscan finds and bowtie2/RIsearch2 miss show measurable repression in
siRNA-transfection RNA-seq.

---

## 10. Naming, license, release — named 2026-11-12, licence changed 2026-11-16

**Name: `oligoscan`.** "bowtie3" was rejected: it reads as the official
successor to Langmead's Bowtie/Bowtie2 (citation confusion, and realistically
a rename request), and `bowtie_oligo` still implies derivation. A web search
found no software collision for `oligoscan`; it does sit in a crowded `Oligo*`
namespace (OligoWalk, OligoScreen, OligoFaktory, OligoShell), which is a
discoverability concern, not a conflict.

**License: changed from MIT to source-available, 2026-11-16 (Jeff).**
Free for academic and non-profit research; commercial use, and in particular
use by large pharmaceutical companies, requires a paid licence. All code is
fresh — nothing derives from bowtie — so no license was inherited, and
**copyright is held by Jeff Du personally** (confirmed 2026-11-12), which is
what makes dual licensing possible at all. `LICENSE` and the `oligoscan.c`
header read `Copyright (c) 2026 Jeff Du`; `CITATION.cff` lists him as sole
author.

**Consequences, and they are not only cosmetic:**

- **It is no longer open source and must never be described as such.** A
  non-commercial restriction is incompatible with the OSI definition by
  construction. The correct phrase is *source-available, free for academic
  and non-profit research*. Audit README, BUILD, CITATION.cff and the
  manuscript for the words "open source" and for "MIT".
- **The journal is fine with this.** OUP's *Bioinformatics* requires only
  that software be "freely available to non-commercial users", with source
  at a stable URL without request at submission time, and kept available two
  years after publication; OSI licensing is encouraged, not required. NAR's
  own policy page could not be retrieved 2026-11-16 — **confirm its exact
  wording before submission.** ISCB Level 0 availability is also met and
  exceeded (they ask only for free binaries to research users).
- **CHOSEN 2026-11-16: PolyForm Noncommercial 1.0.0**, unmodified, taken
  verbatim from `github.com/polyformproject/polyform-licenses`. SPDX id
  `PolyForm-Noncommercial-1.0.0` (verified against the SPDX licence list, so
  `CITATION.cff` validates). BUSL-1.1 was the considered alternative and was
  not taken.
- **Pharma-SPONSORED academic work is FREE, and this is deliberate.** The
  licence's *Noncommercial Organizations* clause permits use by an
  educational institution or public research organisation "regardless of the
  source of funding or obligations resulting from the funding" — the licence
  says so in terms, it is not a gap. A university group running oligoscan
  under industry sponsorship needs no licence. Decided 2026-11-16 (Jeff) in
  preference to a custom rider, which would have slowed institutional legal
  review for little practical gain.
- **The applied files are done:** `LICENSE` (full text with four
  `Required Notice:` lines, the mechanism the licence's *Notices* clause
  specifies), `CITATION.cff` (`license:` field), the `src/oligoscan.c`
  header, and the README's License section. Build re-verified after the
  header edit (`oligoscan 0.9.2`). No "MIT" string remains anywhere in the
  tree outside this handoff's history.
- **STILL OPEN — the commercial contact.** `LICENSE` and the README both
  carry a placeholder. A commercial licence nobody can buy is a refusal with
  extra steps, and the Required Notice line travels with every copy. Set it
  before the repository is public. Jeff's personal email was deliberately NOT
  written into a public file without his say-so.
- **Distribution narrows.** Debian/Ubuntu packaging is ruled out (non-free).
  conda-forge is unlikely; bioconda does carry non-commercial packages
  (ViennaRNA is one), so it may still be possible — check rather than assume.
  GitHub and Zenodo are unaffected.
- **Dual licensing needs a contact route.** A commercial licence nobody can
  buy is a refusal with extra steps. The README needs an address and a named
  holder.
- **Dependency alignment, which happens to help.** ViennaRNA is itself an
  academic-use licence, and `scripts/case_study.py` depends on it, so that
  pipeline was never freely usable commercially under MIT either. The
  SCANNER has no ViennaRNA dependency — keep it that way so oligoscan's own
  terms stay clean and only the case-study script carries the extra
  condition. bowtie1/bowtie2/RIsearch2 sit in the benchmark image as separate
  binaries, not linked, so they impose nothing on oligoscan.

**Public tree** (19 files, 228 KB) builds and passes every test standalone:
`src/oligoscan.c`, `tests/` (selection suite + bowtie1 harness),
`scripts/` (build_indexes, bench_gates, bench_paper), `examples/` (54 kb
self-contained demo with planted sites of known kind), `README.md`,
`BUILD.md`, `LICENSE`, `CITATION.cff`, Dockerfile + compose.

**Deliberately NOT published** (in `.gitignore`): `PROJECT_HANDOFF.md` — this
document references the deployed platform's modules, production measurements
and internal paths; `data/` — reference, indexes, `platform_snapshot/`;
`src/archive/` — the frozen pre-selection build and the superseded Stage 2
suite, both of which still carry sandbox paths. Verified: the public tree
contains no internal path, module name or employer reference.

---

## 11. Prior art — position it, don't get positioned

| tool | overlap |
|---|---|
| **Cas-OFFinder** (2014) | exhaustive mismatch **and bulge** search for short guides, GPU/OpenCL. Strong *algorithmic* precedent; CRISPR/PAM-oriented. |
| **RIsearch2** (2017, NAR) | suffix-array, large-scale siRNA off-target prediction. Direct *application* precedent. |
| **SeedMatchR**, **GESS** | seed-region off-target analysis, mostly from screen/RNA-seq data. |

None is a transcriptome-aware, seed-structured, **uncapped** exhaustive scanner
covering siRNA *and* ASO. That is a real gap — and a narrower, far more
defensible claim than "a better bowtie".

---

## Sources

- [bowtie2 NEWS (2.5.5)](https://codegraph.jelmer.uk/bowtie2/2.5.5-1/NEWS)
- [BenLangmead/bowtie activity](https://summary.ecosyste.ms/projects/356189)
- [RIsearch2](https://pubmed.ncbi.nlm.nih.gov/28108657/)
- [Cas-OFFinder](https://pmc.ncbi.nlm.nih.gov/articles/4016707)
- [Intel Core Ultra 9 285 specifications](https://www.intel.com/content/www/us/en/products/sku/241061/specifications.html) — 24 cores (8P+16E), 24 threads, 36M cache, 5.6 GHz turbo
- Myers, G. "A fast bit-vector algorithm for approximate string matching based
  on dynamic programming." *JACM* 46(3), 1999.

### Prior project docs this consolidates

`claude/OLIGOSCAN_DESIGN_AND_STAGE1.md` and
`claude/OLIGOSCAN_STAGE2_INDELS.md` in the
*siRNA_Off-Target_analysis_withBowtie* project — worth carrying over to the new
project alongside this summary, as they hold the full detail.

## 8o. The Cas-OFFinder-only residual was an oligoscan bug, and it is fixed

**Status: root-caused, fixed, regression-tested. 2026-11-17.**

### What the residual turned out to be

The 2,951 Cas-OFFinder-only loci from the 50-guide bulge-1 benchmark were
**real alignments that oligoscan should have reported**. Every one was
independently re-derived by exhaustive enumeration over all placements with
at most one gap event of at most one gap base (`scripts/verify_cof_only.py`,
which uses neither tool's code):

| verdict | count |
|---|---|
| verified real, cost 4 | 2,929 |
| verified real, cost 3 | 22 |
| outside the cost<=4 gate | 0 |

705 were target bulges, 2,246 oligo bulges. Cas-OFFinder was not
over-crediting; oligoscan was missing them.

### Root cause

`push_aln()` applied the gap budget **after** the traceback had already
chosen a description:

```c
if (gap_bases > c->maxI || a->gap_events > c->maxE) return;
```

The traceback minimised `(edit total, gap bases, gap events)`
lexicographically -- **cost first**. Where the cheapest description at an end
position needed more gap bases or events than the gate allowed, the traceback
returned it, the gate rejected it, and the **whole candidate was discarded** --
including a legal cheaper-in-gaps description that existed at the same end
position.

Minimal real case, oligo `Hu_1943` on `ENST00000529654.5` @ 1721:

| gate | oligoscan 0.9.3 reports |
|---|---|
| `indel 1 events 1`, total<=5 | nothing |
| `indel 1 events 1`, total<=6 | pos 1721, **ungapped, cost 6** |
| `indel 2 events 2`, total<=8 | pos 1722, both_bulges, cost 3 |

The cheapest description there is a 2-event both-bulges at cost 3. Under
`maxIndelBases=1` it was generated, rejected, and the candidate vanished --
even though a 3-mismatch + 1-nt target-bulge alignment at cost 4 with ONE
event exists and is inside the gate. Cas-OFFinder reports exactly that
alignment (3 mismatches + 1 bulge).

This is the same class as the 0.9.2 junction bug: judging a locus by one
chosen description when several exist. The 0.9.3 fix made the generator and
the chooser agree with each other, but neither was constrained by the gate.

### The fix

The gate is now a **constraint on the search**, not a filter on its output.
`traceback()` carries gap bases and gap events in the DP state, bounded by
`maxI`/`maxE`, so the optimum returned is the best *feasible* description and
a feasible description is never lost to an infeasible cheaper one. State is
`(text col, pattern row, gap bases, gap events, last move)`.

Two latent defects were fixed in the same pass:

- the free-text-prefix initialisation seeded all three `last` states, letting
  the first gap move take the *continuation* branch and consume a gap base
  without charging a gap event -- one free event per alignment;
- the DP scratch is now a thread-local buffer. Allocating it per Myers
  candidate took the 50-guide scan from 78 s to over 180 s.

### Result, 50 guides, gate 4 mismatches + 1 nt bulge, compared at cost <= 4

| | oligoscan loci | shared | oligoscan-only | Cas-OFFinder-only | recall |
|---|---|---|---|---|---|
| 0.9.3 | 328,551 | 328,432 | 119 | 2,951 | 99.11% |
| fixed | 331,502 | 331,382 | 120 | **1** | **100.00%** |

Cost: 78.2 s -> 125.5 s on the same machine, same 50 guides, same peak RSS
(931 MB). **1.61x slower**, and on a 5-guide subset 8.31 s -> 10.09 s
(1.21x). The scan recovers 105,770 additional descriptions over 50 guides.

### Regression tests

- `tests/test_bulge_leading_mismatch.py` -- plants a target bulge plus three
  mismatches at six different placements. Mismatches at the leading edge of
  the DP window (window idx 0,1,2) lose the gapped description on 0.9.3 and
  are found after the fix. Exit 1 on the bug.
- `scripts/constrained_oracle.py` -- DP bounding mismatches, gap bases and
  gap events independently. The referee for any question of the form "should
  oligoscan have found this?". A plain edit-distance oracle allows unlimited
  gap events and is NOT a fair referee for a run capped at one.
- `scripts/verify_cof_only.py` -- base-by-base re-derivation and rendering of
  every Cas-OFFinder-only locus.
- `--selftest`, `test_junction_edge.py`, `test_exhaustiveness.py --quick` all
  still pass.

### RESOLVED in 8q -- this was a comparison artifact, not a miss

~~OPEN: one locus still missing~~

`Hu_0864` on `ENST00000879574.1` @ 169, an oligo bulge at cost 4 (three
distinct cost-4 alignments exist in positions 162-187). oligoscan reports
**nothing** in that window even at total<=6, so the candidate is not being
generated at all -- a *different* defect from the one fixed above, likely in
the Myers filter or tile warm-up. One locus in 331,383 (0.0003%). **Section 8q shows this was wrong**: the
locus IS reported, at positions 162 and 169, and the comparison's greedy
locus merge hid it. The Cas-OFFinder-only residual is zero.

### Consequences for the manuscript

Every exhaustiveness and recall number measured with `maxIndelBases >= 1`
predates this fix and must be re-run. The bowtie2 comparison, the
Cas-OFFinder comparison, and the ΔG analysis of "additional sites" are all
affected. The exhaustiveness oracle did **not** catch this, because it
compared minimum edit distance at each end position without constraining gap
events, so it never asked the question the gate actually asks.

## 8p. Scope of the 0.9.4 fix: which measurements it can touch

**Measured 2026-11-18. The re-run changed nothing, and that is provably correct.**

The overnight re-run rescanned all three bowtie2 gates with 0.9.4 and
reproduced the 0.9.3 output **byte for byte** (344,195 / 6,022,546 /
91,438,413 bytes; scans regenerated 04:00-04:01, version confirmed 0.9.4).
The dG re-run likewise reproduced `dg_20261116/dg_sites.tsv` byte for byte.

### Why — the invariant

Every gap base costs one edit, so for any description

    gap_events <= gap_bases <= edit_total <= maxTotal

The 0.9.3 defect fired only when the traceback's cheapest description was
rejected by `gap_bases > maxI` or `gap_events > maxE`. When
`maxI >= maxTotal` and `maxE >= maxTotal`, both caps are slack by the
inequality above and **no description can ever be rejected** -- the defect is
inert. It bites only where the indel or event budget is TIGHTER than the
total budget.

### Audit of every gate the manuscript quotes

| measurement | gate (total:maxI:maxE) | affected |
|---|---|---|
| bowtie2 100 guides | 2:2:2, 3:3:3, 4:4:4 | no -- caps slack |
| bowtie2 1,000 guides | 4:4:4 | no -- caps slack |
| mismatch-only three-way | 2,3,4 (maxI=0) | no -- mismatch path, no traceback |
| tools matrix | gate:0:1 | no -- mismatch path |
| GPU scaling grid (8j) | 3:0:1, 4:0:1, 3:3:3, 4:4:4 | no |
| dG missed-vs-shared (8k) | 4:4:4 | no |
| **Cas-OFFinder bulge (8n)** | **5:1:1** | **YES -- tight, re-run in 8o** |

So §3.6, Figures 2 and 3, the case study, and the GPU scaling table all
**stand unchanged**. The only measurement the defect ever touched is the
Cas-OFFinder bulge comparison, which is the one that exposed it.

This also states the user-facing scope precisely: 0.9.3 and 0.9.4 differ only
for runs where `maxIndelBases < maxTotal` or `maxIndelEvents < maxTotal`.

### Still outstanding

- the Cas-OFFinder bulge timing table: CPU wall clock and the GPU arm
  (`oligoscan-cuda`) both need re-measuring on the benchmark host at 5:1:1.
  The manuscript currently carries a `[TK]` placeholder for the speed
  multiple. Output size is already updated (167.1 MB, machine-independent).
- the single unexplained locus, `Hu_0864` / `ENST00000879574.1` @ 169 (§8o).

## 8q. The residual is zero: oligoscan is a strict superset of Cas-OFFinder

**Measured 2026-11-18. Supersedes the "1 locus outstanding" in 8o.**

The single Cas-OFFinder-only locus left after the 0.9.4 fix was **not an
oligoscan miss**. oligoscan reports it, at `ENST00000879574.1` positions 162
(target bulge, cost 4) and 169 (oligo bulge, cost 4). The comparison lost it.

### Why the comparison lost it

`to_loci()` merges descriptions within `--tol` greedily, left-anchored, and
keys the merged locus at the FIRST position. oligoscan's 162 and 169 are 7 nt
apart, so they collapse to one locus keyed at **162**. Cas-OFFinder reports
loci at 162 and 171. Matching the merged representative, 171 - 162 = 9 > 8,
so Cas-OFFinder's 171 found no partner -- even though oligoscan's 169 sits
2 nt from it.

Judging a locus by one chosen representative position: the same failure mode
as the 0.9.2 junction bug and the 0.9.3 gate-after-traceback bug, now in the
comparison layer rather than the scanner. **Match on raw description
positions, then collapse only for counting.**

### The corrected result, 50 guides, gate 4 mismatches + 1 nt bulge, cost <= 4

| matching | shared | oligoscan-only | Cas-OFFinder-only | CoF recall |
|---|---|---|---|---|
| collapsed representative | 331,382 | 120 | 1 | 99.9997% |
| **raw position (correct)** | **331,383** | **120** | **0** | **100.0000%** |

**oligoscan reports every locus Cas-OFFinder reports.** The relationship is
strict superset, not near-agreement.

### And the 120 are real

All 120 oligoscan-only loci were re-derived independently by exhaustive
enumeration: 117 at cost 4, two at cost 3, one at cost 2. None is a false
positive. Their location explains them:

| where | n |
|---|---|
| flush against the transcript 3' end | 116 |
| beginning at the transcript's first base | 4 |

So the difference is a **boundary effect in Cas-OFFinder**, which does not
report alignments running up against a sequence end, not a disagreement
about transcript interiors. For transcriptome scanning this matters: 3' UTR
ends are exactly where miRNA-like seed pairing is of interest.

### Two bugs fixed in `verify_cof_only.py`

Both produced plausible but false conclusions, and both were mine:

- the three alignment kinds have different text spans (L, L+1, L-1); bounding
  every start by `len(w)-L` never tests a site flush against a 3' end. This
  made all 116 look like oligoscan false positives at cost 5-13.
- the unpaired base may be the first one; starting the gap loops at index 1
  misses a 5'-end oligo bulge. This made the remaining four look like
  cost-5 false positives.

**Standing pattern, now five instances: every error this project has hit has
been in the bookkeeping or verification layer, never in the scanner, and each
produced a believable number rather than an error.** Three flattered
oligoscan, two maligned it. Any claim of the form "tool X missed this" gets
an independent re-derivation before it is believed.

## 8r. Final Cas-OFFinder table, measured on the benchmark host 2026-11-18

oligoscan 0.9.4, 50 guides, Cas-OFFinder DNA/RNA bulge 1 + 4 mismatches
<-> oligoscan total<=5, 1 nt indel, 1 event. Cas-OFFinder reused (it did not
change); both oligoscan arms re-timed.

| | oligoscan CPU | oligoscan GPU | Cas-OFFinder 3.0.0b3 |
|---|---|---|---|
| wall clock | **34.1 s** | 65.0 s | **5,180.3 s** |
| peak memory | 923 MB | 969 MB | — |
| output | 167.1 MB | 167.1 MB | 3,557.9 MB |

**152x faster** (5,180.3 / 34.1), writes **21.3x less** (3,557.9 / 167.1).

Against 0.9.3: CPU 31.5 -> 34.1 s (1.08x), GPU 37.5 -> 65.0 s (1.73x). The
GPU arm pays more because the CUDA backend replaces ONLY the Myers scan --
traceback runs on the host in both arms, and the constrained traceback added
host work, so the GPU's share of the total shrank. The GPU remains slower
than the CPU at fifty guides, now by 1.9x rather than 1.19x.

### Comparison: 149 was wrong, 120 is right

The first run on the host printed oligoscan-only **149**, not the 120 this
handoff reports. Reproduced and diagnosed: a third instance of
representative-anchoring.

`bench_casoffinder_bulge.py` collapsed descriptions at cost <= k (5) and only
then kept loci whose minimum was <= CC (4). `to_loci()` keys each locus at
its FIRST position, so a **cost-5 description can anchor a locus** whose
qualifying member sits up to `--tol` away; matching that anchor then fails
spuriously. Filtering to CC before clustering guarantees every anchor is
itself a qualifying description.

| ordering | oligoscan loci | shared | oligoscan-only | CoF-only |
|---|---|---|---|---|
| collapse at <=5, then keep <=4 (old) | 331,503 | 331,354 | 149 | 0 |
| **filter to <=4, then collapse (correct)** | **331,502** | **331,382** | **120** | **0** |

Verified anchor-independently by matching on ANY member position of a locus
rather than its anchor: identical result, 331,502 / 331,382 / 120 / 0.
Cas-OFFinder reports 331,352 loci and oligoscan reports every one of them.

`bench_casoffinder_bulge.py` now filters before clustering, so a re-run
reproduces 120.

**This is the sixth bookkeeping-layer error and the third of the
representative-anchoring kind specifically** (0.9.2 junction: one chosen
description per locus; 2026-11-18 comparison: merged-locus anchor;
2026-11-18 again: non-qualifying anchor). Whenever a position stands in for
a set, check what happens when the stand-in is not a member of the set.

## 8s. Clean timed bowtie2 comparison — measured 2026-11-18 (`timed_20261118`)

The 2026-11-16 comparison recorded bowtie2's SITES but never its WALL TIME:
`three_way.tsv` was empty and stdout was not captured, so the only evidence
was log-file mtimes. This run fixes that. oligoscan 0.9.4, the same 100
guides, `bench_three.py` so the bowtie2 invocation is byte-identical to the
published comparison.

### Timings

| gate | oligoscan 22t | oligoscan 4t | bowtie2 4t |
|---|---|---|---|
| total<=2 + 2 nt | 15.9 s | 70.6 s | 4,095.4 s |
| total<=3 + 3 nt | 17.2 s | 71.7 s | 3,916.7 s |
| total<=4 + 4 nt | 22.5 s | 87.8 s | 3,962.1 s |

oligoscan peak RSS 931 MB at either thread count; bowtie2 35.3 GB (8l).

| gate | as each can be run (os 22t / bt2 4t) | thread-matched (4t each) |
|---|---|---|
| total<=2 + 2 nt | **258x** | **58x** |
| total<=3 + 3 nt | **228x** | **55x** |
| total<=4 + 4 nt | **176x** | **45x** |

**Report both.** The as-run pair is what a user experiences; the
thread-matched pair preempts the referee question. bowtie2 at 4 threads is
not a handicap: per 8l its `-a` buffers cost ~9 GB/thread nearly
independently of read count, so 22 threads wants ~190 GB against 110
available. The 100-guide run of 2026-11-16 survived `-p 22` only because
each thread saw ~5 reads and none was pathological -- luck, not headroom, and
not a publishable baseline. 8l puts the memory-only ceiling near 11 threads,
so 4 threads understates bowtie2 and the as-run ratio errs in its favour.

### Cross-validation

8l measured bowtie2 alone, 100 guides, 4 threads, `-S /dev/null`: 1:06:18 =
3,978 s. This run, gate 4, inside the harness: 3,962.1 s. **Within 0.4% from
two different harnesses** -- the figure is reproducible.

### Three things this run also establishes

- **Recall reproduced**: 99.32% / 97.17% / 87.07% against the published
  99.32 / 97.17 / 87.05, and site counts identical to 0.9.3 (5,773 /
  101,981 / 1,552,873). Third independent confirmation that the 0.9.4 fix is
  inert where maxI = maxE = maxTotal (8p).
- **Measured cost of the 0.9.4 fix** at these gates, scan time only:
  12.92 -> 14.95 s, 15.37 -> 16.09 s, 16.77 -> 21.45 s, i.e. **1.05x to
  1.28x for byte-identical output**. The honest price of making the gate a
  constraint on the search rather than a filter on its output.
- **bowtie2's output is thread-dependent.** 4 threads: 6,305 bowtie2-only,
  200,799 oligoscan-only. 22 threads (2026-11-16): 7,017 and 201,022. Same
  tool, same data, same parameters, different alignment set. Worth a sentence
  in the paper: it reinforces that the shortfall is architectural, and it is
  a reproducibility caveat for anyone benchmarking against bowtie2.

### oligoscan parallel efficiency, incidentally

22t -> 4t costs 4.44x / 4.17x / 3.90x for a 5.5x reduction in threads, i.e.
~75-80% parallel efficiency, consistent with the map-over-positions claim in
the Discussion.

### Correction recorded

On first reading these results I called the 4-thread bowtie2 configuration
unfair and the as-run ratio unpublishable, having mistaken 8l's memory
finding for one that applied only at 1,000 guides. It applies at any read
count; it is the per-thread buffer that is large. The setup was right.

## 8t. GPU indel arm: a serial host loop, found by a measurement that went the wrong way

**Measured 2026-11-18 (`gpu_arm_20261118`). Supersedes 8j for gapped gates.**

Two defects, one in the harness and one in the scanner's host code.

### 1. The GPU arm was never run on gapped gates

`bench_three.py` guarded it with `if a.gpu and maxI == 0:`, and `oligoscan-gpu`
is the mismatch-only popcount kernel which cannot place a bulge at all. So
every indel run silently reported no GPU column -- `gpu_s` was empty in every
gapped `three_way.tsv` and nobody noticed. Patched: the gate now selects the
binary, `--gpu` for maxI == 0 and `--gpu-cuda` (`oligoscan-cuda --gpu`)
otherwise, with the byte-identity check preserved on both paths.

### 2. Traceback ran on ONE host thread on the GPU path

With the arm finally running, the GPU **lost** at the widest gate. The cause
was not the kernel: the candidate loop in the `use_gpu && maxI > 0` branch had
no `#pragma omp`, so traceback ran single-threaded while 22 cores sat idle.
Harmless on 0.9.3, when traceback was cheap. From 0.9.4 the gate is a
CONSTRAINT on the traceback DP, whose state is (maxI+1)*(maxE+1) times larger
-- 25x at maxIndelBases = maxIndelEvents = 4 -- and the serial loop dominated.

| gate 4+4, 100 guides | CPU scan | GPU scan |
|---|---|---|
| 0.9.3 (8j) | 17.80 s | 6.23 s |
| 0.9.4, serial host loop | 21.83 s | 53.62 s |
| 0.9.4, parallel host loop | 21.98 s | **24.96 s** |

The CPU grew 1.2x while the GPU grew 8.6x; only serial host work explains
that asymmetry. Fixed by mirroring the CPU path: per-thread `CandVec`,
`#pragma omp for schedule(dynamic, 8192)`, merge under `critical`. Safe
because candidates are independent, the traceback scratch is thread-local,
and `all` is sorted before selection, so thread ordering cannot change the
output -- confirmed byte-identical to the CPU build at all three gates.

### Current GPU scan speedup, 0.9.4, 100 guides

| gate | CPU scan | GPU scan | speedup |
|---|---|---|---|
| total<=2 + 2 nt | 15.36 s | 1.43 s | **10.7x** |
| total<=3 + 3 nt | 16.21 s | 3.43 s | **4.7x** |
| total<=4 + 4 nt | 21.98 s | 24.96 s | **0.9x** |

Wall clock: 16.2 -> 2.4 s, 17.0 -> 4.4 s, 22.8 -> 25.9 s.

The GPU wins decisively at narrow indel gates and is level with the CPU at
the widest. This extends 8j's finding rather than contradicting it: 8j
already showed the scan speedup collapsing at the widest gate because
emission goes through a single global `atomicAdd`. On 0.9.4 the larger
traceback state pushes that collapse below parity. Remaining GPU tuning is
possible -- the kernel returns candidates for a whole query across the whole
reference, so host traceback loses the tiled locality the CPU path enjoys --
but it is not needed for the paper.

### Why this one matters

The GPU arm was added because it was expected to improve oligoscan's ratio
against bowtie2. It did the opposite, and that is what exposed the bug. Had
the GPU won at gate 4 by any margin, the serial loop would still be there,
the Discussion's "14.4 times faster" would have looked plausible, and the
paper would have reported a CPU/GPU crossover that does not exist.

**Seventh finding, and the first in the parallel plumbing rather than in a
comparison script.** The scanner's arithmetic has still never been wrong.

### Manuscript consequences

- The Discussion's "running the indel path 14.4 times faster" is a 0.9.3
  figure from the 2,361-oligo total<=3 + 3 nt gate. Replaced with the
  gate-dependent numbers above.
- **8j's GPU table is stale for every gapped gate** and must be re-run on
  0.9.4 before anything from it is quoted. Its mismatch-only rows are
  unaffected (that path does not call the traceback).

## 8u. Where the GPU's gate-4 time actually goes — measured 2026-11-18

Instrumented `osg_scan` and the host traceback loop separately, because from
outside the three candidate explanations were indistinguishable. 100 guides,
total<=4 + 4 nt:

```
GPU 2,166,882 candidates (34.7 MB moved)
kernel 3.33 s | D2H copy 0.03 s | host traceback 14.94 s
   18%        |      0%         |      82%
```

**The constrained traceback DP is the cost, and it is host-side, so both
backends pay it identically.** Three hypotheses died here, all of them
plausible and all of them wrong:

| hypothesis | killed by |
|---|---|
| random candidate order wrecks cache locality | 0.9.3 had the same ordering and ran gate 4 in 6.23 s |
| full INF fill of the DP array (644 GB of writes) | measured 4% on CPU; 281 KB is L2-resident, so the DRAM-bandwidth arithmetic was wrong |
| PCIe transfer / VRAM pressure at scale | 34.7 MB moved, 0.03 s, and device buffers are allocated once (~4.0 GB of 18 GB) regardless of library size |
| contention on the kernel's single global atomicAdd (8j) | kernel is 18%, and candidates are 2.17M not the tens of millions assumed |

Two further facts worth keeping:

- **Candidates are only 1.43x sites** (2,166,882 -> 1,510,546). Gaps do not
  admit the flood of extra end positions that was assumed; the Myers filter
  is already highly selective, so there is little to gain by filtering more.
- **The Myers kernel is only ~2.4x faster than the CPU Myers scan** at this
  gate (3.33 s against ~8.2 s), where the mismatch kernel is 11-15x. The
  indel kernel does far more work per position.

### Consequence

At the widest gate the GPU's ceiling is set by host traceback, not by the
kernel. Optimising traceback helps BOTH backends and is the only lever that
matters there. The 2026-11-18 targeted-fill change is one instance: worth 4%
on CPU but 16% on GPU, because traceback is a quarter of CPU time and nearly
all of GPU time. That alone moved gate 4 at 100 guides from 0.88x to 1.08x.

**Remaining lever, not taken:** events <= gap bases always, so only 15 of the
25 (g,e) planes are reachable; skipping the rest would cut up to 40% of the
DP. It needs care -- a reachable state reads (g-1,e) and so can touch the
e = g+1 boundary, which must then be pre-set to INF -- and the paper's claims
do not depend on it. Recorded as available, deliberately not done.

## 8v. GPU scaling table, 0.9.4 FINAL — measured 2026-11-18 (`gpu_scaling_0.9.4`)

**Replaces 8j wholesale.** Byte-identical CPU/GPU output at all twelve cells,
up to 37.6M sites.

| oligos | gate | sites | CPU scan | GPU scan | scan x |
|---|---|---|---|---|---|
| 100 | total<=3 | 20,068 | 4.97 s | 0.44 s | 11.3 |
| 100 | total<=4 | 192,077 | 5.30 s | 0.43 s | 12.3 |
| 100 | <=3 + 3 nt | 94,102 | 18.14 s | 3.08 s | 5.9 |
| 100 | <=4 + 4 nt | 1,510,546 | 24.46 s | 21.48 s | 1.14 |
| 1,000 | total<=3 | 235,006 | 50.84 s | 4.27 s | 11.9 |
| 1,000 | total<=4 | 2,050,337 | 55.39 s | 4.36 s | 12.7 |
| 1,000 | <=3 + 3 nt | 1,052,765 | 194.37 s | 30.99 s | 6.3 |
| 1,000 | <=4 + 4 nt | 15,666,765 | 268.70 s | 227.27 s | 1.18 |
| 2,361 | total<=3 | 517,754 | 139.07 s | 10.06 s | 13.8 |
| 2,361 | total<=4 | 4,795,910 | 153.67 s | 10.11 s | 15.2 |
| 2,361 | <=3 + 3 nt | 2,488,108 | 499.54 s | 71.16 s | 7.0 |
| 2,361 | <=4 + 4 nt | 37,591,667 | 635.01 s | 603.12 s | 1.05 |

### Where the GPU's time goes, at the indel gates

| oligos | gate | candidates | kernel | D2H copy | host traceback |
|---|---|---|---|---|---|
| 100 | <=3+3 | 125,264 | 39% | 0% | 61% |
| 100 | <=4+4 | 2,166,882 | 18% | 0% | 81% |
| 1,000 | <=3+3 | 1,499,026 | 35% | 0% | 65% |
| 1,000 | <=4+4 | 22,781,180 | 25% | 0% | 74% |
| 2,361 | <=3+3 | 3,412,658 | 34% | 0% | 66% |
| 2,361 | <=4+4 | 54,370,309 | 35% | 0% | 65% |

**D2H copy is 0% everywhere** -- 1.37 s of 603 s at the largest cell, 870 MB
moved. The VRAM/PCIe hypothesis is dead, and so is 8j's attribution of the
gate-4 collapse to the kernel's global atomicAdd: the kernel is 18-39%.
Host traceback is 61-81% at every indel gate, and it is shared code, so
Amdahl bounds what the GPU can do there.

Candidates run a steady **1.43-1.45x the site count** at 4+4 and ~23,000 per
query regardless of library size: the Myers filter is already selective.

### Two cautions that must travel with these numbers

**1. Run-to-run spread is +-15%.** Two consecutive grids of essentially the
same binary differed by -12.6% to +15.0% on CPU scan time. Quoting these to
two decimals implies a precision the machine does not have, and the 4%
targeted-fill change is BELOW the noise floor.

Consequence: **gate 4+4 is PARITY, not a GPU win.** 1.05-1.18x sits inside
the variance band. The GPU's own 24.89 -> 21.48 s improvement does have a
mechanism (the fill fix cuts traceback, which is 81% of GPU time there, so
~13% was predicted and ~14% observed); the CPU's +11% in the same interval
has none and is noise. Report parity.

**2. The Myers kernel slows 2.3x per query at the largest library.**
33.9 -> 50.4 -> 78.0 ms per query at 100 / 1,000 / 2,361 oligos, while
candidates per query stay at ~23,000 -- identical work. The mismatch kernel
shows no such growth (4.3 ms/query at both ends), which rules out anything
structural and points at thermal or power throttling on the much heavier
Myers kernel under sustained load. The 2,361-oligo cells run last, after
~10 minutes of continuous GPU work.

So the large-library GPU figures are probably throttle-limited rather than
algorithm-limited. Anyone quoting them should say so, and 8j's claim that
"scan speedup rises with oligo count as the per-launch cost amortises" holds
for the mismatch gates but is confounded at the indel gates by this.

## 8w. The unreachable-plane lever, taken — measured 2026-11-23

8u recorded this as "available, deliberately not done". Implemented in
`traceback()`: a gap EVENT consumes at least one gap BASE, and g bases span at
most g events, so `e == 0` only alongside `g == 0` and otherwise `1 <= e <= g`.
At maxI = maxE = 4 that is 11 reachable planes of 25 — **56% of the (g,e)
table is unreachable state**. The pre-set loops, the recurrence, the final
best-state search and the walk now all iterate reachable planes only.

Exact, not a heuristic band: an unreachable plane has no predecessor chain
back to a seed, so it evaluated to INF under the old full iteration. Skipping
it cannot change the optimum, only the time spent proving it is INF.

**Two reads genuinely cross the boundary** and are guarded with `REACH`
rather than pre-filled, because a boundary fill would cost about a plane of
stores per g — most of what the pruning saves:

| branch | reads | unreachable when |
|---|---|---|
| gap CONTINUATION | `(g-1, e)` | `e == g` |
| gap OPEN | `(g-1, e-1)` | `e == 1, g > 1` |

The traceback walk's continuation test reads `(g-1, e)` too, and is guarded
the same way.

### Verification

- `--selftest` passes, including *traceback cost vs DP, region split sums to
  total* (0/2000 disagreements against the brute-force DP).
- Byte-identical output against the pre-change build at every gate tried:
  rabbit cDNA / 100 guides at 3:3 indel 1 and 2, 4:4 indel 2; rabbit / 20
  guides at 4:4 and 5:5 with indel 4 nt / 4 events.

### What it is actually worth

Traceback CPU-seconds, measured by timing `on_candidate` directly (human
cDNA, 1 guide, 8 threads):

| gate | indel | traceback old | new | cut |
|---|---|---|---|---|
| 6 | 2 nt / 2 ev | 44.61 s | 39.75 s | 11% |
| 5 | 4 nt / 4 ev | 16.39 s | 12.74 s | **22%** |

End to end on CPU: 0.98–1.08x, because traceback is only about a quarter of
CPU time. At the narrow gates used for real siRNA screening it is 1.00x —
at `total <= 3` only **1,165 candidates** reach traceback in an entire human
cDNA scan, so the path is filter-bound and this changes nothing.

**56% of planes did not become 56% of time.** The pruned cells were already
cheap: their predecessors are INF and both gap branches short-circuit on it,
while the reachable g = 0 plane carries the diagonal, which is the real work.
The plane count is an upper bound on cells, not on cost.

### Consequence for the GPU

8u measured gate 4 + 4 nt at kernel 3.33 s / D2H 0.03 s / host traceback
14.94 s. Applying the 22% cut to the traceback term alone projects
14.94 -> ~11.7 s and 18.3 -> ~15.0 s, i.e. **~1.2x end to end**, moving the
gate-4 GPU/CPU ratio from 1.08x toward ~1.3x. **Projection, not a
measurement** — needs a CUDA box:

```
make gpu && ./oligoscan-gpu <ref>.osr <queries>.tsv 21 4 4 anti 4 4 1 --gpu
```
and compare the `[oligoscan] GPU ... % kernel, % copy, % traceback` line.

**This does not move a headline claim.** It is a modest, verified, free
improvement in the one regime that was the GPU's ceiling; no figure in the
manuscript depends on it, and none should be rewritten around it.

## 8x. The ASO/genome OOM: a 32-bit `long` on Windows — fixed 2026-11-23

**Symptom.** The ASO tab, and `oligoscan-v3.exe` directly, reached 128 GB of
RAM within seconds against the human genome. The siRNA tab was unaffected.
Bisected cleanly by the user: `indel 0` finished in 3 s, `indel 1` exploded.

**Cause.** `long` is 64-bit on Linux (LP64) and **32-bit on Windows**
(LLP64), and the shipped engine is a MinGW cross-build. `on_candidate()`
computed genomic positions in `long`, so every position past 2^31 --
2.147 Gb, the last third of a human genome -- wrapped negative:

| position 2.5 Gb in | Linux | Windows |
|---|---|---|
| `(long)endpos + 1 - L` | 2,499,999,982 | -1,794,967,314 |
| resulting `wn` | 20 | garbage, up to ~1.5e9 |

`wn` sizes the traceback DP buffer, `(wn+1)*(m+1)*GI*GE*3` cells of 8 bytes.
That buffer is `__thread`, grown on demand and never shrunk, so each of 22
threads inflated to multiple GB. Where the garbage `wn` did not allocate, the
`wstart < 0` guard returned instead -- so the indel path **silently dropped
candidates past 2.147 Gb**.

**Why nothing caught it.** Every automated test, every benchmark and the
whole manuscript run on human cDNA (764 Mb), which never reaches 2^31. The
genome is used only by the ASO path, only interactively, and only the indel
path calls `traceback()` at all. All development measurement was on Linux,
where `long` is 64-bit and the bug cannot exist.

**Fix.** Five position computations in `on_candidate()` moved to `int64_t`
(`us`, `wstart`, `tx0`, the `wn` expression, the `wstart + wn` bound check),
plus the unpack loop at load, which would overflow the same way above 8.6 Gb.
`tx_at()` results and `select_loci()` offsets are array indices and
within-sequence offsets, both safe, and were left alone.

Added with it: `if (wn64 > MAXL + BANDMAX) return;`. By construction
`wn <= L + band`, so this never fires -- but it converts any future
arithmetic slip into one dropped candidate instead of a multi-gigabyte
allocation per thread, which is precisely how a wrong number became an
out-of-memory event here.

**Verification.** `--selftest` passes on Windows (0/4000, 0/2000). Linux
output byte-identical before and after on chr1, a 194-sequence split of
chr1, and human cDNA -- expected, since Linux never had the bug. On Windows,
5 oligos vs the whole genome at `9 9 5 anti 1 1`: **41,104 sites in 4.01 s,
peak ~3.3 GB**, against 43,425 predicted by scaling a chr1 measurement
(within 5%). A run still dropping the far third would have returned ~29,000.

**Scope of bad results.** Any *gapped* ASO scan against a *whole genome* run
on *Windows* was missing hits, roughly the last third of the reference, with
no warning. Mismatch-only genome scans are unaffected (`traceback()` is never
called). Everything cDNA-based is unaffected: all siRNA work, the PCSK9 case
study, and every benchmark in the manuscript.

**Eighth finding, and the second in plumbing rather than arithmetic.** The
scanner's numerics have still never been wrong; this was a C type whose width
differs between the platform it was developed on and the platform it ships
on.

### Also fixed in the same pass

`fgets` on the query file now skips a UTF-8 BOM. Excel and Notepad write one
by default, and on a header-less one-column file those three bytes landed in
the first oligo's sequence, which was then rejected as non-ACGTU -- one
warning in a long log, and a silently four-oligo run. The GUI was never
affected (it reads CSVs as `utf-8-sig`); the command line was. Confirmed by
the run above reporting `5 queries` where the previous build reported 4.
