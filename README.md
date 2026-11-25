# oligoscan

**Exhaustive, uncapped off-target enumeration for short oligonucleotides —
siRNA, antisense oligonucleotide and CRISPR guide RNA — against cDNA,
pre-mRNA and whole-genome references.**

oligoscan enumerates *every* alignment of an oligo within a chosen mismatch
and indel budget. There is no mismatch cap, no seeding, no filtering
heuristic and no effort limit, so recall is exhaustive by construction —
and runtime is independent of the error budget.

```
                     ≤3 mismatches   ≤6 mismatches   +2 nt indels
bowtie1 -v              yes            cannot           cannot
bowtie2             heuristic        heuristic        heuristic
oligoscan             exact            exact            exact
```

New here? **[docs/TUTORIAL.md](docs/TUTORIAL.md)** is the walkthrough:
install, the three workflows, how to read the output, and the conventions
that change your answer.

## What it is for

Three applications, one scanner. The gate is defined over oligo positions
rather than over any biological feature of the target, and a reference is
just a packed two-bit file, so nothing about the method is specific to one
of these.

| | what it answers | status |
|---|---|---|
| **siRNA / ASO off-targets** | every transcript an oligo can bind in a stated budget | exact agreement with `bowtie -v 3`; contains every bowtie2 alignment |
| **siRNA duplexes** | both strands, because a RISC-loaded passenger strand silences its own off-targets | one command, `scripts/sirna_duplex_screen.py` |
| **CRISPR guide RNAs** | genome-wide off-target sites, PAM-free then filtered | exact agreement with Cas-OFFinder across 2,355,073 and 12,694,694 sites |
| **ASOs vs pre-mRNA** | intronic sites no cDNA screen can represent | 90.1% of a 2,361-oligo library's transcribed-orientation sites are intronic |

## Why

A 21-mer carries 42 bits of information; locating one position in 560 Mb of
cDNA needs about 29.7. At 4–6 mismatches only 15–17 bases are still
informative, and every filtering scheme loses its selectivity exactly there:

- **q-gram lemma** — a match with ≤k errors shares ≥ L−q+1−kq common
  q-grams, which is *negative* for L=21 at any useful q (q=5, k=4 → −3).
- **pigeonhole partition** — k+1 blocks, one of which must match exactly;
  selectivity against 560 Mb needs blocks ≥ ~13 nt, i.e. reads of
  L ≥ 13(k+1) = 65 nt at k=4. A 21-mer cannot supply it.

With nothing to filter *with*, the honest method is to examine every
position — and for short oligos that is cheap. One byte per reference base,
an inner loop that vectorises to one instruction per 32–64 positions, tiled
so the reference streams once rather than L times. Indels use Myers'
bit-vector edit distance (one 64-bit word per position for patterns ≤64 nt),
which is exact, so it filters candidates with no recall loss.

What that buys structurally rather than by tuning: no error cap, 100% recall
by definition, runtime independent of k *and* of how repetitive the
reference is, and no `-k` truncation, seed length, gap barrier or effort
limit to mis-set.

## Install

**Docker** (also installs bowtie1 and bowtie2 for validation):

```bash
docker compose build          # compiles, then runs the full test suite
docker compose run --rm oligoscan oligoscan --version
```

**Or directly** — one file, no dependencies beyond a C compiler and OpenMP:

```bash
gcc -O3 -march=native -fopenmp -o oligoscan src/oligoscan.c
./oligoscan --selftest
```

`-march=native` targets the building CPU; use `-march=x86-64-v3` (AVX2) for a
portable binary.

Run `--selftest` after every build. It checks the edit-distance core against
a brute-force O(mn) dynamic program and verifies that the per-region
mismatch split sums to the reported total. It needs no input files and takes
seconds. A build that fails it must not be used — the first build of this
code disagreed with the DP on 3,531 of 4,000 random trials, and a hit count
alone would never have shown it.

**CUDA backend** (optional, mismatch gates at 10–32 nt):

```bash
docker compose build oligoscan-gpu
docker compose run --rm oligoscan-gpu oligoscan-gpu --version
```

It reads only packed `.osr` references — run `pack` first. Output is
byte-identical to the CPU build.

## Quickstart

```bash
./oligoscan examples/demo_ref.fa examples/demo_oligos.tsv 1 3 4 anti 1 1
```

A 60 kb bundled reference with planted off-target sites, so this runs in
under a second with nothing to download. See `examples/README.md`.

**Screening a whole siRNA design set, both strands, one command:**

```bash
python3 scripts/sirna_duplex_screen.py \
    --sirna-csv       designs.csv \
    --ref             ref/Homo_sapiens.GRCh38.cdna.all.osr \
    --precheck-fasta  ref/Homo_sapiens.GRCh38.cdna.all.fa \
    --tx2gene         ref/cdna_tx2gene.tsv \
    --target-gene     PCSK9 \
    --outdir          screen_out
```

`designs.csv` needs an ID, a sense and an antisense column; the names are
auto-detected and `U` is read as `T`. Both strands are screened in one scan,
each yielding a mismatch-only and an edit-distance summary, and the result is
one table with a verdict per design naming every criterion it failed.
`--resume` re-filters at new thresholds without rescanning.

## Usage

```
oligoscan <ref.fa|ref.osr> <oligos.tsv> <maxSeedMM> <maxNonseedMM> <maxTotalEdits> \
          [anti|sense|both] [maxIndelBases=0] [maxIndelEvents=1]
oligoscan pack ref.fa ref.osr      # one-time; ref.osr loads ~20x faster
oligoscan --selftest
oligoscan --version
```

`oligos.tsv` is `name<TAB>sequence`, one per line; `U` is accepted and
lengths may be mixed freely between 10 and 64 nt. Lower case and `U` are
accepted. Thread count follows `OMP_NUM_THREADS`.

### The gate

Mismatches are counted in three **disjoint** regions of the antisense oligo
and gated separately:

| region | positions | gate |
|---|---|---|
| 5′ terminal base | 1 | **its own region — counts toward the total only, never toward the body** |
| seed | 2–8 | `maxSeedMM` |
| body | 9..L | `maxNonseedMM` |

So `total_mm = pos1 + seed_mm + nonseed_mm`, and the output's `seed_mm` and
`nonseed_mm` columns exclude position 1 entirely. That separability is not
cosmetic: a position-1 mismatch is conventionally discounted in RISC, and
because the three regions are gated independently you can recover the
discounted distance as `seed_mm + nonseed_mm` without rescanning. Nothing
else in this comparison class exposes that — a `bowtie -v 3` has no notion of
guide position at all.

`maxTotalEdits` caps mismatches *plus* indel bases. The seed region is an
siRNA/RISC concept, so for a chemistry where it does not apply — a gapmer
ASO, say — leave `maxSeedMM` unconstrained (set it to the oligo length) and
gate on the body and total instead. Keeping the regions
separate is the point: a seed-perfect site with five body mismatches and a
seed-mangled site with five body matches share a total and share nothing of
the risk. The seed gate is also what makes a 6-mismatch search usable — at
total ≤6 it cuts sites per oligo 15.6-fold (865 → 55) relative to an
unconstrained seed.

So the deliverable is not "the cap is now 6". It is *the cap is gone, and
the gate is what you choose*.

### Conventions that change your answer

Four choices change the result more than most parameters do, and all four
belong in a methods section. On a complete 3,617-design tiling of PCSK9:

- **Position-1 mismatches.** Discounting them — the field's convention, and
  the conservative direction, since it can only move off-targets closer —
  re-tiers 36% of designs and cuts a "clears three mismatches" shortlist
  from **816 to 189**.
- **Mismatches vs edits.** Run these as separate analyses. A mismatch-only
  screen counts mismatches; an indel-permitting screen should count **edits**,
  an indel and a mismatch each costing one. `edit_total` is already that
  number, and `edit_total = total_mm + gap_events`. Of the 189 designs above,
  admitting one indel base leaves **26**.
- **Low buckets are exact; high ones are floors.** At `maxTotalEdits 4`, a hit
  4 away excluding position 1 *plus* a position-1 mismatch totals 5 and is
  outside the gate. So 0–3 are counts and anything above is a lower bound.
- **Seed matches are a separate axis from total distance.** Across that tiling
  the correlation between nearest-off-target distance and perfect-seed
  off-target count is only **r = −0.36**; within the designs clearing three
  mismatches the perfect-seed count still ranges 0 to 49. Ranking on distance
  alone admits seed-catastrophic designs.

One scan serves both analyses: `class == mismatch_only` (equivalently
`gap_events == 0`) is the mismatch-only subset of an indel-permitting run, so
there is no need to scan twice.

### Output

TSV on stdout, sorted by (oligo, sequence, position) so runs are
reproducible byte for byte at any thread count:

```
oligo  strand  seq_id  pos  seed_mm  nonseed_mm  total_mm  gap_events
       oligo_bulge  target_bulge  class  edit_total
```

`pos` is 0-based within the sequence. `class` is `mismatch_only`,
`oligo_bulge` (a base in the oligo with no partner), `target_bulge` (a base
in the transcript with no partner) or `both_bulges`. Progress and timing go
to stderr.

With `maxIndelBases > 0` exactly **one row per locus** is reported: the
best-scoring description, ranked by edit total, then gap bases, gap events,
seed / position-1 / body mismatches, then position. Equal cost prefers the
ungapped description — near the 3′ terminus a 1-nt bulge and a single
mismatch describe the same alignment equally well, and the edit total is
identical either way.

## Validation

Everything below is reproducible from this repository.

**Against bowtie1 `-v 3`, full Ensembl human cDNA.** bowtie1's `-v` mode is
exhaustive backtracking within its 3-mismatch cap, so it is a true oracle
there — and the only one available, since nothing else is exhaustive above
3 MM. 328,868 transcripts / 743 Mb, 40 real 21-mers, seed≤3 body≤3 total≤3:

```
bowtie1   : 1,917 sites
oligoscan : 1,917 sites
shared 1,917 · only bowtie1 0 · only oligoscan 0 · count disagreements 0
```

`tests/validate_oligoscan.py` — a set comparison, not a count.

**Against brute force.** `--selftest` checks the Myers recurrence against an
O(mn) DP before anything else runs (0/4000 disagreements), and the banded
traceback against the same DP (0/2000). A transcription error in the
recurrence is invisible in a hit count: the first build of this code
disagreed with the DP on 3,531 of 4,000 random trials.

**Selection and packing.** `tests/test_oligoscan_select.py` is
self-contained — no bowtie, no external data — and checks the mismatch path
against a numpy Hamming scan, the indel path against a gap-constrained
brute-force DP (every within-budget locus reported, each at the optimum
cost), one row per locus, transcript-junction safety, and byte-identical
output from a packed reference and across thread counts. The Docker build
fails if any of this fails.

### Against Cas-OFFinder

Whole GRCh38 primary assembly. Cas-OFFinder
searches mismatches and bulges exhaustively for short guides and is the
closest algorithmic precedent. The configuration deliberately favours it: an
NGG pattern is its native published case and restricts its candidate set
roughly sixteen-fold, while oligoscan applies no PAM constraint and must
enumerate the entire PAM-free neighbourhood and filter afterwards — at four
mismatches it reports and discards 6,958,129 sites to keep 2,355,073.

500 SpCas9 protospacers, stratified across genomic context and
error-neighbourhood size, drawn from 1,233 genes, both strands:

```
3 mismatches :   600,632 of   600,632 sites   agreement exact, both directions
4 mismatches : 2,355,073 of 2,355,073 sites   agreement exact, both directions
repeat-derived guides : 12,694,694 sites      agreement exact
```

Nothing in either tool's exclusive column, and exact agreement within each of
sixteen context-by-difficulty strata taken separately. Peak memory 3.69 GB
against 5.67 GB.

**CPU–GPU byte-identity.** The CUDA backend reproduces the CPU output byte
for byte at every gate and library size tested, up to 37.6 million sites.
This is evidence about the port, not about the algorithm: the two backends
share the candidate-handling path, so a defect there passes unnoticed — one
did, and was caught instead by cross-checking against another tool. It has
earned its place anyway, having found a four-row deficit in thirty-seven
million arising from ambiguous-base handling that no other test detected.

## Performance

Measured on the full Ensembl human cDNA (743 Mb), 40 21-mers, antisense,
22 threads, reference loaded from `.osr`. Machine: Intel Core Ultra 9 285
(24 cores / 24 threads, no SMT, 36 MB L3), 128 GB RAM, Docker on Windows.
This CPU has **no AVX-512**, so the scan vectorised to 32 bytes per
instruction rather than 64 — these figures are a floor, and a part with
AVX-512 should be faster on the same code.

| gate (seed / body / total, indel) | sites or loci | scan |
|---|---|---|
| 3 / 3 / 3 | 1,917 | 1.76 s |
| 1 / 3 / 4 | 7,607 | 1.75 s |
| 0 / 6 / 6 | 47,218 | 1.77 s |
| 1 / 5 / 6 | 207,262 | 1.81 s |
| unconstrained / 6 | 681,388 | 2.08 s |
| 1 / 3 / 4, 1 nt indel | 28,996 | 5.24 s |
| 1 / 3 / 4, 2 nt indel | 33,923 | 5.44 s |

**Flat across the gate** — 1.76 s at 3 mismatches, 2.08 s at 6 — and flat
across the indel budget: 1 nt and 2 nt cost the same. This is the property a
seeded aligner structurally cannot have. End to end at ≤3 MM, including
loading, oligoscan takes ~2.8 s against bowtie1's 7.6 s, so it is faster
even in bowtie1's best case; at 4–6 MM and with indels bowtie1 cannot run at
all.

Scaling linearly in oligo count, 1,750 oligos is ~76 s mismatch-only and
~230 s with indels on the same 22 threads.

`oligoscan pack` converts a reference to 2 bits per base (838 MB → 207 MB
for human cDNA), loading in ~1 s instead of ~20 s, with byte-identical
output. Run `scripts/bench_gates.sh` to reproduce the table.

### At genome scale, and on GPU

Runtime decomposes into two terms, which is why no single speed multiple is
quoted: a per-query scan cost and a per-site reporting cost.

```
oligoscan CPU  : 0.219 s/query + 1.49 us/site
oligoscan CUDA : 4.07 s fixed + 0.0171 s/query + 1.41 us/site
Cas-OFFinder   : ~51 s fixed  + 95-128 us/site
```

The CUDA kernel replaces the Myers scan only; traceback, the region gate and
selection stay on the host as shared code. So what the GPU is worth depends
on how much of the work it replaces: **11-15x** at mismatch-only gates where
the scan is nearly everything, 6-7x with a 3-nt indel, and level with the CPU
at 4 nt, where traceback is two thirds to four fifths of the work and no
kernel touches it. The claim is not a single speedup but that the method
admits acceleration at all, in proportion to the share of the work that is a
regular map over positions.

Cas-OFFinder's ~51 s of fixed cost amortises across guides, so the margin
depends entirely on scale: an eight-fold advantage at ten guides is a
statement about startup latency, while at 500 guides the honest figures are
1.56x at three mismatches and 2.49x at four. What does generalise is flatness
in the error budget — raising the gate from three to four mismatches
multiplied reported sites by 6.14 and cost oligoscan **3.0%** more time,
against Cas-OFFinder's **64.5%**.

### Against bowtie2

Two comparisons, in order of how much weight they carry.

#### At a verified-equivalent gate

The strongest form of the comparison: both tools search **one and the same
edit sphere**, both are allowed gaps, and the two site sets are compared
reciprocally at **exact transcript coordinates** rather than by counting.
bowtie2 is scored with `--score-min C,-k,0` and `--rdg 0,1 --rfg 0,1`, so a
gap of length N costs exactly N — matching oligoscan's edit model instead of
bowtie2's default affine penalty, which silently forbids a bulge plus two
substitutions at budget 3. One hundred published 19-nt guides, whole human
cDNA, 22 threads (oligoscan 0.9.3):

| search allowance | oligoscan sites | bowtie2 alignments | bowtie2 recall | bowtie2 precision | oligoscan | bowtie2 |
|---|---|---|---|---|---|---|
| ≤2 edits, ≤2 indel bases | 5,773 | 5,734 | **99.32%** | 100% | **14.1 s** | 2,460.1 s |
| ≤3 edits, ≤3 indel bases | 101,981 | 99,189 | **97.17%** | 99.90% | **15.3 s** | 2,381.2 s |
| ≤4 edits, ≤4 indel bases | 1,552,873 | 1,358,868 | **87.05%** | 99.48% | **17.2 s** | 2,411.8 s |

Sites found only by oligoscan: 39, 2,888, 201,022.

**Every bowtie2 alignment is contained in oligoscan's output.** Matching at
locus level — same guide and transcript, start positions within 8 nt —
oligoscan accounts for **1,358,868 of 1,358,868** bowtie2 alignments at the
widest gate, and for all of them at the other two: nothing bowtie2 finds is
absent from oligoscan's set. The 0, 96 and 7,017 alignments that do not match
at an *exact* start position are the same loci with the boundary drawn one or
two bases differently, not sites oligoscan missed.

Times are the scan, which is flat: a **269× rise in sites costs 1.23× more
time**, while bowtie2 runs 140–175× slower and is flat for its own reason
(`-a` with flat gap costs dominates).

This is the comparison to quote. It is also where the runtime gap is honest:
~520× at matched settings, because `-a` with flat gap costs is what exhaustive
enumeration actually costs bowtie2.

#### Earlier, wider-gate comparison

> **Provenance warning.** The figures below come from `example_oligos.tsv`,
> which is 40 consecutive 21-nt windows sliding one base at a time across
> ~60 nt of a single transcript — one sequence sampled 40 times, not 40
> independent queries. Treat the site counts as indicative only; the ≤6 row
> in particular runs well above the per-oligo average for that reason. They
> are kept for continuity and will be replaced by a library-scale re-run.
> The configuration claim below is also known to be wrong: a seed length one
> base longer than the "provable" bound scores 22 points higher recall and
> runs 23% faster, so `-L` needs sweeping rather than deriving.

Same reference, 10 21-mers. bowtie2 is given its best honest configuration,
not its defaults — which cannot express this search at all: scoring
flattened to edit distance, `--gbar 1`, `-a` rather than `-k`, `-N 1`, and
the longest seed length that can provably still nucleate a k-error
alignment. `scripts/bench_paper.py` reproduces it.

| edit budget | oligoscan | bowtie2 | bowtie2 recall |
|---|---|---|---|
| ≤4 | 5,600 sites, 1.5 s | 5,967 alignments, 1,444.5 s | **93.8%** |
| ≤6 | 252,831 sites, 2.5 s | 213,691 alignments, 1,709.6 s | **58.0%** |

Three qualifications, all of which make the gap wider rather than narrower:

- oligoscan ran **mismatch-only** while bowtie2 was allowed gaps (12% and
  33% of its alignments carry an indel), so bowtie2 searched a strictly
  larger space and still missed this much. These are upper bounds on its
  completeness.
- Some of what it misses lies **inside bowtie1's exhaustive range** — 4 and
  18 of the misses are 3-edit sites, which `bowtie -v 3` finds by
  construction and oligoscan matches exactly.
- Its recall is **non-monotonic in the budget**: of the 5,600 sites within
  ≤4 edits, bowtie2 covered 5,254 at budget 4 but only 4,671 at budget 6.
  Widening the search lost 735 sites it had already found and gained 152.
  The deficit is not a fixed, characterisable gap that a wider search
  repairs — it moves with the parameters.

This is not a knock on bowtie2, which is an excellent read mapper doing the
job it was designed for. Its multiseed filter is documented as heuristic,
and 21-mers at 4–6 errors are outside the regime it targets.

### Against RIsearch2

RIsearch2 (Alkan et al., *NAR* 2017) is the closest tool by purpose: a
suffix-array seed-and-extend engine built for large-scale siRNA off-target
prediction. It answers a different question — where a duplex is energetically
favourable, not where an oligo aligns within an error budget — so three
adjustments are needed before the numbers mean anything, and
`scripts/bench_risearch.py` applies all three explicitly:

- **G:U wobble.** RIsearch2 pairs G:U by default, so a "perfect seed" can be
  almost entirely wobbles. `--noGUseed` is always set. oligoscan counts G:U
  as a mismatch, so the two still are not identical.
- **Local vs full-length.** RIsearch2 reports the energy-optimal *local*
  duplex. In a planted-site control, a 19/19 alignment with 2 mismatches came
  back as a 16-nt duplex with 1 unpaired base at shifted coordinates.
  Matching is therefore by **locus** (same transcript, ≥ half the guide
  overlapping), which is *generous to RIsearch2*.
- **Energy, not edits.** The `-e` threshold has no principled value for
  enumeration, so it is swept rather than chosen.

100 library 19-mers, whole human cDNA, 22 threads, gate seed-perfect /
≤4 edits / 2 bulge bases (196,433 sites). Recall is of oligoscan's output,
so no recall figure is meaningful for oligoscan itself:

| `-e` | 0 edits | 1 | 2 | 3 | 4 | all | unmatched interactions | RIsearch2 | oligoscan |
|---|---|---|---|---|---|---|---|---|---|
| −30 | 60.8% | 38.9% | 27.0% | 6.3% | 1.2% | **2.5%** | 3,681 | 44.6 s | **14.6 s** |
| −25 | 99.4% | 62.8% | 49.5% | 21.1% | 9.2% | **11.5%** | 76,836 | 46.0 s | **14.6 s** |
| −20 (default) | 100% | 73.5% | 64.4% | 37.9% | 24.9% | **27.2%** | 652,043 | 46.3 s | **14.6 s** |
| −15 | 100% | 86.7% | 70.7% | 53.5% | 40.2% | **42.4%** | 2,521,758 | 46.2 s | **14.6 s** |
| −10 | 100% | 86.7% | 71.1% | 57.3% | 48.2% | **49.8%** | 5,178,169 | 47.0 s | **14.6 s** |

**RIsearch2 is fast, and the case against it is not speed.** Three to nine
times slower here, not three orders of magnitude — and that excludes a
one-time 292 s index build producing an 11.9 GB suffix array, where oligoscan
builds no index. Its runtime is also flat in `-e`, so flatness in a search
parameter is not unique to oligoscan; what differs is that oligoscan is flat
in the **error budget**.

**The sharpest result is the leftmost column, and it is not a recall
figure.** At `-e -30` RIsearch2 reports 60.8% of *perfectly complementary*
sites — 505 of 1,287 absent. Its search is fine; those are AU-rich 19-mers
whose duplexes never reach −30 kcal/mol. **An energy cutoff discriminates by
base composition, not by complementarity**, so a screen built on one is
biased toward GC-rich targets, and lowering the threshold to admit AU-rich
true sites admits everything else in proportion.

**Completeness is not reachable by relaxing the threshold.** From −20 to −10
recovers ~44,300 sites and adds 4,526,126 interactions matching nothing in
the budget — about **100 extra interactions per additional true site** —
while the recall curve flattens (+15.2 then +7.4 points).

Those unmatched interactions are **not errors**: partial duplexes and
wobble-containing pairings are exactly what RIsearch2 is designed to find and
oligoscan is not. They are what a user must triage to reach the sites an
exhaustive gate returns directly.

*Caveat.* The seed specification used here (`-s 2:8/7`) mirrors oligoscan's
seed region rather than the parameters the RIsearch2 authors recommend in
their own siRNA off-target pipeline. Those should be measured before these
figures are called that tool's best configuration.

## Limitations

- **Oligos 10–64 nt.** The upper bound is Myers' bit-vector needing the
  pattern in one 64-bit word; the lower is information-theoretic (a 10-mer
  carries 20 bits where a transcriptome needs ~30, so shorter queries cannot
  constrain the search). Lengths may be mixed freely in one file — each oligo
  is scanned at its own length, and the result is identical to running it
  alone. Anything out of range, or carrying a degenerate base, is **skipped
  with a reason on stderr** rather than silently truncated or altered.
- **Edit distance, not thermodynamics.** No energy model, no G:U wobble, no
  target accessibility: a G:U pair counts as a mismatch and every mismatch
  costs the same. oligoscan tells you *where* the near-complementary sites
  are, exhaustively; it does not predict how much knockdown each causes.
  Its output is a natural input to an energy-based scorer.
- **Counts are per transcript, not per gene.** Ensembl cDNA lists every
  isoform, so one genomic site appears once per isoform containing it.
  Collapse to gene level in reporting.
- **Memory is one byte per reference base during the CPU scan** (~0.8 GB for
  human cDNA, ~3.1 GB for GRCh38). The CUDA build scans the packed 2-bit form
  directly and needs 856 MB for the genome, so that saving is realised there
  but not yet in the CPU build.
- **The CUDA build covers mismatch gates only, at 10–32 nt.** Gapped gates on
  GPU use a separate binary in which the kernel replaces the Myers scan alone,
  leaving traceback on the host.
- **Counts are per sequence.** On an isoform-rich annotation one genomic site
  is counted once per isoform; the genomic runs use the primary assembly,
  which avoids the same multiple counting across alt haplotypes and patches.
- **Off-target counts are not comparable across species.** Annotation depth
  differs several-fold — human 328,868 transcripts, mouse 116,358, cyno
  51,376, rabbit 41,711 — so the same design shows fewer off-target genes in
  a less-annotated species for reasons unrelated to safety. Read each species
  against its own transcriptome. Cross-species *activity* is a different
  question, answered from the orthologue alone by
  `scripts/crossreact_scan.py`.
- **Non-ACGT reference bases count as mismatches**, so a window containing
  `N` is unlikely to pass a tight gate. bowtie1 instead declines to align
  across ambiguous bases; the two agree everywhere else.
- **No index, so no sub-second single-oligo lookup.** A full-reference scan
  for one oligo is ~2 s.
- Version 0.9.4: the interface may still change.

## Citation

A manuscript is in preparation. Until then please cite the software; see
`CITATION.cff`.

## Related work

| tool | relation |
|---|---|
| [Cas-OFFinder](https://pmc.ncbi.nlm.nih.gov/articles/4016707) | exhaustive mismatch *and* bulge search for short guides, GPU/OpenCL — the closest algorithmic precedent. [Measured above](#against-cas-offinder): exact agreement across 2,355,073 and 12,694,694 sites. Differs in being PAM-bound, in gating on one total rather than three regions, and in requiring an OpenCL runtime |
| [RIsearch2](https://pubmed.ncbi.nlm.nih.gov/28108657/) | suffix-array, energy-based large-scale siRNA off-target prediction; seed-and-extend, so heuristic with respect to the full interaction space — [measured above](#against-risearch2) |
| bowtie1 / bowtie2 | general read mappers; `-v` is exhaustive but capped at 3 mismatches, and bowtie2's multiseed is heuristic by design |

To our knowledge no existing tool is a transcriptome-aware, seed-structured,
**uncapped** exhaustive scanner covering siRNA, ASO *and* guide RNA against
both transcriptomic and genomic references.

## License

**PolyForm Noncommercial License 1.0.0** — see `LICENSE`.

Repository: <https://github.com/breezedu/OligoScan>

oligoscan is **free of charge for academic and non-profit research**. The
licence permits use by educational institutions, public research
organisations, charitable organisations and government institutions
*regardless of the source of their funding* — so a university group working
under industry sponsorship is covered.

**Commercial use requires a separate licence.** This includes use by
pharmaceutical and biotechnology companies, contract research organisations,
and any for-profit entity. Enquiries: Jeff Du, jeff.g.du@gmail.com.

This is a **source-available** licence, not an open-source one: the
noncommercial restriction is incompatible with the OSI definition. Source,
tests and benchmarks are public and reproducible, and academic users need ask
no permission.

`scripts/case_study.py` additionally depends on ViennaRNA, which carries its
own academic-use licence. The scanner itself has no third-party dependencies.
