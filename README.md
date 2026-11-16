# oligoscan

**Exhaustive, uncapped off-target scanning for short oligonucleotides
(siRNA / ASO) against cDNA and whole-genome references.**

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

## Quickstart

```bash
./oligoscan examples/demo_ref.fa examples/demo_oligos.tsv 1 3 4 anti 1 1
```

A 60 kb bundled reference with planted off-target sites, so this runs in
under a second with nothing to download. See `examples/README.md`.

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
| 5′ terminal base | 1 | counts toward non-seed and total |
| seed | 2–8 | `maxSeedMM` |
| body | 9..L | `maxNonseedMM` |

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

### Against bowtie2

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
- **Memory is one byte per reference base during the scan** (~0.8 GB for
  human cDNA, ~3.1 GB for GRCh38). Scanning the packed 2-bit form directly
  would quarter it and is not yet implemented.
- **Non-ACGT reference bases count as mismatches**, so a window containing
  `N` is unlikely to pass a tight gate. bowtie1 instead declines to align
  across ambiguous bases; the two agree everywhere else.
- **No index, so no sub-second single-oligo lookup.** A full-reference scan
  for one oligo is ~2 s.
- Version 0.9: the interface may still change.

## Citation

A manuscript is in preparation. Until then please cite the software; see
`CITATION.cff`.

## Related work

| tool | relation |
|---|---|
| [Cas-OFFinder](https://pmc.ncbi.nlm.nih.gov/articles/4016707) | exhaustive mismatch *and* bulge search for short guides, GPU/OpenCL — the closest algorithmic precedent, but CRISPR/PAM-oriented |
| [RIsearch2](https://pubmed.ncbi.nlm.nih.gov/28108657/) | suffix-array, energy-based large-scale siRNA off-target prediction; seed-and-extend, so heuristic with respect to the full interaction space |
| bowtie1 / bowtie2 | general read mappers; `-v` is exhaustive but capped at 3 mismatches, and bowtie2's multiseed is heuristic by design |

To our knowledge no existing tool is a transcriptome-aware, seed-structured,
**uncapped** exhaustive scanner covering both siRNA and ASO.

## License

**PolyForm Noncommercial License 1.0.0** — see `LICENSE`.

oligoscan is **free of charge for academic and non-profit research**. The
licence permits use by educational institutions, public research
organisations, charitable organisations and government institutions
*regardless of the source of their funding* — so a university group working
under industry sponsorship is covered.

**Commercial use requires a separate licence.** This includes use by
pharmaceutical and biotechnology companies, contract research organisations,
and any for-profit entity. Enquiries: _[contact — open item]_.

This is a **source-available** licence, not an open-source one: the
noncommercial restriction is incompatible with the OSI definition. Source,
tests and benchmarks are public and reproducible, and academic users need ask
no permission.

`scripts/case_study.py` additionally depends on ViennaRNA, which carries its
own academic-use licence. The scanner itself has no third-party dependencies.
