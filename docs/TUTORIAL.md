# OligoScan tutorial

OligoScan enumerates **every** site in a reference where a short oligo could
bind within a stated edit budget. Not the best site, not a sample — all of them.
It is built for the three places that matters: siRNA, antisense oligonucleotide
(ASO), and CRISPR guide-RNA design.

Version documented here: **oligoscan 0.9.4**. Companion to `README.md`, which
covers the design and the benchmarks; this file is about using it.

---

## 1. What it is, and what it is not

**It is** an exhaustive scanner. Given an oligo and a gate (how many
mismatches, where, and how many indels), it reports every locus that satisfies
the gate. Nothing is capped, sampled, or heuristically pruned. If a site is
inside your budget, it is in the output.

**It is not** a design tool. It does not pick siRNAs, score potency, model
secondary structure, or predict knockdown. It answers one question — *where
else could this bind?* — and answers it completely.

Three things follow from that:

- **Output can be large.** A 19-mer at ≤3 mismatches plus one indel against the
  human transcriptome produces tens of millions of rows. Plan disk accordingly;
  a full siRNA tiling run is multiple GB.
- **The gate is the result.** "No off-targets" is meaningless without the gate
  that produced it. Every number OligoScan reports is complete *within a stated
  budget*, and the budget belongs in your methods section.
- **Conventions change the answer more than you expect.** See §8. This is not a
  footnote — on a real dataset the position-1 convention alone moved 36% of
  designs into a worse tier.

---

## 2. Install

### Build from source

```bash
gcc -O3 -march=native -fopenmp -o build/oligoscan src/oligoscan.c -lm
build/oligoscan --version      # oligoscan 0.9.4
build/oligoscan --selftest     # Myers bit-vector vs brute-force DP
```

Run `--selftest` after every build. It checks the edit-distance core against a
brute-force O(mn) dynamic-programming implementation on random data, and checks
that the per-region mismatch split sums to the reported total. It needs no input
files and takes seconds. A build that fails it must not be used.

### Docker

```bash
docker compose build oligoscan-casoffinder
docker compose run --rm oligoscan-casoffinder oligoscan --version
```

The compose services mount `./data` at `/data`, so reference and input paths
inside the container are `/data/...`. The binary is on `PATH`.

### Hardware

CPU is the normal path and scales with cores. A CUDA build
(`src/oligoscan_gpu.cu`) speeds the *scan* substantially but not the
*reporting*, so the benefit depends on how many sites your gate admits: a tight
gate is scan-dominated and gains a lot, a loose gate is report-dominated and
gains little. The GPU build reads only packed `.osr` references.

---

## 3. One-time: pack your reference

```bash
build/oligoscan pack ref.fa ref.osr
```

`.osr` is a 2-bit-per-base packed form. It loads far faster than FASTA and is
required by the GPU build. Use it for everything except the inverted-repeat
precheck in §4, which needs sequence text.

A caution on sequence names: **`seq_id` in the output is the first
whitespace-delimited token of the FASTA header**, nothing more. Ensembl cDNA
headers begin `>ENST00000302118.5 cdna chromosome:...`, so `seq_id` is
`ENST00000302118.5`. If your headers put something unhelpful first, fix them
before packing.

---

## 4. Workflow A — screen an siRNA duplex set (start here)

This is the main entry point and it is one command.

**Why both strands.** An siRNA is a duplex. If the sense (passenger) strand is
loaded into RISC it acts as a guide and silences *its own* off-targets.
Screening only the antisense strand covers half the molecule.

### Input

A CSV with an ID column, a sense column and an antisense column. Names are
auto-detected case-insensitively and a UTF-8 BOM is tolerated; override with
`--id-col` / `--sense-col` / `--antisense-col`. Sequences may be RNA or DNA — `U`
is read as `T`. Extra columns are ignored.

```csv
siRNA Name,Target (sense),Antisense
siRNA_0001,GCUGUAAGAGCUGAUAUGA,UCAUAUCAGCUCUUACAGCUU
```

### Run

```bash
python3 scripts/sirna_duplex_screen.py \
    --sirna-csv       designs.csv \
    --ref             data/ref/Homo_sapiens.GRCh38.cdna.all.osr \
    --precheck-fasta  data/ref/Homo_sapiens.GRCh38.cdna.all.fa \
    --tx2gene         data/ref/cdna_tx2gene.tsv \
    --target-gene     PCSK9 \
    --oligoscan       build/oligoscan \
    --outdir          screen_out
```

`--tx2gene` is a TSV of `transcript, gene_id, gene_symbol, ...` — it is what
turns transcript hits into gene-level counts. `--target-gene` takes a symbol or
an ENSG id and is excluded from off-target totals for both strands.

### What it does

1. Parses and validates the CSV, checking that
   `antisense[:19] == revcomp(sense[:19])` for every design and reporting any
   that fail — a mis-paired duplex or a different overhang convention.
2. Deduplicates identical cores, so a shared sequence is scanned once.
3. Runs the **inverted-repeat precheck** (§6).
4. Runs **one** oligoscan pass with both strands as `as_*` / `se_*` queries, so
   the reference is loaded and traversed once.
5. Aggregates to gene level, splitting by strand, and verifies every query is
   accounted for — a truncated scan aborts rather than producing quiet
   undercounts.
6. Writes `duplex_screen_results.csv` and prints a per-criterion funnel.

### Resuming

The scan is the slow step. `--resume` reuses an existing complete `scan.tsv`,
so you can re-filter at new thresholds in seconds. Threshold changes never
require a rescan.

---

## 5. Workflow B — a single oligo set (ASO, or one strand)

For ASOs, or when you want raw loci rather than a gene-level verdict, call
oligoscan directly.

```bash
build/oligoscan ref.osr oligos.tsv 3 3 4 anti 1 1 > hits.tsv
```

`oligos.tsv` is two columns, no header: `id<TAB>sequence`.

```
oligoscan ref oligos.tsv maxSeedMM maxNonseedMM maxTotal [anti|sense|both] [maxIndelBases] [maxIndelEvents]
```

| argument | meaning |
|---|---|
| `maxSeedMM` | mismatches allowed in the seed (positions 2–8) |
| `maxNonseedMM` | mismatches allowed in the body (position 9 onward) |
| `maxTotal` | total edits allowed; the binding constraint |
| mode | `anti`, `sense`, or `both` — see below |
| `maxIndelBases` | `0` for a mismatch-only scan |
| `maxIndelEvents` | distinct gap events allowed |

**Mode is the thing people get wrong.** The three regions (position 1, seed 2–8,
body 9+) are numbered in the **input oligo's** own coordinates, and mode decides
what gets searched:

- **`anti`** — searches for the **reverse complement** of your query and numbers
  positions from the query's 5′ end. This is RISC numbering, and it is correct
  for *any* strand acting as a guide: an siRNA antisense strand, an siRNA sense
  strand, a gapmer ASO.
- **`sense`** — searches for the **literal** sequence as given.
- **`both`** — runs each query both ways.

Use `anti` for anything that will be loaded into RISC. Use `sense` only when you
genuinely want direct sequence occurrences. Feeding an siRNA sense core to
`sense` mode just re-finds the target mRNA and tells you nothing about passenger
risk.

Under a flat gate (`maxSeedMM == maxNonseedMM`) the region split does not affect
which hits pass, so mode affects only the search direction. Under an uneven gate
it affects both.

---

## 6. The inverted-repeat precheck

The sense strand normally has **no** on-target: it carries the target mRNA's own
sequence, so as a guide it seeks sites *complementary* to itself, which the
target does not contain.

The exception is real. If the transcript contains an **inverted repeat**
spanning a design site — the sequence property that forms a hairpin — then the
antisense sequence *does* occur in the target, and the sense strand has a
genuine on-target that must not be scored as an off-target.

The precheck establishes that count independently, before the scan:

```bash
python3 scripts/case_study_selfcomp_precheck.py \
    --cdna data/ref/Homo_sapiens.GRCh38.cdna.all.fa \
    --tx2gene data/ref/cdna_tx2gene.tsv \
    --oligos screen_out/oligos_both.tsv \
    --gene-symbol PCSK9
```

This makes the orientation check decisive rather than suggestive. The two
failure modes are separable by **prevalence and position**:

- **Inverted orientation (a bug):** nearly every design hits the target gene at
  0MM, each at its own design site.
- **Inverted repeat (biology):** a handful of designs hit, at positions other
  than their own site.

The driver runs this automatically when given `--precheck-fasta` and compares
the scan's result against it.

Two scope limits worth knowing: it tests full-length self-complementarity
only, ungapped. mRNA hairpin stems are typically 5–15 bp, so secondary structure
is common while a full-length inverted repeat is rare — and a stem covering only
the seed (≥7 bp) is a seed-match question, which a full-length screen does not
address.

---

## 7. Workflow C — CRISPR guide off-targets

OligoScan has no notion of a PAM; it is a sequence scanner. Guide screening is
therefore two steps: scan, then filter by PAM.

```bash
build/oligoscan genome.osr guides.tsv 4 4 4 both 0 > guide_hits.tsv
python3 -c "
import sys; sys.path.insert(0,'scripts')
from grna_pam import Faidx, pam_at, pam_ok
fx = Faidx('genome.fa')          # needs genome.fa.fai alongside
for line in open('guide_hits.tsv'):
    p = line.rstrip().split('\t')
    if p[0] == 'oligo': continue
    pam = pam_at(fx, p[2], int(p[3]), p[1], guide_len=20)
    if pam and pam_ok(pam, 'NGG'): print(line, end='')
"
```

The PAM convention, which is easy to invert:

- strand `+` → PAM at `seq[pos+L : pos+L+3]`, must match `NGG`
- strand `−` → PAM at `seq[pos-3 : pos]`, must be `CC` (reverse complement of
  `NGG` is `CCN`)

If you are comparing against Cas-OFFinder, its reported location differs from
OligoScan's `pos` by **0 on the plus strand and −3 on the minus strand**. Getting
that sign wrong produces a clean-looking disagreement that is purely
bookkeeping.

`scripts/bench_grna_casoffinder.py` is a **benchmark harness**, not a user tool —
it exists to establish agreement with Cas-OFFinder, and it runs Cas-OFFinder
too. There is no one-command guide screen yet; this section is the supported
route.

---

## 8. Five conventions that change your answer

These are not stylistic. Each materially moves results, and each belongs in a
methods section.

**1. Position-1 mismatches.** The guide's 5′ nucleotide is anchored in Ago2's
MID domain and conventionally excluded from mismatch counts. OligoScan's gate
keeps position 1 in its own region precisely so you can address it:
`total_mm = pos1 + seed_mm + nonseed_mm`, so the pos-1-discounted distance is
`seed_mm + nonseed_mm`. On a 3,617-design tiling this moved **36%** of designs
to a closer nearest-off-target and cut a "≥3 mismatches" shortlist from 816 to
189. Discounting is the conservative direction: it can only move off-targets
closer. The driver discounts by default; `--count-pos1` turns it off.

**2. Mismatches versus edits.** Run these as separate analyses. A mismatch-only
screen counts mismatches. An indel-permitting screen should count **edits**,
where an indel and a mismatch are each one edit, rather than reporting a
mismatch count plus a separate "needed an indel" flag — the latter is two
incommensurable dimensions and invites double-counting. `edit_total` in the
output is already this number, and `edit_total = total_mm + gap_events`.

**3. Indels at position 1.** Discounting applies to *mismatches*. An indel
shifts the register and is not a position-1 mismatch, so OligoScan counts indels
wherever they fall.

**4. Buckets 0–3 are exact; higher buckets are lower bounds.** At `maxTotal 4`,
a hit that is 4 away excluding position 1 *plus* a position-1 mismatch totals 5
and is outside the gate. So 0/1/2/3 are counts and anything above is a floor.
The driver names those columns `gt3*_LOWER_BOUND` so they cannot be misread.
Raise `--total` to 5 to make the 4 bucket exact, at a much larger scan.

**5. Seed matches are a separate axis from total distance.** Across 3,617
designs, the correlation between nearest-off-target distance and perfect-seed
off-target count was only **r = −0.36** — directionally right, but roughly 13%
of shared variance. Designs in the "≥3 mismatches" tier ranged from 0 to 736
perfect-seed off-target genes. Ranking on distance alone admits
seed-catastrophic designs; the driver reports both and filters on both.

---

## 9. Choosing thresholds

```
--min-mm-as 3 --min-edit-as 3     antisense: nearest off-target gene
--min-mm-ss 2 --min-edit-ss 2     sense: nearest off-target gene
--max-seed 3                      perfect-seed off-target genes, per strand
--min-ontarget-tx 1               on-target isoform coverage
```

**The strands deserve different thresholds.** A well-designed duplex loads the
antisense strand preferentially via thermodynamic asymmetry, so a sense-strand
off-target is discounted by the loading ratio. Scoring both at `≥3` treats them
as equally dangerous when they are not — hence the asymmetric defaults.

**`--max-seed` is the most aggressive knob.** Perfect-seed counts have a median
around 22 and a long tail, so `≤3` sits near the 10th percentile; applied to
both strands independently it dominates every other criterion. On a 30-design
test it passed 2/30 and 6/30 while the distance criteria passed 3/30 and 25/30.
Loosen it first if nothing survives.

Also note that the perfect-seed count includes seed matches **anywhere in the
transcript**, while seed-mediated silencing is a 3′UTR phenomenon — so it
overstates canonical risk, and demanding zero is stricter than the biology
warrants.

**Expect the conjunction to be brutal, and treat that as information.** Seven
criteria across two strands can leave very few survivors from thousands of
candidates. If a complete tiling yields nothing clean on both strands, that is a
finding about the target — and an argument for solving passenger loading
chemically, through modifications that block it, rather than by sequence
selection alone.

Because `--resume` re-filters without rescanning, explore thresholds freely.

---

## 10. Output columns

Raw scan (`scan.tsv`):

| column | meaning |
|---|---|
| `oligo` | your query id |
| `strand` | `+` or `−`, reference strand of the hit |
| `seq_id` | first token of the FASTA header |
| `pos` | **0-based**, within-sequence, end-exclusive |
| `seed_mm` | mismatches in positions 2–8 |
| `nonseed_mm` | mismatches in position 9 onward |
| `total_mm` | all mismatches, position 1 included |
| `gap_events` | distinct gap events |
| `oligo_bulge` / `target_bulge` | bulge bases on each side |
| `class` | `mismatch_only`, `oligo_bulge`, `target_bulge` |
| `edit_total` | `total_mm + gap_events` |

Position-1 mismatches are recovered as `total_mm − seed_mm − nonseed_mm`, which
is 0 or 1. `class == mismatch_only` is equivalent to `gap_events == 0`, so a
mismatch-only analysis is a subset of an indel-permitting scan — you do not need
to scan twice.

`duplex_screen_results.csv` carries, per design, both strands × both analyses:
`AS_*` / `SS_*` prefixes, `minMM` and `minEDIT`, exact 0–3 buckets, lower-bound
`gt3` columns, `SeedPerfect_OffGenes`, on-target transcript count,
`SS_target_gene_bestMM` (the orientation check), and a `verdict` naming every
criterion a design failed.

---

## 11. Troubleshooting

| symptom | cause |
|---|---|
| `Critical error! The length of target sequences should match...` | Cas-OFFinder needs queries padded to the pattern length — append `NNN` for a 3-nt PAM |
| GPU build: `No OpenCL devices found` | no usable OpenCL device; `clinfo -l` to check. A CPU-only OpenCL runtime is not a GPU — do not report it as one |
| GPU build cannot read the reference | the CUDA build reads only packed `.osr`; run `pack` first |
| Aggregation covers fewer queries than expected | truncated scan. The driver aborts on this rather than undercounting. Rerun without `--resume` |
| Everything fails `--max-seed` | expected; see §9. Loosen it first |
| Output is enormous | the gate is loose. Lower `maxTotal`, or drop indels with `maxIndelBases 0` |
| Hits on the wrong strand | mode confusion. Re-read §5; use `anti` for anything RISC-loaded |
| Sense strand hits the target gene | not necessarily a bug — see §6. Run the precheck |

---

## 12. Reproducing the paper's case study

A complete 19-nt tiling of PCSK9 (`ENST00000302118.5`), both strands:

```bash
python3 scripts/sirna_duplex_screen.py \
    --sirna-csv       data/case_study/PCSK9_ENST00000302118.5_siRNAs.csv \
    --ref             data/ref/Homo_sapiens.GRCh38.cdna.all.osr \
    --precheck-fasta  data/ref/Homo_sapiens.GRCh38.cdna.all.fa \
    --tx2gene         data/ref/cdna_tx2gene.tsv \
    --target-gene     PCSK9 \
    --oligoscan       build/oligoscan \
    --outdir          data/case_study/screen
```

Reference: Ensembl GRCh38 cDNA. Gate: `seed ≤3, body ≤3, total ≤4`, 1 indel
base / 1 indel event. Position-1 mismatches discounted.

---

## Citing and contact

See `README.md` for citation. Issues and questions: `[CONTACT — OPEN ITEM]`.
