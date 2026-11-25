# oligoscan — mismatches, indels, per-locus selection

Exhaustive off-target scanner for short oligos (siRNA / ASO) against a
transcriptome FASTA. No mismatch cap, no seeding, no heuristic: every position
of the reference is examined, so recall is 100% by construction.

## Build

    gcc -O3 -march=native -fopenmp -o oligoscan src/oligoscan.c

### Or in Docker (oligoscan + bowtie1 + bowtie2)

    docker compose build        # compiles, runs --selftest and the selection suite
    docker compose run --rm oligoscan build_indexes.sh /data/ref/<cdna>.fa
    docker compose run --rm oligoscan            # shell in /data

Put the reference FASTA in `data/ref/` on the host; it appears as `/data/ref/`
in the container. Indexes land in `data/index/bowtie1/` and
`data/index/bowtie2/`, outside the image, so they survive rebuilds. oligoscan
itself needs no index -- bowtie1's is the validation oracle, bowtie2's is for
comparison with the current platform.

## Run

    ./oligoscan ref.fa oligos.tsv <maxSeedMM> <maxNonseedMM> <maxTotalEdits> \
                [anti|sense|both] [maxIndelBases=0] [maxIndelEvents=1] [--all-descriptions]
    ./oligoscan --selftest
    ./oligoscan pack ref.fa ref.osr      # one-time; then pass ref.osr instead of ref.fa

`ref.osr` is the reference packed 2 bits per base (about a quarter of the
FASTA) plus the sequence table and the positions of padding and non-ACGT
bases. It loads in about a second instead of re-parsing text, and gives
byte-identical output. It is tied to the build that wrote it (a version byte
is checked), so re-run `pack` after upgrading oligoscan if told to.
Output is sorted (oligo, sequence, position) in every mode, so it is
identical across runs and thread counts.

With `maxIndelBases = 0` the byte-compare mismatch path runs and every passing
offset is reported (this is the path that matches bowtie1 `-v`). With
`maxIndelBases > 0` the Myers indel path runs and reports **one row per locus**:
the best description, ranked by edit total, then gap bases, gap events, seed /
pos-1 / body mismatches, then position. Equal-cost ties go to the ungapped
description. Two descriptions are the same locus when their start or end is
within R = max(1, maxIndelBases) bases. Output in this mode is sorted and
byte-identical across thread counts. `--all-descriptions` turns selection off
(debugging only).

`oligos.tsv` is `name<TAB>sequence` per line; U and lower case are accepted,
and lengths may be mixed between 10 and 64 nt. Out-of-range oligos and any
carrying a degenerate base are skipped with a reason on stderr.
Output is TSV on stdout: oligo, strand, seq_id, pos (0-based within the
sequence), seed_mm, nonseed_mm, total_mm, gap_events, oligo_bulge,
target_bulge, class, edit_total. Progress goes to stderr.

Mismatches are counted in three disjoint regions of the antisense oligo:
position 1 (5' terminal base), positions 2-8 (seed), positions 9..L (body).

## Validate

    ./oligoscan --selftest
    python3 tests/test_oligoscan_select.py ./oligoscan      # no external data needed

The selection test checks the mismatch path against a numpy brute force, and
the indel path against a gap-constrained brute-force DP: one row per locus, at
the optimum cost, nothing within budget missed, junctions respected, output
identical at 1 and 2 threads.

    python3 tests/validate_oligoscan.py ref.fa oligos.tsv ./oligoscan [bowtie_index]

Asserts byte-identical output against `bowtie -v 3` at the matching gate.
bowtie1's `-v` mode is exhaustive backtracking within its 3-MM cap, so it is a
true oracle there — and the only one available, since nothing else is
exhaustive above 3 MM.

## Measured (80.5 Mb reference with repeat and paralog structure, 40 real
## 21-mers, 2 cores)

    gate                              sites   per oligo   scan
    seed<=0 nonseed<=3 total<=3           2         0.1   2.17s
    seed<=1 nonseed<=3 total<=4         175         4.4   2.2s
    seed<=0 nonseed<=6 total<=6       2,216        55.4   2.2s
    seed<=1 nonseed<=5 total<=6       9,819       245.5   2.26s
    unconstrained seed, total<=6     34,583       864.6   2.40s

Runtime is independent of k and of the gate — 2.17s at 3 MM, 2.40s at 6 MM.
The small growth is output writing, not search.

## Known gaps
- Reference held as one byte per base; 2-bit packing would cut it to a quarter.
- No index, so there is no sub-second single-oligo lookup. A 560 Mb scan for
  one oligo is ~0.2s on 20 cores, which is interactive enough.
