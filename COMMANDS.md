# oligoscan — benchmark command reference

Every command below is **a single line**. PowerShell does not accept the `\`
continuations used in bash; if a command wraps in your editor it is still one
line. Run them from `C:\WebUI\oligoscan`.

Shared paths, assumed throughout:

| what | path |
|---|---|
| cDNA FASTA | `/data/ref/Homo_sapiens.GRCh38.cdna.all.fa` |
| packed reference | `/data/ref/Homo_sapiens.GRCh38.cdna.all.osr` |
| RIsearch2 suffix array | `/data/ref/Homo_sapiens.GRCh38.cdna.all.suf` (11.9 GB, built) |
| bowtie1 index | `/data/index/bowtie1/Homo_sapiens.GRCh38.cdna.all` |
| bowtie2 index | `/data/index/bowtie2/Homo_sapiens.GRCh38.cdna.all` |
| library | `/opt/oligoscan/tests/huesken_guides.tsv` (2,361 × 19 nt) |
| outputs | `data/validation/<prefix>_<timestamp>/` |

---

## 1. Gapped bowtie2 comparison, 100 guides — **RUN THIS ONE**

Replaces the ten-guide table in §3.6 and the warned-off table in the README.
Three gates, both tools gapped at matched cost.

**Launcher:**

```
.\scripts\run_gapped100.cmd
```

**Or directly (one line):**

```
docker compose run --rm -d oligoscan python3 /opt/oligoscan/scripts/bench_three.py --ref /data/ref/Homo_sapiens.GRCh38.cdna.all.fa --scan-ref /data/ref/Homo_sapiens.GRCh38.cdna.all.osr --oligos /opt/oligoscan/tests/huesken_guides.tsv --bt1 /data/index/bowtie1/Homo_sapiens.GRCh38.cdna.all --bt2 /data/index/bowtie2/Homo_sapiens.GRCh38.cdna.all --n-oligos 100 --gates none --indel-gates 2:2:2,3:3:3,4:4:4
```

**Expect ~10 hours.** Basis: §3.6 measured bowtie2 at 1,196.6 / 1,198.0 /
1,205.6 s for *ten* guides at these same three gates — essentially flat,
because `-a` with flat gap costs dominates. Index load is negligible against
1,200 s, so the scaling is linear in guide count: ~3.3 h per gate, three
gates. This is an extrapolation from a measured per-oligo rate at the same
gates, not a guess.

`-d` detaches it, so closing the terminal will not kill it — that is what
killed the 2026-11-13 attempt. Follow it with:

```
docker ps
docker logs -f <container-name>
```

`runtime.tsv` and `A_recall.tsv` are rewritten after **every** measurement, so
a crash costs the remaining gates, never the finished ones.

---

## 2. Measurement B — bulge recall  ✅ done 2026-11-16

Corrected gap scoring (`--rdg 0,1`); replaces the voided §3.3 figures.

```
.\scripts\run_measB.cmd
```

```
docker compose run --rm oligoscan python3 /opt/oligoscan/scripts/bench_paper.py --ref /data/ref/Homo_sapiens.GRCh38.cdna.all.fa --scan-ref /data/ref/Homo_sapiens.GRCh38.cdna.all.osr --oligos /opt/oligoscan/tests/huesken_guides.tsv --bt1 /data/index/bowtie1/Homo_sapiens.GRCh38.cdna.all --bt2 /data/index/bowtie2/Homo_sapiens.GRCh38.cdna.all --only B --n-guides 10
```

~30 s. Result: oligoscan 588/588 (100%), bowtie2 `--gbar 1` 98.6%, bowtie2
default `--gbar 4` 91.3%. The old "11% collapse" was a scoring artifact.

---

## 3. RIsearch2 comparison, 100 guides  ✅ done 2026-11-16

```
.\scripts\run_risearch100.cmd
```

```
docker compose run --rm oligoscan python3 /opt/oligoscan/scripts/bench_risearch.py --ref /data/ref/Homo_sapiens.GRCh38.cdna.all.fa --scan-ref /data/ref/Homo_sapiens.GRCh38.cdna.all.osr --oligos /opt/oligoscan/tests/huesken_guides.tsv --suf /data/ref/Homo_sapiens.GRCh38.cdna.all.suf --n-oligos 100 --threads 22 --gates 0:3:3:0:0,0:4:4:2:1 --energies=-30,-25,-20,-15,-10
```

~4.5 min now the index exists (first run spent 291.6 s building it).

> **`--energies` must use `=`.** Written as `--energies -30,...` argparse
> reads the leading minus as another flag and the run dies at startup.

**Still open:** this used `-s 2:8/7` to mirror oligoscan's seed region, not
the parameters the RIsearch2 authors recommend in their own siRNA off-target
pipeline. Measure theirs before calling these figures RIsearch2's best
configuration.

---

## 4. Case study — enumeration + thermodynamic scoring  ⏳ needs a gene

Needs an image rebuild first (adds ViennaRNA and pyfaidx):

```
docker compose build oligoscan
```

Then, with a guides TSV (`name<TAB>sequence`, U or T both fine — export from
the designer app or hand-write it):

```
docker compose run --rm oligoscan python3 /opt/oligoscan/scripts/case_study.py --guides /data/<your_guides>.tsv --ref /data/ref/Homo_sapiens.GRCh38.cdna.all.fa --scan-ref /data/ref/Homo_sapiens.GRCh38.cdna.all.osr --gate 1:3:4:2:1 --bt2 /data/index/bowtie2/Homo_sapiens.GRCh38.cdna.all --on-target <ENST000000000> --threads 22 --top 25
```

Ranks every site by ViennaRNA `duplexfold` ΔG and marks which ones bowtie2
also found. The intended output is a sentence of the form *"of the N
strongest predicted off-target duplexes, M were absent from bowtie2's
output."*

Drop `--bt2` to skip the cross-check and halve the runtime.

---

## 5. GPU

Build (needs `Dockerfile.gpu`; `CUDA_ARCH=89` for the RTX 4000 Ada):

```
docker compose build oligoscan-gpu
```

Equivalence check, CPU vs GPU, whole library:

```
docker compose run --rm oligoscan-gpu python3 /opt/oligoscan/tests/test_gpu_equivalence.py --ref /data/ref/Homo_sapiens.GRCh38.cdna.all.osr --oligos /opt/oligoscan/tests/huesken_guides.tsv
```

Last measured: gapped path, full library, gate `19 19 3 anti 3 3` —
2,487,939 sites both sides, **byte-identical**, CPU 357.69 s vs GPU 25.17 s
(**14.2×**).

---

## Things that have bitten us

| symptom | cause | fix |
|---|---|---|
| `Missing expression after unary operator '--'` | bash `\` continuations pasted into PowerShell | one line, or use the `.cmd` launchers |
| `error: argument --energies: expected one argument` | argparse reads a leading `-` as a flag | `--energies=-30,-20` |
| `invalid literal for int(): 'pos'` | oligoscan writes a header row | parsers use `csv.DictReader` now |
| Dockerfile step fails right after printing RIsearch2 help | `RIsearch2 -h` exits **1** on success | build checks `grep -q 'RIsearch v2.1'` |
| A run dies silently overnight | terminal closed, or the machine slept | `-d` to detach; check power settings |
| An edited file commits with its old contents | `device_commit_files` caches by staged path | stage under a fresh name, verify on disk |
| Docker daemon wedges mid-run | a gate at the noise floor wrote a 14.2 GB TSV | `--max-gb` now deletes the file and exits |

---

*Generated 2026-11-16. Paths and timings reflect the state of
`C:\WebUI\oligoscan` on that date.*

---

## 6. Re-run after the 0.9.3 junction-edge fix — **DO THIS NEXT**

§8i of the handoff: 0.9.2 lost 71 sites at `total<=4` and 9 at `total<=3`,
every one flush against a transcript 5′ end. The scan must be redone; **the
bowtie2 side must not be.** The SAM files in
`data/validation/three_20261116_142320/` are still valid and cost ten hours,
so they are reused as they stand.

**Launcher — this is the whole job, and it is the one to run:**

```
.\scripts\rerun_after_fix.cmd
```

It rebuilds the image, checks the version, moves the 0.9.2 scans into
`pre_0.9.3\`, re-runs the three scans, recounts, and prints what the result
has to say. The steps below are the same thing by hand.

Rebuild the image first — it now refuses to build if the regression test fails:

```
docker compose build oligoscan
```

Then the three scans, into the existing run directory (one line each; ~18 s per
gate, so under a minute in total):

```
docker compose run --rm oligoscan oligoscan /data/ref/Homo_sapiens.GRCh38.cdna.all.osr /data/validation/three_20261116_142320/oligos.tsv 19 19 2 anti 2 2 > data/validation/three_20261116_142320/scan_total_le2_p_2_nt_indel.tsv
```

```
docker compose run --rm oligoscan oligoscan /data/ref/Homo_sapiens.GRCh38.cdna.all.osr /data/validation/three_20261116_142320/oligos.tsv 19 19 3 anti 3 3 > data/validation/three_20261116_142320/scan_total_le3_p_3_nt_indel.tsv
```

```
docker compose run --rm oligoscan oligoscan /data/ref/Homo_sapiens.GRCh38.cdna.all.osr /data/validation/three_20261116_142320/oligos.tsv 19 19 4 anti 4 4 > data/validation/three_20261116_142320/scan_total_le4_p_4_nt_indel.tsv
```

Keep the old scans first if you want the before/after on record:

```
powershell -c "Get-ChildItem data\validation\three_20261116_142320\scan_*.tsv | Rename-Item -NewName { $_.Name -replace '\.tsv$','.0.9.2.tsv' }"
```

Then recount. This replaces `three_way.tsv`'s own figures, which are not
self-consistent at gate 3:

```
docker compose run --rm oligoscan python3 /opt/oligoscan/scripts/recount_three_way.py /data/validation/three_20261116_142320
```

**Expect:** the `(oligo,transcript) oligoscan never reported` line to read **0
at all three gates**. Anything else is a second bug. oligoscan's site counts
rise by ~71 and ~9, so bowtie2's recall falls very slightly below the
published 87.06% / 97.18% — §3.6 and Figure 2 both need the new numbers.

### 6b. Regression test alone, no reference needed

```
docker compose run --rm oligoscan python3 /opt/oligoscan/tests/test_junction_edge.py --bin /usr/local/bin/oligoscan
```

Expect `30/30 planted sites reported`. On 0.9.2 it fails with exactly one
missing case, `BULGE_p1_off0`.

---

## 7. CPU vs GPU across oligo count and gate

The GPU evidence is currently spread over three runs at different oligo
counts and gates, so no single table shows how the speedup moves with either
axis. This measures one grid: **100 / 1,000 / 2,361 oligos × 4 gates ×
2 backends**, every cell on one machine with one reference.

**Launcher:**

```
.\scripts\run_gpu_scaling.cmd
```

**Expect 1.5–2 hours.** Rows are written as each cell finishes, so a crash
costs only the remaining cells.

### Two things this gets right that a hand-run grid would not

**There are two GPU binaries and they are not interchangeable.**
`oligoscan-gpu` is the mismatch-only kernel. `oligoscan-cuda --gpu` is the
Myers kernel for the indel path — and in `oligoscan.c` the GPU is taken only
when `use_gpu && maxI > 0`, so **at a mismatch-only gate `oligoscan-cuda
--gpu` silently runs on the CPU** and would put a CPU time in a GPU column
with nothing to warn you. The harness picks the binary from the gate.

**Three times are reported, not one** — reference load, scan, and wall clock.
At wide gates most of the wall clock is writing the output, which no kernel
touches: at ≤5 edits the kernel was 27 s of a 102.6 s run, so the same
measurement is 7.2× or 1.9× depending on which number you quote. A
wall-clock-only table understates the kernel; a kernel-only table overstates
what a user experiences. The paper must say which it means.

### Reading the result

`data\validation\gpu_scaling\gpu_scaling.tsv`. **Check the `identical`
column before reading any timing.** Every GPU row must read `yes`: outputs
are md5'd and compared whole, which is how the 4-row ambiguous-base deficit
in 37 million was caught. Outputs are deleted after hashing, because the full
library at the widest gate writes a few GB per run.

### Narrower slices

```
docker compose run --rm oligoscan-gpu python3 /opt/oligoscan/scripts/bench_gpu_scaling.py --ref /data/ref/Homo_sapiens.GRCh38.cdna.all.osr --oligos /opt/oligoscan/tests/huesken_guides.tsv --out /data/validation/gpu_scaling --sizes 100 --gates 4:4:4
```

`--sizes 0` means the whole library. `--cpu-only` skips the GPU arm.
`--keep` retains the TSVs instead of deleting them after hashing.

---

## 8. oligoscan vs bowtie2 at 1,000 guides — overnight

Extends §3.6 tenfold. Everything else about the comparison is unchanged:
same matched scoring, same gates, same exact-and-locus reporting.

**Launcher:**

```
.\scripts\run_gapped1000.cmd
```

**Budget 6.5–7 hours for the first gate.** bowtie2 is flat in the gate and
took ~2,412 s for 100 guides, so it scales with guide count, not with the
error budget.

**Gate order is deliberate: `4:4:4` runs first**, because that is the gate
carrying the argument — at ≤2 and ≤3 recall is already 99% and nothing is in
dispute. `3:3:3` follows as a bonus. `runtime.tsv` and `three_way.tsv` are
rewritten after every measurement, so stopping early costs only the
unfinished gate.

Runs **detached** (`-d`). Closing the terminal will not kill it; that is what
killed the 2026-11-13 attempt.

### In the morning

```
dir /od data\validation
```

```
docker compose run --rm oligoscan python3 /opt/oligoscan/scripts/recount_three_way.py /data/validation/<new-dir>
```

Use the recount rather than `three_way.tsv`: the harness's own counts are not
self-consistent at every gate (§8i), and the recount recomputes them from the
scans and SAMs under both matching criteria.

**Expected, if nothing is wrong:** `on an (oligo,transcript) oligoscan never
reported at all` reads **0**, and `bt2_only` under the *locus* criterion also
reads **0** — that is the containment result from §3.6 holding at ten times
the scale. Anything else is worth stopping on.

Disk: roughly 0.9 GB of oligoscan TSV and ~2.4 GB of SAM per gate.
