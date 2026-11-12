# Benchmarks

The scripts that generated the paper's simulated data and scored it with BICAMERA,
a Fisher exact test and ASMS. Every simulated region is drawn from its own fixed
seed, so each table below regenerates row for row.

## What is simulated

A region is one 10-CpG BICAMERA window. Depth is the number of reads per region;
each read covers 3–10 consecutive CpGs. A read is tagged H1 or H2 with probability
(1 − *u*)/2 each and left untagged with probability *u*, the untagged fraction. H1
and H2 reads come from their own allele, untagged reads from either at 50/50. The
generator is `em/synthetic/generate_synthetic.py`.

| Truth | Methylation per CpG |
|---|---|
| Shared (negative) | one value, uniform 0–1, the same on both alleles |
| ASM (positive) | allele 1 uniform 0.65–1.00, allele 2 uniform 0.00–0.35 |
| Partial ASM | ASM at *n* of the 10 CpGs (*n* = 0–10), Shared at the rest |

The same reads are scored under four tagging conditions:

| Condition | Key in the tables | Reads scored |
|---|---|---|
| Standard | `as_is` | tags as generated |
| Hap Only | `drop_noh` | untagged reads removed |
| Untagged | `all_noh` | every tag removed |
| Tag Randomized | `randomized` | tags shuffled among the reads |

by three detectors:

| Detector | Called ASM when | Script |
|---|---|---|
| BICAMERA | `BICresult = 1` | `run_sweep.py`, using `em/` with 5 EM restarts |
| Fisher exact test | `fisher_p < 0.001` (the paper's α); an empty p is no call | `fisher/run_fisher_sweep.py` |
| ASMS | `asms_call = 1` | `asms/run_asms_sweep.py` |

Fisher pools a region's calls into one 2 × 2 table (H1/H2 × methylated/unmethylated),
so it cannot use untagged reads and makes no call under Untagged. ASMS ignores tags;
only Hap Only changes what it sees.

Recall is the fraction of ASM regions called ASM in a (condition, depth, untagged
fraction) cell. Precision is TP / (TP + FP), with the false positives taken from the
Shared regions of the same cell.

## Paper figures

| Paper | Sweep | Output table |
|---|---|---|
| Fig. 2B–C; Supp. Figs S1–S6 | `bic`, `fisher`, `asms` | `bic_paper.tsv`, `fisher_paper.tsv`, `asms_paper.tsv` |
| Fig. 3; Supp. Figs S7–S8 | `ncpg_bic`, `ncpg_fisher`, `ncpg_asms` | `ncpg_bic.tsv`, `ncpg_fisher.tsv`, `ncpg_asms.tsv` |
| Methods: five random EM initializations | `restarts` | `restarts.tsv` |

The paper grid is depth 10–100 in steps of 10 × untagged fraction 0.0–1.0 in steps
of 0.1, with 100 regions per truth class per cell. The partial-ASM sweep uses the
same grid for every *n*, with its own seed base, so its regions are not the full-ASM
regions.

## Run

```bash
pip install -r requirements.txt          # numpy, pandas, scipy

bash benchmarks/run_benchmarks.sh bic    /scratch/bench/bic
bash benchmarks/run_benchmarks.sh fisher /scratch/bench/fisher
bash benchmarks/run_benchmarks.sh asms   /scratch/bench/asms
bash benchmarks/run_benchmarks.sh ncpg_bic    /scratch/bench/ncpg_bic
bash benchmarks/run_benchmarks.sh ncpg_fisher /scratch/bench/ncpg_fisher
bash benchmarks/run_benchmarks.sh ncpg_asms   /scratch/bench/ncpg_asms
bash benchmarks/run_benchmarks.sh restarts    /scratch/bench/restarts
```

Each call queues one task per depth (or per *n*, or per restart-sweep depth) and a
merge held on them. The merged table lands in `<OUT_DIR>/<label>.tsv`, shards in
`<OUT_DIR>/shards/`, logs in `<OUT_DIR>/log/`. As with the pipeline, the SGE backend
is used when `qsub` exists; `BICAMERA_SCHEDULER=local` runs the tasks one after
another instead.

Anything after `<OUT_DIR>` goes to the Python script of every task. For a quick
check, `--reps 10` scores the first 10 regions per cell. A region's seed depends only
on its own class (or *n*), depth, untagged fraction and replicate number, so those
rows are identical to the first 10 replicates of a full run.

### The fine grid

The paper's tables were cut from a finer sweep: depth 1–100 × untagged fraction
0.00–1.00 in steps of 0.01. The default `GRID=paper` draws the same regions directly,
because it seeds each untagged fraction with its index on the fine grid. To run the
fine sweep as the paper did:

```bash
GRID=fine bash benchmarks/run_benchmarks.sh bic /scratch/bench/bic_fine
bash benchmarks/subset_paper_grid.sh /scratch/bench/bic_fine/bic_fine.tsv bic_paper.tsv
```

The subset has the same rows as `GRID=paper`, in a different order. `GRID=smoke` is a
five-depth grid for testing; its seeds don't match the other two.

## Output

| Table | Columns |
|---|---|
| `bic_*` | `truth_class depth noh_frac rep approach total_reads BICsinglecomp BICmiddlecomp BICresult mean_theta_diff` |
| `fisher_*` | `truth_class depth noh_frac rep approach fisher_n_reads fisher_p` |
| `asms_*` | `truth_class depth noh_frac rep approach asms_n_reads asms_call asms_absmdiff asms_pval` |
| `ncpg_*` | as above, with `n_asm` in place of `truth_class` |
| `restarts` | `truth_class depth noh_frac rep approach total_reads restart_idx BIC_M0 BIC_M1` |

`approach` holds the condition keys above. In the partial-ASM tables `n_asm = 0` is
the Shared baseline. Regions with fewer than 2 reads get BICAMERA's low-coverage
sentinel (`BICsinglecomp = 0`, `BICmiddlecomp = 1`, Shared). In `restarts`, a region
is ASM at restart count *R* when the smallest of its first *R* `BIC_M1` values is
below `BIC_M0`; its untagged fractions run 0.00–1.00 in steps of 0.05 at depths 30,
50 and 100.

## Compute

Measured with one core per task:

| Sweep | Run | Tasks | Longest task | All tasks | Peak memory per task |
|---|---|---|---|---|---|
| `bic` | paper grid | 10 | 4.5 min | 32 min | 0.9 GB |
| `fisher` | paper grid | 10 | 1.1 min | 8 min | 0.9 GB |
| `asms` | paper grid | 10 | 3.0 min | 19 min | 0.9 GB |
| `ncpg_bic` | `--reps 10` | 11 | 2.8 min | 20 min | 0.9 GB |
| `ncpg_fisher` | `--reps 10` | 11 | 0.8 min | 8 min | 0.9 GB |
| `ncpg_asms` | `--reps 10` | 11 | 1.5 min | 13 min | 0.9 GB |
| `restarts` | defaults | 3 | 3.7 min | 7 min | 0.9 GB |

Run time grows about linearly with the number of regions, so a full partial-ASM run
(100 replicates) takes roughly ten times the `--reps 10` figures. The fine grid has
92 times as many cells as the paper grid.

## Reproducibility

These scripts were checked against the tables behind the paper's figures by
re-running them and comparing every field as text. All rows were identical:

| Sweep | Rows compared |
|---|---|
| `bic`, `fisher`, `asms` (paper grid) | 88,000 each, the whole table |
| `ncpg_bic`, `ncpg_fisher`, `ncpg_asms` | 48,400 each, the first 10 replicates |
| `restarts` | 126,000, the whole table |
| `fisher` fine grid | 242,400, the first 3 replicates |

The cut-down fine-grid table from `subset_paper_grid.sh` also matched the paper-grid
run row for row. The check ran with Python 3.12.0, numpy 1.26.2, pandas 2.2.2 and
scipy 1.11.4.

On one machine the output is the same on every run. On a different CPU, the BIC
values of the ASM model can differ in the last bits, without changing calls; see
[docs/how_it_works.md](../docs/how_it_works.md#reproducibility).

## ASMS

`asms/asms_scorer.py` is a Python port of the `cluster` and `filter` steps of ASMS,
so that ASMS scores the same simulated reads as the other two detectors. It covers
the EM clustering, the filter thresholds and the 5,000-shuffle permutation test; the
`scan` step and BAM reading are not ported. Its shuffles come from numpy, so its
p-values match the Rust tool only up to Monte-Carlo error.

`asms/validate_asms_scorer.py` checks the port against ASMS's own unit-test values
and, given a built `asms` binary, against `asms pval` on simulated regions:

```bash
python benchmarks/asms/validate_asms_scorer.py /path/to/asms   # or set ASMS_BIN
```

The port and the copied test values are adapted from ASMS by E. Raineri,
M. A. Esteve Marco, E. Casals and A. Esteve Codina, "asms: finding allele-specific
methylation in human genomes without phasing", *NAR Genomics and Bioinformatics*
8(3):lqag072 (2026), [doi:10.1093/nargab/lqag072](https://doi.org/10.1093/nargab/lqag072).
Source: the ecmra/asms Zenodo release,
[doi:10.5281/zenodo.20323001](https://doi.org/10.5281/zenodo.20323001)
(commit `c7178db`), used under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The changes are listed above.
