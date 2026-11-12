# BICAMERA

Finds allele-specific methylation (ASM) in a single sample from one phased modBAM.

For each short window of CpGs, BICAMERA fits two models — one methylation pattern
shared by both alleles, or a separate pattern per allele — and reports which one BIC
prefers. Reads carrying a haplotype tag are used as such; reads without one still
contribute, so unphased coverage isn't thrown away.

Output is one TSV per chromosome: a row per window, both BIC values, and the winner.

## Requirements

| | |
|---|---|
| [modkit](https://github.com/nanoporetech/modkit/releases) | ≥ 0.4.1 |
| samtools, bedtools | any recent version |
| Python | 3.9+ with numpy and pandas |

A scheduler is optional. If `qsub` is present BICAMERA submits per-chromosome array
jobs; if not, it runs the same stages locally, one chromosome at a time.

## Install

```bash
git clone https://github.com/rmateosr/BICAMERA.git
cd BICAMERA
pip install -r requirements.txt
```

modkit is not bundled — its licence doesn't allow redistribution. Put it on your
`PATH`, or point at it directly:

```bash
export MODKIT_BIN=/path/to/modkit
```

## Input

A modBAM with:

- **MM/ML tags** — per-read methylation calls from Dorado. Guppy-era basecalls
  won't work; the pipeline checks and refuses them.
- **HP tags** — haplotype assignments, e.g. from Clair3 + WhatsHap. Not strictly
  required, but they're the point of the tool.
- **an index** (`.bai`) beside it.

ONT data needs one prep step first to strip 5hmC. See
[docs/input_requirements.md](docs/input_requirements.md).

## Run

Start with a single chromosome to check your setup:

```bash
bash examples/run_example.sh sample.bam ref.fa /scratch/asm_test chr22
```

Then the whole genome:

```bash
bash examples/run_example.sh sample.bam ref.fa /scratch/asm_full
```

Results land in `<OUT_DIR>/ASM/<chr>_ASM.tsv`, logs in `<OUT_DIR>/log/`, and a
record of the run in `<OUT_DIR>/RUN_LOG.md`.

Put `OUT_DIR` somewhere with room — intermediates reach roughly 130 GB per
chromosome on a 48× genome. Once a run finishes, these subdirectories are
regenerable and safe to delete:

```
pileup/  read_format/  read_format_split/
modkit_splitbam_referenced_bychrandhap/
modkit_splitbam_referenced_bychrandhap_filtered/
modkit_referenced_splitbam_mergedbychr/
```

Keep `ASM/`, `hpmap/`, `potential_asm_regions/` and `log/`.

For scale: chr22 of a 48× HiFi genome takes about 8 minutes and peaks near 4 GB
with `WORKERS=4`, on one machine with no scheduler.

## Output

| Column | |
|---|---|
| `coord_evaluated` | first CpG in the window |
| `window_last_coord` | last CpG in the window |
| `BICsinglecomp` | BIC of the shared model |
| `BICmiddlecomp` | BIC of the per-allele (ASM) model |
| `BICresult` | **1 = ASM wins**, 0 = shared. Equivalently `deltaBIC > 0` |
| `BIC_winner` | same value as `BICresult`; kept for backward compatibility. Prefer `BICresult` |
| `windows_size` | CpGs in the window |
| `total_reads` | reads covering it |
| `em_iterations` | iterations the ASM fit needed |
| `mean_theta_diff` | distance between the two alleles' methylation — effect size |
| `window_cpgs` | the window's CpG coordinates, comma-separated |
| `deltaBIC` | `BICsinglecomp − BICmiddlecomp`, i.e. BIC<sub>Shared</sub> − BIC<sub>ASM</sub>. **Positive favours ASM**, larger is stronger |

Lower BIC wins, so `BICresult = 1` means `BICmiddlecomp < BICsinglecomp` — equivalently
`deltaBIC > 0`. `deltaBIC` uses the same sign convention as the paper, so you can read
it directly without deriving or negating anything.

`deltaBIC` sits last rather than beside the two BIC columns because some downstream
tooling addresses this file by column number. New columns are appended for the same
reason — read by name, and the order stays stable.

`em/asm_bic_to_bed.py` turns these TSVs into a BED for a genome browser.

## Tuning

Pass any of these as environment variables:

| | Default | |
|---|---|---|
| `CHROMS` | chr1–22, X, Y | chromosomes to process |
| `WINDOWSIZE` | 10 | CpGs per window |
| `WINDOW_STEP` | 5 | CpGs between window starts (1 = every CpG) |
| `WORKERS` | 10 | parallel EM processes |
| `MIN_COV` | 10 | minimum coverage for a candidate CpG |
| `FRAC_LOW` / `FRAC_HIGH` | 0.325 / 0.675 | methylation band a candidate CpG must fall in |
| `MERGE_DIST` | 1000 | bp within which candidate CpGs merge into a region |
| `MIN_REGION_CPGS` | 3 | drop regions with fewer CpGs than this |
| `RUN_LABEL` | `bam` | prefix for job names and logs |

Every default lives in `config/defaults.sh`.

BICAMERA only scores windows inside pre-selected candidate regions, not the whole
genome — so the fraction of windows called ASM is not a genome-wide rate. The
`FRAC` band is the main knob on what gets considered;
[docs/how_it_works.md](docs/how_it_works.md) explains where 0.325–0.675 comes from.

## Environment

| | |
|---|---|
| `BICAMERA_SCHEDULER` | `sge` or `local`. Default: `sge` if `qsub` exists. |
| `BICAMERA_SITE` | loads `config/site.<name>.sh` for clusters where tools come from modules. `none` to use `PATH` only. |
| `MODKIT_BIN`, `SAMTOOLS`, `BEDTOOLS`, `PYTHON` | override any tool path |

On a cluster whose tools come from environment modules, copy
`config/site.shirokane.sh` to `config/site.<yoursite>.sh`, edit the module names,
and run with `BICAMERA_SITE=<yoursite>`.

## Benchmarks

[benchmarks/](benchmarks/README.md) holds the scripts that generated the paper's
simulated datasets and scored them with BICAMERA, a Fisher exact test and ASMS.
They need `scipy` in addition to the requirements above.

## Tests

```bash
pip install pytest
pytest
```

## Citation

[Paper reference to be added upon publication]
