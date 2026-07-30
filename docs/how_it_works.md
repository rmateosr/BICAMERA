# How it works

## The two models

For a window of J CpGs, BICAMERA compares:

| | Column in output | Parameters |
|---|---|---|
| **Shared** | `BICsinglecomp` | J — one methylation probability per CpG, common to both alleles |
| **ASM** | `BICmiddlecomp` | 2J — a separate probability per CpG per allele |

The ASM model is diploid-aware: the per-allele mixture weight is pinned at 0.5
rather than fitted, because a diploid genome has two alleles in equal proportion.
That constraint is what stops the extra J parameters from being spent explaining
coverage imbalance instead of methylation difference. It lives in `compute_M1_BIC`
in `em/EMfunctions_ASM.py`.

Both models are fitted by EM. Lower BIC wins, so `BICresult = 1` means the ASM
model paid for its extra parameters.

Haplotype tags are used semi-supervised. A read tagged `HP:i:1` or `HP:i:2` is
**conditioned on** that tag — its allele is treated as observed, and the mixture
prior drops out of its likelihood. An untagged read keeps its allele latent, and the
EM marginalises over both possibilities using the 0.5 prior. So phasing information
is taken as fact where it exists, and inferred where it doesn't, which is why a
partially-phased sample still yields calls. See `_observed_loglikelihood` in
`em/EMfunctions_ASM.py`.

## Candidate regions, not the whole genome

BICAMERA doesn't run EM everywhere. Doing so would be slow and would test millions
of windows that cannot show ASM. Instead it pre-selects.

```
pileup      modkit pileup, one row per CpG      population methylation
   ↓
regions     keep CpGs with coverage ≥ MIN_COV and
            FRAC_LOW ≤ methylation ≤ FRAC_HIGH,
            merge survivors within MERGE_DIST bp     candidate regions
   ↓
hpmap       read_id → H1/H2/noH from the HP tags
   ↓
modkitref   read-level extract, restricted to those regions
   ↓
filter      keep valid m/- calls
   ↓
coordfix    join haplotype labels, fix the −strand off-by-one
   ↓
m2r         one row per read
   ↓
bucket      split reads into per-region files
   ↓
ASM EM      windows of WINDOWSIZE CpGs → BIC per window
```

**Consequence worth understanding:** the fraction of windows called ASM is a rate
over pre-selected candidates, not a genome-wide rate. A high proportion is expected
and is not evidence of a problem.

## Where 0.325–0.675 comes from

The `FRAC` band selects CpGs whose *pooled* methylation — both alleles together,
phase ignored — is intermediate. A locus that is genuinely allele-specific has one
allele methylated and one not, so pooled it looks half-methylated. A locus that is
uniformly methylated or uniformly unmethylated cannot be allele-specific, so
scoring it wastes time.

The bounds come from the per-fragment cutoffs used in the DNA methylation atlas
(Loyfer et al., *Nat Commun* 2025): hypomethylated ≤ 0.35, hypermethylated ≥ 0.65.
For a balanced mix of one hypo and one hyper allele, the pooled mean falls between:

```
lower = 0.5 × 0.65        = 0.325
upper = 0.5 × (1 + 0.35)  = 0.675
```

Widening the band scores more windows and costs runtime; narrowing it risks missing
partially-methylated ASM. Both bounds are environment variables.

## Reproducibility

The ASM model's EM starts from a seed derived from a hash of the window's data
(`_seed_from_X`), so on a given machine the same input gives the same BIC every time —
no run-to-run drift from random restarts, and no dependence on `WORKERS`.

**Across different hardware, expect agreement to about 1e-10, not to the bit.** The
shared model (`BICsinglecomp`) is order-invariant and reproduces exactly. The ASM
model goes through a BLAS matrix product whose last bits depend on the CPU and the
kernel the library dispatches to, and the EM iteration amplifies that slightly.
Measured on chr22 of a 48× HiFi genome, same input, two different Intel CPUs:

| | |
|---|---|
| `BICsinglecomp` | identical on all 94,423 windows |
| `BICmiddlecomp` | differs on 11,621, max 3.9e-10 |
| `em_iterations` | differs on 4 |
| `BICresult` | **0 flips** |

So calls are stable and results are unaffected; only the low-order digits of the ASM
BIC move. If you need bit-identical output for a regression test, pin the machine.

This is also why `modkit extract calls` runs single-threaded. Threading it doesn't
corrupt anything, but it reorders reads between threads, the read order reaches the
EM's input matrix, and a different matrix order changes the seed and therefore the
EM trajectory. Values stay valid, but they stop matching a previous run. The pileup
stage *is* threaded, because per-site counts don't depend on read order.

One caveat: the driver collects results with `imap_unordered`, so the *row order*
within an output TSV can vary between runs even though every value is
deterministic. Sort before diffing two runs.

## Layout

```
run_pipeline.sh          submits or runs the DAG
config/defaults.sh       every tunable default
config/site.sh           tool resolution (+ optional site profiles)
scheduler/sge.sh         one backend per file
scheduler/local.sh
scripts/                 the nine stages above
em/                      EM + BIC engine, tests, synthetic data generator
examples/                worked launcher
```

`em/synthetic/generate_synthetic.py` produces truth-labelled regions (shared and
ASM), which is how the test suite checks the models recover what they should.
