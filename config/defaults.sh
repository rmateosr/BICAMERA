# ABOUTME: Single source of truth for every BICAMERA tunable default.
# ABOUTME: Sourced by run_pipeline.sh and by every per-chromosome inner job.
#
# Every value here is `${VAR:-default}`, so an exported override always wins.
# The launcher exports these and the scheduler backend propagates them to the
# inner jobs; sourcing this file in the inners as well means a stage run by hand
# gets the same values instead of a second, independently-maintained copy.
#
# NOTE: CHROMS is deliberately NOT defaulted here. It is set once by
# run_pipeline.sh, which uses it to size the job array (-t 1-N); each inner then
# derives its own chromosome by indexing the SAME list. If an inner could fall
# back to its own CHROMS default, a lost environment would make it index a
# different list and silently process the wrong chromosome. Inners therefore
# keep a hard `: "${CHROMS:?...}"` guard and fail loudly instead.

# --- EM windowing ---
export WINDOWSIZE="${WINDOWSIZE:-10}"        # CpG sites per EM window
export WINDOW_STEP="${WINDOW_STEP:-5}"       # CpGs advanced between windows (1 = every CpG)
export WORKERS="${WORKERS:-10}"              # parallel processes in the EM driver
export ASM_RESTARTS="${ASM_RESTARTS:-}"      # empty -> driver default (5)
export EM_TOL="${EM_TOL:-}"                  # empty -> driver default (1e-4)

# --- Potential-ASM region pre-selection (pileup -> select -> merge) ---
# A CpG is a candidate iff valid coverage >= MIN_COV and
# FRAC_LOW <= pooled methylation fraction <= FRAC_HIGH; candidates within
# MERGE_DIST bp merge into one potential-ASM region.
export MIN_COV="${MIN_COV:-10}"

# FRAC_LOW/FRAC_HIGH bound the pooled (phase-ignored) methylation fraction. The
# defaults are the population-mean envelope of a balanced hypo/hyper site under
# the Atlas (Nat Commun 2025) per-fragment cutoffs (hypo <= 0.35, hyper >= 0.65):
#   0.5 * 0.65    = 0.325  (lower)
#   0.5 * (1+0.35) = 0.675 (upper)
# See docs/region_preselection.md for the full rationale.
export FRAC_LOW="${FRAC_LOW:-0.325}"
export FRAC_HIGH="${FRAC_HIGH:-0.675}"
export MERGE_DIST="${MERGE_DIST:-1000}"

# Drop potential-ASM regions with fewer than this many distinct in-region CpGs.
# Regions with >= this but < WINDOWSIZE get one whole-region EM window.
export MIN_REGION_CPGS="${MIN_REGION_CPGS:-3}"

# --- SGE backend only: per-slot memory for the ASM EM array ---
# Raise to 4G for ONT / unfiltered inputs. Ignored by the local backend.
export ASM_EM_VMEM="${ASM_EM_VMEM:-2G}"
