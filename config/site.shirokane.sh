# ABOUTME: Site profile for SHIROKANE (SGE) — puts module-provided tools on PATH.
# ABOUTME: Selected with BICAMERA_SITE=shirokane; must never fail its caller.
#
# Sourced by config/site.sh before it resolves tools from PATH. Every command is
# guarded with `|| true` because this file is also sourced inside SGE array tasks,
# which run under `set -o errexit`: a module that has been renamed or retired must
# degrade to "tool not on PATH" (caught later by bicamera_require, with a useful
# message) rather than killing the job with no explanation.
#
# `module` is a shell function, so it reaches an SGE job only via `qsub -V`
# (which the SGE backend always passes) or a login shell that defines it.
#
# Copy this file to config/site.<yoursite>.sh and edit the module names to adapt
# BICAMERA to another cluster.

module use /usr/local/package/modulefiles/ >/dev/null 2>&1 || true
module load samtools/1.19                 >/dev/null 2>&1 || true
module load bedtools/2.31.1               >/dev/null 2>&1 || true
module load python/3.12.0                 >/dev/null 2>&1 || true
