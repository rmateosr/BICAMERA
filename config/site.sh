# ABOUTME: Resolves BICAMERA's external tools — PATH first, optional site profile.
# ABOUTME: Sourced by run_pipeline.sh and every inner job; safe to source twice.
#
# BICAMERA needs four external tools: modkit (>= 0.4.1), samtools, bedtools and
# python3. Each is resolved from PATH and can be overridden by exporting
# MODKIT_BIN / SAMTOOLS / BEDTOOLS / PYTHON directly.
#
# On clusters where tools come from environment modules rather than the login
# PATH, put the `module load` calls in a site profile — config/site.<name>.sh —
# and select it with BICAMERA_SITE=<name>. config/site.shirokane.sh is a worked
# example. A profile must never fail its caller: a missing module system is not
# an error anywhere except on that one site.

: "${BICAMERA_ROOT:?BICAMERA_ROOT not set — source this from run_pipeline.sh or an inner job}"

# Auto-detect a profile only when BICAMERA_SITE is genuinely UNSET, so that on a
# site with a module tree the pipeline works with no configuration. Set
# BICAMERA_SITE=none to opt out explicitly — needed both for testing the
# PATH-only path and for any cluster that happens to have this directory but is
# not this site. The `-d` test is only ever a probe: if it fails, tool resolution
# falls through to PATH. This is the one site-specific path in the portable code.
if [ -z "${BICAMERA_SITE+set}" ] && [ -d /usr/local/package/modulefiles ]; then
    BICAMERA_SITE=shirokane
fi
export BICAMERA_SITE="${BICAMERA_SITE:-none}"

if [ "${BICAMERA_SITE}" != "none" ]; then
    _bicamera_profile="${BICAMERA_ROOT}/config/site.${BICAMERA_SITE}.sh"
    if [ -f "${_bicamera_profile}" ]; then
        . "${_bicamera_profile}"
    else
        echo "WARN: BICAMERA_SITE=${BICAMERA_SITE} but ${_bicamera_profile} not found; resolving tools from PATH only." >&2
    fi
    unset _bicamera_profile
fi

# modkit is not bundled (Oxford Nanopore Public License v1.0 — we do not
# redistribute it). Fall back to an optional in-repo bin/modkit for users who
# drop the release binary there.
MODKIT_BIN="${MODKIT_BIN:-$(command -v modkit || true)}"
[ -n "${MODKIT_BIN}" ] || MODKIT_BIN="${BICAMERA_ROOT}/bin/modkit"
SAMTOOLS="${SAMTOOLS:-$(command -v samtools || true)}"
BEDTOOLS="${BEDTOOLS:-$(command -v bedtools || true)}"
PYTHON="${PYTHON:-$(command -v python3 || true)}"
export MODKIT_BIN SAMTOOLS BEDTOOLS PYTHON

# Minimum modkit version. 0.4.1 moved read-level extraction under the `extract
# calls` subcommand, and the pipeline reads that subcommand's TSV by fixed column
# position — $3 ref_position and $12 call_code in the filter and coordfix stages,
# $10/$11 of the pileup bedMethyl in region selection. An older modkit has no
# `calls` subcommand at all; a layout change would shift those columns and yield
# silently wrong methylation calls rather than an error, which is why the version
# is enforced here instead of merely noted in a comment.
BICAMERA_MODKIT_MIN="0.4.1"
export BICAMERA_MODKIT_MIN

# Verify $MODKIT_BIN is at least $BICAMERA_MODKIT_MIN. Returns non-zero on a
# too-old version. `modkit --version` prints e.g. "mod_kit 0.4.1"; an unparseable
# version warns rather than blocks, so an unexpected build string cannot make the
# pipeline unrunnable.
bicamera_check_modkit_version() {
    _mk_ver=$("${MODKIT_BIN}" --version 2>/dev/null | awk 'NR==1{print $NF}')
    # Only compare something that actually looks like a dotted version. Without
    # this guard a build string such as "some weird build" reaches sort -V, whose
    # ordering against a numeric version is arbitrary — it could just as easily
    # report "modkit build is too old" as accept it.
    case "${_mk_ver}" in
        *[!0-9.]*|''|.*|*.)
            echo "WARN: could not read a version from '${MODKIT_BIN} --version'" >&2
            echo "      (got '${_mk_ver:-<empty>}'); skipping the modkit >= ${BICAMERA_MODKIT_MIN} check." >&2
            echo "      Column positions are assumed to match the 0.4.x layout." >&2
            unset _mk_ver
            return 0 ;;
    esac
    # sort -V puts the lower version first; if that is still the minimum, we are at
    # or above it (equal versions also leave the minimum first).
    if [ "$(printf '%s\n%s\n' "${BICAMERA_MODKIT_MIN}" "${_mk_ver}" | sort -V | head -1)" != "${BICAMERA_MODKIT_MIN}" ]; then
        echo "ERROR: modkit ${_mk_ver} is too old (need >= ${BICAMERA_MODKIT_MIN})." >&2
        echo "       binary: ${MODKIT_BIN}" >&2
        echo "       Read-level extraction moved under 'extract calls' in ${BICAMERA_MODKIT_MIN}; older" >&2
        echo "       builds lack the subcommand. Get a release binary from" >&2
        echo "       https://github.com/nanoporetech/modkit/releases" >&2
        unset _mk_ver
        return 1
    fi
    unset _mk_ver
    return 0
}

# bicamera_require <tool> [<tool> ...]
# Abort with an actionable message if a needed tool is unresolved. Each stage
# declares only the tools it actually runs, so a missing bedtools cannot fail a
# stage that never calls bedtools.
bicamera_require() {
    _req_missing=0
    for _req_tool in "$@"; do
        case "${_req_tool}" in
            modkit)
                if [ ! -x "${MODKIT_BIN}" ]; then
                    echo "ERROR: modkit not found (searched PATH, then ${MODKIT_BIN})." >&2
                    echo "       BICAMERA requires modkit >= ${BICAMERA_MODKIT_MIN}; get a release binary from" >&2
                    echo "       https://github.com/nanoporetech/modkit/releases" >&2
                    echo "       then re-run with MODKIT_BIN=/path/to/modkit" >&2
                    _req_missing=1
                elif ! bicamera_check_modkit_version; then
                    _req_missing=1
                fi ;;
            samtools)
                [ -n "${SAMTOOLS}" ] || { echo "ERROR: samtools not on PATH; re-run with SAMTOOLS=/path/to/samtools" >&2; _req_missing=1; } ;;
            bedtools)
                [ -n "${BEDTOOLS}" ] || { echo "ERROR: bedtools not on PATH; re-run with BEDTOOLS=/path/to/bedtools" >&2; _req_missing=1; } ;;
            python)
                [ -n "${PYTHON}" ] || { echo "ERROR: python3 not on PATH; re-run with PYTHON=/path/to/python3" >&2; _req_missing=1; } ;;
            *)
                echo "ERROR: bicamera_require: unknown tool '${_req_tool}'" >&2; _req_missing=1 ;;
        esac
    done
    unset _req_tool
    if [ "${_req_missing}" -ne 0 ]; then
        unset _req_missing
        exit 1
    fi
    unset _req_missing
}
