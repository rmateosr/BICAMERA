# ABOUTME: Integration test — run BIC_algorithm_ASM on a synthetic dataset and
# ABOUTME: check detection rate per region class (NULL low, STATIC_ASM high).
import os
import subprocess
import sys
import tempfile

import pandas as pd
import pytest

THIS = os.path.abspath(os.path.dirname(__file__))
EM_DIR = os.path.abspath(os.path.join(THIS, ".."))
SYNTHETIC_SCRIPT = os.path.join(EM_DIR, "synthetic", "generate_synthetic.py")
BIC_DRIVER = os.path.join(EM_DIR, "BIC_algorithm_ASM.py")

# Detection-rate thresholds are loose because compute_M1_BIC reseeds via
# np.random.seed() inside the EM. With seed=42 and reads_per_region=120,
# detection on STATIC_ASM should hit 80%+ and false-positive on NULL stays
# below 25%.
DETECT_HI = 0.80
DETECT_LO = 0.25


@pytest.fixture(scope="module")
def synthetic_run():
    with tempfile.TemporaryDirectory() as workdir:
        data_dir = os.path.join(workdir, "synth")
        out_tsv = os.path.join(workdir, "asm.tsv")

        subprocess.check_call([
            sys.executable, SYNTHETIC_SCRIPT,
            "--n-regions", "40",
            "--reads-per-region", "120",
            "--cpgs-per-region", "30",
            "--seed", "42",
            "--output-dir", data_dir,
            "--class-mix", "1,1",
        ])

        subprocess.check_call([
            sys.executable, BIC_DRIVER,
            data_dir,
            "--output", out_tsv,
            "--windowsize", "10",
            "--workers", "2",
            "--asm-restarts", "1",
        ])

        manifest = pd.read_csv(os.path.join(data_dir, "manifest.tsv"), sep="\t")
        result = pd.read_csv(out_tsv, sep="\t")
        yield manifest, result


def _merge(manifest, result):
    return manifest.merge(
        result[["coord_evaluated", "BICresult"]],
        on="coord_evaluated",
        how="inner",
        suffixes=("_truth", "_pred"),
    )


def test_shared_region_low_detection(synthetic_run):
    manifest, result = synthetic_run
    merged = _merge(manifest, result)
    shared = merged[merged["true_class"] == "SHARED"]
    assert len(shared) > 0, "no SHARED regions in manifest"
    rate = shared["BICresult"].mean()
    assert rate < DETECT_LO, f"SHARED false-positive rate too high: {rate:.3f}"


def test_static_asm_high_detection(synthetic_run):
    manifest, result = synthetic_run
    merged = _merge(manifest, result)
    static = merged[merged["true_class"] == "STATIC_ASM"]
    assert len(static) > 0, "no STATIC_ASM regions in manifest"
    rate = static["BICresult"].mean()
    assert rate > DETECT_HI, f"STATIC_ASM detection rate too low: {rate:.3f}"


def test_BIC_winner_matches_BICresult(synthetic_run):
    """BIC_winner column should be 0 iff BICsinglecomp < BICmiddlecomp."""
    _, result = synthetic_run
    expected = (result["BICmiddlecomp"] < result["BICsinglecomp"]).astype(int)
    assert (result["BIC_winner"] == expected).all()
    assert (result["BICresult"] == expected).all()


# Column ORDER, not just membership. Several analysis harnesses address this TSV by
# column number — BICresult at 5 and window_cpgs at 11, plus raw $3/$4/$7/$8/$11 in
# two awk scripts, three of which have no header guard and would silently read the
# wrong field. Any new column therefore goes at the END. This test is what stops a
# well-meaning tidy-up from reordering the schema.
EXPECTED_COLUMN_ORDER = [
    "coord_evaluated",    # 1
    "window_last_coord",  # 2
    "BICsinglecomp",      # 3
    "BICmiddlecomp",      # 4
    "BICresult",          # 5
    "BIC_winner",         # 6
    "windows_size",       # 7
    "total_reads",        # 8
    "em_iterations",      # 9
    "mean_theta_diff",    # 10
    "window_cpgs",        # 11
    "deltaBIC",           # 12
]


def test_column_order_is_stable(synthetic_run):
    _, result = synthetic_run
    assert list(result.columns) == EXPECTED_COLUMN_ORDER


def test_empty_and_populated_schemas_agree():
    """_BIC_SUMMARY_COLUMNS and the pd.DataFrame dict are two independent
    declarations of the same order, and nothing else keeps them in sync.

    A region with no scoreable windows returns _empty_bic_summary(), built from the
    tuple; every other region is built from the dict. The driver concatenates them
    and writes the header from whichever frame comes first, so if the two ever
    diverge the output header order becomes dependent on region ordering. Anchoring
    both to EXPECTED_COLUMN_ORDER is what stops that.
    """
    import EMfunctions_ASM as em

    assert list(em._BIC_SUMMARY_COLUMNS) == EXPECTED_COLUMN_ORDER
    assert list(em._empty_bic_summary().columns) == EXPECTED_COLUMN_ORDER


def test_deltaBIC_uses_the_manuscript_sign(synthetic_run):
    """deltaBIC = BIC_Shared - BIC_ASM, so positive favours ASM."""
    _, result = synthetic_run
    expected = result["BICsinglecomp"] - result["BICmiddlecomp"]
    # Not exact equality: this fixture reads the TSV back from disk, so deltaBIC
    # round-trips through decimal text while `expected` is recomputed from the
    # re-parsed operands. That disagrees in the last bits (~1e-13 on values ~300,
    # i.e. float64 epsilon), which says nothing about the arithmetic being right.
    assert (result["deltaBIC"] - expected).abs().max() < 1e-9
    # Equivalent to the boolean call, on every row.
    assert ((result["deltaBIC"] > 0).astype(int) == result["BICresult"]).all()
    # A real run must contain both signs, or the assertion above proves nothing.
    assert (result["deltaBIC"] > 0).any() and (result["deltaBIC"] < 0).any()
