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
