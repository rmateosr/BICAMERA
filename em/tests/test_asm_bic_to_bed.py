# ABOUTME: Lock the dBIC sign convention of asm_bic_to_bed.py: dBIC = BIC_Shared - BIC_ASM, >0 = ASM.
# ABOUTME: The sign silently disagreed with the manuscript once; these tests stop it drifting again.
import os
import subprocess
import sys

import pytest

THIS = os.path.abspath(os.path.dirname(__file__))
EM_DIR = os.path.abspath(os.path.join(THIS, ".."))
BED_SCRIPT = os.path.join(EM_DIR, "asm_bic_to_bed.py")

HEADER = (
    "coord_evaluated\twindow_last_coord\tBICsinglecomp\tBICmiddlecomp\tBICresult\t"
    "BIC_winner\twindows_size\ttotal_reads\tem_iterations\tmean_theta_diff\twindow_cpgs"
)

# BICresult is int(BICmiddlecomp < BICsinglecomp), i.e. exactly dBIC > 0, so the flag
# and the sign agree by construction. What these tests pin down is that the BED
# reports that sign the way the manuscript defines it.
#
#   coord 1000  ASM      single=120.0 middle=100.0  -> dBIC = +20.0
#   coord 2000  shared   single=100.0 middle=115.5  -> dBIC = -15.5
#   coord 3000  ASM      near-tie                   -> dBIC = +0.000331
ROW_ASM = "1000\t1040\t120.0\t100.0\t1\t1\t10\t30\t12\t0.48\t1000,1010,1040"
ROW_SHARED = "2000\t2040\t100.0\t115.5\t0\t0\t10\t30\t9\t0.03\t2000,2010,2040"

# A near-tie ASM window, taken from a real chr22 window (coord 20146759). At 2 decimal
# places this printed as "0.00", making a red ASM window look like a tie; the emitted
# text has to keep the sign. 3 of chr22's 8,251 ASM windows are in this regime.
ROW_ASM_NEAR_TIE = (
    "3000\t3040\t258.1997970473942\t258.19946603769466\t1\t1\t10\t30\t40\t0.02\t3000,3010,3040"
)

ASM_RGB = "228,26,28"
NON_ASM_RGB = "55,126,184"


@pytest.fixture(scope="module")
def bed_run(tmp_path_factory):
    """Run asm_bic_to_bed.py over a three-window TSV and return the parsed outputs."""
    work = tmp_path_factory.mktemp("bed")
    asm_dir = work / "ASM"
    asm_dir.mkdir()
    (asm_dir / "chr9_ASM.tsv").write_text(
        "\n".join([HEADER, ROW_ASM, ROW_SHARED, ROW_ASM_NEAR_TIE]) + "\n"
    )
    prefix = str(work / "out")
    subprocess.check_call(
        [sys.executable, BED_SCRIPT, str(asm_dir), "-o", prefix],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    def rows(path, skip_track=False):
        lines = [l for l in open(path).read().splitlines() if l.strip()]
        if skip_track:
            lines = [l for l in lines if not l.startswith("track ")]
        return [l.split("\t") for l in lines]

    return {
        "scored": rows(prefix + "_scored.bed"),
        "merged": rows(prefix + "_merged.bed"),
        "colored": rows(prefix + "_all_colored.bed", skip_track=True),
        "colored_raw": open(prefix + "_all_colored.bed").read(),
    }


def by_start(rows, start):
    hits = [r for r in rows if r[1] == start]
    assert len(hits) == 1, f"expected exactly one row starting at {start}"
    return hits[0]


def test_scored_bed_holds_only_asm_windows(bed_run):
    scored = bed_run["scored"]
    assert len(scored) == 2, "only BICresult==1 windows belong in _scored.bed"
    assert sorted(r[1] for r in scored) == ["1000", "3000"]
    # half-open end covers the last CpG base: window_last_coord + 1
    assert by_start(scored, "1000")[2] == "1041"
    assert all(r[0] == "chr9" for r in scored)


def test_asm_window_delta_bic_is_positive(bed_run):
    """dBIC = BIC_Shared - BIC_ASM, so an ASM window must score POSITIVE.

    This is the manuscript's convention. If someone reinstates
    BICmiddlecomp - BICsinglecomp, this is the test that fails.
    """
    delta = float(by_start(bed_run["scored"], "1000")[3])
    assert delta > 0, f"ASM window must have dBIC > 0, got {delta}"
    assert delta == pytest.approx(20.0)


def test_shared_window_delta_bic_is_negative(bed_run):
    colored = bed_run["colored"]
    assert len(colored) == 3, "every analyzed window belongs in _all_colored.bed"
    shared = by_start(colored, "2000")
    assert float(shared[3]) == pytest.approx(-15.5)
    assert float(shared[3]) < 0, "a shared window must have dBIC < 0"


def test_near_tie_asm_window_keeps_its_sign_in_the_text(bed_run):
    """A window that only just wins must not render as a signless 0.

    Real chr22 data has ASM windows with dBIC as small as 0.00033. Printed at 2
    decimals those became "0.00" — a red window that reads as a tie. The emitted
    string itself, not just the float behind it, has to stay positive.
    """
    for source in ("scored", "colored"):
        text = by_start(bed_run[source], "3000")[3]
        assert float(text) > 0, f"{source}: near-tie ASM window rendered as {text!r}"
        assert float(text) == pytest.approx(0.000331, abs=1e-5)


def test_sign_of_delta_bic_agrees_with_the_colour(bed_run):
    """dBIC > 0 and the ASM colour must never disagree.

    BICresult drives the colour while the arithmetic drives the sign; they are
    equivalent by construction, and this catches it if that ever stops being true.
    """
    for row in bed_run["colored"]:
        delta, rgb = float(row[3]), row[8]
        if rgb == ASM_RGB:
            assert delta > 0, f"red (ASM) window with dBIC {delta} <= 0"
        else:
            assert rgb == NON_ASM_RGB
            assert delta <= 0, f"blue (non-ASM) window with dBIC {delta} > 0"


def test_col5_score_is_non_negative_for_every_window(bed_run):
    """BED column 5 must be a non-negative int in [0, 1000].

    Guards the abs() that becomes tempting to drop once dBIC > 0 means ASM:
    non-ASM windows are negative, so without abs() they would clamp to 0 and lose
    all display intensity.
    """
    for row in bed_run["colored"]:
        score = int(row[4])
        assert 0 <= score <= 1000
        assert score == round(abs(float(row[3])))
    shared = by_start(bed_run["colored"], "2000")
    assert int(shared[4]) == 16, "|-15.5| must survive as magnitude, not clamp to 0"


def test_track_line_documents_the_convention(bed_run):
    raw = bed_run["colored_raw"]
    assert raw.startswith("track "), "_all_colored.bed needs the itemRgb track header"
    assert 'itemRgb="On"' in raw
    assert "BIC_Shared - BIC_ASM" in raw, "the track description must state the sign convention"


def test_merged_bed_is_bed3_over_asm_windows_only(bed_run):
    merged = bed_run["merged"]
    # The two ASM windows are far apart, so they stay separate blocks; the shared
    # window between them must not appear at all.
    assert merged == [["chr9", "1000", "1041"], ["chr9", "3000", "3041"]]
