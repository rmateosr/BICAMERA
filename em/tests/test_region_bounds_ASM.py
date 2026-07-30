# ABOUTME: Tests the filename-bound window-start rule in EMBIC_bin_path_ASM.
# ABOUTME: Regions are standalone/non-overlapping; every in-region CpG that fits a full window is evaluated.
import os
import uuid

import numpy as np
import pandas as pd
import pytest

import EMfunctions_ASM as em


def _synth_coords(J, offset, rng):
    spacings = rng.exponential(200, size=J).astype(int).clip(min=10)
    return (np.cumsum(spacings) + offset).astype(int)


def _synth_reads(coords, n_reads, rng, p_meth=0.5):
    """Minimal read generator: each read covers a random contiguous CpG span."""
    J = len(coords)
    rows = []
    for _ in range(n_reads):
        span = rng.randint(3, J + 1)
        start = rng.randint(0, J - span + 1)
        sel = slice(start, start + span)
        status = (rng.rand(span) < p_meth).astype(int)
        hap = rng.choice(["H1", "H2", "noH"], p=(0.3, 0.3, 0.4))
        rows.append({
            "readlabel": uuid.uuid4().hex[:16],
            "chrom": "chrSYN",
            "strand": "+",
            "startcoord": ",".join(str(c) for c in coords[sel]),
            "status": ",".join(str(s) for s in status),
            "haplotype": hap,
        })
    return rows


def _write_region(rows, path):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_parse_region_bounds_round_trip(tmp_path):
    p = tmp_path / "methylationfraction_12345_67890_.tsv"
    p.write_text("")
    assert em.parse_region_bounds(str(p)) == (12345, 67890)


def test_parse_region_bounds_rejects_bad_filename(tmp_path):
    bad = tmp_path / "not_a_region.tsv"
    bad.write_text("")
    with pytest.raises(ValueError):
        em.parse_region_bounds(str(bad))


def _full_span_read(coords_slice):
    """A single read covering every CpG in coords_slice, so J_full is deterministic."""
    return {
        "readlabel": "fullspan",
        "chrom": "chrSYN",
        "strand": "+",
        "startcoord": ",".join(str(int(c)) for c in coords_slice),
        "status": ",".join("1" for _ in coords_slice),
        "haplotype": "noH",
    }


def test_no_cede_evaluates_through_region_right_edge(tmp_path):
    """Standalone, non-overlapping region: every in-region CpG that can start a
    full W-window is evaluated, including the rightmost one. Region = the whole
    CpG ladder, so with J_full = J the window-starts are indices [0, J-W]
    (J-W+1 windows); coords[J-W] (the right-edge window) must be present.
    """
    W = 10
    rng = np.random.RandomState(123)

    J = 30
    coords = _synth_coords(J, offset=3_000_000, rng=rng)
    s, e = int(coords[0]), int(coords[J - 1])

    # Full-span read guarantees every CpG is present (J_full == J), plus random
    # reads for depth.
    rows = _synth_reads(coords, 50, rng)
    rows.append(_full_span_read(coords))
    p = tmp_path / f"methylationfraction_{s}_{e}_.tsv"
    _write_region(rows, p)

    out = em.EMBIC_bin_path_ASM(str(p), windowsize=W, asm_restarts=1)
    coord_set = set(out["coord_evaluated"].astype(int).tolist())

    expected = {int(coords[i]) for i in range(0, J - W + 1)}  # [0, J-W]
    assert coord_set == expected, (
        f"mismatch: extra={sorted(coord_set - expected)}, "
        f"missing={sorted(expected - coord_set)}"
    )
    assert int(coords[J - W]) in coord_set       # right-edge window evaluated
    assert all(s <= c <= e for c in coord_set)   # never past region end


def test_window_step_strides_window_starts(tmp_path):
    """window_step keeps every k-th window-start. With J_full=J=30, W=10 the
    step-1 starts are indices [0, J-W]=[0,20] (21 windows); step 5 keeps indices
    {0,5,10,15,20} (5 windows), a strict subset; step 1 reproduces the default.
    """
    W = 10
    rng = np.random.RandomState(321)
    J = 30
    coords = _synth_coords(J, offset=4_000_000, rng=rng)
    s, e = int(coords[0]), int(coords[J - 1])
    rows = _synth_reads(coords, 50, rng)
    rows.append(_full_span_read(coords))
    p = tmp_path / f"methylationfraction_{s}_{e}_.tsv"
    _write_region(rows, p)

    out1 = em.EMBIC_bin_path_ASM(str(p), windowsize=W, asm_restarts=1, window_step=1)
    out5 = em.EMBIC_bin_path_ASM(str(p), windowsize=W, asm_restarts=1, window_step=5)
    set1 = set(out1["coord_evaluated"].astype(int).tolist())
    set5 = set(out5["coord_evaluated"].astype(int).tolist())

    assert set1 == {int(coords[i]) for i in range(0, J - W + 1)}      # step 1 = every start
    assert set5 == {int(coords[i]) for i in range(0, J - W + 1, 5)}   # step 5 = strided
    assert set5 < set1                                               # strict subset
    assert len(set5) == 5

    with pytest.raises(ValueError):
        em.EMBIC_bin_path_ASM(str(p), windowsize=W, asm_restarts=1, window_step=0)


def test_small_region_uses_adaptive_window(tmp_path):
    """A region with 3 <= J < windowsize in-region CpGs is evaluated as ONE
    whole-region window (adaptive), not skipped: exactly one coord_evaluated row,
    windows_size == J, anchored at the first in-region CpG."""
    W = 10
    rng = np.random.RandomState(99)
    J = 6  # 3 <= 6 < W
    coords = _synth_coords(J, offset=5_000_000, rng=rng)
    s, e = int(coords[0]), int(coords[J - 1])
    rows = _synth_reads(coords, 40, rng)
    rows.append(_full_span_read(coords))
    p = tmp_path / f"methylationfraction_{s}_{e}_.tsv"
    _write_region(rows, p)

    out = em.EMBIC_bin_path_ASM(str(p), windowsize=W, asm_restarts=1)
    assert len(out) == 1                                       # one whole-region window
    assert int(out["windows_size"].iloc[0]) == J              # window shrunk to J
    assert int(out["coord_evaluated"].iloc[0]) == int(coords[0])


def test_tail_rescue_anchors_final_window(tmp_path):
    """Right-edge tail rescue (always on): when striding leaves orphan trailing
    CpGs, a final full-W window anchored on the last in-region CpG is appended.
    J=17, W=10, step=5: strided starts are indices {0,5}; the last 2 CpGs
    (indices 15,16) would be orphaned, so the anchored start at index J-W=7 is
    added, whose window covers indices 7..16 — every window stays W=10 CpGs wide.
    """
    W = 10
    rng = np.random.RandomState(2026)
    J = 17
    coords = _synth_coords(J, offset=6_000_000, rng=rng)
    s, e = int(coords[0]), int(coords[J - 1])
    rows = _synth_reads(coords, 50, rng)
    rows.append(_full_span_read(coords))
    p = tmp_path / f"methylationfraction_{s}_{e}_.tsv"
    _write_region(rows, p)

    out = em.EMBIC_bin_path_ASM(str(p), windowsize=W, asm_restarts=1, window_step=5)
    coord_set = set(out["coord_evaluated"].astype(int).tolist())

    # {0,5} strided + {7} anchored; window at index 7 covers indices 7..16.
    assert coord_set == {int(coords[0]), int(coords[5]), int(coords[J - W])}
    assert all(int(w) == W for w in out["windows_size"])         # every window is W CpGs
    # The anchored window ends on the last in-region CpG (the orphan tail is scored).
    anchored = out[out["coord_evaluated"].astype(int) == int(coords[J - W])]
    assert int(anchored["window_last_coord"].iloc[0]) == int(coords[J - 1]) == e


def test_tail_rescue_noop_when_no_orphans(tmp_path):
    """When (n_in - W) % window_step == 0 the last strided window already reaches
    the last CpG, so no anchored window is added. J=15, W=10, step=5: (15-10)%5==0
    -> strided starts {0,5} already end on the last CpG."""
    W = 10
    rng = np.random.RandomState(2027)
    J = 15
    coords = _synth_coords(J, offset=7_000_000, rng=rng)
    s, e = int(coords[0]), int(coords[J - 1])
    rows = _synth_reads(coords, 50, rng)
    rows.append(_full_span_read(coords))
    p = tmp_path / f"methylationfraction_{s}_{e}_.tsv"
    _write_region(rows, p)

    out = em.EMBIC_bin_path_ASM(str(p), windowsize=W, asm_restarts=1, window_step=5)
    assert set(out["coord_evaluated"].astype(int)) == {int(coords[0]), int(coords[5])}
    assert len(out) == 2


def test_hanging_cpgs_outside_filename_are_not_evaluated(tmp_path):
    """A read whose CpGs extend past regionend supports boundary windows but must
    NOT generate coord_evaluated rows past regionend."""
    W = 10
    rng = np.random.RandomState(7)

    J_chrom = 40
    coords = _synth_coords(J_chrom, offset=2_000_000, rng=rng)

    r_left, r_right = 0, 19
    s, e = int(coords[r_left]), int(coords[r_right])

    # Reads spill to index 29 — past regionend (index 19). A full-span read over
    # 0..29 fixes J_full = 30 so the in-region window-starts are deterministic.
    rows = _synth_reads(coords[r_left:30], 30, rng)
    rows.append(_full_span_read(coords[r_left:30]))
    p = tmp_path / f"methylationfraction_{s}_{e}_.tsv"
    _write_region(rows, p)

    out = em.EMBIC_bin_path_ASM(str(p), windowsize=W, asm_restarts=1)
    coord_set = set(out["coord_evaluated"].astype(int).tolist())

    leaked = {c for c in coord_set if c > e}
    assert leaked == set(), (
        f"coord_evaluated leaked past regionend ({e}): {sorted(leaked)}"
    )

    # No-cede: every in-region CpG (indices 0..19) fits a full window
    # (19 + W = 29 <= J_full = 30), so all are evaluated.
    expected = {int(coords[i]) for i in range(r_left, r_right + 1)}
    assert coord_set == expected, (
        f"coord_set mismatch: extra={sorted(coord_set - expected)}, "
        f"missing={sorted(expected - coord_set)}"
    )
