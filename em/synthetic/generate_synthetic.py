# ABOUTME: Generates synthetic methylation data for evaluating the static ASM pipeline.
# ABOUTME: Produces one per-region TSV per region (no time periods) plus an answer-key manifest.
"""Synthetic data generator for static (Shared vs ASM) BIC pipeline evaluation.

Generates region files in the format consumed by BIC_algorithm_ASM.py:
  - Tab-separated with header: readlabel, chrom, strand, startcoord, status, haplotype
  - One file per region: methylationfraction_<start>_<end>_.tsv

Two region classes:
  SHARED      -- identical methylation on both alleles (no ASM signal)
  STATIC_ASM  -- alleles differ; the allele-specific signal the ASM model detects

The diploid per-allele marginal constraint pi[h, 0] = 0.5 is enforced in the EM;
this generator places every read into one of the two alleles (H1/H2/noH) and
samples methylation from a per-allele theta. No altered/time-emergent component
exists.

Also writes manifest.tsv (answer key): one row per sliding window, with the
expected BICresult (1 if STATIC_ASM, 0 if SHARED).

Note: the SHARED class is named "SHARED" rather than "NULL" because pandas'
read_csv treats the literal "NULL" as a missing value by default.
"""

import argparse
import os
import uuid

import numpy as np
import pandas as pd

SHARED = "SHARED"
STATIC_ASM = "STATIC_ASM"
CLASSES = [SHARED, STATIC_ASM]


def generate_coords(J, chrom_offset, rng):
    """Generate J sorted CpG coordinates with jittered spacing (mean 200 bp)."""
    spacings = rng.exponential(200, size=J).astype(int).clip(min=10)
    coords = np.cumsum(spacings) + chrom_offset
    return coords.astype(int)


def generate_region_theta(
    region_class, J, rng,
    shared_range=(0.0, 1.0),
    asm_hyper_range=(0.65, 1.00),
    asm_hypo_range=(0.00, 0.35),
):
    """Generate true per-allele theta for a region.

    Each theta[j, h, k] is drawn independently per CpG from the (h, k) uniform
    range; there is no per-region base value and no per-CpG noise.

    Parameters
    ----------
    region_class : str
        One of CLASSES.
    J : int
        Number of CpGs.
    rng : np.random.RandomState
    shared_range : tuple of 2 floats
        (lo, hi) uniform range for SHARED regions; the same per-CpG value is
        placed on both alleles (no ASM signal). Default (0.0, 1.0) — matches
        the SHARED_RANGE used by the predictive_power_2026-06 sweep.
    asm_hyper_range, asm_hypo_range : tuple of 2 floats
        (lo, hi) uniform ranges for the two STATIC_ASM alleles at k=0 (the ASM
        signal). Defaults (0.65, 1.00) and (0.00, 0.35) — the ASM methylation
        bands reported in Nature s41467-025-57433-1.

    Returns
    -------
    theta : ndarray, shape (J, 2, 2)
        theta[j, h, k] = P(methylated at CpG j | allele h, component k).
        Only k=0 is observed at sampling time (no altered component); k=1 is
        filled with free per-CpG values so the EM has to learn "k=1 unused"
        from pi alone, not from a degenerate theta similarity.
    """
    theta = np.full((J, 2, 2), 0.5)

    if region_class == SHARED:
        vals = rng.uniform(shared_range[0], shared_range[1], J)
        theta[:, :, :] = vals[:, None, None]

    elif region_class == STATIC_ASM:
        # k=0: strong allele difference (the ASM signal).
        theta[:, 0, 0] = rng.uniform(asm_hyper_range[0], asm_hyper_range[1], J)
        theta[:, 1, 0] = rng.uniform(asm_hypo_range[0], asm_hypo_range[1], J)
        # k=1: free per-CpG values (unused by sampler — pi[:,1]=0).
        theta[:, 0, 1] = rng.uniform(0.01, 0.99, J)
        theta[:, 1, 1] = rng.uniform(0.01, 0.99, J)

    else:
        raise ValueError(f"Unknown region_class: {region_class}. Expected one of {CLASSES}.")

    return theta


def generate_reads_for_region(
    coords, theta, n_reads, chrom, rng,
    hap_probs=(0.3, 0.3, 0.4), min_cpg_span=3,
):
    """Generate reads for one region under the static (no-altered-component) model.

    Each read is assigned a true allele (h=0 or h=1) and a true k=0 (altered
    component is disallowed). Methylation status per CpG is sampled from
    Bernoulli(theta[j, h, 0]). H1/H2 reads carry their true allele in the
    haplotype tag; noH reads have h drawn 50/50 from the two alleles and the
    tag set to "noH".

    Parameters
    ----------
    coords : ndarray of int, shape (J,)
    theta : ndarray, shape (J, 2, 2)
    n_reads : int
    chrom : str
    rng : np.random.RandomState
    hap_probs : tuple of 3 floats
        (P(H1), P(H2), P(noH)). Must sum to 1.
    min_cpg_span : int
        Minimum CpGs a read must cover.

    Returns
    -------
    rows : list of dict
    """
    J = len(coords)
    rows = []

    for _ in range(n_reads):
        hap = rng.choice(["H1", "H2", "noH"], p=hap_probs)
        if hap == "H1":
            true_h = 0
        elif hap == "H2":
            true_h = 1
        else:
            true_h = int(rng.choice([0, 1]))
        true_k = 0

        span = rng.randint(min_cpg_span, J + 1)
        start = rng.randint(0, J - span + 1)
        sel = slice(start, start + span)
        sel_coords = coords[sel]

        probs = theta[sel, true_h, true_k]
        status = (rng.rand(span) < probs).astype(int)

        rows.append({
            "readlabel": uuid.uuid4().hex[:16],
            "chrom": chrom,
            "strand": "+",
            "startcoord": ",".join(str(c) for c in sel_coords),
            "status": ",".join(str(s) for s in status),
            "haplotype": hap,
        })

    return rows


def generate_synthetic_dataset(
    n_regions=50,
    reads_per_region=80,
    J=50,
    class_mix=None,
    windowsize=10,
    chrom="chrSYN",
    output_dir=None,
    seed=None,
    hap_probs=None,
):
    """Generate a complete synthetic dataset for static ASM evaluation.

    Writes per-region TSVs into output_dir/ directly (no t1/.../tT/ subdirs)
    and a manifest.tsv answer key.

    Parameters
    ----------
    n_regions : int
    reads_per_region : int
    J : int
        CpG sites per region.
    class_mix : dict, optional
        {class_name: weight} over CLASSES. Default: equal weights.
    windowsize : int
        Sliding-window size used to enumerate manifest rows.
    chrom : str
    output_dir : str
    seed : int, optional
    hap_probs : tuple of 3 floats, optional

    Returns
    -------
    manifest : DataFrame
    """
    if class_mix is None:
        class_mix = {c: 1 for c in CLASSES}
    if output_dir is None:
        output_dir = "synthetic_output"
    if hap_probs is None:
        hap_probs = (0.3, 0.3, 0.4)
    hap_probs = tuple(float(x) for x in hap_probs)
    if len(hap_probs) != 3 or abs(sum(hap_probs) - 1.0) > 1e-6:
        raise ValueError(
            f"hap_probs must be a length-3 tuple summing to 1; got {hap_probs}"
        )

    rng = np.random.RandomState(seed)

    weights = np.array([class_mix.get(c, 0) for c in CLASSES], dtype=float)
    weights /= weights.sum()
    region_classes = rng.choice(CLASSES, size=n_regions, p=weights)

    os.makedirs(output_dir, exist_ok=True)

    manifest_rows = []
    offset = 1_000_000

    for ridx in range(n_regions):
        rc = region_classes[ridx]
        coords = generate_coords(J, offset, rng)
        rstart = int(coords[0])
        rend = int(coords[-1])
        offset = rend + 100_000

        theta = generate_region_theta(rc, J, rng)

        reads = generate_reads_for_region(
            coords, theta, reads_per_region, chrom, rng, hap_probs=hap_probs,
        )
        fname = f"methylationfraction_{rstart}_{rend}_.tsv"
        df = pd.DataFrame(reads)
        df.to_csv(os.path.join(output_dir, fname), sep="\t", index=False)

        n_windows = max(1, J - windowsize + 1)
        for w in range(n_windows):
            coord_eval = int(coords[w]) if J > windowsize else int(coords[0])
            manifest_rows.append({
                "region_idx": ridx,
                "chrom": chrom,
                "region_start": rstart,
                "region_end": rend,
                "coord_evaluated": coord_eval,
                "true_class": rc,
                "true_theta_h0_k0": float(np.mean(theta[:, 0, 0])),
                "true_theta_h1_k0": float(np.mean(theta[:, 1, 0])),
                "expected_BICresult": 0 if rc == SHARED else 1,
            })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(os.path.join(output_dir, "manifest.tsv"), sep="\t", index=False)
    return manifest


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic methylation data for static ASM evaluation."
    )
    parser.add_argument("--n-regions", type=int, default=50)
    parser.add_argument("--reads-per-region", type=int, default=80)
    parser.add_argument("--cpgs-per-region", type=int, default=50, dest="J")
    parser.add_argument("--windowsize", type=int, default=10)
    parser.add_argument("--chrom", type=str, default="chrSYN")
    parser.add_argument("--output-dir", type=str, default="synthetic_output")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--class-mix", type=str, default=None,
        help="Comma-separated weights for SHARED,STATIC_ASM (e.g. 1,1)."
    )
    parser.add_argument(
        "--hap-probs", type=str, default=None,
        help="Comma-separated P(H1),P(H2),P(noH) summing to 1. Default 0.3,0.3,0.4."
    )
    parser.add_argument(
        "--noh-frac", type=float, default=None,
        help="Fraction of reads tagged noH; H1/H2 split the remainder evenly. "
             "Mutually exclusive with --hap-probs."
    )
    args = parser.parse_args()

    class_mix = None
    if args.class_mix:
        vals = [float(x) for x in args.class_mix.split(",")]
        class_mix = dict(zip(CLASSES, vals))

    if args.hap_probs is not None and args.noh_frac is not None:
        parser.error("--hap-probs and --noh-frac are mutually exclusive")
    hap_probs = None
    if args.hap_probs is not None:
        hap_probs = tuple(float(x) for x in args.hap_probs.split(","))
    elif args.noh_frac is not None:
        f = args.noh_frac
        hap_probs = ((1 - f) / 2, (1 - f) / 2, f)

    manifest = generate_synthetic_dataset(
        n_regions=args.n_regions,
        reads_per_region=args.reads_per_region,
        J=args.J,
        class_mix=class_mix,
        windowsize=args.windowsize,
        chrom=args.chrom,
        output_dir=args.output_dir,
        seed=args.seed,
        hap_probs=hap_probs,
    )
    print(f"Wrote {len(manifest)} manifest rows to {args.output_dir}/manifest.tsv")
    print(f"Wrote {args.n_regions} region TSVs to {args.output_dir}/")


if __name__ == "__main__":
    main()
