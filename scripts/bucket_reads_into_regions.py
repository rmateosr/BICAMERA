# ABOUTME: Bucket a chromosome's per-read methylation rows into per-region TSVs.
# ABOUTME: Replaces the divide step; region boundaries come from potential_asm_regions/<chr>.bed.
import argparse
import bisect
import os
import sys

# Usage: python bucket_reads_into_regions.py <chr> <output_folder>
#                                            [--min-region-cpgs N]
#
# Inputs
#   <out>/read_format/<chr>_reads.tsv           (header: readlabel haplotype chrom startcoord status)
#   <out>/potential_asm_regions/<chr>.bed       (chrom  start  end ; 0-based, non-overlapping, sorted)
#
# Output (one file per region, the input layout the ASM EM driver consumes):
#   <out>/read_format_split/<chr>/methylationfraction_<start>_<end>_.tsv
#
# A read is written to every region it overlaps: read_min_cpg <= region_end AND
# read_max_cpg >= region_start (the same overlap test the divide R script used).
# Regions are standalone and non-overlapping, so run the EM driver with
# --standalone-regions (no rightmost-W ceding).

parser = argparse.ArgumentParser(
    description="Bucket a chromosome's per-read methylation rows into per-region TSVs."
)
parser.add_argument("chrom", help="chromosome name, e.g. chr15")
parser.add_argument("output_folder", help="pipeline output dir for this sample")
# Drop regions with fewer than this many distinct in-region CpGs: too small to
# carry a meaningful Shared-vs-ASM test. The EM evaluates a single whole-region
# window for surviving regions below windowsize.
#
# Passed explicitly by asmregion_bucket_inner.sh from config/defaults.sh. It used
# to be read from $MIN_REGION_CPGS here, which meant a lost environment silently
# filtered regions at a different threshold than the run was launched with.
parser.add_argument(
    "--min-region-cpgs",
    type=int,
    default=3,
    help="drop regions with fewer distinct in-region CpGs (default: 3)",
)
args = parser.parse_args()

chrom = args.chrom
OUTPUT_FOLDER = args.output_folder
MIN_REGION_CPGS = args.min_region_cpgs

reads_path = os.path.join(OUTPUT_FOLDER, "read_format", f"{chrom}_reads.tsv")
bed_path = os.path.join(OUTPUT_FOLDER, "potential_asm_regions", f"{chrom}.bed")
out_dir = os.path.join(OUTPUT_FOLDER, "read_format_split", chrom)
os.makedirs(out_dir, exist_ok=True)

# ---- load regions (sorted, non-overlapping) ----
region_starts = []
region_ends = []
if os.path.exists(bed_path):
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            region_starts.append(int(f[1]))
            region_ends.append(int(f[2]))

if not region_starts:
    print(f"{chrom}: 0 potential-ASM regions; nothing to bucket.")
    sys.exit(0)

# Regions are emitted by bedtools merge already sorted by start and disjoint,
# so ends are sorted too. Keep an explicit guard in case that ever changes.
order = sorted(range(len(region_starts)), key=lambda i: region_starts[i])
region_starts = [region_starts[i] for i in order]
region_ends = [region_ends[i] for i in order]

HEADER = "readlabel\thaplotype\tchrom\tstartcoord\tstatus\n"

# Accumulate clipped output lines per region, then write one file at a time.
# Do NOT open one handle per region: a chromosome can have tens of thousands of
# regions (genome-wide, >100k), and holding that many open file objects blows up
# memory / the open-file ulimit (observed: MemoryError on chr22's 8623 regions).
# buckets: region index -> list of output line strings.
buckets = {}

n_assignments = 0
with open(reads_path) as fh:
    fh.readline()  # skip header
    for line in fh:
        row = line.rstrip("\n")
        if not row:
            continue
        cols = row.split("\t")
        if len(cols) < 5:
            continue
        coords = cols[3].split(",")
        statuses = cols[4].split(",")
        try:
            cpgs = [int(float(c)) for c in coords]
        except ValueError:
            continue
        if not cpgs or len(statuses) != len(cpgs):
            continue
        rmin = min(cpgs)
        rmax = max(cpgs)
        # Overlapping regions form a contiguous index range because regions are
        # sorted and disjoint:
        #   lo = first region whose end >= rmin
        #   hi = last region whose start <= rmax
        lo = bisect.bisect_left(region_ends, rmin)
        hi = bisect.bisect_right(region_starts, rmax) - 1
        prefix = f"{cols[0]}\t{cols[1]}\t{cols[2]}\t"
        for idx in range(lo, hi + 1):
            s = region_starts[idx]
            e = region_ends[idx]
            # CLIP the read's CpGs to this region's interval. ONT reads are long
            # and span several regions; writing the full CpG list would let the
            # EM's coords_array (and thus a 10-CpG window near the region edge)
            # pull in CpGs from neighbouring regions across the inter-region gap.
            # Each region file must be self-contained and local. EM's in-region
            # test is coords >= region_start & <= region_end, so clip to [s, e].
            keep = [j for j, c in enumerate(cpgs) if s <= c <= e]
            if not keep:
                continue
            kc = ",".join(str(cpgs[j]) for j in keep)
            ks = ",".join(statuses[j] for j in keep)
            buckets.setdefault(idx, []).append(f"{prefix}{kc}\t{ks}\n")
            n_assignments += 1

# Write one region file at a time (single open handle), atomic .tmp -> rename.
# Regions with no overlapping reads get no file; the EM driver simply does not
# evaluate them (nothing to test there).
n_written = 0
n_dropped_small = 0
for idx, lines in buckets.items():
    s = region_starts[idx]
    e = region_ends[idx]
    # Distinct in-region CpG coords across all reads in this region (clipped col 3).
    coords = set()
    for ln in lines:
        parts = ln.split("\t")
        if len(parts) >= 4:
            coords.update(parts[3].split(","))
    coords.discard("")
    if len(coords) < MIN_REGION_CPGS:
        n_dropped_small += 1
        continue
    final = os.path.join(out_dir, f"methylationfraction_{s}_{e}_.tsv")
    tmp = final + ".tmp"
    with open(tmp, "w") as out:
        out.write(HEADER)
        out.writelines(lines)
    os.replace(tmp, final)
    n_written += 1

print(f"{chrom}: {len(region_starts)} regions, "
      f"{len(buckets)} with reads, {n_assignments} read-region assignments; "
      f"{n_written} written, {n_dropped_small} dropped (<{MIN_REGION_CPGS} CpGs).")
