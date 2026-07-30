# ABOUTME: Join the per-chr read_id->HP map onto the single modkit extract for one chromosome.
# ABOUTME: Streams per haplotype (H1->H2->noH order); applies the negative-strand off-by-one fix.
import sys
import os

# Input chromosome (e.g., "chr1") from command line argument

if len(sys.argv) != 3:
    print("Usage: python script.py <chromosome> <output_folder>")
    sys.exit(1)

chromosome = sys.argv[1]
OUTPUT_FOLDER = sys.argv[2]
# Define the directory paths and filenames based on the input chromosome
directory_path = OUTPUT_FOLDER + '/modkit_splitbam_referenced_bychrandhap_filtered/'
directory_path_output = OUTPUT_FOLDER + '/modkit_referenced_splitbam_mergedbychr/'
print(directory_path)
print(chromosome)

# Ensure the output directory exists
os.makedirs(directory_path_output, exist_ok=True)

output_file_name = os.path.join(directory_path_output, f'modkit_extract_{chromosome}_merged_filtered_coordmod.tsv')
tmp_file_name = output_file_name + '.tmp'

# Single filtered extract for this chromosome (one modkit pass over the whole
# per-chr BAM). The haplotype label is carried by a separate read_id->HP map built
# from primary alignments (see hpmap_inner.sh) rather than by file provenance.
filtered_path = os.path.join(directory_path, f'modkit_extract_{chromosome}_filtered.tsv')
hpmap_path = os.path.join(OUTPUT_FOLDER, 'hpmap', f'{chromosome}.tsv')

# Load the read_id -> HP map. modkit extract only emits primary alignments, so every
# read_id in the filtered extract MUST be present exactly once in this primary-only
# map; a read with no HP tag maps to 'noH'.
#
# That invariant is load-bearing, so it is verified rather than assumed. A read
# absent from the map would fall through to 'noH' below — turning a phased read
# unphased, weakening or erasing the ASM signal, and reporting no error at all. A
# truncated map (upstream samtools dying mid-pipe) is exactly how that happens, so
# any unmatched read_id aborts this stage instead of being silently relabelled.
read_to_hap = {}
if os.path.exists(hpmap_path):
    with open(hpmap_path) as fh:
        for line in fh:
            rid, hap = line.rstrip('\n').split('\t')
            read_to_hap[rid] = hap

# Columns (0-indexed) in the modkit extract we keep: read_id (0), ref_position (2),
# chrom (3), ref_strand (5, used only for the negative-strand off-by-one fix),
# call_code (11). call_code carries the methylation decision (m | -); we drop
# call_prob (10) because it is the confidence of the winning class, not
# P(methylated) (see CHANGES.md 2026-05-25 modkit filter fix).
#
# Output structure (no header, tab-separated): read_id, ref_position, chrom, call_code, Haplotag
# Emitted in H1 -> H2 -> noH order (one streaming pass per haplotype over the single
# filtered extract) so the row order is identical to the old per-haplotype concatenation.
# An empty haplotype bin is a valid state (e.g. an untagged BAM has 0 H1/H2 reads;
# a haplotagged-only BAM has 0 noH reads). A missing / 0-byte extract is treated as empty.
min_coord = None
max_coord = None
have_input = os.path.exists(filtered_path) and os.path.getsize(filtered_path) > 0
# Unmatched read_ids, counted on the first pass only (all three passes see the same
# rows). A few examples are kept for the error message; the count is what matters.
missing_count = 0
missing_examples = []
with open(tmp_file_name, 'w') as out:
    for hap in ('H1', 'H2', 'noH'):
        if not have_input:
            continue
        with open(filtered_path) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                hap_of_read = read_to_hap.get(f[0])
                if hap_of_read is None:
                    if hap == 'H1':
                        missing_count += 1
                        if len(missing_examples) < 3:
                            missing_examples.append(f[0])
                    hap_of_read = 'noH'
                if hap_of_read != hap:
                    continue
                pos = int(f[2])
                # Adjust coordinate for the negative strand (off-by-one fix).
                if f[5] == '-':
                    pos -= 1
                out.write(f'{f[0]}\t{pos}\t{f[3]}\t{f[11]}\t{hap}\n')
                if min_coord is None or pos < min_coord:
                    min_coord = pos
                if max_coord is None or pos > max_coord:
                    max_coord = pos

if missing_count:
    os.remove(tmp_file_name)
    sys.exit(
        f"ERROR: {chromosome}: {missing_count} extract rows reference read_ids absent "
        f"from the haplotype map.\n"
        f"       map     : {hpmap_path} ({len(read_to_hap)} reads)\n"
        f"       extract : {filtered_path}\n"
        f"       examples: {', '.join(missing_examples)}\n"
        f"       modkit extract emits primary alignments only, so every read_id must "
        f"appear in the map. A shortfall means the map is truncated (see the hpmap "
        f"stage log). Continuing would relabel those reads noH and quietly weaken the "
        f"ASM signal, so this stage aborts. Re-run hpmap for {chromosome}."
    )

os.replace(tmp_file_name, output_file_name)

# Print the min and max coordinates for the final merged file
print(f"{chromosome} Merged Min and Max:", (min_coord, max_coord))
