# Preparing the input BAM

BICAMERA needs one BAM per sample, aligned, indexed, and carrying two kinds of tag.

## What the tags are

**MM/ML** — per-read base modification calls, written by Dorado during basecalling.
These are the methylation data; BICAMERA reads them with modkit and cannot work
without them. The launcher samples the first 5000 reads and aborts if it finds no
MM tag, so a wrongly-basecalled BAM fails immediately rather than producing empty
output.

**HP** — haplotype assignment, `HP:i:1` or `HP:i:2`, written by a phasing tool.
Reads with no HP tag are labelled `noH` and still used: they inform the shared
methylation pattern even though they can't be assigned to an allele.

## Getting there

Basecall with modification calling on, align, call variants, phase, then haplotag:

```bash
# 1. basecall with 5mC calling (produces MM/ML)
dorado basecaller <model> pod5/ --modified-bases 5mCG_5hmCG > calls.bam

# 2. align to your reference
dorado aligner ref.fa calls.bam > aligned.bam

# 3. variants + phasing, then haplotag (writes HP)
#    Clair3 -> WhatsHap phase -> WhatsHap haplotag is one common route.
```

Any pipeline that ends in an indexed BAM with MM/ML and HP tags will do. The
reference you pass to BICAMERA must be the one the BAM was aligned to.

## ONT: strip 5hmC first

If you basecalled with a 5mC+5hmC model, run the prep step before the pipeline:

```bash
bash prep_ont_5mc_mapq10.sh -q 10 sample.bam
# prints the path of the BAM to use
```

It does two things:

1. `modkit adjust-mods --ignore h` — drops 5hmC calls.
2. `samtools view -q 10` — keeps MAPQ ≥ 10.

The first matters more than it looks. With 5hmC left in, `modkit extract calls`
emits `h` calls that BICAMERA's filter silently discards, *and* the
canonical-vs-5mC decision becomes a three-way call instead of two-way. Both make
ONT results incomparable to 5mC-only HiFi data. Stripping 5hmC keeps the two
platforms on the same footing.

The step is idempotent and writes atomically, so re-running it is harmless.

## Reference builds

BICAMERA doesn't verify that your reference matches the BAM — it can't. On the
primary chromosomes the common GRCh38 flavours agree, but they differ at masked,
decoy and alt loci. If you mix builds, coordinates in those regions will be wrong
and nothing will warn you. Check the BAM header against the FASTA you intend to
pass.
