#!/bin/bash
# Example usage and roundtrip test for extraction from FLARE outputs.
CHROM=21
OUTD="test_flare/recon_${CHROM}"
STEM="agd_250k_chr${CHROM}_flare"
INFILE="./test_flare/$STEM.anc.vcf"


#### File processing options ##################################################
## 1) Decompress the input file (can use any decompress cmd) (fastest)
#if [[ ! -f "$INFILE" ]]; then
#    pigz -d --keep "$INFILE.gz";
#fi
#python -m tractor_prep_flare "$INFILE" \
#    --outdir="./$OUTD" \
#    --stem="$STEM" \
#    --num-ancestors=6

## 2) Pipe a decompressed stream into the python script stdin
#pigz -dc --keep "$INFILE.gz" |\
#python -m tractor_prep_flare "stdin" \
#    --outdir="./$OUTD" \
#    --stem="$STEM" \
#    --num-ancestors=6

## 3) Read a gzip compressed file directly (only supports gzip)
#python -m tractor_prep_flare "$INFILE.gz" --vcf-gzip \
#    --outdir="./$OUTD" \
#    --stem="$STEM" \
#    --num-ancestors=6

## Do something like this to compress the output files if the python script
## writes ascii (patterns 1-3 above). Any compression program can be used.
#for file in `readlink -f ./$OUTD/*`; do
#    pigz --keep --best $file
#done

# 4) Read *and write* gzip compressed files directly (only supports gzip).
#    Threading only supported if ISAL installed, otherwise ignored. "0" is
#    single threaded mode, "-1" uses as many cores as ISAL wants.
#    XXX: Threading does not seem faster - does not seem input throttled. IPC
#    issues? Buffer size issues? Anyway, THREADS=0 is faster by far for this
#    set of test data. Or omit --gzip-threads.
THREADS=0
python -m tractor_prep_flare "$INFILE.gz" \
    --vcf-gzip \
    --out-gzip --gzip-threads=${THREADS} \
    --outdir="./$OUTD" \
    --stem="$STEM" \
    --num-ancestors=6
# Now need to decompress files before they can be tested with diff (see below)
for file in `readlink -f ./$OUTD/*`; do
    pigz --keep --decompress $file
done


#### Check the output files ###################################################
# Check the tracts for diff
for idx in `seq 0 5`; do
    diff -q "./test_flare/tract_${CHROM}/${STEM}.anc.anc${idx}.dosage.txt" \
            "./test_flare/recon_${CHROM}/${STEM}.anc${idx}.dosage.txt"

    diff -q "./test_flare/tract_${CHROM}/${STEM}.anc.anc${idx}.hapcount.txt" \
            "./test_flare/recon_${CHROM}/${STEM}.anc${idx}.hapcount.txt"
done
