#!/bin/bash
# Example usage and roundtrip testing of extraction from RFMix outputs.
set -euo pipefail

# IID cohort list
IIDS="test_rfmix/sample_id_list.txt"

# Parts and Merged directories (PTSDIR for pt1 / pt2 outputs, etc)
PTSDIR="test_rfmix/tract_chr1_parts"
MERGED="test_rfmix/tract_chr1"

mkdir -p "$PTSDIR"
mkdir -p "$MERGED"

function run_extraction() {
    SOURCE_MSP="test_rfmix/AGD_$1_RFmix_dose.vcf.gz.msp.tsv.gz"
    SOURCE_VCF="test_rfmix/phased_$1_agd35k_test5.vcf.gz"

    time python -m tractor_prep extract_rfmix "$SOURCE_VCF" \
        --outdir "$PTSDIR" \
        --stem "$1" \
        --vcf-gzip \
        --msp-gzip \
        --out-gzip \
        --num-ancestors 6 \
        --msp "$SOURCE_MSP" \
        --unkink \
        --emit-msp \
        --emit-vcf \
        --cohort "$IIDS"
}

function run_merge() {
    time python -m tractor_prep merge_tracts \
        --vcf-gzip \
        --out-gzip \
        --outdir="$MERGED" \
        --output="chr1.$1.txt.gz" \
        "$PTSDIR/chr1_pt1.$1.txt.gz" "$PTSDIR/chr1_pt2.$1.txt.gz"
}

# Run each part's extraction in tandem
run_extraction "chr1_pt1" &
run_extraction "chr1_pt2" &
wait

# Finally, merge the output dosage / hapcount files (for split inputs)
for idx in `seq 0 5`; do
    run_merge "anc${idx}.dosage" &
    run_merge "anc${idx}.hapcount" &
done
wait
