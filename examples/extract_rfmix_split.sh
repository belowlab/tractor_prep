#!/bin/bash
set -euo pipefail
# Example using RFMix outputs (GT formatted VCFs and MSP files). Six of the
# chromosomes have been split into two parts each before running RFMix and need
# to be recombined after tract extraction.

CHR_SPLIT=(1 2 3 4 5 9)
CHR_NOSPLIT=(6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22)

IIDS="./sample_ids.txt"
mkdir -p results/{splits,tracts,unkinked}
mkdir -p timings

function run_extraction_single() {
    { time python -m tractor_prep extract_rfmix "./src_vcf/chrom_$1.vcf.gz" \
        --msp="./src_msp/chrom_$1.msp" \
        --outdir "./results/tracts/" \
        --iids "./sample_ids.txt" \
        --stem "chrom_$1" \
        --vcf-gzip \
        --out-gzip \
        --num-ancestors 6 \
        --unkink \
        --emit-msp \
        --emit-vcf ;
    } 2>"timings/chrom_$1.txt"
    # Move the unkinked files out of the tracts directory
    mv "./results/tracts/chrom_$1_switches.msp.gz" results/unkinked/
    mv "./results/tracts/chrom_$1_unkinked.msp.gz" results/unkinked/
    mv "./results/tracts/chrom_$1_unkinked.vcf.gz" results/unkinked/
}

function run_merge() {
    { time python -m tractor_prep merge_tracts \
        --vcf-gzip \
        --out-gzip \
        --outdir="./results/tracts/" \
        --output="chrom_$1.anc$2.$3.txt.gz" \
        "./results/splits/chrom_$1_pt1.anc$2.$3.txt.gz" \
        "./results/splits/chrom_$1_pt2.anc$2.$3.txt.gz" ;
    } 2>"timings/merge_chrom_$1.txt"
}

function run_extraction_split() {
    for part in `seq 1 2`; do
        STEM="chrom_$1_pt${part}"
        { time python -m tractor_prep extract_rfmix "./src_vcf/${STEM}.vcf.gz" \
            --msp="./src_msp/${STEM}.msp" \
            --outdir "./results/splits/" \
            --iids "./sample_ids.txt" \
            --stem "${STEM}" \
            --vcf-gzip \
            --out-gzip \
            --num-ancestors 6 \
            --unkink \
            --emit-msp \
            --emit-vcf ;
        } 2>"timings/${STEM}.txt"
        # Move the unkinked files out of the tracts directory
        mv "./results/splits/${STEM}_switches.msp.gz" results/unkinked/
        mv "./results/splits/${STEM}_unkinked.msp.gz" results/unkinked/
        mv "./results/splits/${STEM}_unkinked.vcf.gz" results/unkinked/
    done
    for anc in `seq 0 5`; do
        run_merge $1 $anc "dosage"
        run_merge $1 $anc "hapcount"
    done
}
# run strictly sequentially (no multi-threaded ISAL either). It would be easy
# to run all in parallel, or to run the merging step only in parallel, etc.
for chrom in "${CHR_NOSPLIT[@]}"; do
    run_extraction_single $chrom
done
for chrom in "${CHR_SPLIT[@]}"; do
    run_extraction_split $chrom
done
