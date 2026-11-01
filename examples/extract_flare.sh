#!/bin/bash
set -euo pipefail
#### TIMING RESULTS AND SUMMARY DATA ##########################################
# Timings are single-cpu; process never required a RES size above 200 Mb. All
# 22 chromosome source VCFs used the same Sample ID set of ~250k samples.
#
# CPU: Intel(R) Xeon(R) CPU E5-2683 v4 @ 2.10GHz
#
# If we assume that local ancestry inference has been run for some number of
# supopulations, e.g. 6: ##ANCESTRY=<CSA=0,AFR=1,OCE=2,EUR=3,MID=4,EAS=5>, then
# to extract hapcount and dosages for all six ancestries one would use the flag
# `--ancestors=6`. If, however, one were only interested in using two
# subpopulations (say AFR and EUR) in downstream analysis, one could use the
# flag: `--ancestors 1 3` to specify a list of ancestries to extract. In the
# following stats, the first time is to extract all six ancestries, the second
# time is to extract only AFR and EUR (the two most populous in this dataset).
#
# N Samples: 249,927 (AGD 250k)
# Chromosome | Filesize | n positions | `time` all | `time` subset
# ------------------------------------------------------------------------
# chr1  | 7.2G  | 49106     | 29m31.680s    | 14m2.195s
# chr2  | 8.1G  | 55504     | 33m22.149s    | 15m4.164s
# chr3  | 7.0G  | 48002     | 29m4.948s     | 13m20.501s
# chr4  | 6.4G  | 43412     | 27m21.410s    | 11m35.564s
# chr5  | 5.7G  | 38700     | 23m55.450s    | 10m21.016s
# chr6  | 7.0G  | 46088     | 27m53.850s    | 12m24.881s
# chr7  | 5.4G  | 36746     | 22m10.307s    | 9m50.844s
# chr8  | 5.0G  | 34499     | 20m45.196s    | 9m16.088s
# chr9  | 4.2G  | 28541     | 17m16.853s    | 7m52.254s
# chr10 | 4.7G  | 32081     | 19m24.479s    | 8m37.584s
# chr11 | 4.5G  | 30858     | 18m34.022s    | 8m17.105s
# chr12 | 4.4G  | 29854     | 17m54.579s    | 7m59.287s
# chr13 | 3.4G  | 23067     | 13m56.816s    | 6m9.456s
# chr14 | 3.1G  | 20810     | 12m32.033s    | 5m43.000s
# chr15 | 2.9G  | 19468     | 11m45.769s    | 5m14.426s
# chr16 | 3.2G  | 21730     | 13m6.545s     | 5m49.334s
# chr17 | 2.8G  | 18851     | 11m25.120s    | 5m6.644s
# chr18 | 2.9G  | 19228     | 11m41.839s    | 5m13.570s
# chr19 | 2.0G  | 13830     | 8m37.575s     | 3m43.396s
# chr20 | 2.5G  | 16685     | 10m7.753s     | 4m30.807s
# chr21 | 1.4G  | 9037      | 5m30.842s     | 2m26.375s
# chr22 | 1.5G  | 9972      | 6m1.281s      | 2m46.517s
# -----------------------------------------------------
# TOTAL | 95G   | 646069    | 392m0.639s    | 175m25.297s
#
# Total output directory size: all ancestries 171G; AFR/EUR only 122G


#### RUNNING THE SCRIPT ON FLARE INPUTS #######################################
RESULTS="./results"
TIMINGS="./timings"
mkdir -p "$RESULTS" "$TIMINGS"

function run_extraction() {
    STEM="agd_250k_chr$1_flare"
    { time python -m tractor_prep extract_flare "./src_vcf/$STEM.anc.vcf.gz" \
        --outdir "./$RESULTS/" \
        --stem "$STEM" \
        --vcf-gzip \
        --out-gzip \
        --num-ancestors 6 ;
    } 2>"$TIMINGS/chrom_$1.txt"
}
function run_sequential() {
    for chrom in `seq 1 22`; do
        run_extraction $chrom
    done
}
# Run all the chromosomes sequentially
{ time run_sequential ;} 2> "$TIMINGS/total.txt"


#### TESTING THE RESULTS ######################################################
function test_output() {
    STEM="agd_250k_chr$1_flare"
    TRACT="./src_tracts/agd_250k_chr$1_flare.anc.anc$2.$3.txt.gz"
    RECON="./results/agd_250k_chr$1_flare.anc$2.$3.txt.gz"
    # Need to decompress the output files before they can be tested with diff
    diff -q <(zcat "$TRACT") <(zcat "$RECON")
}
# check each chromosome's ancs outputs in parallel against the original
for chrom in `seq 1 22`; do
    for idx in `seq 0 5`; do
        test_output "$chrom" "$idx" "dosage" &
        test_output "$chrom" "$idx" "hapcount" &
    done
    wait
done


#### GENERATE SUMMARY DATA (see start of this file) ###########################
mkdir -p vcf_info
SAMPLES='./samples.txt'
SUMMARY="./summary.txt"

# Check that all the files have the same samples
bcftools query -l "./src_vcf/agd_250k_chr22_flare.anc.vcf.gz" >${SAMPLES}
for i in `seq 1 22`; do
    VCF="./src_vcf/agd_250k_chr${i}_flare.anc.vcf.gz"
    diff -q ${SAMPLES} <(bcftools query -l "$VCF") &
done
wait

# Write number of samples to the summary file
wc -l ${SAMPLES} >${SUMMARY}
echo >>${SUMMARY}

# Extract information about each VCF file and how long it took to run
for i in `seq 1 22`; do
    VCF="./src_vcf/agd_250k_chr${i}_flare.anc.vcf.gz"
    SIZE=`du -sh "$VCF" | cut -f1 | tr -d '\n'`
    NPOS=`bcftools query -f '%POS\n' "$VCF" | wc -l | tr -d '\n'`
    TIME=`sed -n '2 p' timings/chrom_${i}.txt | cut -f2 | tr -d '\n'`
    printf "${i}\t${SIZE}\t${NPOS}\t${TIME}\n" >>$SUMMARY
    echo "Processed ${i}"
done
