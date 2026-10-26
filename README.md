Tractor Prep
============
A Python tool for creating `*.dosage.txt[.gz]` and `*.hapcount.txt[.gz]` files,
from FLARE or RFMix outputs (VCF and MSP files), for downstream Tractor
analysis. Three commands are exposed: `extract_rfmix`, `extract_flare`, and
`merge_tracts`. Example usage can be found in `test_tractor_prep_flare.sh` and
`test_tractor_prep_rfmix.sh`.

If extracting tracts from RFMix outputs, the input VCF **must** be formatted as
Phased Genotypes (`GT`) only; i.e. the main body of each row of the VCF is a
repeated 4-character pattern of `N|N\t` for `N` in `{0,1}`. Similary, the MSP
file **must** be formatted as a repeating 2-character pattern of `N\t`. I.e.,
the subpopulation codes are tab-delimited single integers.

If extracting tracts from FLARE outputs, the input VCF file **must** have
format `GT:AN1:AN2`. I.e., an 8-character repeating pattern of `N|N:A:A\t`,
with `N` in `{0,1}` and `A` a single-integer subpopulation code.

If the inputs have, for whatever reason, been split by sample IIDs across
files, then each split should be extracted from independently (2-way unkinking
and tract extraction are per-sample). The resulting files can then be merged by
`merge_tracts`.

See `python -m tractor_prep --help` and the example scripts for usage.


Performance
-----------
All commands are implemented as line-by-line stream filters. Therefore the
memory overhead should be low.

Under the assumption that input VCF files have fixed-width formats for the
sample information (i.e. `GT` for RFMix and `GT:AN1:AN2` for FLARE), line
processing is very efficient.

Running the example test scripts using single-core processing, without gzip
multi-threading, and ISAL installed:

`test_tractor_prep_flare.sh` extracts from a VCF file for chrom 21 with approx.
9k rows and 250k samples in 5m46s. The directory of compressed outputs requires
2.5G disk space. It extracts from a VCF file for chrom 1 with 49k rows and 250k
samples in 29m25s.

`test_tractor_prep_rfmix.sh` extracts from two VCF files, each with 1,117,575
rows and a total of 35k samples between the two files. An IID file with 5015
IIDs was used. The pt1 file is filtered to 3006 samples, and the pt2 is
filtered down to 2409. The MSP and VCF files are unkinked and written to disk.
Pt1 was processed in 15m23s, pt2 in 12m48s. The resulting six sets of pt1/pt2
output files were merged in between 2 and 3.5 minutes each. The compressed
outputs require 7.5G total.
