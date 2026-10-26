Tractor Prep
============

A Python tool for creating `*.dosage.txt[.gz]` and `*.hapcount.txt[.gz]` files,
from FLARE or RFMix outputs (VCF and MSP files), for downstream Tractor
analysis. Three commands are exposed: `extract_rfmix`, `extract_flare`, and
`merge_tracts`. Example usage can be found in `test_tractor_prep_flare.sh` and
`test_tractor_prep_rfmix.sh`.

If extracting tracts from RFMix outputs, the input VCF **must** be formatted as
Phased Genotypes (`GT`) only; i.e. the main body of each row of the VCF is a
repeated 4-character pattern of `N|N\t` for `N` in `{0, 1}`. Similary, the MSP
file **must** be formatted as a repeating 2-character pattern of `N\t`. I.e.,
the subpopulation codes are tab-delimited single integers.

If extracting tracts from FLARE outputs, the input VCF file **must** have
format `GT:AN0:AN1`. I.e., an 8-character repeating pattern of `N|N:A:A\t`,
with `N` in `{0,1}` and `A` a single-integer subpopulation code.

If the inputs have, for whatever reason, been split by sample IIDs across
files, then each split should be extracted from independently (2-way unkinking
and tract extraction are per-sample). The resulting files can then be merged by
`merge_tracts`.

See `python -m tractor_prep --help` and the example scripts for usage.
