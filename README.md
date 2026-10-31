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


Attribution
-----------
The functionality implemented here was developed and first implemented by
Atkinson Lab. Please see [https://github.com/Atkinson-Lab/Tractor].

If you use Tractor in your research, please cite the following article:
> Atkinson, E.G., Maihofer, A.X., Kanai, M. et al. Tractor uses local ancestry
> to enable the inclusion of admixed individuals in GWAS and to boost power.
> Nat Genet 53, 195–204 (2021).
> [Link](https://doi.org/10.1038/s41588-020-00766-y)
