"""
NOTE: These are primarily not unit tests but end-to-end tests, designed to be
run interactively in a specific environment with access to large data files.
For example:

    $ ipython -i test_tractor_prep.py
    In [1]: %load_ext line_profiler
    In [2]: %lprun -f unkink_msp test_unkink_msp()
    ... # etc

"""
import csv
import filecmp
import pathlib

from tractor_prep import *


def files_eq(f1, f2):
    filecmp.clear_cache()
    return filecmp.cmp(f1, f2, shallow=False)


def cleanup_test_dir(skip=None, glob=None):
    test_dir = pathlib.Path('./test_data')
    assert test_dir.exists()

    keep = [
        # MSP and related files
        'test.msp.tsv',                 # 25M
        'test.msp.switches.tsv',        # 36M
        'test.msp.unkinked.tsv',        # 25M
        # PhasedVCF and related files
        'phased_chr1_pt1.vcf',          # 12.6G
        'phased_chr1_pt1.unkinked.vcf', # 12.6G
    ]
    if skip is not None:
        keep.extend(skip)
    if glob is not None:
        keep.extend([f.name for f in test_dir.glob(glob)])

    for path in test_dir.iterdir():
        if path.name not in keep:
            path.unlink()


def _test_roundtrip(cls, test_file, cleanup=True):
    tsv_path = pathlib.Path(f'./test_data/{test_file}')
    assert tsv_path.exists(), tsv_path

    alt_path = tsv_path.with_name(f'{tsv_path.name}.2')
    npz_path = tsv_path.with_suffix('.npz')

    txt = cls.from_tsv(tsv_path)
    txt.to_npz(npz_path)
    npz = cls.from_npz(npz_path)
    npz.to_tsv(alt_path)

    assert txt == npz
    assert files_eq(tsv_path, alt_path)

    if cleanup:
        npz_path.unlink()
        alt_path.unlink()


def test_msp_roundtrip(cleanup=True):
    _test_roundtrip(MSP, 'test.msp.tsv', cleanup)

def test_vcf_roundtrip(cleanup=True):
    _test_roundtrip(PhasedVCF, 'phased_chr1_pt1.vcf', cleanup)

def test_msp_switches_roundtrip(cleanup=True):
    _test_roundtrip(MSPSwitches, 'test.msp.switches.tsv', cleanup)


def test_unkink_msp():
    tsv_path = pathlib.Path(f'./test_data/test.msp.tsv')
    assert tsv_path.exists()

    switches_path = tsv_path.with_suffix('.switches.tsv')
    unkinked_path = tsv_path.with_suffix('.unkinked.tsv')

    msp = MSP.from_tsv(tsv_path)
    switches_orig = MSPSwitches.from_tsv(switches_path)
    unkinked_orig = MSP.from_tsv(unkinked_path)

    unkinked, switches = unkink_msp(msp, inplace=False)
    assert not unkinked == msp
    assert unkinked == unkinked_orig
    assert switches == switches_orig

    unkinked, switches = unkink_msp(msp, inplace=True)
    assert unkinked == msp
    assert unkinked == unkinked_orig
    assert switches == switches_orig


def test_unkink_geno(setup_npz=False):
    if setup_npz:
        vcf = PhasedVCF.from_tsv('./test_data/phased_chr1_pt1.vcf')
        vcf.to_npz('./test_data/phased_chr1_pt1.npz')
        unk = PhasedVCF.from_tsv('./test_data/phased_chr1_pt1.unkinked.vcf')
        unk.to_npz('./test_data/phased_chr1_pt1.unkinked.npz')
    else:
        vcf = PhasedVCF.from_npz('./test_data/phased_chr1_pt1.npz')
        unk = PhasedVCF.from_npz('./test_data/phased_chr1_pt1.unkinked.npz')
    switches = MSPSwitches.from_tsv('./test_data/test.msp.switches.tsv')
    unkinked = unkink_geno(vcf, switches, inplace=False)
    assert unk == unkinked
