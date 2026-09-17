"""
NOTE: These are not unit tests but end-to-end tests, designed to be run
interactively in a specific environment with access to large data files. E.g.:

    $ ipython -i test_tractor_prep.py
    In [1]: %load_ext line_profiler
    In [2]: %lprun -f unkink_msp test_unkink_msp()
    ... # etc
"""
import filecmp
import functools
import pathlib

import numpy as np
from tractor_prep import *


test_files = [
    # MSP unkinking
    'test.msp.tsv',                 # 25M
    'test.msp.switches.tsv',        # 36M
    'test.msp.unkinked.tsv',        # 25M
    # VCF unkinking
    'phased_chr1_pt1.vcf',          # 12.6G
    'phased_chr1_pt1.unkinked.vcf', # 12.6G
    # TODO Tract Extraction
]


def files_eq(f1, f2):
    filecmp.clear_cache()
    return filecmp.cmp(f1, f2, shallow=False)

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

def test_anc_roundtrip(cleanup=True):
    _test_roundtrip(AncVCF, 'tracts/phased_chr1.anc0.dosage.tsv', cleanup)


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


def setup_test_extract_tracts():
    vcf = PhasedVCF.from_tsv('./test_data/phased_chr1.unkinked_stripped.vcf')
    vcf.to_npz('./test_data/phased_chr1.unkinked_stripped.npz')
    msp = MSP.from_tsv('./test_data/chr1.msp.tsv')
    msp.to_npz('./test_data/chr1.msp.npz')


def test_extract_tracts(trunc=None):
    msp = MSP.from_npz('./test_data/chr1.msp.npz')
    vcf = PhasedVCF.from_npz('./test_data/phased_chr1.unkinked_stripped.npz')
    if trunc is not None:
        vcf.index = vcf.index[:trunc]
        vcf.X0 = vcf.X0[:trunc]
        vcf.X1 = vcf.X1[:trunc]
    multi_anc = MultiAncVCF(*extract_tracts(vcf, msp, 6))
    multi_anc.to_npz('./test_data/chr1_all_tracts.npz')


def rewrite_tracts():
    data = MultiAncVCF.from_npz('./test_data/chr1_all_tracts.npz')
    for idx in data.ancestors:
        hapfile_old = f'./test_data/tracts/phased_chr1.anc{idx}.hapcount.tsv'
        hapfile_new = hapfile_old + '.2'
        data.get_anc(idx, 'hapcount').to_tsv(hapfile_new)
        assert files_eq(hapfile_old, hapfile_new)

        dosfile_old = f'./test_data/tracts/phased_chr1.anc{idx}.dosage.tsv'
        dosfile_new = dosfile_old + '.2'
        data.get_anc(idx, 'dosage').to_tsv(dosfile_new)
        assert files_eq(dosfile_old, dosfile_new)
