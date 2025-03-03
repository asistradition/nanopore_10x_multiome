"""
Generated with cursor/claude-3.5-sonnet and then fixed to actually work
"""


import pytest
import regex
from nanopore_10x_multiome.utils import get_barcode_parasail, RC

@pytest.fixture
def setup_regex_and_comparison():
    # Example pattern for testing - adjust based on your actual use case
    pattern = r'(?:ACGT){3}(.{16})(?:TGCA){3}'
    compiled_regex = regex.compile(pattern, regex.BESTMATCH)
    comparison_seq = 'ACGTACGTACGTNNNNNNNNNNNNNNNNTGCATGCATGCA'
    return compiled_regex, comparison_seq


def test_successful_barcode_match(setup_regex_and_comparison):
    compiled_regex, comparison_seq = setup_regex_and_comparison
    # Perfect match scenario
    seq = 'ACGTACGTACGTATCGATCGATCGATCGTGCATGCATGCA'
    qual = 'F' * len(seq)  # High quality scores
    bc_len = 16

    barcode, quality, position = get_barcode_parasail(
        seq, qual, compiled_regex, comparison_seq, bc_len
    )

    assert barcode == 'ATCGATCGATCGATCG'
    assert quality == 'F' * bc_len
    assert position == 12  # Position where barcode starts


def test_failed_barcode_match_because_RC(setup_regex_and_comparison):
    compiled_regex, comparison_seq = setup_regex_and_comparison
    # Perfect match scenario
    seq = RC('ACGTACGTACGTATCGATCGATCGATCGTGCATGCATGCA')
    qual = 'F' * len(seq)  # High quality scores
    bc_len = 16

    barcode, quality, position = get_barcode_parasail(
        seq, qual, compiled_regex, comparison_seq, bc_len
    )

    assert barcode is None
    assert quality is None
    assert position is None

def test_no_regex_match(setup_regex_and_comparison):
    compiled_regex, comparison_seq = setup_regex_and_comparison
    # Sequence that won't match the regex pattern
    seq = 'GGGGCCCCAAAATTTT'
    qual = 'F' * len(seq)
    bc_len = 16

    barcode, quality, position = get_barcode_parasail(
        seq, qual, compiled_regex, comparison_seq, bc_len
    )

    assert barcode is None
    assert quality is None
    assert position is None

def test_barcode_with_mismatches(setup_regex_and_comparison):
    compiled_regex, comparison_seq = setup_regex_and_comparison
    # Sequence with one mismatch in the barcode region
    seq = 'ACGTACGTACGTANCGATCGATCGATCGTGCATGCATGCA'
    qual = 'F' * len(seq)
    bc_len = 16

    barcode, quality, position = get_barcode_parasail(
        seq, qual, compiled_regex, comparison_seq, bc_len
    )

    assert barcode == 'ANCGATCGATCGATCG'
    assert quality == 'F' * bc_len
    assert position == 12

def test_short_barcode(setup_regex_and_comparison):
    compiled_regex, comparison_seq = setup_regex_and_comparison
    # Sequence that would result in a barcode shorter than bc_len
    seq = 'ACGTACGTACGTATCGATGCATGCATGCA'  # Missing some barcode bases
    qual = 'F' * len(seq)
    bc_len = 16

    barcode, quality, position = get_barcode_parasail(
        seq, qual, compiled_regex, comparison_seq, bc_len
    )

    assert barcode is None
    assert quality is None
    assert position is None

def test_barcode_with_extra_base(setup_regex_and_comparison):
    compiled_regex, comparison_seq = setup_regex_and_comparison
    # Sequence that might introduce gaps in alignment
    seq = 'ACGTACGTACGTATCGAATCGATCGATCGTGCATGCATGCA'
    qual = 'F' * len(seq)
    bc_len = 16

    barcode, quality, position = get_barcode_parasail(
        seq, qual, compiled_regex, comparison_seq, bc_len
    )

    assert barcode is None
    assert quality is None
    assert position is None


def test_quality_string_slicing(setup_regex_and_comparison):
    compiled_regex, comparison_seq = setup_regex_and_comparison
    seq = 'ACGTACGTACGTATCGATCGATCGATCGTGCATGCATGCA'
    qual = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' + '!' * 15  # Varying quality scores
    bc_len = 16

    barcode, quality, position = get_barcode_parasail(
        seq, qual, compiled_regex, comparison_seq, bc_len
    )

    assert len(quality) == bc_len
    assert quality == qual[position:position + bc_len]
