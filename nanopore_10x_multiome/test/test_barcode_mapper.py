"""
Generated with cursor/claude-3.5-sonnet and then fixed to actually work
"""

import pytest
from nanopore_10x_multiome.barcodes._correct_barcodes import barcode_correction_table

def test_empty_barcode_list():
    """Test behavior with empty input list"""
    result = barcode_correction_table([])
    assert result == {}

def test_single_barcode():
    """Test with a single barcode"""
    result = barcode_correction_table(["ACGT"])
    # Should include original and all 1-edit distance variants that aren't ambiguous
    assert result["ACGT"] == "ACGT"  # Original maps to itself
    assert result["CCGT"] == "ACGT"  # Substitution
    assert result["CGT"] == "ACGT"   # Deletion
    assert result["AACGT"] == "ACGT" # Insertion

def test_ambiguous_corrections():
    """Test handling of ambiguous corrections"""
    # Two barcodes that differ by 2 edits
    result = barcode_correction_table(["AAAA", "CCAA"])
    assert result["AAAA"] == "AAAA"  # Original maps to itself
    assert result["CCAA"] == "CCAA"  # Original maps to itself
    assert "CAAA" not in result      # Ambiguous correction (could be either)
    assert "ACAA" not in result      # Ambiguous correction (could be either)

def test_similar_barcodes():
    """Test with similar barcodes that could create ambiguity"""
    # Barcodes that differ by only one base
    result = barcode_correction_table(["ACGT", "ACTT"])
    assert result["ACGT"] == "ACGT"  # Original maps to itself
    assert result["ACTT"] == "ACTT"  # Original maps to itself

    assert "CCTT" in result
    assert "ACAT" not in result
    assert "ACCT" not in result

def test_progress_bar():
    """Test that progress bar option doesn't affect results"""
    barcodes = ["AAAA", "CCCC"]
    result_with_pbar = barcode_correction_table(barcodes, pbar=True)
    result_without_pbar = barcode_correction_table(barcodes, pbar=False)
    assert result_with_pbar == result_without_pbar

def test_special_characters():
    """Test handling of N in barcodes"""
    result = barcode_correction_table(["ACGT"])
    assert "ANGT" in result  # N substitution should be correctable
    assert result["ANGT"] == "ACGT"

def test_identical_barcodes():
    """Test with duplicate barcodes in input"""
    result = barcode_correction_table(["AAAA", "AAAA"])
    assert result["AAAA"] == "AAAA"
    # Should behave same as if there was only one copy

def test_different_length_barcodes():
    """Test with barcodes of different lengths"""
    result = barcode_correction_table(["AAA", "AAAA"])
    assert result["AAA"] == "AAA"
    assert result["AAAA"] == "AAAA"
    # Insertions/deletions that could map to either should be excluded
    assert "AAAT" not in result
