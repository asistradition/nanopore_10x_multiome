"""
Generated with cursor/claude-3.5-sonnet and then fixed to actually work
"""

import os
import pytest
import pysam
import tempfile

from nanopore_10x_multiome.utils._bam import split_bam_by_barcode, write_bam_record

@pytest.fixture
def temp_bam(tmp_path):
    """Create a temporary BAM file with test data."""
    bam_path = tmp_path / "test.bam"
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'chr1'}]
    }
    
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        # Write test records with different barcodes
        write_bam_record(out, "read1", "ACGT", "FFFF", flag=0, CB="AAAA")
        write_bam_record(out, "read2", "ACGT", "FFFF", flag=0, CB="BBBB")
        write_bam_record(out, "read3", "ACGT", "FFFF", flag=0, CB="AAAA")
        write_bam_record(out, "read4", "ACGT", "FFFF", flag=4, CB="CCCC")  # unmapped
        write_bam_record(out, "read5", "ACGT", "FFFF", flag=0, CB="DDDD")  # not in lookup
    
    # Index the BAM file
    pysam.index(str(bam_path))
    return bam_path


def test_written_bam(temp_bam, tmp_path):
    with pysam.AlignmentFile(temp_bam, "rb") as infile:
        records = [r for r in infile]

        assert len(records) == 5
        assert [r.get_tag('CB') for r in records] == ['AAAA', 'BBBB', 'AAAA', 'CCCC', 'DDDD']
        assert [r.query_sequence for r in records] == ['ACGT'] * 5


def test_standard_split(temp_bam, tmp_path):
    """Test standard BAM splitting with valid inputs."""
    lookup_table = {
        "AAAA": "group1",
        "BBBB": "group2",
        "CCCC": "group3"
    }
    
    result = split_bam_by_barcode(
        temp_bam,
        lookup_table,
        out_path=str(tmp_path),
        out_prefix="test_"
    )
    
    assert result == {"group1": 2, "group2": 1, "group3": 0}
    assert os.path.exists(tmp_path / "test_group1_bam.bam")
    assert os.path.exists(tmp_path / "test_group2_bam.bam")
    assert os.path.exists(tmp_path / "test_group3_bam.bam")

def test_custom_output_files(temp_bam, tmp_path):
    """Test splitting with custom output file paths."""
    lookup_table = {"AAAA": "group1", "BBBB": "group2"}
    output_files = {
        "group1": str(tmp_path / "custom1.bam"),
        "group2": str(tmp_path / "custom2.bam")
    }
    
    result = split_bam_by_barcode(
        temp_bam,
        lookup_table,
        output_files=output_files
    )
    
    assert result == {"group1": 2, "group2": 1}
    assert os.path.exists(tmp_path / "custom1.bam")
    assert os.path.exists(tmp_path / "custom2.bam")

def test_mismatched_output_files(temp_bam, tmp_path):
    """Test error handling for mismatched output files."""
    lookup_table = {"AAAA": "group1", "BBBB": "group2"}
    output_files = {
        "group1": str(tmp_path / "custom1.bam"),
        "wrong_group": str(tmp_path / "custom2.bam")
    }
    
    with pytest.raises(ValueError, match="Non-overlapping output_files"):
        split_bam_by_barcode(
            temp_bam,
            lookup_table,
            output_files=output_files
        )

def test_empty_bam(tmp_path):
    """Test handling of empty BAM file."""
    # Create empty BAM file
    empty_bam = tmp_path / "empty.bam"
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'chr1'}]
    }
    with pysam.AlignmentFile(empty_bam, "wb", header=header) as out:
        pass
    
    pysam.index(str(empty_bam))
    
    lookup_table = {"AAAA": "group1"}
    result = split_bam_by_barcode(
        empty_bam,
        lookup_table,
        out_path=str(tmp_path)
    )
    
    assert result == {"group1": 0}

def test_progress_bar(temp_bam, tmp_path):
    """Test progress bar functionality."""
    lookup_table = {"AAAA": "group1"}
    result = split_bam_by_barcode(
        temp_bam,
        lookup_table,
        out_path=str(tmp_path),
        pbar=True
    )
    
    assert result == {"group1": 2}

def test_custom_barcode_tag(temp_bam, tmp_path):
    """Test using a different barcode tag."""
    # Create BAM with different tag
    bam_path = tmp_path / "custom_tag.bam"
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1000, 'SN': 'chr1'}]
    }
    
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        write_bam_record(out, "read1", "ACGT", "FFFF", flag=0, BC="AAAA")
    
    pysam.index(str(bam_path))
    
    lookup_table = {"AAAA": "group1"}
    result = split_bam_by_barcode(
        bam_path,
        lookup_table,
        out_path=str(tmp_path),
        barcode_tag='BC'
    )
    
    assert result == {"group1": 1}