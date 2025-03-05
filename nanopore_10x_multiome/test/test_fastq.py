"""
Generated with cursor/claude-3.5-sonnet and then fixed to actually work
"""

import pytest
from io import StringIO
from nanopore_10x_multiome.utils._fastq import fastqProcessor, convert_qual_illumina, fastq_gen
from nanopore_10x_multiome.utils import RC

def test_convert_qual_illumina():
    # Test basic quality conversion
    assert convert_qual_illumina("!") == [0]  # Lowest quality (33)
    assert convert_qual_illumina("I") == [40]  # Common quality value
    assert convert_qual_illumina("~") == [93]  # Highest quality (126)
    assert convert_qual_illumina("!!~") == [0, 0, 93]

def create_fastq_file(sequences):
    """Helper function to create StringIO objects with FASTQ content"""
    content = ""
    for name, seq, qual in sequences:
        content += f"@{name}\n{seq}\n+\n{qual}\n"
    return StringIO(content)

def test_single_fastq_processing():
    # Test processing a single FASTQ file
    data = [("read1", "ACGT", "IIII")]
    fh = create_fastq_file(data)
    
    processor = fastqProcessor()
    result = next(processor.fastq_gen(fh))
    
    assert len(result) == 1
    control, sequence, quality = result[0]
    assert control == "@read1"
    assert sequence == "ACGT"
    assert quality == [40, 40, 40, 40]

def test_paired_fastq_processing():
    # Test processing paired FASTQ files
    data1 = [("read1", "ACGT", "IIII")]
    data2 = [("read1", "TGCA", "IIII")]
    
    fh1 = create_fastq_file(data1)
    fh2 = create_fastq_file(data2)
    
    processor = fastqProcessor()
    result = next(processor.fastq_gen(fh1, fh2))
    
    assert len(result) == 2
    assert result[0][1] == "ACGT"
    assert result[1][1] == "TGCA"

def test_id_verification():
    # Test ID verification functionality
    data1 = [("read1", "ACGT", "IIII")]
    data2 = [("read2", "TGCA", "IIII")]  # Mismatched ID
    
    fh1 = create_fastq_file(data1)
    fh2 = create_fastq_file(data2)
    
    processor = fastqProcessor(verify_ids=True)
    
    with pytest.raises(AssertionError):
        next(processor.fastq_gen(fh1, fh2))

def test_n_records_limit():
    # Test n_records limit functionality
    data = [
        ("read1", "ACGT", "IIII"),
        ("read2", "ACGT", "IIII"),
        ("read3", "ACGT", "IIII")
    ]
    fh = create_fastq_file(data)
    
    processor = fastqProcessor(n_records=2)
    results = list(processor.fastq_gen(fh))
    assert len(results) == 2

def test_fastq_gen_wrapper():
    # Test the fastq_gen wrapper function
    data = [("read1", "ACGT", "IIII")]
    fh = create_fastq_file(data)
    
    result = next(fastq_gen(fh))
    assert result[0] == "@read1"
    assert result[1] == "ACGT"

def test_invalid_phred_type():
    # Test invalid phred type handling
    with pytest.raises(ValueError):
        fastqProcessor(phred_type="invalid")

def test_empty_fastq():
    # Test handling of empty FASTQ file
    fh = StringIO("")
    processor = fastqProcessor()
    results = list(processor.fastq_gen(fh))
    assert len(results) == 0

def test_malformed_fastq():
    # Test handling of malformed FASTQ
    malformed = StringIO("@read1\nACGT\n+\n")  # Missing quality scores
    processor = fastqProcessor()
    with pytest.raises(StopIteration):
        next(processor.fastq_gen(malformed))

def test_extract_control_id():
    # Test control ID extraction
    processor = fastqProcessor()
    assert processor.extract_control_id("@read1 description") == "@read1"
    assert processor.extract_control_id("@read1") == "@read1"
    assert processor.extract_control_id("") is None

def test_rc():

    seq = "ATGCNCGTA"
    rc_seq = "TACGNGCAT"

    assert RC(seq) == rc_seq
