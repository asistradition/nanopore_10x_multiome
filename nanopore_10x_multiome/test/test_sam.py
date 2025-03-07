"""
Generated with cursor/claude-3.5-sonnet and then fixed to actually work
"""

import pytest
import os
from nanopore_10x_multiome.utils._sam import sam_comment_to_tag
from tempfile import NamedTemporaryFile

def test_sam_comment_to_tag_basic():
    # Create temporary input and output files
    with NamedTemporaryFile(mode='w', delete=False) as in_file, \
         NamedTemporaryFile(mode='w', delete=False) as out_file:
        
        in_file.write("@HD\tVN:1.0\n")  # Header line
        in_file.write("read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\t####\tCB=ACGT\n")  # Normal line
        in_file.flush()
        
        sam_comment_to_tag(in_file.name, out_file.name)
        
        with open(out_file.name, 'r') as f:
            lines = f.readlines()
            
        assert len(lines) == 2
        assert lines[0].strip() == "@HD\tVN:1.0"
        assert "CB:Z:ACGT" in lines[1]
        
    # Cleanup
    os.unlink(in_file.name)
    os.unlink(out_file.name)

def test_sam_comment_to_tag_empty_file():
    with NamedTemporaryFile(mode='w', delete=False) as in_file, \
         NamedTemporaryFile(mode='w', delete=False) as out_file:
        
        in_file.write("")
        in_file.flush()
        
        sam_comment_to_tag(in_file.name, out_file.name)
        
        with open(out_file.name, 'r') as f:
            lines = f.readlines()
            
        assert len(lines) == 0
        
    os.unlink(in_file.name)
    os.unlink(out_file.name)

def test_sam_comment_to_tag_multiple_tags():
    with NamedTemporaryFile(mode='w', delete=False) as in_file, \
         NamedTemporaryFile(mode='w', delete=False) as out_file:
        
        in_file.write("read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\t####\tCB=ACGT UB=TGCA\n")
        in_file.flush()
        
        sam_comment_to_tag(in_file.name, out_file.name)
        
        with open(out_file.name, 'r') as f:
            line = f.readline().strip()
            
        assert "CB:Z:ACGT" in line
        assert "UB:Z:TGCA" in line
        
    os.unlink(in_file.name)
    os.unlink(out_file.name)

def test_sam_comment_to_tag_custom_prefixes():
    with NamedTemporaryFile(mode='w', delete=False) as in_file, \
         NamedTemporaryFile(mode='w', delete=False) as out_file:
        
        in_file.write("read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\t####\tXX=ACGT\n")
        in_file.flush()
        
        sam_comment_to_tag(
            in_file.name, 
            out_file.name,
            comment_prefix='XX=',
            tag_prefix='XX:Z:'
        )
        
        with open(out_file.name, 'r') as f:
            line = f.readline().strip()
            
        assert "XX:Z:ACGT" in line
        
    os.unlink(in_file.name)
    os.unlink(out_file.name)

def test_sam_comment_to_tag_no_tags():
    with NamedTemporaryFile(mode='w', delete=False) as in_file, \
         NamedTemporaryFile(mode='w', delete=False) as out_file:
        
        in_file.write("read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\t####\tNO_TAGS_HERE\n")
        in_file.flush()
        
        sam_comment_to_tag(in_file.name, out_file.name)
        
        with open(out_file.name, 'r') as f:
            line = f.readline().strip()
            
        assert line == "read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\t####\tNO_TAGS_HERE"
        
    os.unlink(in_file.name)
    os.unlink(out_file.name)

def test_sam_comment_to_tag_invalid_input():
    with pytest.raises(AssertionError):
        sam_comment_to_tag(
            "dummy.sam",
            "output.sam",
            comment_prefix=['CB='],
            tag_prefix=['CB:Z:', 'extra']  # Mismatched lengths
        )