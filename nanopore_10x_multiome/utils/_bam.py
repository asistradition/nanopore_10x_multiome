import os
from collections import Counter

import pysam
import tqdm
import pandas as pd

def write_bam_record(
    handle,
    header,
    sequence,
    qual,
    flag=0,
    **tags
):
    
    a = pysam.AlignedSegment()
    a.query_name = header
    a.query_sequence = sequence
    a.flag = flag
    a.query_qualities = pysam.qualitystring_to_array(qual)
    a.tags = [
        (k, v)
        for k, v in tags.items()
        if v is not None
    ]

    handle.write(a)


def bam_summarize_barcodes(
    bam_file,
    pbar=False,
    barcode_tag='CB'
):
    
    # Set up progress bar if requested
    if pbar:
        iterer = tqdm.tqdm
    else:
        def iterer(x, **kwargs):
            return iter(x)

    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        bamlen = bamfile.mapped + bamfile.unmapped

        _data = [
            (
                r.get_tag(barcode_tag),
                len(r.query_sequence),
                r.is_mapped
            )
            for r in iterer(bamfile, total=bamlen)
        ]

    return pd.DataFrame(
        _data,
        columns=['barcode', 'length', 'is_mapped']
    )



def split_bam_by_barcode(
    bam_file,
    lookup_table,
    out_path=None,
    output_files=None,
    out_prefix='',
    barcode_tag='CB',
    pbar=False
):
    """Split a BAM file into multiple files based on barcode assignments.

    :param bam_file: Path to input BAM file
    :type bam_file: str
    :param lookup_table: Dictionary mapping barcodes to their assigned groups/classes
    :type lookup_table: dict
    :param out_path: Directory path for output files if output_files not provided
    :type out_path: str or None
    :param output_files: Dictionary mapping groups to output file paths
    :type output_files: dict or None
    :param out_prefix: Prefix to add to output filenames if using out_path
    :type out_prefix: str
    :param barcode_tag: BAM tag containing the barcode
    :type barcode_tag: str
    :param pbar: Whether to show progress bar
    :type pbar: bool
    :return: Dictionary with number of reads written to each output file
    :rtype: dict
    """
    # Open input BAM file
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    bamlen = bamfile.mapped + bamfile.unmapped

    # Count number of entries for each class in lookup table
    n_classes = Counter([x for x in lookup_table.values()])

    # Generate output filenames if not provided
    if output_files is None:
        output_files = {
            x: os.path.join(out_path, f'{out_prefix}{x}_bam.bam')
            for x in n_classes.keys()
        }
    else:
        # Verify output_files matches lookup_table classes
        _mismatches = set(list(output_files.keys())).symmetric_difference(
            set(list(n_classes.keys()))
        )

        if len(_mismatches) > 0:
            raise ValueError(
                f"Non-overlapping output_files and lookup_table entries: {_mismatches}"
            )
            
    # Open output BAM files
    output_handles = {
        x: pysam.AlignmentFile(output_files[x], "wb", template=bamfile)
        for x in n_classes.keys()
    }

    # Set up progress bar if requested
    if pbar:
        iterer = tqdm.tqdm
    else:
        def iterer(x, **kwargs):
            return iter(x)

    try:
        # Track number of reads written to each file
        n_written = {x: 0 for x in n_classes.keys()}

        # Iterate through BAM file
        for r in iterer(bamfile, total=bamlen):

            if r.is_unmapped:
                continue

            try:
                # Look up which output file this barcode belongs to
                _mapped = lookup_table[r.get_tag(barcode_tag)]
            except KeyError:
                continue

            # Write read to appropriate output file and increment counter
            output_handles[_mapped].write(r)
            n_written[_mapped] = n_written[_mapped] + 1
        
    except:
        raise

    finally:
        # Ensure all files are closed properly
        bamfile.close()
        for x in output_handles.values():
            x.close()

    return n_written
