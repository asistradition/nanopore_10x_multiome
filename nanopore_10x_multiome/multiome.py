import itertools

import numpy as np
import joblib

from nanopore_10x_multiome.utils import (
    fastqProcessor,
    get_file_writer,
    file_opener
)
from nanopore_10x_multiome.atac import (
    get_atac_anchors,
    process_atac_tags
)
from nanopore_10x_multiome.gex import (
    get_gex_anchors,
    process_gex_tags
)
from nanopore_10x_multiome.barcodes import (
    load_missing_multiome_barcode_info,
    BarcodeHolder
)

###############################################################################
# 10x multiome ATAC tags
# CB - Corrected barcode and translated to match GEX sequence
# CR - Raw barcode off instrument
# CY - Barcode quality scores from instrument
#
# 10x multiome GEX tags
# CB - Corrected barcode
# CR - Raw barcode off instrument
# CY - Barcode quality scores from instrument
# UB - Corrected UMI
# CR - Raw UMI off instrument
# CY - UMI quality scores from instrument
#
# ATAC and GEX barcodes are not the same - use translation table!
###############################################################################


def split_multiome_preamp_fastq(
    in_file_name,
    atac_file_name,
    gex_file_name,
    other_file_name,
    atac_technical_file_name=None,
    n_jobs=None,
    write_only_valid_barcodes=False,
    keep_runoff_fragments=False,
    verbose=0
):
    """
    Split multiome pre-amplification FASTQ file(s) into ATAC, GEX and other reads.

    :param in_file_name: Input FASTQ file path
    :type in_file_name: str
    :param atac_file_name: Output FASTQ file path for ATAC reads
    :type atac_file_name: str
    :param gex_file_name: Output FASTQ file path for GEX reads
    :type gex_file_name: str
    :param other_file_name: Output FASTQ file path for unidentified reads
        (genomic, too many errors near barcode, etc)
    :type other_file_name: str
    :param atac_technical_file_name: Optional output file for ATAC technical sequences
    :type atac_technical_file_name: str or None
    :param n_jobs: Number of parallel processes for joblib, defaults to None
    :type n_jobs: int or None
    :param keep_runoff_fragments: Keep ATAC fragments where the barcode end is intact,
        but no Tn5 site is located on the other end. Defaults to False.
    :type keep_runoff_fragments: bool
    :param write_only_valid_barcodes: Only write reads with valid barcodes
    :type write_only_valid_barcodes: bool
    :param verbose: Verbose parameter for joblib.Parallel
    :type verbose: int

    :return: Array of counts [ATAC reads, GEX reads, other reads]
    :rtype: numpy.ndarray
    """

    
    if not isinstance(in_file_name, (tuple, list)):
        return _split_multiome_preamp_fastq(
            in_file_name,
            atac_file_name,
            gex_file_name,
            other_file_name,
            atac_technical_file_name,
            write_only_valid_barcodes=write_only_valid_barcodes,
            keep_runoff_fragments=keep_runoff_fragments,
        )
    
    load_missing_multiome_barcode_info(pbar=verbose > 0)

    return np.stack([
        r
        for r in joblib.Parallel(
            n_jobs=n_jobs,
            batch_size=1,
            verbose=verbose
        )(
            joblib.delayed(_split_multiome_preamp_fastq)(
                *files,
                write_only_valid_barcodes=write_only_valid_barcodes,
                keep_runoff_fragments=keep_runoff_fragments
            )
            for files in zip(
                in_file_name,
                atac_file_name,
                gex_file_name,
                other_file_name,
                atac_technical_file_name
            )
        )
    ])


def _split_multiome_preamp_fastq(
    in_file_name,
    atac_file_name,
    gex_file_name,
    other_file_name,
    atac_technical_file_name=None,
    n_records=None,
    write_only_valid_barcodes=False,
    keep_runoff_fragments=False
):
    """
    Split a multiome pre-amplification FASTQ file into ATAC, GEX and other reads.

    :param in_file_name: Input FASTQ file path
    :type in_file_name: str
    :param atac_file_name: Output FASTQ file path for ATAC reads
    :type atac_file_name: str
    :param gex_file_name: Output FASTQ file path for GEX reads
    :type gex_file_name: str
    :param other_file_name: Output FASTQ file path for unidentified reads
        (genomic, too many errors near barcode, etc)
    :type other_file_name: str
    :param atac_technical_file_name: Optional output file for ATAC technical sequences
    :type atac_technical_file_name: str or None
    :param n_records: Number of records to process (None for all)
    :type n_records: int or None
    :param gex_barcodes: List of valid GEX barcodes
    :type gex_barcodes: list or None
    :param atac_barcodes: List of valid ATAC barcodes
    :type atac_barcodes: list or None
    :param gex_correction_table: Correction table for GEX barcodes
    :type gex_correction_table: dict or None
    :param atac_correction_table: Correction table for ATAC barcodes
    :type atac_correction_table: dict or None
    :param atac_gex_translation_table: Translation table between ATAC and GEX barcodes
    :type atac_gex_translation_table: dict or None
    :param write_only_valid_barcodes: Only write reads with valid barcodes
    :type write_only_valid_barcodes: bool
    :param keep_runoff_fragments: Keep ATAC fragments where the barcode end is intact,
        but no Tn5 site is located on the other end. Defaults to False.

    :return: Array of counts [ATAC reads, GEX reads, other reads]
    :rtype: numpy.ndarray
    """

    # Initialize counters for ATAC, GEX and other reads
    result_counts = np.zeros(3, dtype=int)

    # Load any missing barcode information
    load_missing_multiome_barcode_info(pbar=False)

    # Initialize FASTQ processor
    processor = fastqProcessor(
        verify_ids=False,
        phred_type='raw',
        n_records=n_records
    )

    # Open input and output files
    with (
        file_opener(in_file_name, mode='r') as fh,
        file_opener(atac_file_name, mode='w') as atac_fh,
        file_opener(gex_file_name, mode='w') as gex_fh,
        file_opener(other_file_name, mode='w') as other_fh
    ):
        
        # Get file writers for each output
        atac_writer = get_file_writer(atac_file_name)
        gex_writer = get_file_writer(gex_file_name)
        other_writer = get_file_writer(other_file_name)

        # Handle optional ATAC technical file
        if atac_technical_file_name is not None:
            atac_tech_fh = file_opener(atac_technical_file_name, mode='w')
            atac_tech_writer = get_file_writer(atac_technical_file_name)
        else:
            atac_tech_fh = None
            atac_tech_writer = None

        try:
            # Process each FASTQ record
            for x in processor.fastq_gen(fh):

                c, s, q = x[0]  # header, sequence, quality scores

                # First try to identify as ATAC read
                _bc, _bc_qual, tn5_locs = get_atac_anchors(
                    s,
                    q,
                    keep_runoff_fragments=keep_runoff_fragments
                )

                if _bc is not None:
                    # Process ATAC barcode and check validity
                    _tags, _valid = process_atac_tags(
                        _bc,
                        _bc_qual,
                        BarcodeHolder.atac_correction_table,
                        BarcodeHolder.atac_gex_translation_table
                    )

                    if write_only_valid_barcodes and not _valid:
                        continue

                    # Write ATAC read
                    atac_writer(
                        atac_fh,
                        c,
                        s[tn5_locs[1]:tn5_locs[2]],
                        q[tn5_locs[1]:tn5_locs[2]],
                        **_tags
                    )

                    # Write technical sequence if requested
                    if atac_tech_fh is not None:
                        atac_tech_writer(
                            atac_tech_fh,
                            c,
                            s[:tn5_locs[1]] + '----' + s[tn5_locs[2]:],
                            q[:tn5_locs[1]] + '----' + q[tn5_locs[2]:],
                           **_tags
                        )

                    result_counts[0] += 1
                    continue

                # If not ATAC, try to identify as GEX read
                _bc, _umi, gex_locs = get_gex_anchors(s, q)

                if _bc is not None:
                    # Process GEX barcode and UMI, check validity
                    _tags, _valid = process_gex_tags(
                        _bc[0],
                        _bc[1],
                        _umi[0],
                        _umi[1],
                        BarcodeHolder.gex_correction_table
                    )

                    if write_only_valid_barcodes and not _valid:
                        continue

                    # Write GEX read
                    gex_writer(
                        gex_fh,
                        c,
                        s[gex_locs[0]:gex_locs[1]],
                        q[gex_locs[0]:gex_locs[1]],
                        **_tags
                    )
                    result_counts[1] += 1
                    continue

                # If neither ATAC nor GEX, write to other file
                other_writer(
                    other_fh,
                    c,
                    s,
                    q
                )
                result_counts[2] += 1

        finally:
            # Ensure technical file is closed if it was opened
            if atac_tech_fh is not None:
                atac_tech_fh.close()

    return result_counts
