import itertools

import numpy as np
import joblib

from nanopore_10x_multiome.utils import (
    fastqProcessor,
    write_record
)
from nanopore_10x_multiome.atac import (
    get_atac_anchors,
    process_atac_header
)
from nanopore_10x_multiome.gex import (
    get_gex_anchors,
    process_gex_header
)
from nanopore_10x_multiome.barcodes import (
    load_gex_barcodes,
    load_atac_barcodes,
    load_translations,
    barcode_correction_table
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
    n_records=None,
    n_jobs=None
):
    """
    Split a 10x multiome pre-amp FASTQ file into ATAC, GEX and other reads.

    If lists of files are provided, process all of them in parallel based on
    the number of jobs (n_jobs).

    :param in_file_name: Input FASTQ file path
    :type in_file_name: str
    :param atac_file_name: Output FASTQ file path for ATAC reads
    :type atac_file_name: str 
    :param gex_file_name: Output FASTQ file path for GEX reads
    :type gex_file_name: str
    :param other_file_name: Output FASTQ file path for unclassified reads
    :type other_file_name: str
    :param atac_technical_file_name: Optional output FASTQ file path for ATAC technical sequences
    :type atac_technical_file_name: str or None
    :param n_records: Number of records to process (None for all)
    :type n_records: int or None
    :return: Array of counts [ATAC reads, GEX reads, other reads]
    :rtype: numpy.ndarray
    """
    
    if n_jobs is None or not isinstance(in_file_name, (tuple, list)):
        return _split_multiome_preamp_fastq(
            in_file_name,
            atac_file_name,
            gex_file_name,
            other_file_name,
            atac_technical_file_name,
            n_records
        )
    
    if atac_technical_file_name is None:
        atac_technical_file_name = itertools.repeat(None)

    return np.stack([
        r
        for r in joblib.Parallel(n_jobs=n_jobs)(
            joblib.delayed(_split_multiome_preamp_fastq)(
                *files,
                n_records=n_records
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
    n_records=None
):
    """
    Split a 10x multiome pre-amp FASTQ file into ATAC, GEX and other reads.

    :param in_file_name: Input FASTQ file path
    :type in_file_name: str
    :param atac_file_name: Output FASTQ file path for ATAC reads
    :type atac_file_name: str 
    :param gex_file_name: Output FASTQ file path for GEX reads
    :type gex_file_name: str
    :param other_file_name: Output FASTQ file path for unclassified reads
    :type other_file_name: str
    :param atac_technical_file_name: Optional output FASTQ file path for ATAC technical sequences
    :type atac_technical_file_name: str or None
    :param n_records: Number of records to process (None for all)
    :type n_records: int or None
    :return: Array of counts [ATAC reads, GEX reads, other reads]
    :rtype: numpy.ndarray
    """
    
    result_counts = np.zeros(3, dtype=int)

    gex_barcodes = load_gex_barcodes()
    atac_barcodes = load_atac_barcodes()

    gex_correction_table = barcode_correction_table(gex_barcodes)
    atac_correction_table = barcode_correction_table(atac_barcodes)
    atac_gex_translation_table = load_translations(
        atac_barcodes,
        gex_barcodes
    )

    processor = fastqProcessor(
        verify_ids=False,
        phred_type='raw',
        n_records=n_records
    )

    with (
        open(in_file_name, mode='r') as fh,
        open(atac_file_name, mode='w') as atac_fh,
        open(gex_file_name, mode='w') as gex_fh,
        open(other_file_name, mode='w') as other_fh
    ):
        
        if atac_technical_file_name is not None:
            atac_tech_fh = open(atac_technical_file_name, mode='w')
        else:
            atac_tech_fh = None

        try:
            for x in processor.fastq_gen(fh):

                c, s, q = x[0]

                # Find the anchor and extract the
                # cell-specific barcode
                _bc, _bc_qual, tn5_locs = get_atac_anchors(s, q)

                if _bc is not None:

                    _header = process_atac_header(
                        c,
                        _bc,
                        _bc_qual,
                        atac_correction_table,
                        atac_gex_translation_table
                    )

                    write_record(
                        _header,
                        s[tn5_locs[1]:tn5_locs[2]],
                        q[tn5_locs[1]:tn5_locs[2]],
                        atac_fh
                    )

                    if atac_tech_fh is not None:
                        write_record(
                            _header,
                            s[:tn5_locs[1]] + '----' + s[tn5_locs[2]:],
                            q[:tn5_locs[1]] + '----' + q[tn5_locs[2]:],
                            atac_tech_fh
                        )

                    result_counts[0] += 1
                    continue

                # Check for GEX anchors
                _bc, _umi, gex_locs = get_gex_anchors(s)

                if _bc is not None:
                    write_record(
                        process_gex_header(
                            c,
                            _bc[0],
                            _bc[1],
                            _umi[0],
                            _umi[1],
                            gex_correction_table
                        ),
                        s[gex_locs[0]:gex_locs[1]],
                        q[gex_locs[0]:gex_locs[1]],
                        gex_fh
                    )
                    result_counts[1] += 1
                    continue

                # Put the remains into other
                write_record(
                    c,
                    s,
                    q,
                    other_fh
                )
                result_counts[2] += 1

        finally:
            if atac_tech_fh is not None:
                atac_tech_fh.close()

    return result_counts
