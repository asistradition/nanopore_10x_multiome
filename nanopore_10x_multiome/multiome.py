import itertools

import numpy as np
import joblib
import tqdm

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
    n_jobs=None,
    progress_bar=False,
    write_only_valid_barcodes=False,
    keep_runoff_fragments=False,
    verbose=0
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
            n_records=n_records,
            write_only_valid_barcodes=write_only_valid_barcodes,
            keep_runoff_fragments=keep_runoff_fragments
        )
    
    if atac_technical_file_name is None:
        atac_technical_file_name = itertools.repeat(None)

    (
        gex_barcodes,
        atac_barcodes,
        gex_correction_table,
        atac_correction_table,
        atac_gex_translation_table
    ) = _load_missing_multiome_barcode_info()

    if progress_bar:
        iterer = tqdm.tqdm
    else:
        def iterer(x, **kwargs):
            return iter(x)

    return np.stack([
        r
        for r in iterer(joblib.Parallel(
            n_jobs=n_jobs,
            batch_size=1,
            verbose=verbose
        )(
            joblib.delayed(_split_multiome_preamp_fastq)(
                *files,
                n_records=n_records,
                gex_barcodes=gex_barcodes,
                atac_barcodes=atac_barcodes,
                gex_correction_table=gex_correction_table,
                atac_correction_table=atac_correction_table,
                atac_gex_translation_table=atac_gex_translation_table,
                write_only_valid_barcodes=write_only_valid_barcodes
        )
            for files in zip(
                in_file_name,
                atac_file_name,
                gex_file_name,
                other_file_name,
                atac_technical_file_name
            )
        ),
        total=len(in_file_name)
    )])


def _split_multiome_preamp_fastq(
    in_file_name,
    atac_file_name,
    gex_file_name,
    other_file_name,
    atac_technical_file_name=None,
    n_records=None,
    gex_barcodes=None,
    atac_barcodes=None,
    gex_correction_table=None,
    atac_correction_table=None,
    atac_gex_translation_table=None,
    write_only_valid_barcodes=False
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

    (
        gex_barcodes,
        atac_barcodes,
        gex_correction_table,
        atac_correction_table,
        atac_gex_translation_table
    ) = _load_missing_multiome_barcode_info(
        gex_barcodes,
        atac_barcodes,
        gex_correction_table,
        atac_correction_table,
        atac_gex_translation_table
    )

    processor = fastqProcessor(
        verify_ids=False,
        phred_type='raw',
        n_records=n_records
    )

    with (
        file_opener(in_file_name, mode='r') as fh,
        file_opener(atac_file_name, mode='w') as atac_fh,
        file_opener(gex_file_name, mode='w') as gex_fh,
        file_opener(other_file_name, mode='w') as other_fh
    ):
        
        atac_writer = get_file_writer(atac_file_name)
        gex_writer = get_file_writer(gex_file_name)
        other_writer = get_file_writer(other_file_name)

        if atac_technical_file_name is not None:
            atac_tech_fh = file_opener(atac_technical_file_name, mode='w')
            atac_tech_writer = get_file_writer(atac_technical_file_name)
        else:
            atac_tech_fh = None
            atac_tech_writer = None

        try:
            for x in processor.fastq_gen(fh):

                c, s, q = x[0]

                # Find the anchor and extract the
                # cell-specific barcode
                _bc, _bc_qual, tn5_locs = get_atac_anchors(s, q)

                if _bc is not None:

                    _tags, _valid = process_atac_tags(
                        _bc,
                        _bc_qual,
                        atac_correction_table,
                        atac_gex_translation_table
                    )

                    if write_only_valid_barcodes and not _valid:
                        continue

                    atac_writer(
                        atac_fh,
                        c,
                        s[tn5_locs[1]:tn5_locs[2]],
                        q[tn5_locs[1]:tn5_locs[2]],
                        **_tags
                    )

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

                # If the ATAC search failed, check for GEX anchors
                _bc, _umi, gex_locs = get_gex_anchors(s, q)

                if _bc is not None:
                    _tags, _valid = process_gex_tags(
                        _bc[0],
                        _bc[1],
                        _umi[0],
                        _umi[1],
                        gex_correction_table
                    )

                    if write_only_valid_barcodes and not _valid:
                        continue

                    gex_writer(
                        gex_fh,
                        c,
                        s[gex_locs[0]:gex_locs[1]],
                        q[gex_locs[0]:gex_locs[1]],
                        **_tags
                    )
                    result_counts[1] += 1
                    continue

                # Put the remains into other
                other_writer(
                    other_fh,
                    c,
                    s,
                    q
                )
                result_counts[2] += 1

        finally:
            if atac_tech_fh is not None:
                atac_tech_fh.close()

    return result_counts


def _load_missing_multiome_barcode_info(
    gex_barcodes=None,
    atac_barcodes=None,
    gex_correction_table=None,
    atac_correction_table=None,
    atac_gex_translation_table=None
):
    
    if gex_barcodes is None:
        gex_barcodes = load_gex_barcodes()

    if atac_barcodes is None:
        atac_barcodes = load_atac_barcodes()

    if gex_correction_table is None:
        gex_correction_table = barcode_correction_table(gex_barcodes)

    if atac_correction_table is None:
        atac_correction_table = barcode_correction_table(atac_barcodes)

    if atac_gex_translation_table is None:
        atac_gex_translation_table = load_translations(
            atac_barcodes,
            gex_barcodes
        )

    return (
        gex_barcodes,
        atac_barcodes,
        gex_correction_table,
        atac_correction_table,
        atac_gex_translation_table
    )
