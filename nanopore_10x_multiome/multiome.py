import numpy as np

from nanopore_10x_multiome.utils import fastqProcessor
from nanopore_10x_multiome.atac import get_atac_anchors
from nanopore_10x_multiome.gex import get_gex_anchors

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
    n_records=None
):

    result_counts = np.zeros(3, dtype=int)

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
                    print_record(
                        c + " barcode=" + _bc + " bcqual=" + _bc_qual,
                        s[tn5_locs[1]:tn5_locs[2]],
                        q[tn5_locs[1]:tn5_locs[2]],
                        atac_fh
                    )

                    if atac_tech_fh is not None:
                        print_record(
                            c + " barcode=" + _bc + " bcqual=" + _bc_qual,
                            s[:tn5_locs[1]] + '----' + s[tn5_locs[2]:],
                            q[:tn5_locs[1]] + '----' + q[tn5_locs[2]:],
                            atac_tech_fh
                        )

                    result_counts[0] += 1
                    continue

                # Check for GEX anchors
                _bc_umi, gex_locs = get_gex_anchors(s)

                if _bc_umi is not None:
                    print_record(
                        c + " barcode=" + _bc_umi[0] + " umi=" + _bc_umi[1],
                        s[gex_locs[0]:gex_locs[1]],
                        q[gex_locs[0]:gex_locs[1]],
                        gex_fh
                    )
                    result_counts[1] += 1
                    continue

                # Put the remains into other
                print_record(
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


def print_record(
    header,
    seq,
    qual,
    out_fh
):

    print(header, file=out_fh)
    print(seq, file=out_fh)
    print("+", file=out_fh)
    print(qual, file=out_fh)
