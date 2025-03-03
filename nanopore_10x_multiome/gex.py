import regex

from nanopore_10x_multiome.utils import RC, REV, get_barcode_parasail
from nanopore_10x_multiome.barcodes import correct_barcode

###############################################################################
# 10x multiome ATAC sequence
# AATGATACGGCGACCACCGAGATCTACAC-N16-CGCGTCTG-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
# Barcode read RC
#
# 10x multiome GEX sequence
# ACACTCTTTCCCTACACGACGCTCTTCCGATCT-N16-N12
# Barcode (N16) and UMI (N12) read directly
#
# ATAC and GEX barcodes are not the same - use translation table!
###############################################################################

TENX_GEX_ADAPTER = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTT'


gex_re = regex.compile(
    '(CTACACGACGCTCTTCCGA){e<=3}TCT([ATGCN]{28})(TTT)',
    regex.IGNORECASE
)

def get_gex_anchors(
    seq,
    qual,
    min_len=25,
    bc_len=16,
    umi_len=12
):

    n = len(seq)
    bc_umi_len = bc_len + umi_len
    seq = seq.upper()

    # Find 10x GEX Barcode on the forward strand
    _bc, _bc_qual, _bc_pos = get_barcode_parasail(
        seq,
        qual,
        gex_re,
        TENX_GEX_ADAPTER,
        bc_len=28
    )
    
    # If not on the forward strand, look on the reverse strand
    if _bc is None:
        _bc, _bc_qual, _bc_pos = get_barcode_parasail(
            RC(seq),
            REV(qual),
            gex_re,
            TENX_GEX_ADAPTER,
            bc_len=bc_umi_len
        )

        if _bc is None:
            return None, None, None

        _seq_loc = 0, max(_bc_pos - bc_umi_len, 0)
    else:
        _seq_loc = min(_bc_pos - bc_umi_len, n), n

    if (_seq_loc[1] - _seq_loc[0]) < min_len:
        return None, None, None
    
    barcode = _bc[0:bc_len], _bc_qual[0:bc_len]
    umi = _bc[bc_len:], _bc_qual[bc_len:]

    return barcode, umi, _seq_loc

def process_gex_header(
    header,
    barcode,
    barcode_quality,
    umi,
    umi_quality,
    gex_correction_table
):
    
    corrected_barcode = correct_barcode(
        barcode,
        barcode_quality,
        gex_correction_table
    )

    if corrected_barcode is not None:
        bc_tags = f"CB={corrected_barcode} CR={barcode} CY={barcode_quality}"
    else:
        bc_tags = f"CR={barcode} CY={barcode_quality}"

    umi_tags = f"UB={umi} UR={umi} UY={umi_quality}"

    return f"{header} {bc_tags} {umi_tags}"
