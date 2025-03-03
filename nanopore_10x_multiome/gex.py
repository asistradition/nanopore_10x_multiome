import regex

from nanopore_10x_multiome.utils import RC, REV, get_barcode_parasail

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
    min_len=25
):

    n = len(seq)
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
            bc_len=28
        )
        _seq_loc = 0, max(_bc_pos - 28, 0)
    else:
        _seq_loc = min(_bc_pos - 28, n), n

    if _bc is None:
        return None, None, None
    
    if (_seq_loc[1] - _seq_loc[0]) < min_len:
        return None, None, None
    
    return _bc, _bc_qual, _seq_loc
