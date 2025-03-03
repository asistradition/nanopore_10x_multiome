import regex

from nanopore_10x_multiome.utils import RC, REV, get_barcode_parasail
from nanopore_10x_multiome.barcodes import translate_barcode, correct_barcode

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


TENX_ATAC_ADAPTER = 'AATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNNNNNNNCGCGTCTGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'


tenx_re = regex.compile(
    '(AGATCTACAC){e<=2}([ATGCN]{16})(CGCGTCTGTCGTCGGCAGCG){e<=3}',
    regex.IGNORECASE
)

tn5_re = regex.compile(
    '(AGATGTGTATAAGAGACAG){e<=3}',
    regex.IGNORECASE
)
tn5_rev_re = regex.compile(
    '(CTGTCTCTTATACACATCT){e<=3}',
    regex.IGNORECASE
)


def get_atac_anchors(
    seq,
    qual,
    keep_runoff_fragments=False,
    min_len=10
):
    """
    Search for ATAC technical sequences from a sequence string

    :param seq: Sequence to search
    :type seq: str
    :param qual: Quality string to slice
    :type qual: str
    :param keep_runoff_fragments: Keep fragments that have a barcode and single
        Tn5 insertion, if there is no matching Tn5 on the other side, defaults to False
    :type keep_runoff_fragments: bool, optional
    :param min_len: Minimum genomic insertion to retain, defaults to 10
    :type min_len: int, optional

    :return: Tuple of (
        barcode sequence,
        barcode quality,
        tn5 insert locations as start, stop, start, stop
    )
    :rtype: (str, str, (int, int, int, int))
    """

    n = len(seq)
    seq = seq.upper()

    # Find 10x ATAC Barcode on the forward strand
    _bc, _bc_qual, _bc_pos = get_atac_barcode_parasail(seq, qual)
    _fwd = True
    
    # If not on the forward strand, look on the reverse strand
    if _bc is None:
        _bc, _bc_qual, _bc_pos = get_atac_barcode_parasail(RC(seq), REV(qual))
        _fwd = False

    # If no atac anchors, return Nones
    if _bc is None:
        return None, None, None

    # Count number of tn5 MEs
    tn5_searches = [
        y
        for y in tn5_re.finditer(seq)
    ] + [
        y
        for y in tn5_rev_re.finditer(seq)
    ]

    # IF there's two Tn5 insertions, find the spot between them
    if len(tn5_searches) == 2:
        tn5_locs = sorted(tn5_searches[0].span() + tn5_searches[1].span())

        # Check for overlapping/no genomic Tn5 insertions
        if (
            ((tn5_locs[2] - tn5_locs[1]) < min_len) or
            ((tn5_locs[3] - tn5_locs[0]) < (38 + min_len))
        ):
            return None, None, None

    # IF there's one Tn5 insertion and the runoff flag is set,
    # check that the barcode is on the correct side of the Tn5
    # and then go to the end of the sequence
    elif keep_runoff_fragments and (len(tn5_searches) == 1):

        _single_tn5 = sorted(tn5_searches[0].span())

        if _fwd and _bc_pos < _single_tn5[1]:
            tn5_locs = [0, _single_tn5[1]] + [n, n]

        elif _bc_pos > _single_tn5[0]:
            tn5_locs = [0, 0] + [_single_tn5[0], n]

        # Check for overlapping/no genomic Tn5 insertions
        if (tn5_locs[2] - tn5_locs[1]) < min_len:
            return None, None, None

    # IF there's 3+ or 0 Tn5 insertions return nothing
    else:
        return None, None, None

    return _bc, _bc_qual, tn5_locs


def get_atac_barcode_parasail(seq, qual):
    """
    Find an ATAC barcode by

    1. fuzzy regex
    2. smith-waterman local sequence alignment against regex match

    :param seq: Sequence to search
    :type seq: str
    :param qual: Quality string to slice
    :type qual: str

    :return: Tuple of (
        Barcode sequence string,
        Barcode quality string,
        Barcode start position
    )
    :rtype: (str, str, int)
    """
    
    return get_barcode_parasail(
        seq,
        qual,
        tenx_re,
        TENX_ATAC_ADAPTER,
        16
    )

def process_atac_header(
    header,
    barcode,
    barcode_quality,
    atac_correction_table,
    atac_gex_translation_table
):
    
    corrected_barcode = translate_barcode(
        correct_barcode(
            barcode,
            barcode_quality,
            atac_correction_table
        ),
        atac_gex_translation_table
    )

    if corrected_barcode is not None:
        tags = f"CB={corrected_barcode} CR={barcode} CY={barcode_quality}"
    else:
        tags = f"CR={barcode} CY={barcode_quality}"

    return f"{header} {tags}"
