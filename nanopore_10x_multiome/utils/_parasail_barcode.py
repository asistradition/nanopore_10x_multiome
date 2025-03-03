import parasail

PARASAIL_MATRIX = parasail.matrix_create("ACGTN", 5, -1)
for i in [4, 10, 16, 22, 24, 25, 26, 27]:
    PARASAIL_MATRIX.pointer[0].matrix[i] = 1

def get_barcode_parasail(
    seq,
    qual,
    compiled_regex,
    comparison_sequence,
    bc_len
):
    """
    Find an barcode by

    1. fuzzy regex using a precompiled regular expression
    2. smith-waterman local sequence alignment against regex match

    :param seq: Sequence to search
    :type seq: str
    :param qual: Quality string to slice
    :type qual: str
    :param compiled_regex: Precompiled fuzzy regular expression
    :type compiled_regex: regex.Regex
    :comparison_sequence: Sequence for S-W alignment, with barcode
        masked using Ns
    :comparison_sequecne: str
    :param bc_len: Barcode length
    :type bc_len: int

    :return: Tuple of (
        Barcode sequence string,
        Barcode quality string,
        Barcode start position
    )
    :rtype: (str, str, int)
    """
    
    _bc = compiled_regex.search(seq)

    if _bc is None:
        return None, None, None

    _span = _bc.span()
    seq = seq[_span[0]:_span[1]]
    qual = qual[_span[0]:_span[1]]

    result = parasail.sw_trace(
        s1=seq,
        s2=comparison_sequence,
        open=2,
        extend=4,
        matrix=PARASAIL_MATRIX,
    )
    _position = result.traceback.ref.find('N')

    if _position == -1:
        return None, None, None
    
    _barcode = result.traceback.query[
        _position:_position+bc_len
    ].replace('-', 'N')

    if len(_barcode) < bc_len:
        return None, None, None

    _loc = seq.find(_barcode)
    _bcq = qual[_loc:_loc + len(_barcode)]

    return _barcode, _bcq, _span[0] + _position
