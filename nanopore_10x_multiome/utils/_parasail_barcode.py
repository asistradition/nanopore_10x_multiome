import parasail

# Substitution matrix where N isn't penalized for any match
PARASAIL_MATRIX = parasail.matrix_create("ACGTN", 5, -1)
for i in [4, 10, 16, 22, 24, 25, 26, 27]:
    PARASAIL_MATRIX.pointer[0].matrix[i] = 4

def get_barcode_parasail(
    seq,
    qual,
    compiled_regex,
    comparison_sequence,
    bc_len,
    split_barcode=None
):
    """
    Find an barcode by:

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
        extend=1,
        matrix=PARASAIL_MATRIX,
    )

    # Get the barcode start position
    _position = result.traceback.ref.find('N')

    if _position == -1:
        return None, None, None
    
    # Fix gaps right next to barcode as they're probably insertions
    if (_position > 0) and (result.traceback.ref[_position - 1] == '-'):
        _position = _position - 1
        bc_len = bc_len + 1
    elif (_position < len(result.traceback.ref)) and (result.traceback.ref[_position + 1] == '-'):
        bc_len = bc_len + 1

    _barcode = result.traceback.query[
        _position:_position+bc_len
    ].replace('-', '')

    # Find any padding in the query before the barcode
    # so the output position can be fixed
    _extra_offset = result.traceback.query.count(
        '-', 0, _position
    )

    # Allow at most one deletion
    if len(_barcode) < (bc_len - 1):
        return None, None, None

    _loc = seq.find(_barcode)
    _bcq = qual[_loc:_loc + len(_barcode)]

    return _barcode, _bcq, _span[0] + _position - _extra_offset
