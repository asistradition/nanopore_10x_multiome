import numpy as np

from nanopore_10x_multiome.utils import convert_qual_illumina

def barcode_correction_table(barcodes):

    def _oneoffs(x):

        return [
            x[:s] + c + x[s+1:]
            for s in range(len(x))
            for c in ["A", "T", "G", "C", "N"]
        ]

    return {
        mm: b
        for b in barcodes
        for mm in _oneoffs(b)
    }


def correct_barcode(
    barcode,
    qual,
    correction_lookup_table,
    max_dist=1,
    valid_barcodes=None,
    valid_barcodes_char_table=None,
    min_weight_dist=None
):
    
    """
    Assign barcodes to the closest valid barcode.

    :param barcode: Barcode sequence
    :type barcode: str
    :param qual: Barcode quality scores
    :type qual: str
    :param correction_lookup_table: Correction lookup table
    :type correction_lookup_table: dict
    :param max_dist: Maximum distance for assignment, if this is greater than 1
        the other kwargs in this function must be provided. Defaults to 1.
    :type max_dist: int, optional
    :param valid_barcodes: Valid barcodes
    :type valid_barcodes: list[str]
    :param valid_barcodes_char_table: Valid barcodes as character table
    :type valid_barcodes_char_table: np.ndarray
    :param correction_lookup_table: Correction lookup table
    :type correction_lookup_table: dict
    :param min_weight_dist: Minimum weight distance for assignment,
        defaults to 0.5
    :type min_weight_dist: float, optional

    :return: Assigned barcode or None
    :rtype: str or None
    """

    try:
        return correction_lookup_table[barcode]
    except KeyError:
        pass

    if max_dist <= 1:
        return None
    
    if valid_barcodes is None:
        raise RuntimeError(
            'If max_dist > 1, pass valid_barcodes, '
            'valid_barcodes_char_table, and min_weight_dist'
        )

    weights = (convert_qual_illumina(qual) - 15) / 15
    weights = np.maximum(weights, 0) + 1

    # Encode quality scores
    quality_scores = np.array([ord(x) for x in barcode])
    
    # Hamming distance (+1 for each mismatch)
    distance = valid_barcodes_char_table != quality_scores[None, :]

    # Weight by PHRED score
    # So high scoring mismatches are more distant than
    # low scoring mismatches
    wdistance = np.sum(distance * weights[None, :], axis=1)
    distance = np.sum(distance, axis=1)

    sort_order = np.argsort(wdistance)[0:4]

    if distance[sort_order[0]] > max_dist:
        correction_lookup_table[barcode] = None
        return None

    if wdistance[sort_order[0]] < (wdistance[sort_order[1]] - min_weight_dist):
        return valid_barcodes[sort_order[0]]

    else:
        return None
