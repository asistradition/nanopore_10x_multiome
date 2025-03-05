import pytest
from nanopore_10x_multiome.atac import (
    get_atac_anchors,
    TENX_ATAC_ADAPTER,
    tenx_re
)
from nanopore_10x_multiome.utils import (
    RC,
    get_barcode_parasail
)
from nanopore_10x_multiome.utils.test import (
    create_sequence,
    create_qual,
    BASE_SEQ
)
ATAC_BARCODE = 'AATGATACGGCGACCACCGAGATCTACACATGCATGCATGCATGCCGCGTCTGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
ABN = len(ATAC_BARCODE)

TN5_SEQ = 'CTGTCTCTTATACACATCT'
TN5N = len(TN5_SEQ)

BSN = len(BASE_SEQ)


def create_atac_sequence(
    barcode_seq,
    barcode_seq_loc,
    tn_seq,
    tn_seq_loc,
    rc=False
):
    
    return create_sequence(
        barcode_seq,
        barcode_seq_loc,
        tn_seq,
        tn_seq_loc,
        rc=rc
    )


def test_create_seq():
    seq = ATAC_BARCODE + BASE_SEQ + TN5_SEQ

    seq2 = create_sequence(
        ATAC_BARCODE,
        0,
        TN5_SEQ,
        1311
    )

    assert seq == seq2

def test_perfect_atac_parasail():
    seq = create_atac_sequence(
        ATAC_BARCODE,
        0,
        TN5_SEQ,
        1311
    )

    psbc, psqual, psloc = get_barcode_parasail(
        seq,
        create_qual(len(seq), 29, 16),
        tenx_re,
        TENX_ATAC_ADAPTER,
        16
    )

    assert len(psbc) == 16
    assert psbc == "ATGCATGCATGCATGC"
    assert psqual == "A" * 16
    assert psloc == 29


def test_extra_bases_atac_parasail():
    
    seq = create_atac_sequence(
        'AATGATACGNGCGACCACCNGAGATCTACACATGCATGCATGCATGCCGCGTCTGTCGTCGGCAGNCGTCAGATGTGTATANAGAGACAG',
        0,
        'CTGTCTCTTNATACACNATCT',
        1311
    )

    psbc, psqual, psloc = get_barcode_parasail(
        seq,
        create_qual(len(seq), 31, 16),
        tenx_re,
        TENX_ATAC_ADAPTER,
        16
    )

    assert len(psbc) == 16
    assert psbc == "ATGCATGCATGCATGC"
    assert psqual == "A" * 16
    assert psloc == 31


def test_missing_bases_atac_parasail():
    seq = create_atac_sequence(
        'ATACGGCGACCACCGATCTACACATGCATGCATGCATGCCGCGTCTTCGTCGGCGCGTCAGATGTGTATAAGAGACAG',
        0,
        'CTGTTCTTATACACTCT',
        1311
    )

    psbc, psqual, psloc = get_barcode_parasail(
        seq,
        create_qual(len(seq), 23, 16),
        tenx_re,
        TENX_ATAC_ADAPTER,
        16
    )

    assert len(psbc) == 16
    assert psbc == "ATGCATGCATGCATGC"
    assert psqual == "A" * 16
    assert psloc == 23


def test_missing_bases_moved_atac_parasail():
    seq = create_atac_sequence(
        'ATACGGCGACCACCGATCTACACATGCATGCATGCATGCCGCGTCTTCGTCGGCGCGTCAGATGTGTATAAGAGACAG',
        1310,
        'CTGTTCTTATACACTCT',
        1311
    )

    psbc, psqual, psloc = get_barcode_parasail(
        seq,
        create_qual(len(seq), 1310 + 23, 16),
        tenx_re,
        TENX_ATAC_ADAPTER,
        16
    )

    assert len(psbc) == 16
    assert psbc == "ATGCATGCATGCATGC"
    assert psqual == "A" * 16
    assert psloc == 1310 + 23



def test_perfect_atac_full():
    seq = create_atac_sequence(
        ATAC_BARCODE,
        0,
        TN5_SEQ,
        1311
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 29, 16)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == BASE_SEQ

    assert len(bc) == 16
    assert bc == "ATGCATGCATGCATGC"
    assert bcq == "A" * 16
    assert locs == [67, 67 + TN5N, ABN + BSN, ABN + BSN + TN5N]


def test_perfect_atac_full_rc():
    seq = create_atac_sequence(
        ATAC_BARCODE,
        0,
        TN5_SEQ,
        1311,
        rc=True
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 29, 16, rev=True)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == RC(BASE_SEQ)

    assert len(bc) == 16
    assert bc == "ATGCATGCATGCATGC"
    assert bcq == "A" * 16
    assert locs == [0, 0 + TN5N, BSN + TN5N, BSN + 2 * TN5N]


def test_extra_bases_atac_full():
    seq = create_atac_sequence(
        'AATGATACGNGCGACCACCNGAGATCTACACATGCATGCATGCATGCCGCGTCTGTCGTCGGCAGNCGTCAGATGTGTATANAGAGACAG',
        0,
        'CTGTCTCTTNATACACNATCT',
        1311
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 31, 16)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == BASE_SEQ

    assert len(bc) == 16
    assert bc == "ATGCATGCATGCATGC"
    assert bcq == "A" * 16
    assert locs == [67 + 3, 67 + 3 + 1 + TN5N, ABN + BSN + 4, ABN + BSN + TN5N + 4 + 2]


def test_extra_bases_atac_full_rc():
    seq = create_atac_sequence(
        'AATGATACGNGCGACCACCNGAGATCTACACATGCATGCATGCATGCCGCGTCTGTCGTCGGCAGNCGTCAGATGTGTATANAGAGACAG',
        0,
        'CTGTCTCTTNATACACNATCT',
        1311,
        rc=True
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 31, 16, rev=True)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == RC(BASE_SEQ)

    assert len(bc) == 16
    assert bc == "ATGCATGCATGCATGC"
    assert bcq == "A" * 16
    assert locs == [0, 0 + TN5N + 2, 2 + BSN + TN5N, BSN + 2 * TN5N + 2 + 1]


def test_missing_bases_atac_full():
    seq = create_atac_sequence(
        'ATACGGCGACCACCGATCTACACATGCATGCATGCATGCCGCGTCTTCGTCGGCGCGTCAGATGTGTATAAGAGACAG',
        0,
        'CTGTCTCTTNATACACNATCT',
        1311
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 23, 16)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == BASE_SEQ

    assert len(bc) == 16
    assert bc == "ATGCATGCATGCATGC"
    assert bcq == "A" * 16
    assert locs == [67 - 8, 67 - 8 + TN5N, ABN + BSN - 8, ABN + BSN + TN5N - 6]


def test_missing_bases_atac_full_rc():
    seq = create_atac_sequence(
        'ATACGGCGACCACCGATCTACACATGCATGCATGCATGCCGCGTCTTCGTCGGCGCGTCAGATGTGTATAAGAGACAG',
        0,
        'CTGTCTCTTNATACACNATCT',
        1311,
        rc=True
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 23, 16, rev=True)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == RC(BASE_SEQ)

    assert len(bc) == 16
    assert bc == "ATGCATGCATGCATGC"
    assert bcq == "A" * 16
    assert locs == [0, 0 + TN5N + 2, 2 + BSN + TN5N, BSN + 2 * TN5N + 2]


def test_near_barcode_errors_atac_full():
    seq = create_atac_sequence(
        'AATGATACGGCGACCACCGAGATCCACATGCATGCATGCATGCCGTCTGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',
        0,
        TN5_SEQ,
        1311
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 27, 16)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == BASE_SEQ

    assert len(bc) == 16
    assert bc == "ATGCATGCATGCATGC"
    assert bcq == "A" * 16
    assert locs == [67 - 4, 67 - 4 + TN5N, ABN + BSN - 4, ABN + BSN + TN5N - 4]


def test_within_barcode_insertion_atac_full():
    seq = create_atac_sequence(
        'AATGATACGGCGACCACCGAGATCTACACATGCATGCGATGCATGCCGCGTCTGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',
        0,
        TN5_SEQ,
        1311
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 29, 17)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == BASE_SEQ

    assert len(bc) == 17
    assert bc == "ATGCATGCGATGCATGC"
    assert bcq == "A" * 17
    assert locs == [67 + 1, 67 + TN5N + 1, ABN + BSN + 1, ABN + BSN + TN5N + 1]


def test_within_barcode_loss_atac_full():
    seq = create_atac_sequence(
        'AATGATACGGCGACCACCGAGATCTACACATGCATGATGCATGCCGCGTCTGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',
        0,
        TN5_SEQ,
        1311
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 29, 15)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == BASE_SEQ

    assert len(bc) == 15
    assert bc == "ATGCATGATGCATGC"
    assert bcq == "A" * 15
    assert locs == [67 - 1, 67 + TN5N - 1, ABN + BSN - 1, ABN + BSN + TN5N - 1]


def test_within_barcode_loss_atac_full_rc():
    seq = create_atac_sequence(
        'AATGATACGGCGACCACCGAGATCTACACATGCATGATGCATGCCGCGTCTGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',
        0,
        TN5_SEQ,
        1311,
        rc=True
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        create_qual(len(seq), 29, 15, rev=True)
    )

    genomic = seq[locs[1]:locs[2]]
    assert genomic == RC(BASE_SEQ)

    assert len(bc) == 15
    assert bc == "ATGCATGATGCATGC"
    assert bcq == "A" * 15
    assert locs == [0, TN5N, TN5N + BSN, BSN + 2 * TN5N]



def test_no_barcode_atac_full():
    seq = create_atac_sequence(
        RC(TN5_SEQ),
        0,
        TN5_SEQ,
        1311
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        "I" * len(seq)
    )

    assert bc is None
    assert bcq is None
    assert locs is None


def test_no_second_tn5_atac_full():
    seq = create_atac_sequence(
        ATAC_BARCODE,
        0,
        None,
        1311
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        "I" * len(seq)
    )

    assert bc is None
    assert bcq is None
    assert locs is None


def test_no_second_tn5_runoff_atac_full():
    seq = create_atac_sequence(
        ATAC_BARCODE,
        0,
        None,
        1311
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        "I" * len(seq),
        keep_runoff_fragments=True
    )

    assert bc == "ATGCATGCATGCATGC"
    assert bcq == "IIIIIIIIIIIIIIII"
    assert locs == [0, ABN, ABN + BSN, ABN + BSN]


def test_no_second_tn5_runoff_atac_full_rc():
    seq = create_atac_sequence(
        ATAC_BARCODE,
        0,
        None,
        1311,
        rc=True
    )

    bc, bcq, locs = get_atac_anchors(
        seq,
        "I" * len(seq),
        keep_runoff_fragments=True
    )

    assert bc == "ATGCATGCATGCATGC"
    assert bcq == "IIIIIIIIIIIIIIII"
    assert locs == [0, 0, BSN, ABN + BSN]
