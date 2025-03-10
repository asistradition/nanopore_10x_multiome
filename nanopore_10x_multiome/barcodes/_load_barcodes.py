import os
from pathlib import Path
import pandas as pd


def load_atac_barcodes(test=False):
    return pd.read_csv(
        os.path.join(Path(__file__).parent.absolute(), '737K-arc-v1_atac.txt.gz'),
        sep='\t',
        header=None,
        nrows=100 if test else None
    ).values.ravel()


def load_gex_barcodes(test=False):
    return pd.read_csv(
        os.path.join(Path(__file__).parent.absolute(), '737K-arc-v1_rna.txt.gz'),
        sep='\t',
        header=None,
        nrows=100 if test else None
    ).values.ravel()


def load_translations(
    atac_barcodes=None,
    gex_barcodes=None
):
    

    return dict(zip(
        atac_barcodes if atac_barcodes is not None else load_atac_barcodes(),
        gex_barcodes if gex_barcodes is not None else load_gex_barcodes()
    ))


def translate_barcode(
    barcode,
    translation_table
):
    
    if barcode is None:
        return None

    try:
        return translation_table[barcode]
    except KeyError:
        return barcode
