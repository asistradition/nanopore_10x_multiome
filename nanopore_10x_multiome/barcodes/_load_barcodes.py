import os
from pathlib import Path
import pandas as pd


def load_atac_barcodes():
    return pd.read_csv(
        os.path.join(Path(__file__).parent.absolute(), '737K-arc-v1_atac.txt.gz'),
        sep='\t',
        header=None
    ).values.ravel()


def load_gex_barcodes():
    return pd.read_csv(
        os.path.join(Path(__file__).parent.absolute(), '737K-arc-v1_rna.txt.gz'),
        sep='\t',
        header=None
    ).values.ravel()


def load_translations():
    return dict(zip(load_atac_barcodes(), load_gex_barcodes()))

