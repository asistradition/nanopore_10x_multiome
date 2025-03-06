from ._load_barcodes import (
    load_atac_barcodes,
    load_gex_barcodes,
    load_translations,
    translate_barcode
)

from ._correct_barcodes import (
    correct_barcode,
    barcode_correction_table
)

def load_missing_multiome_barcode_info(
    gex_barcodes=None,
    atac_barcodes=None,
    gex_correction_table=None,
    atac_correction_table=None,
    atac_gex_translation_table=None
):
    
    if gex_barcodes is None:
        gex_barcodes = load_gex_barcodes()

    if atac_barcodes is None:
        atac_barcodes = load_atac_barcodes()

    if gex_correction_table is None:
        gex_correction_table = barcode_correction_table(gex_barcodes)

    if atac_correction_table is None:
        atac_correction_table = barcode_correction_table(atac_barcodes)

    if atac_gex_translation_table is None:
        atac_gex_translation_table = load_translations(
            atac_barcodes,
            gex_barcodes
        )

    return (
        gex_barcodes,
        atac_barcodes,
        gex_correction_table,
        atac_correction_table,
        atac_gex_translation_table
    )
