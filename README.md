# nanopore_10x_multiome
Processing for nanopore sequenced 10x multiome preamp libraries


```
split_multiome_preamp_fastq(
    in_file_name,
    atac_file_name,
    gex_file_name,
    other_file_name,
    atac_technical_file_name=None,
    n_jobs=None,
    write_only_valid_barcodes=False,
    keep_runoff_fragments=False,
    verbose=0
)

Split multiome pre-amplification FASTQ file(s) into ATAC, GEX and other reads.

:param in_file_name: Input FASTQ file path(s)
:type in_file_name: str, list
:param atac_file_name: Output FASTQ file path(s) for ATAC reads
:type atac_file_name: str, list
:param gex_file_name: Output FASTQ file path(s) for GEX reads
:type gex_file_name: str, list
:param other_file_name: Output FASTQ file path(s) for unidentified reads
    (genomic, too many errors near barcode, etc)
:type other_file_name: str, list
:param atac_technical_file_name: Optional output file(s) for ATAC technical sequences
:type atac_technical_file_name: str, list, or None
:param n_jobs: Number of parallel processes for joblib, defaults to None
:type n_jobs: int or None
:param keep_runoff_fragments: Keep ATAC fragments where the barcode end is intact,
    but no Tn5 site is located on the other end. Defaults to False.
:type keep_runoff_fragments: bool
:param write_only_valid_barcodes: Only write reads with valid barcodes
:type write_only_valid_barcodes: bool
:param verbose: Verbose parameter for joblib.Parallel
:type verbose: int

:return: Array of counts with n_files x [ATAC reads, GEX reads, other reads]
:rtype: numpy.ndarray
```



```
split_bam_by_barcode(
    bam_file,
    lookup_table,
    out_path=None,
    output_files=None,
    out_prefix='',
    barcode_tag='CB',
    pbar=False
)

Split a BAM file into multiple files based on barcode assignments.

:param bam_file: Path to input BAM file
:type bam_file: str
:param lookup_table: Dictionary mapping barcodes to their assigned groups/classes
:type lookup_table: dict
:param out_path: Directory path for output files if output_files not provided
:type out_path: str or None
:param output_files: Dictionary mapping groups to output file paths
:type output_files: dict or None
:param out_prefix: Prefix to add to output filenames if using out_path
:type out_prefix: str
:param barcode_tag: BAM tag containing the barcode
:type barcode_tag: str
:param pbar: Whether to show progress bar
:type pbar: bool
:return: Dictionary with number of reads written to each output file
:rtype: dict
```