import pysam
import gzip as gz


from ._fastq import (
    fastq_gen,
    fastqProcessor,
    convert_qual_illumina,
    write_fastq_record
)

from ._bam import (
    write_bam_record
)

from ._sam import (
    sam_comment_to_tag
)

from ._sequence import (
    RC,
    REV
)

from ._parasail_barcode import (
    get_barcode_parasail
)


def file_opener(file_name, mode='r', file_format=None, gzip=False, header=None):

    if file_format is None:
        if file_name.endswith('.bam'):
            file_format = 'bam'
        elif file_name.endswith('.fastq.gz'):
            file_format = 'fastq'
            gzip = True
        elif file_name.endswith('.fastq'):
            file_format = 'fastq'
        else:
            raise ValueError(f"Unknown file format: {file_name}")

    if file_format == 'bam':
        if 'b' not in mode:
            mode = mode + 'b'

        return pysam.AlignmentFile(file_name, mode, header=header)
    
    elif gzip:
        return gz.open(file_name, mode=mode)
    
    else:
        return open(file_name, mode)


def get_file_writer(file_name=None, file_format=None):

    if file_format is None:
        if file_name.endswith('.bam'):
            file_format = 'bam'
        elif file_name.endswith('.fastq.gz'):
            file_format = 'fastq'
        elif file_name.endswith('.fastq'):
            file_format = 'fastq'
        else:
            raise ValueError(f"Unknown file format: {file_name}")

    if file_format == 'fastq':
        return write_fastq_record
    elif file_format == 'bam':
        return write_bam_record
    else:
        raise ValueError(f"Unknown file format: {file_format}")
