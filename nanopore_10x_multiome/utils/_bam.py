import pysam

def write_bam_record(
    handle,
    header,
    sequence,
    qual,
    flag=0,
    **tags
):
    
    a = pysam.AlignedSegment()
    a.query_name = header
    a.query_sequence = sequence
    a.flag = flag
    a.query_qualities = pysam.qualitystring_to_array(qual)
    a.tags = [
        (k, v)
        for k, v in tags.items()
        if v is not None
    ]

    handle.write(a)
