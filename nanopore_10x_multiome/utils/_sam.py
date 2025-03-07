def sam_comment_to_tag(
    sam_file,
    output_file,
    comment_prefix=['CB=', 'CR=', 'CY=', 'UB=', 'UR=', 'UY='],
    tag_prefix=['CB:Z:', 'CR:Z:', 'CY:Z:', 'UB:Z:', 'UR:Z:', 'UY:Z:']
):
    """Convert SAM comment fields to standard SAM tags.

    Takes comment entries like 'CB=ACGT' and converts them to SAM format tags in the correct
    positions.

    :param sam_file: Path to input SAM file
    :type sam_file: str
    :param output_file: Path to output SAM file
    :type output_file: str 
    :param comment_prefix: List of comment prefixes to convert (e.g. ['CB=', 'CR='])
    :type comment_prefix: list[str] or str
    :param tag_prefix: List of SAM tag formats to convert to (e.g. ['CB:Z:', 'CR:Z:'])
    :type tag_prefix: list[str] or str
    """
    
    # Convert single strings to lists for consistent handling
    if not isinstance(comment_prefix, (tuple, list)):
        comment_prefix = [comment_prefix]
    if not isinstance(tag_prefix, (tuple, list)):
        tag_prefix = [tag_prefix]

    # Validate matching prefix lengths
    assert len(comment_prefix) == len(tag_prefix)
    
    # Pre-calculate prefix lengths for efficiency
    _prelen = [len(x) for x in comment_prefix]

    with open(sam_file, mode='r') as sam_fh:
        with open(output_file, mode='w') as out_fh:
            for line in sam_fh:
                line = line.strip()
                # Pass through header lines unchanged
                if line.startswith("@"):
                    print(line, file=out_fh)
                    continue
                            
                # Split line into fields and get comments section
                _header = line.split('\t')
                _comments = _header[-1]

                # Store converted tags
                _new_tag = []

                # Process each prefix pair
                for (
                    c_pref,
                    tag_pref,
                    c_len
                ) in zip(
                    comment_prefix,
                    tag_prefix,
                    _prelen
                ):
                    # Look for comment prefix in comments section
                    _tag_loc = _comments.find(c_pref)

                    if _tag_loc == -1:
                        continue

                    trailing_comment = _comments[_tag_loc + c_len:]

                    if len(trailing_comment) == 0:
                        extracted_tag = ''
                    elif trailing_comment[0] == ' ':
                        extracted_tag = ''
                    else:
                        extracted_tag = trailing_comment.split()[0]

                    # Extract value and convert to SAM tag format
                    _new_tag.append(tag_pref + extracted_tag)

                # Output original line if no tags converted
                if len(_new_tag) == 0:
                    print(line, file=out_fh)
                else:
                    # Insert new tags before comments field
                    _header = _header[:-1] + _new_tag + [_comments]
                    print("\t".join(_header), file=out_fh)