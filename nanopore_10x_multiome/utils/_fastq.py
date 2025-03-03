### Pure python FASTQ parser ###

# Converts a quality ASCII string to a list of qualities
# This is the 33-offset illumina quality scoring
def convert_qual_illumina(qstr):
    return [ord(ch) - 33 for ch in qstr]

# Dict to map strings to quality score functions
PHRED_SCORE = {
    "illumina": convert_qual_illumina,
    "raw": lambda x: x
}


# Wrapper for fastqProcessor that takes a single file handle
# Returns a tuple of (control, sequence, quality) instead of a list of tuples
def fastq_gen(fh, phred_type='illumina'):
    for rec in fastqProcessor(phred_type=phred_type).fastq_gen(fh):
        yield rec[0]


# The fastqProcessor class takes an arbitrary number of linked reads and
# yields a list of tuples for each sequence
# The list is in the same order as the files that were given as arguments
# to fastq_gen
class fastqProcessor:
    phred = None
    verify_ids = True
    n_records = None

    def __init__(
        self,
        phred_type='illumina',
        verify_ids=True,
        n_records=None
    ):

        try:
            phred = PHRED_SCORE[phred_type]
        except KeyError:
            raise ValueError("Score type {} unknown".format(phred_type))

        self.phred = phred
        self.verify = verify_ids
        self.n_records = n_records

    def fastq_gen(self, *fhs):

        i = 0
        while True:

            try:
                record = []
                cid = None
                for fh in fhs:
                    (c, s, q) = self.fastq_process_file(fh, self.phred)

                    # If verify_ids was set, assert that the sequence IDs
                    # all match
                    if self.verify and cid is not None:
                        try:
                            assert c.startswith(cid)
                        except AssertionError:
                            print(f"ID mismatch (Record {i}): {cid} != {c}")
                            raise

                    if self.verify:
                        cid = self.extract_control_id(c)

                    record.append((c, s, q))

                i += 1
                yield record

                if self.n_records is not None and (i >= self.n_records):
                    break

            except StopIteration:
                break

    @staticmethod
    def fastq_process_file(fh, phred):
        cont, seq, qual = None, None, None
        line_id = -1
        for line in fh:
            line = line.strip()
            if line.startswith("+"):
                continue
            elif line.startswith("@"):
                line_id = 0
                cont = line
            elif line_id == 0:
                line_id = 1
                seq = line
            elif line_id == 1:
                qual = phred(line)
                return cont, seq, qual
        raise StopIteration

    @staticmethod
    def extract_control_id(con):
        try:
            return con.strip().split()[0]
        except IndexError:
            return None
