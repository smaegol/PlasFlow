
import re


def get_kmer_counts(seq, K):
    """Return a list of kmers in a sequence."""
    out = {}
    for start in range(len(seq) - K + 1):
        kmer = seq[start:start + K]
        out[kmer] = 1 + out.get(kmer, 0)
    return out


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size."""
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def re_select_columns(row, tbl, reg):
    return row[[col for col in tbl.columns if re.match(reg, col)]]
