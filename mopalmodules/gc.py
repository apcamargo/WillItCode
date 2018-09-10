import numpy as np
from Bio import SeqUtils


def get_gc_content(orf_record):
    sequence_str = str(orf_record.seq.upper())
    if sequence_str:
        gc_content = SeqUtils.GC123(sequence_str)
        total_gc_content = gc_content[0]
        codon_gc_content = gc_content[1:3]
        gc_content_bias = np.max(codon_gc_content) - np.mean(codon_gc_content)
        return total_gc_content, gc_content_bias
    else:
        return 0, 0
