from collections import Counter

import numpy as np
import scipy.stats as ss

from willitcode.orf import get_codons


def get_entropy(record):
    sequence_str = str(record.seq.upper())
    orf_counter = Counter(codon[0] for codon in orf.get_codons(sequence_str))
    total_codon_count = sum(orf_counter.values())
    if total_codon_count <= 1:
        return 1
    else:
        orf_frequencies = []
        for codon in orf_counter.keys():
            orf_frequencies.append(orf_counter[codon]/total_codon_count)
        entropy = ss.entropy(orf_frequencies)
        norm_entropy = entropy/np.log(len(orf_frequencies))
        return norm_entropy
