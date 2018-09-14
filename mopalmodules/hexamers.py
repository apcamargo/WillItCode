from collections import Counter, OrderedDict
from itertools import product

import numpy as np


def sliding_window(sequence, win_size, step):
    # https://scipher.wordpress.com/2010/12/02/simple-sliding-window-iterator-in-python/
    if win_size > len(sequence):
        raise Exception('Error: window size must not be larger than sequence length')
    n_chunks = int(((len(sequence)-win_size)/step)+1)
    for i in range(0, n_chunks*step, step):
        yield sequence[i:i+win_size]


def all_possible_kmer(win_size=6):
    for kmer in product(['A', 'T', 'C', 'G'], repeat=win_size):
        yield ''.join(kmer)


def get_hexamer_score(sequence_str, hex_table, win_size=6, step=3):
    if len(sequence_str) < win_size:
        return 0
    log_ratio_sum = 0
    hexamer_count = 0
    coding = hex_table['coding']
    noncoding = hex_table['noncoding']
    for hexamer in sliding_window(sequence_str, win_size, step):
        if (hexamer not in coding) or (hexamer not in noncoding):
            continue
        elif coding[hexamer] > 0 and noncoding[hexamer] > 0:
            log_ratio_sum += np.log(coding[hexamer]/noncoding[hexamer])
        elif coding[hexamer] > 0 and noncoding[hexamer] == 0:
            log_ratio_sum += 1
        elif coding[hexamer] == 0 and noncoding[hexamer] == 0:
            continue
        elif coding[hexamer] == 0 and noncoding[hexamer] > 0:
            log_ratio_sum -= 1
        hexamer_count += 1
    return log_ratio_sum/hexamer_count


def get_hexamer_bias(orf_record, hex_table):
    sequence_str = str(orf_record.seq.upper())
    frames_hexamer_score = []
    for frame in range(0, 3):
        frame_sequence_str = sequence_str[frame:]
        frames_hexamer_score.append(get_hexamer_score(frame_sequence_str, hex_table))
    hexamer_score = max(frames_hexamer_score)
    hexamer_score_distance = 0
    for frame in range(0, 3):
        hexamer_score_distance += hexamer_score-frames_hexamer_score[frame]
    hexamer_score_distance = np.log1p(hexamer_score_distance/2)
    return hexamer_score, hexamer_score_distance


def get_hex_frequency_table(cds_records, noncoding_records, win_size=6, step=3):
    hex_counter_coding = Counter()
    hex_counter_noncoding = Counter()
    hex_freq = OrderedDict({'coding': OrderedDict(), 'noncoding': OrderedDict()})
    for sequence in cds_records:
        if len(sequence) < win_size:
            continue
        hex_counter_coding.update(sliding_window(str(sequence.seq), win_size, step))
    total_hex_coding = sum(hex_counter_coding.values())
    for sequence in noncoding_records:
        if len(sequence) < win_size:
            continue
        hex_counter_noncoding.update(sliding_window(str(sequence.seq), win_size, step))
    total_hex_noncoding = sum(hex_counter_noncoding.values())
    for kmer in all_possible_kmer(win_size):
        if kmer in hex_counter_coding:
            hex_freq['coding'][kmer] = hex_counter_coding[kmer]/total_hex_coding
        else:
            hex_freq['coding'][kmer] = 0
        if kmer in hex_counter_noncoding:
            hex_freq['noncoding'][kmer] = hex_counter_noncoding[kmer]/total_hex_noncoding
        else:
            hex_freq['noncoding'][kmer] = 0
    return hex_freq
