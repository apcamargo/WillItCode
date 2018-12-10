import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_codons(sequence_str, frame_number=0):
    coordinate = frame_number
    while coordinate + 3 <= len(sequence_str):
        yield (sequence_str[coordinate:coordinate+3], coordinate)
        coordinate += 3


def find_longest_in_frame(sequence_str, frame_number, start_codon, stop_codon):
    result = [0, 0, 0, 0]
    longest = 0
    triplet_got = get_codons(sequence_str, frame_number)
    while True:
        try:
            codon, index = next(triplet_got)
        except StopIteration:
            break
        if codon in start_codon and codon not in stop_codon:
            orf_start = index
            end_extension = False
            while True:
                try:
                    codon, index = next(triplet_got)
                except StopIteration:
                    end_extension = True
                    integrity = 0
                if codon in stop_codon:
                    end_extension = True
                    integrity = 1
                if end_extension:
                    orf_end = index + 3
                    length = (orf_end - orf_start)
                    if length > longest:
                        longest = length
                        result = [orf_start, orf_end, length, integrity]
                    break
    return result


def get_orf_record(record, start_codon=None, stop_codon=None):
    if start_codon is None:
        start_codon = {'ATG'}
    if stop_codon is None:
        stop_codon = {'TAG', 'TAA', 'TGA'}
    sequence_str = str(record.seq.upper())
    orf_sequence = ''
    longest_orf_result = [0, 0, 0, 0]
    for frame_number in range(3):
        result = find_longest_in_frame(sequence_str, frame_number, start_codon, stop_codon)
        if result[2] > longest_orf_result[2]:
            longest_orf_result = result
    orf_sequence = sequence_str[longest_orf_result[0]:longest_orf_result[1]]
    orf_record = SeqRecord(Seq(orf_sequence), id=record.id, name=record.id, description=record.id)
    orf_integrity = longest_orf_result[3]
    return orf_record, orf_integrity


def get_lengths(record, orf_record):
    sequence_length = len(record)
    orf_length = len(orf_record)
    log_sequence_length = np.log1p(sequence_length)
    log_orf_length = np.log1p(orf_length)
    orf_ratio = orf_length/sequence_length
    return log_sequence_length, log_orf_length, orf_ratio
