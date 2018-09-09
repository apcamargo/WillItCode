import argparse
import pickle
import sys

import numpy as np

from mopalmodules.classifier import get_feature_matrix


def main(fasta_file, hex_table_file, output_file, hmmer_cpu):
    """Compute the feature matrix of a input FASTA and saves it into a file."""
    hex_table = pickle.load(open(hex_table_file, 'rb'))
    print('* Computing feature matrix.')
    sequence_id_list, feature_matrix = get_feature_matrix(fasta_file, hex_table,
                                                          hmmer_cpu=hmmer_cpu)
    columns = ('sequence_id\torf_integrity\trecord_log_sequence_length\tlog_orf_length\t'
               'orf_ratio\tfickett_score\tgc_content\tgc_skew\thexamer_bias\t'
               'hexamer_bias_distance\tprotein_pi\tsnr\thmmer_score')
    matrix = np.column_stack((sequence_id_list, feature_matrix))
    with open(output_file, 'w') as output:
        output.write(columns)
        output.write('\n')
        for row in matrix:
            output.write('\t'.join(row))
            output.write('\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Classify sequences from a input FASTA file.')
    parser.add_argument('fasta_file',
                        help='FASTA file containing transcript sequences.')
    parser.add_argument('hex_table_file',
                        help='Hexamer frequency table file.')
    parser.add_argument('output_file',
                        help='Name of the output file containing the feature matrix.')
    parser.add_argument('--hmmer_cpu',
                        default=1,
                        help='Number of parallel CPU to use for multithreads in HMMER.')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    main(**vars(args))
