import argparse
import pickle
import sys

import numpy as np

from mopalmodules.classifier import classify_fasta


def main(fasta_file, classification_model_file, hex_table_file, output_features,
         output_file, hmmer_cpu):
    """Classify sequences from a input FASTA file."""
    classification_model = pickle.load(open(classification_model_file, 'rb'))
    hex_table = pickle.load(open(hex_table_file, 'rb'))
    (sequence_id_list, feature_matrix,
     prediction_proba, prediction_label) = classify_fasta(fasta_file, hex_table,
                                                          classification_model,
                                                          hmmer_cpu=hmmer_cpu)
    if output_features == 'no':
        columns = 'sequence_id\tprediction'
        matrix = np.column_stack((sequence_id_list, prediction_label))
    else:
        columns = ('sequence_id\torf_integrity\trecord_log_sequence_length\tlog_orf_length\t'
                   'orf_ratio\tfickett_score\tgc_content\tgc_skew\thexamer_bias\t'
                   'hexamer_bias_distance\tprotein_pi\tsnr\thmmer_score\tcoding_probability\t'
                   'prediction')
        matrix = np.column_stack((sequence_id_list, feature_matrix, prediction_proba,
                                  prediction_label))
    if output_file is None:
        print(columns)
        for row in matrix:
            print('\t'.join(row))
    else:
        with open(output_file, 'w') as output:
            output.write(columns)
            output.write('\n')
            for row in matrix:
                output.write('\t'.join(row))
                output.write('\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Classify sequences from a input FASTA file.'
    )
    parser.add_argument('fasta_file', help='FASTA file containing complete sequences of protein-coding transcripts.')
    parser.add_argument('classification_model_file', help='FASTA file containing complete sequences of noncoding transcripts.')
    parser.add_argument('hex_table_file', help='Hexamer frequency table file.')
    parser.add_argument('--output_features', choices=['yes', 'no'], default='no', help='Output computed features and the coding probability.')
    parser.add_argument('--output_file', help='Save output to a file. If not set, the output will be printed on the screen.')
    parser.add_argument('--hmmer_cpu', default=1, help='Number of parallel CPU to use for multithreads in HMMER.')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    main(**vars(args))
