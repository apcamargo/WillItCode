import argparse
import pickle
import sys

import numpy as np

from Bio import SeqIO
from willitcode.classifier import classify_fasta


def main(fasta_file, classification_model_file, hex_table_file, output_features,
         output_file, output_fasta, hmmer_cpu):
    """Classify sequences from a input FASTA file."""
    classification_model = pickle.load(open(classification_model_file, 'rb'))
    hex_table = pickle.load(open(hex_table_file, 'rb'))
    (sequence_id_list, feature_matrix,
     prediction_proba, prediction_label,
     coding_protein_record_list) = classify_fasta(fasta_file, hex_table,
                                                          classification_model,
                                                          hmmer_cpu=hmmer_cpu)
    if output_features:
        columns = ('sequence_id\torf_integrity\tsequence_length\torf_length\t'
                   'orf_ratio\tfickett_score\tgc_content\tgc_bias\thexamer_bias\t'
                   'hexamer_bias_distance\tprotein_pi\tcodon_entropy\tsnr\thmmer_score\t'
                   'coding_probability\tprediction')
        matrix = np.column_stack((sequence_id_list, feature_matrix, prediction_proba,
                                  prediction_label))
    else:
        columns = 'sequence_id\tcoding_probability\tprediction'
        matrix = np.column_stack((sequence_id_list, prediction_proba, prediction_label))
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
    if output_fasta is not None:
        with open(output_fasta, 'w') as fasta_file:
            SeqIO.write(coding_protein_record_list, fasta_file, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Classify sequences from a input FASTA file.')
    parser.add_argument('fasta_file',
                        help='FASTA file containing transcripts that will be classified.')
    parser.add_argument('classification_model_file',
                        help='File containing a pre-trained WillItCode model.')
    parser.add_argument('hex_table_file',
                        help='Hexamer frequency table file.')
    parser.add_argument('--output_features',
                        action='store_true',
                        help='Output computed features probability.')
    parser.add_argument('--output_file',
                        help='Save output to a file. If not set, the output will be printed on the screen.')
    parser.add_argument('--output_fasta',
                        help='Save protein sequences of the predicted coding ORFs into a FASTA file.')
    parser.add_argument('--hmmer_cpu',
                        default=1,
                        help='Number of parallel CPU to use for multithreads in HMMER.')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    main(**vars(args))
