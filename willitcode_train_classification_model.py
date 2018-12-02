import argparse
import pickle
import sys

from willitcode.classifier import train_classifier


def main(coding_file, noncoding_file, hex_table_file, hmmer_cpu, output_file):
    """Train a classification model from training data and saves it into a file."""
    hex_table = pickle.load(open(hex_table_file, 'rb'))
    classification_model = train_classifier(coding_file, noncoding_file, hex_table, hmmer_cpu=hmmer_cpu)
    with open(output_file, 'wb') as output:
        pickle.dump(classification_model, output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Train a classification model from training data and saves it into a file.')
    parser.add_argument('coding_file',
                        help='FASTA file containing complete sequences of protein-coding transcripts.')
    parser.add_argument('noncoding_file',
                        help='FASTA file containing complete sequences of noncoding transcripts.')
    parser.add_argument('hex_table_file',
                        help='Hexamer frequency table file.')
    parser.add_argument('output_file',
                        help='Output file containing the classification model.')
    parser.add_argument('--hmmer_cpu',
                        default=1, help='Number of parallel CPU to use for multithreads in HMMER.')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    main(**vars(args))
