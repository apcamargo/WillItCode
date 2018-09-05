import argparse
import pickle
import sys

from mopalmodules.hexamers import get_hex_frequency_table


def main(coding_file, noncoding_file, output_file):
    """Calculate hexamer frequences in coding and noncoding transcripts and save them into a file."""
    hex_table = get_hex_frequency_table(coding_file, noncoding_file)
    with open(output_file, 'wb') as output:
        pickle.dump(hex_table, output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate hexamer frequences in coding and noncoding transcripts and save them into a file.'
    )
    parser.add_argument('coding_file', help='FASTA file containing CDS sequences of protein-coding transcripts.')
    parser.add_argument('noncoding_file', help='FASTA file containing complete sequences of noncoding transcripts.')
    parser.add_argument('output_file', help='Output file containing the hexamer frequency matrix.')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    main(**vars(args))
