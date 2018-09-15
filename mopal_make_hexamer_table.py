import argparse
import pickle
import sys

from Bio import SeqIO

from mopalmodules.hexamers import get_hex_frequency_table
from mopalmodules.orf import get_orf_record


def main(protein_fasta, noncoding_fasta, output_file):
    """Calculate hexamer frequences in coding and noncoding transcripts and save them into a file."""
    cds_records = []
    noncoding_records = []
    for record in SeqIO.parse(protein_fasta, 'fasta'):
        orf_record = get_orf_record(record)[0]
        if len(orf_record) >= 6:
            cds_records.append(orf_record)
    for record in SeqIO.parse(noncoding_fasta, 'fasta'):
        orf_record = get_orf_record(record)[0]
        if len(orf_record) >= 6:
            noncoding_records.append(record)
    hex_table = get_hex_frequency_table(cds_records, noncoding_records)
    with open(output_file, 'wb') as output:
        pickle.dump(hex_table, output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate hexamer frequences in coding and noncoding ORFs and save them into a file.')
    parser.add_argument('protein_fasta',
                        help='FASTA file containing complete sequences of protein-coding transcripts.')
    parser.add_argument('noncoding_fasta',
                        help='FASTA file containing complete sequences of noncoding transcripts.')
    parser.add_argument('output_file',
                        help='Output file containing the hexamer frequency matrix.')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    main(**vars(args))
