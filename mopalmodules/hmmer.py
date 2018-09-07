import os
import shutil
import subprocess
import sys

import numpy as np
from Bio import SearchIO, SeqIO


def get_complete_path():
    dirname = os.path.dirname(sys.argv[0])
    complete_path = os.path.abspath(dirname)
    return complete_path


def save_protein_fasta(protein_record_list, protein_fasta_path):
    with open(protein_fasta_path, 'w') as handle:
        SeqIO.write(protein_record_list, handle, 'fasta')


def get_hmmer(protein_record_list, hmmer_directory_name='tmp.HMMER',
              hmmer_results='hmmscan.tblout', protein_fasta='translated_orfs.fa', n_cpu=1):
    complete_path = get_complete_path()
    hmmer_directory_path = os.path.join(complete_path, hmmer_directory_name)
    hmmer_output_path = os.path.join(complete_path, hmmer_directory_path, 'hmmscan.out')
    hmmer_tbl_output_path = os.path.join(complete_path, hmmer_directory_path, hmmer_results)
    protein_fasta_path = os.path.join(hmmer_directory_path, protein_fasta)
    pfam_path = os.path.join(complete_path, 'Pfam/Pfam-A.hmm')
    if os.path.isdir(hmmer_directory_path):
        shutil.rmtree(hmmer_directory_path)
    os.makedirs(hmmer_directory_path)
    save_protein_fasta(protein_record_list, protein_fasta_path)
    print('* Executing HMMER.')
    subprocess.Popen('hmmscan --cpu {} -o {} --tblout {} {} {}'.format(n_cpu, hmmer_output_path, hmmer_tbl_output_path, pfam_path, protein_fasta_path), shell=True, cwd=complete_path).wait()
    hmmer_score_dict = dict()
    for query in SearchIO.parse(hmmer_tbl_output_path, format='hmmer3-tab'):
        hmmer_score_dict[query.id] = max([hit.bitscore for hit in query.hits])
    hmmer_score_list = []
    for record in protein_record_list:
        if record.id in hmmer_score_dict:
            hmmer_score_list.append(hmmer_score_dict[record.id])
        else:
            hmmer_score_list.append(0)
    os.remove(protein_fasta_path)
    os.remove(hmmer_output_path)
    os.remove(hmmer_tbl_output_path)
    os.rmdir(hmmer_directory_path)
    return np.array(hmmer_score_list)
