import os

import numpy as np
import xgboost as xgb
from Bio import SeqIO

from willitcode.fickett import get_fickett_score
from willitcode.gc import get_gc_content
from willitcode.hexamers import get_hexamer_bias
from willitcode.hmmer import get_hmmer
from willitcode.orf import get_lengths, get_orf_record
from willitcode.protein import get_protein_pi, get_protein_record
from willitcode.snr import get_snr


def get_feature_matrix(fasta_file, hex_table, hmmer_cpu=1):
    fasta_file_basename = os.path.basename(fasta_file)
    protein_record_list = []
    sequence_features_list = []
    sequence_id_list = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequence_id_list.append(record.id)
        (orf_record, orf_integrity) = get_orf_record(record)
        protein_record = get_protein_record(orf_record)
        (record_log_sequence_length, record_log_orf_length,
         record_orf_ratio) = get_lengths(record, orf_record)
        record_fickett_score = get_fickett_score(orf_record)
        record_gc_content, record_gc_bias = get_gc_content(orf_record)
        record_hexamer_bias, record_hexamer_bias_distance = get_hexamer_bias(orf_record, hex_table)
        record_protein_pi = get_protein_pi(protein_record)
        record_snr = get_snr(orf_record)
        protein_record_list.append(protein_record)
        sequence_features_list.append([orf_integrity, record_log_sequence_length,
                                       record_log_orf_length, record_orf_ratio,
                                       record_fickett_score, record_gc_content, record_gc_bias,
                                       record_hexamer_bias, record_hexamer_bias_distance,
                                       record_protein_pi, record_snr])
    sequence_features_list = np.array(sequence_features_list)
    hmmer_feature_array = get_hmmer(protein_record_list,
                                    hmmer_directory_name=fasta_file_basename+'.HMMER',
                                    n_cpu=hmmer_cpu)
    feature_matrix = np.column_stack((sequence_features_list, hmmer_feature_array))
    return np.array(sequence_id_list), feature_matrix


def train_classifier(coding_file, noncoding_file, hex_table, hmmer_cpu=1):
    print('* Computing feature matrix of the coding sequences.')
    pc_matrix = get_feature_matrix(coding_file, hex_table, hmmer_cpu=hmmer_cpu)[1]
    print('* Computing feature matrix of the non-coding sequences.')
    nc_matrix = get_feature_matrix(noncoding_file, hex_table, hmmer_cpu=hmmer_cpu)[1]
    x_train = np.concatenate([pc_matrix, nc_matrix])
    y_train = np.repeat([1, 0], [pc_matrix.shape[0], nc_matrix.shape[0]])
    classification_model = xgb.XGBClassifier()
    print('* Training the classification model.')
    classification_model.fit(x_train, y_train)
    return classification_model


def classify_fasta(fasta_file, hex_table, classification_model, hmmer_cpu=1):
    print('* Computing feature matrix of the input sequences.')
    sequence_id_list, feature_matrix = get_feature_matrix(fasta_file, hex_table, hmmer_cpu=hmmer_cpu)
    print('* Classifying the input sequences.')
    prediction_proba = classification_model.predict_proba(feature_matrix)[:, 1]
    prediction_label = np.array(['Coding' if i > 0.5 else 'Non-coding' for i in prediction_proba])
    return sequence_id_list, feature_matrix, prediction_proba, prediction_label
