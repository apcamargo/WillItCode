from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import ProtParam


def get_protein_record(orf_record, codon_table='Standard'):
    orf_length = len(orf_record)
    if orf_length != 0:
        if orf_length%3 != 0:
            orf_record_trim = orf_record[:orf_length-orf_length%3]
        else:
            orf_record_trim = orf_record
        protein_sequence = str(orf_record_trim.translate(table=codon_table, to_stop=True).seq)
    else:
        protein_sequence = ''
    protein_record = SeqRecord(Seq(protein_sequence), id=orf_record.id, name=orf_record.id, description=orf_record.id)
    return protein_record

def get_protein_pi(protein_record):
    protein_sequence = str(protein_record.seq)
    if protein_sequence:
        protparam_obj = ProtParam.ProteinAnalysis(protein_sequence)
        protein_pi = protparam_obj.isoelectric_point()
    else:
        protein_pi = 0
    return protein_pi
