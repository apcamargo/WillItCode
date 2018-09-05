from Bio import SeqUtils


def get_gc_content(orf_record):
    sequence_str = str(orf_record.seq.upper())
    if sequence_str:
        gc_content = SeqUtils.GC123(sequence_str)
        total_gc_content = gc_content[0]
        codon_gc_content = gc_content[1:3]
        try:
            gc_content_skew = max(codon_gc_content)/sum(codon_gc_content)
        except ZeroDivisionError:
            gc_content_skew = 0
        return total_gc_content, gc_content_skew
    else:
        return 0, 0
