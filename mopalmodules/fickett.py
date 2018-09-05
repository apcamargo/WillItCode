def look_up_position_prob(value, base):
    """Look up positional probability by base and value"""
    position_prob_table = {
        'A': [0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22],
        'C': [0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23],
        'G': [0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08],
        'T': [0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09]
    }
    position_weight = {'A': 0.26, 'C': 0.18, 'G': 0.31, 'T': 0.33}
    position_para = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
    if value < 0:
        return None
    for idx, val in enumerate(position_para):
        if value >= val:
            pos_prob = position_prob_table[base][idx] * position_weight[base]
            return pos_prob


def look_up_content_prob(value, base):
    """Look up content probability by base and value"""
    content_prob_table = {
        'A': [0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.67, 0.65, 0.81, 0.21],
        'C': [0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.43, 0.44, 0.39, 0.31],
        'G': [0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29],
        'T': [0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58]
    }
    content_weight = {'A': 0.11, 'C': 0.12, 'G': 0.15, 'T': 0.14}
    content_para = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0]
    if value < 0:
        return None
    for idx, val in enumerate(content_para):
        if value >= val:
            con_prob = content_prob_table[base][idx] * content_weight[base]
            return con_prob


def get_fickett_score(orf_record):
    """Calculate the Fickett score for the whole sequence"""
    sequence_str = str(orf_record.seq.upper())
    if len(sequence_str) < 2:
        return 0
    fickett_score = 0
    total_base = len(sequence_str)
    a_content = sequence_str.count('A')/total_base
    c_content = sequence_str.count('C')/total_base
    g_content = sequence_str.count('G')/total_base
    t_content = sequence_str.count('T')/total_base
    phase_0 = [sequence_str[i] for i in range(0, len(sequence_str)) if i % 3 == 0]
    phase_1 = [sequence_str[i] for i in range(0, len(sequence_str)) if i % 3 == 1]
    phase_2 = [sequence_str[i] for i in range(0, len(sequence_str)) if i % 3 == 2]
    a_position = (max(phase_0.count('A'), phase_1.count('A'), phase_2.count('A'))
                  / (min(phase_0.count('A'), phase_1.count('A'), phase_2.count('A')) + 1))
    c_position = (max(phase_0.count('C'), phase_1.count('C'), phase_2.count('C'))
                  / (min(phase_0.count('C'), phase_1.count('C'), phase_2.count('C')) + 1))
    g_position = (max(phase_0.count('G'), phase_1.count('G'), phase_2.count('G'))
                  / (min(phase_0.count('G'), phase_1.count('G'), phase_2.count('G')) + 1))
    t_position = (max(phase_0.count('T'), phase_1.count('T'), phase_2.count('T'))
                  / (min(phase_0.count('T'), phase_1.count('T'), phase_2.count('T')) + 1))
    fickett_score += look_up_content_prob(a_content, 'A')
    fickett_score += look_up_content_prob(c_content, 'C')
    fickett_score += look_up_content_prob(g_content, 'G')
    fickett_score += look_up_content_prob(t_content, 'T')
    fickett_score += look_up_position_prob(a_position, 'A')
    fickett_score += look_up_position_prob(c_position, 'C')
    fickett_score += look_up_position_prob(g_position, 'G')
    fickett_score += look_up_position_prob(t_position, 'T')
    return fickett_score
