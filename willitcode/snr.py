import numpy as np


def get_snr(record):
    sequence_str = str(record.seq.upper())
    if len(sequence_str) < 3:
        return 0
    a_fft = np.fft.fft([1 if i == 'A' else 0 for i in sequence_str])
    t_fft = np.fft.fft([1 if i == 'T' else 0 for i in sequence_str])
    c_fft = np.fft.fft([1 if i == 'C' else 0 for i in sequence_str])
    g_fft = np.fft.fft([1 if i == 'G' else 0 for i in sequence_str])
    power_spectrum = np.abs(a_fft)**2 + np.abs(t_fft)**2 + np.abs(c_fft)**2 + np.abs(g_fft)**2
    signal = power_spectrum[round((len(power_spectrum)-1)/3)]
    snr = np.log(signal/np.average(power_spectrum[1:]))
    return snr
