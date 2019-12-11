def butter_lowpass(cutoff, Fs, order=5):
    nyquist = 0.5 * Fs
    normal_cutoff = cutoff/nyquist
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, Fs, order=5):
    b, a = butter_lowpass(cutoff, Fs, order=5)
    filtered_data = signal.filtfilt(b, a, data)
    return filtered_data


''' Butterworth Lowpass Filter '''
order = 5
cutoff = Fb
b, a = butter_lowpass(cutoff, Fs, order)
butter_filtered_data = butter_lowpass_filter(data, cutoff, Fs, order)

w1, h1 = signal.freqz(b, a, worN=8000)
