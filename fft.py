def FFT(Ts, filtered_data):
    NFFT = len(filtered_data)
    f_values = np.linspace(0.0, 1.0/(2.0*Ts), NFFT//2)
    fft_values = fftpack.fft(filtered_data)
    fft_values = 2.0/NFFT * np.abs(fft_values[0:NFFT//2])
    return f_values, fft_values

''' FFT '''
f_values1, fft_values1 = FFT(Ts, butter_filtered_data)
f_values2, fft_values2 = FFT(Ts, kaiser_filtered_data) 
