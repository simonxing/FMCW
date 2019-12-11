def kaiser(data, cutoff, Fs, ripple_db):
    nyquist = 0.5 * Fs
    width = 4.e6/nyquist
    normal_cutoff = cutoff/nyquist
    M, beta = signal.kaiserord(ripple_db, width)
    taps = signal.firwin(M, normal_cutoff, window=('kaiser', beta))
    filtered_data = signal.lfilter(taps, 1.0, data)
    return M, taps, filtered_data

''' Kaiser Window '''
ripple_db = 25.0
cutoff = Fb
M, taps, kaiser_filtered_data = kaiser(data, cutoff, Fs, ripple_db)

w2, h2 = signal.freqz(taps, worN=8000)
