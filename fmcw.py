#!/usr/bin/env python3
import numpy as np
import scipy.signal as signal
import scipy.fftpack as fftpack 
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')

def tx(F, Fs, alpha, tx_time):
    freq = F + alpha * tx_time
    signal = np.sin(2*np.pi*freq*time)
    return freq, signal

def rx(F, Fs, alpha, rx_time, time, tau):
    freq = 0 if time < tau else F + alpha * (rx_time - tau)
    signal = np.sin(2*np.pi*freq*time)
    return freq, signal

def butter_lowpass(cutoff, Fs, order=5):
    nyquist = 0.5 * Fs
    normal_cutoff = cutoff/nyquist
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, Fs, order=5):
    b, a = butter_lowpass(cutoff, Fs, order=5)
    filtered_data = signal.filtfilt(b, a, data)
    return filtered_data

def kaiser(data, cutoff, Fs, ripple_db):
    nyquist = 0.5 * Fs
    width = 4.e6/nyquist
    normal_cutoff = cutoff/nyquist
    M, beta = signal.kaiserord(ripple_db, width)
    taps = signal.firwin(M, normal_cutoff, window=('kaiser', beta))
    filtered_data = signal.lfilter(taps, 1.0, data)
    return M, taps, filtered_data

def FFT(Ts, filtered_data):
    NFFT = len(filtered_data)
    f_values = np.linspace(0.0, 1.0/(2.0*Ts), NFFT//2)
    fft_values = fftpack.fft(filtered_data)
    fft_values = 2.0/NFFT * np.abs(fft_values[0:NFFT//2])
    return f_values, fft_values

F = 77.e9               # Start frequency (77 GHz)
B = 4.e9                # Chirp bandwidth (4 GHz)
T = 40.e-6              # Chirp period (40 us)
alpha = B/T             # Chirp rate (100 MHz/us)

R = 300                 # Target range (300 m)
c = 3.e8                # Speed of light (in m/s)
tau = 2.*R/c            # Two-way transit time (2 us) 
Fb = alpha*tau          # Beat frequency (200 MHz)
Fs = 2*B                # Sampling frequency (8 GHz)
Ts = 1./Fs              # Sampling period (0.125 ns)
N = round(T/Ts)         # Number of samples per sweep period (320000) 

file1 = open('datafile1.txt', 'w')
file2 = open('datafile2.txt', 'w')

IF_signal = []

samples, spacing = np.linspace(0, T, N, retstep=True)
time_sweep = spacing

initial, final = 0, T
tx_time, rx_time = 0, tau

for n in range(3):
    for time in np.linspace(initial, final, N):
        tx_freq, tx_signal = tx(F, Fs, alpha, tx_time)
        rx_freq, rx_signal = rx(F, Fs, alpha, rx_time, time, tau)
        IFsignal = tx_signal * rx_signal
        IF_signal.append(IFsignal)
        file1.write("%e \t %e \t %e\n" %(time, tx_freq, tx_signal))
        file2.write("%e \t %e \t %e\n" %(time, rx_freq, rx_signal))
        tx_time = tx_time + time_sweep
        rx_time = tau if np.all(time < tau) else rx_time + time_sweep
        tx_freq_ = "%.6e" %(tx_freq)
        if float(tx_freq_) == F+B: tx_time = 0
        rx_freq_ = "%.6e" %(rx_freq)
        if float(rx_freq_) == F+B: rx_time = tau
    initial += T; final += T

file1.close()
file2.close()

data = IF_signal
data1 = np.loadtxt('datafile1.txt')
data2 = np.loadtxt('datafile2.txt')

''' Butterworth Lowpass Filter '''
order = 5
cutoff = Fb
b, a = butter_lowpass(cutoff, Fs, order)
butter_filtered_data = butter_lowpass_filter(data, cutoff, Fs, order)

''' Kaiser Window '''
ripple_db = 25.0
M, taps, kaiser_filtered_data = kaiser(data, cutoff, Fs, ripple_db)

''' FFT '''
f_values1, fft_values1 = FFT(Ts, butter_filtered_data)
f_values2, fft_values2 = FFT(Ts, kaiser_filtered_data) 

''' Plot results '''
w1, h1 = signal.freqz(b, a, worN=8000)
w2, h2 = signal.freqz(taps, worN=8000)

''' chirp plot '''
plt.plot(data1[:,0], data1[:,1], label='TX chirp')
plt.plot(data2[:,0], data2[:,1], label='RX chirp')
plt.legend()
plt.xlim(0, 3*T)
plt.ylim(75.e9, 83.e9)
plt.xlabel('Time (S)')
plt.ylabel('Frequency (Hz)')

''' frequency response of butterworth lowpass filter '''
plt.plot(0.5*Fs*w1/np.pi, np.abs(h1), 'b', label='order = 5')
plt.plot(cutoff, 0.5*np.sqrt(2), 'ko')
plt.axvline(cutoff, color='k')
plt.xlim(0, 0.6e9)
plt.ylim(0, 1.2)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.legend()
plt.grid()

''' data and butter filtered data '''
plt.plot(data1[:,0], data, label='Data')
plt.plot(data1[:,0], butter_filtered_data, 'r-', label='Filtered Data')
plt.legend()
plt.xlim(2.017e-5, 2.020e-5)
plt.ylim(-1.25, 1.25)
plt.xlabel('Time (S)')
plt.ylabel('Amplitude')
plt.grid()

''' fft plot '''
plt.plot(f_values1, fft_values1, 'g')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.xlim(0, 8.e8)

''' impulse response '''
plt.plot(taps, 'b-', linewidth=2)
plt.xlim(1.e3, 1.4e3)
plt.ylim(-1.5e-2, 5.5e-2)
plt.xlabel('n')
plt.ylabel('h(n)')
plt.grid()

''' magnitude response '''
plt.plot((w2/np.pi)*(0.5*Fs), np.absolute(h2), linewidth=2)
plt.xlim(1.75e8, 2.25e8)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.grid()

''' data and kaiser window filtered data '''
plt.plot(data1[:,0], data, label = 'Data')
delay = 0.5*(M-1)/Fs
plt.plot(data1[:,0]-delay, kaiser_filtered_data, 'r-', label='Filtered Data')
plt.xlim(2.017e-5, 2.020e-5)
plt.ylim(-1.25, 1.25)
plt.xlabel('Time (S)')
plt.ylabel('Amplitude')
plt.grid()

''' fft plot '''
plt.plot(f_values2, fft_values2, 'g')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.xlim(0, 8.e8)

plt.tight_layout()
plt.show()
