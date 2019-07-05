#!/usr/bin/env python3

import numpy as np
import scipy.fftpack
import matplotlib.pyplot as plt

def tx(F, Fs, alpha, time, tx_time):
    freq = F + alpha * tx_time
    signal = np.sin(2*np.pi*freq*time/Fs)
    return freq, signal

def rx(F, Fs, alpha, time, tau, rx_time):
    if (time < tau):
        freq = 0; signal = 0
    else:
        freq = F + alpha * (rx_time-tau)
        signal = np.sin(2*np.pi*freq*time/Fs)
    return freq, signal

def if_(tx_freq, rx_freq):
    return tx_freq - rx_freq

def FFT(Ts, signal):
    M = len(signal)
    f_values = np.linspace(0.0, 1.0/(2.0*Ts), M//2)
    fft_values_ = scipy.fftpack.fft(signal)
    fft_values = 2.0/M * np.abs(fft_values_[0:M//2])
    return f_values, fft_values

F = 77.e9
B = 4.e9
T = 40.e-6
alpha = B/T

c = 3.e8
tau = 2.e-6
R = c*tau/2
print(R)

Fs = 4.e4
Ts = 1./Fs
N = int(T/Ts)
print(N)

signal = []

file = open('datafile.txt', 'w')

samples, spacing = np.linspace(0, T, Fs, retstep=True)
time_sweep = spacing

initial = 0; final = T
tx_time = 0; rx_time = tau

for n in range(8):
    for time in np.linspace(initial, final, Fs):

        tx_freq, tx_signal = tx(F, Fs, alpha, time, tx_time)
        rx_freq, rx_signal = rx(F, Fs, alpha, time, tau, rx_time)

        if_freq = if_(tx_freq, rx_freq)
        #if_signal = np.sin(2*np.pi*(if_freq)*time/4.e3)
        #if_signal = np.sin(2*np.pi*(tx_freq*rx_freq)*time/4.e4) 
        if_signal = tx_signal * rx_signal

        signal.append(if_signal)

        file.write("%e \t %e \t %e\n" %(time, tx_signal, if_signal))

        tx_time = tx_time + time_sweep
    
        if (time < tau): rx_time = tau
        else: rx_time = rx_time + time_sweep

        tx_freq_ = "%.6e" %(tx_freq)
        if float(tx_freq_) == F+B: tx_time = 0

        rx_freq_ = "%.6e" %(rx_freq)
        if float(rx_freq_) == F+B: rx_time = tau

    initial = initial + T
    final = final + T

file.close()

# FFT
f_values, fft_values = FFT(Ts, signal)

data = np.loadtxt('datafile.txt')

plt.plot(data[:,0], data[:,2], color='red')
plt.plot(f_values, fft_values)
plt.grid(True)
plt.show()
