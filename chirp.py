#!/usr/bin/env python3
import numpy as np
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
data2 = np.loadtxt('datafile2.txt', skiprows=16000)


''' Plot results '''

plt.plot(data1[:,0], data1[:,1], label='TX chirp')
plt.plot(data2[:,0], data2[:,1], label='RX chirp')
plt.legend()
plt.xlim(0, 3*T)
plt.ylim(75.e9, 83.e9)
plt.xlabel('Time (S)')
plt.ylabel('Frequency (Hz)')
plt.show()
