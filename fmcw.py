#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('seaborn-ticks')

def tx(F, Fs, alpha, tx_time):
    freq = F + alpha * tx_time
    signal = np.sin(2*np.pi*freq*tx_time/Fs)
    return freq, signal

def rx(F, Fs, alpha, rx_time, time, tau):
    freq = 0 if time < tau else F + alpha * (rx_time - tau)
    signal = 0 if time < tau else np.sin(2*np.pi*freq*(rx_time - tau)/Fs)
    return freq, signal


F = 77.e9
B = 4.e9
T = 40.e-6
alpha = B/T

R = 300
c = 3.e8
tau = 2.*R/c
Fb = alpha*tau
Fs = 2*(F+B)
Ts = 1./Fs
N = int(T/Ts)

file1 = open('datafile1.txt', 'w')
file2 = open('datafile2.txt', 'w')

signal = []

samples, spacing = np.linspace(0, T, N, retstep=True)
time_sweep = spacing

initial, final = 0, T
tx_time, rx_time = 0, tau

for n in range(2):
    for time in np.linspace(initial, final, N):
        tx_freq, tx_signal = tx(F, Fs, alpha, tx_time)
        rx_freq, rx_signal = rx(F, Fs, alpha, rx_time, time, tau)

        IF_signal = tx_signal * rx_signal
        signal.append(IF_signal)

        file1.write("%e \t %e \t %e \t %e\n" %(time, tx_freq, tx_signal, IF_signal))
        file2.write("%e \t %e \t %e \t %e\n" %(time, rx_freq, rx_signal, IF_signal))

        tx_time = tx_time + time_sweep
        rx_time = tau if np.all(time < tau) else rx_time + time_sweep

        tx_freq_ = "%.6e" %(tx_freq)
        if float(tx_freq_) == F+B: tx_time = 0

        rx_freq_ = "%.6e" %(rx_freq)
        if float(rx_freq_) == F+B: rx_time = tau

    initial += T; final += T

file1.close()
file2.close()

data = signal
data1 = np.loadtxt('datafile1.txt')
data2 = np.loadtxt('datafile2.txt')

plt.plot(data1[:, 0], data1[:, 2], label='TX chirp')
plt.plot(data2[:, 0], data2[:, 2], label='RX chirp')
plt.gca().xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
plt.gca().yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', scilimits=(0, 0), axis='both')
plt.show()
