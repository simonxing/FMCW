#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')

def tx(F, Fs, alpha, tx_time):
    freq = F + alpha * tx_time
    return freq

F = 77.e9
B = 4.e9
T = 40.e-6
alpha = B/T

Fs = 2*F
Ts = 1./Fs
N = round(T/Ts)

initial, final = 0, T
tx_time = 0

file = open('datafile.txt', 'w')

for n in range(5):
    for time in np.linspace(initial, final, N):
        tx_freq = tx(F, Fs, alpha, tx_time)
        file.write("%e \t %e\n" %(time, tx_freq))
        tx_time += Ts
        tx_freq_ = "%.6e" %(tx_freq)
        if float(tx_freq_) == F+B: tx_time = 0
    initial += T; final += T

file.close()

data = np.loadtxt('datafile.txt')

plt.plot(data[:, 0], data[:, 1])
plt.show()
