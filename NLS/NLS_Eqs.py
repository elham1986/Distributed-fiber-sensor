# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 09:28:18 2021

@author: elham.rahmati
"""
import numpy as np
import scipy
from mpmath import sech
import math
from matplotlib import pyplot as plt

# SSFM code for solving the normalized NLS equation

fiblen = 5    # fiber length in units of the dispersion length L_D
beta2 = -1    # Sign of GVD parameter beta_2
N = 1         # soliton order

############### set simulation parameters ###############
nt = 1024     # FFT points and window size
Tmax = 32     # FFT points and window size
step_num = round(20*fiblen*N**2)  # No. of Z steps
deltaZ = fiblen/step_num          # Step size in Z
dtau = (2*Tmax)/nt                # Step size in tau

tau = np.arange(-nt/2,nt/2 + 1)*dtau   # Time array
omega = scipy.fft.fftshift(np.arange(-nt/2,nt/2 + 1))*(np.pi/Tmax) # Omega array
uu = 1/np.cosh(tau) # Sech pulse shape (can be modified)

############### plot input pulse shape and spectrum ###############
temp = scipy.fft.fftshift(scipy.fft.ifft(uu))   # Fourier transform
spect = abs(temp)**2                       # Input spectrum
spect = spect/np.max(spect)                # Normalize
freq = scipy.fft.fftshift(omega)/(2*np.pi) # Freq. Array

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(tau,abs(uu)**2, '--r')
ax1.set(xlabel='Normalized Time', ylabel='Normalized Power')
ax1.set_xlim(-5, 5)
ax1.set_ylim(0, 1)

ax2.plot(freq,spect, '--r')
ax2.set(xlabel='Normalized Frequency', ylabel='Spectral power')
ax2.set_xlim(-.5, .5)
ax2.set_ylim(0, 1)


############### store dispersive phase shifts to speedup code ###############
dispersion = np.exp(0.5j*beta2*omega**2*deltaZ) # Phase factor
hhz = 1j*N**2*deltaZ    # Nonlinear phase factor

############### begining of main loop ###############
# Scheme : 1/2N -> D -> 1/2N    first half step nonlinear
temp = uu*np.exp(abs(uu)**2*hhz/2)   # note hhz/2
for i in np.arange(0,step_num):
    f_temp = scipy.fft.ifft(temp)*dispersion
    uu = scipy.fft.fft(f_temp)
    temp = uu*np.exp(abs(uu)**2*hhz)
###############   End of main loop  ###############

###############   Plot output pulse shape and spectrum  ###############
temp = scipy.fft.fftshift(scipy.fft.ifft(uu))   # Fourier transform
spect = abs(temp)**2
spect = spect/np.max(spect)
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(tau,abs(uu)**2, '-b')
ax1.set(xlabel='Normalized Time', ylabel='Normalized Power')
ax1.set_xlim(-5, 5)
ax1.set_ylim(0, 1)

ax2.plot(freq,spect, '-b')
ax2.set(xlabel='Normalized Frequency', ylabel='Spectral power')
ax2.set_xlim(-.5, .5)
ax2.set_ylim(0, 1)
