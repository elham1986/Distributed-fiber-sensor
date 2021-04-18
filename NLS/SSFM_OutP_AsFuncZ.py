# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 21:24:28 2021

@author: Elham
"""


import numpy as np
import scipy
from mpmath import sech
import math
from matplotlib import pyplot as plt
plt.close('all')

# SSFM code for solving the normalized NLS equation

fiblen = 50000        # fiber length in units of the dispersion length L_D
N = 1                 # soliton order
lambda_p = 1450e-9    # in m
lambda_s = 1550e-9    # in m
lambda0 = 1e-6        # in m
coeff1 = 0.5
alpha_dBs = 0.18
alpha_dBp = 0.19
c = 3e8               # Light speed in m/s
vg_s = c/1.457        # Group velocity for laser pulse
vg_p = c/1.457        # Group velocity for Pump
Aeff = 100e-12        # in m**2
gR = (1e-13)*lambda0*coeff1/(Aeff * lambda_p)    # Raman-gain coefficient in m/W
wp = 2 * np.pi * c / lambda_p     # Pump frequency in rad/s
ws = 2 * np.pi * c / lambda_s     # Laser frequency in rad/s
g_s = gR
g_p = (wp / ws) * gR
t0 = 5e-9             # in s

# function that returns alpha (dB) to alpha (1/m)
def dB_To_mInv(dB):
    return (dB/10)*(10**-3)*np.log(10)

# Convert alpha (dB) to alpha (1/m) using function dB_To_mInv
alpha_s = dB_To_mInv(alpha_dBs)
alpha_p = dB_To_mInv(alpha_dBp)

############### set simulation parameters ###############
nt = 1024                         # FFT points and window size
Tmax = 32                         # FFT points and window size
step_num = round(5*fiblen*N**2)  # No. of Z steps
deltaZ = fiblen/step_num          # Step size in Z
dtau = (2*Tmax)/nt                # Step size in tau

def pulse(t, width, gamma, amp):
    if t > width/2 :
        return amp*np.exp(-((t - width/2)/t0)**2*gamma)
    elif t < -width/2 :
        return amp*np.exp(-((t + width/2)/t0)**2*gamma)
    else:
        return amp
    
Ps0 = 0.2
Pp0 = 1
tau = np.arange(-nt/2,nt/2 + 1)*dtau                 # Time array
omega = scipy.fft.fftshift(np.arange(-nt/2,nt/2 + 1))*(np.pi/Tmax) # Omega array
uu_s = np.array([pulse(t,25,0.1*t0**2,np.sqrt(Ps0)) for t in tau])     # Tunable pulse shape (can be modified)
uu_p = np.array([pulse(t,25,0.1*t0**2,np.sqrt(Pp0)) for t in tau])   # Tunable pulse shape (can be modified)

plt.plot(tau,uu_s)
plt.plot(tau,uu_p)

############### plot input pulse shape and spectrum ###############
temp_s = scipy.fft.fftshift(scipy.fft.ifft(uu_s))    # Fourier transform
temp_p = scipy.fft.fftshift(scipy.fft.ifft(uu_p))    # Fourier transform
spect_s = abs(temp_s)**2                             # Input spectrum
spect_p = abs(temp_p)**2                             # Input spectrum
spect_s = spect_s/np.max(spect_s)                    # Normalize
spect_p = spect_p/np.max(spect_p)                    # Normalize
freq = scipy.fft.fftshift(omega)/(2*np.pi)    
 # Freq. Array

############### store dispersive phase shifts to speedup code ###############
dispersion_s = np.exp((1j*omega/vg_s - alpha_s/2)*deltaZ) # Phase factor
dispersion_p = np.exp((1j*omega/vg_p - alpha_p/2)*deltaZ) # Phase factor
hhz_s = g_s/2*deltaZ    # Nonlinear phase factor
hhz_p = -g_p/2*deltaZ    # Nonlinear phase factor

############### begining of main loop ###############
# Scheme : 1/2N -> D -> 1/2N    first half step nonlinear
temp_s = uu_s*np.exp(abs(uu_p)**2*hhz_s/2)   # note hhz/2
temp_p = uu_p*np.exp(abs(uu_s)**2*hhz_p/2)   # note hhz/2

Zarray = np.zeros((step_num,1))
Psarray = np.zeros((step_num,1))
Pparray = np.zeros((step_num,1))
P0_s = max(abs(uu_s)**2)
P0_p = max(abs(uu_p)**2)
Zarray[0] = 0
Psarray[0] = P0_s
Pparray[0] = P0_p
for i in np.arange(0,step_num-1):
    print(i/(step_num - 2)*100)
    f_temp_s = scipy.fft.ifft(temp_s)*dispersion_s
    f_temp_p = scipy.fft.ifft(temp_p)*dispersion_p
    uu_s = scipy.fft.fft(f_temp_s)
    uu_p = scipy.fft.fft(f_temp_p)
    temp_s = uu_s*np.exp(abs(uu_p)**2*hhz_s)
    temp_p = uu_p*np.exp(abs(uu_s)**2*hhz_p)
    P_s = max(abs(uu_s)**2)
    P_p = max(abs(uu_p)**2)
    Zarray[i + 1] = (i + 1) * deltaZ
    Psarray[i + 1] = P_s
    Pparray[i + 1] = P_p
    

plt.figure()
plt.plot(Zarray, Psarray,'--r', label='Ps')
plt.plot(Zarray, Pparray, '--b', label='Pp')
plt.legend()
plt.show()
