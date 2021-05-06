# -*- coding: utf-8 -*-
"""
Created on Sun May  2 22:34:02 2021

@author: Elham
"""

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from numpy import savetxt
from matplotlib import pyplot as plt
plt.close('all')

# SSFM code for solving the normalized NLS equation

# Parameters
n = 10000
fiblen = 50000           # fiber length in units of the dispersion length L_D
L_D = 1                # km
N = 1                 # soliton order
lambda_p = 1450e-12    # in km
lambda_s = 1550e-12    # in km
lambda0 = 1e-9        # in km
coeff1 = 0.5
alpha_dBs = 0.18
c = 3 * 1e-7                 # Light speed in km/ps
vg_s_D = c/1.457        # Group velocity for laser pulse
ws = 2 * np.pi * c / lambda_s     # Laser frequency in rad/s
t0 = 5                         # in ps
vg_s = vg_s_D/L_D        # Group velocity for laser pulse Dimensionless
BetaS_D = 1/vg_s
BetaS2_D = 0.02         # ps/sqrt(km)
# BetaS2 =BetaS2_D /c

# function that returns alpha (dB) to alpha (1/m)
def dB_To_mInv(dB):
    return (dB/20)*(10**-3)*np.log(10)

# Convert alpha (dB) to alpha (1/m) using function dB_To_mInv
alpha_s_D = dB_To_mInv(alpha_dBs)    # D in the name of parameter means the parameter with original unit
# alpha_s = alpha_s_D * L_D

############### set simulation parameters ###############
nt = 1024                              # FFT points
Taumax = 150                        # window size
step_num = round(fiblen * N**2)  # No. of Z steps
deltaξ = fiblen /step_num            # Step size in Z
dtau = (2*Taumax)/nt                  # Step size in tau

def pulse(t, width, tau0, amp):
    if t > width/2 :
        return amp*np.exp(-((t - width/2)/(np.sqrt(2)*tau0))**2)
    elif t < -width/2 :
        return amp*np.exp(-((t + width/2)/(np.sqrt(2)*tau0))**2)
    else:
        return amp
    
Ps0 = 1.0
tau = np.arange(-nt/2,nt/2 - 1)*dtau                              # Time array dimensionless
omega = scipy.fft.fftshift(np.arange(-nt/2,nt/2 - 1))*(np.pi/Taumax)   # Omega array
pulse_s0 = np.array([pulse(t,0.0, t0, np.sqrt(Ps0)) for t in tau])        # Tunable pulse shape (can be modified)

uu_s = pulse_s0 / np.sqrt(Ps0)     # Field

############### plot input pulse shape and spectrum ###############
temp_s = scipy.fft.fftshift(scipy.fft.ifft(uu_s))    # Fourier transform
spect_s = abs(temp_s)**2                             # Input spectrum
spect_s = spect_s/np.max(spect_s)                    # Normalize
freq = scipy.fft.fftshift(omega)/(2*np.pi)           # Freq. Array

############### store dispersive phase shifts to speedup code ###############
#dispersion_s = np.exp((-1j*omega**2*BetaS2/2- alpha_s/2) * deltaξ) # Phase factor
dispersion_s = np.exp((-1j*omega**2*BetaS2_D/2) * deltaξ) # Phase factor

############### begining of main loop ###############
# Scheme : 1/2N -> D -> 1/2N    first half step nonlinear
hhz_s = 0
temp_s = uu_s * np.exp(hhz_s/2)   # note hhz/2

Zarray = np.zeros((step_num,1))
Psarray = np.zeros((step_num,1))
P0_s = max(abs(uu_s)**2)
Zarray[0] = 0
Psarray[0] = P0_s
for i in np.arange(0,step_num-1):
    print(i/(step_num - 2)*100)
    f_temp_s = scipy.fft.ifft(temp_s)*dispersion_s
    uu_s = scipy.fft.fft(f_temp_s)
    temp_s = uu_s * np.exp(hhz_s)
    P_s = max(abs(uu_s)**2)
    pulse_s = abs(uu_s) * np.sqrt(Ps0)
    Zarray[i + 1] = (i + 1) * deltaξ * L_D
    Psarray[i + 1] = P_s

z = np.array(np.linspace(0.0,L_D*fiblen,1025)) 

Z,Tau = np.meshgrid(z,tau)
pulse0_Analyt = np.sqrt(Ps0)*np.exp(-(tau)**2/(2*t0**2))
US_Analyt = t0*np.exp(-(Tau**2)/(2*(t0**2 - 1j*BetaS2_D*Z)))/np.sqrt(t0**2 - 1j*BetaS2_D*Z)
pulse_Analyt = abs(US_Analyt) * np.sqrt(Ps0)

# plot usingwith matplotlib package using pyplot
# setup figures
fig = plt.figure()
plt.plot(np.array(tau),pulse_s0,'-b', label='PulseS0_SSFM')
plt.plot(np.array(tau), pulse_s,'-m', label='PulseS_SSFM')
plt.plot(tau, pulse0_Analyt, '--r', label='PulseS0_Analyt')
plt.plot(tau, pulse_Analyt[:,-1], '--g', label='PulseS_Analyt')
plt.legend(loc='best')
plt.suptitle("""Power of pump and signal laser in the Fiber """, fontweight ="bold")  
plt.xlabel("Time (s)")
plt.ylabel("Pulse")
plt.grid()
#plt.savefig('Pulse_shape_with_dispersion.jpg')
plt.show()
