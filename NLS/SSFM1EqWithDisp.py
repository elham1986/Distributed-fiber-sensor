# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 21:40:47 2021

@author: Elham
"""


# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 09:44:41 2021

@author: elham.rahmati
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
fiblen = 1            # fiber length in units of the dispersion length L_D
L_D = 50000
N = 1                 # soliton order
lambda_s = 1550e-9    # in m
lambda0 = 1e-6        # in m
coeff1 = 0.5
alpha_dBs = 0.18
alpha_dB = 0.18
c = 3e8               # Light speed in m/s
vg_s_D = c/1.457        # Group velocity for laser pulse
Aeff = 100e-12        # Effective aera in m**2        
ws = 2 * np.pi * c / lambda_s     # Laser frequency in rad/s
t0 = 5e-9                         # in s

# function that returns alpha (dB) to alpha (1/m)
def dB_To_mInv(dB):
    return (dB/20)*(10**-3)*np.log(10)

# Convert alpha (dB) to alpha (1/m) using function dB_To_mInv
alpha_s_D = dB_To_mInv(alpha_dBs)    # D in the name of parameter means the parameter with original unit
alpha_D = dB_To_mInv(alpha_dB)       # D in the name of parameter means the parameter with original unit

#Dimensionless parameters
alpha_s = alpha_s_D * L_D
vg_s = vg_s_D / L_D

############### set simulation parameters ###############
nt = 1024                              # FFT points
Taumax = 150e-9                        # window size
step_num = round(10001 * fiblen * N**2)  # No. of Z steps
deltaξ = fiblen/step_num               # Step size in Z
dtau = (2*Taumax)/nt                   # Step size in tau

def pulse(t, width, tau0, amp):
    if t > width/2 :
        return amp*np.exp(-((t - width/2)/(2*tau0))**2)
    elif t < -width/2 :
        return amp*np.exp(-((t + width/2)/(2*tau0))**2)
    else:
        return amp
    
Ps0 = 1.0
tau = np.arange(-nt/2,nt/2 + 1)*dtau /t0                               # Time array dimensionless
omega = scipy.fft.fftshift(np.arange(-nt/2,nt/2 + 1))*(np.pi/Taumax)   # Omega array
pulse_s = np.array([pulse(t,10, t0, np.sqrt(Ps0)) for t in tau])       # Tunable pulse shape (can be modified)

list_s = pulse_s

middle_pulse_s = int(len(pulse_s) / 2) # = int(2.5) = 2
print(list_s[middle_pulse_s]) # 3


uu_s = pulse_s / np.sqrt(Ps0)     # Field

plt.figure()
plt.plot(tau,pulse_s,'-r', label='Ps')
plt.legend()
plt.show()

# save to csv file
#savetxt('Pulse_Narrow.csv', AllDataPulse, delimiter=',')
#savetxt('Pulse_Thick.csv', AllDataPulse, delimiter=',')
#savetxt('Pulse_NewScale.csv', AllDataPulse, delimiter=',')

############### plot input pulse shape and spectrum ###############
temp_s = scipy.fft.fftshift(scipy.fft.ifft(uu_s))    # Fourier transform
spect_s = abs(temp_s)**2                             # Input spectrum
spect_s = spect_s/np.max(spect_s)                    # Normalize
freq = scipy.fft.fftshift(omega)/(2*np.pi)           # Freq. Array

############### store dispersive phase shifts to speedup code ###############
dispersion_s = np.exp((1j*omega/vg_s - alpha_s/2)*deltaξ) # Phase factor


############### begining of main loop ###############
# Scheme : 1/2N -> D -> 1/2N    first half step nonlinear
hhz_s = 0
temp_s = uu_s * np.exp( hhz_s/2)   # note hhz/2

Zarray = np.zeros((step_num,1))
Psarray = np.zeros((step_num,1))
P0_s = max(abs(uu_s)**2) * Ps0
pulse_s = abs(uu_s)**2 * Ps0 
Zarray[0] = 0
Psarray[0] = P0_s
for i in np.arange(0,step_num-1):
    print(i/(step_num - 2)*100)
    f_temp_s = scipy.fft.ifft(temp_s)*dispersion_s
    uu_s = scipy.fft.fft(f_temp_s)
    temp_s = uu_s * np.exp(hhz_s)
    P_s = max(abs(uu_s)**2) * Ps0
    pulse_s = abs(uu_s)**2 * Ps0 
    Zarray[i + 1] = (i + 1) * deltaξ * L_D
    Psarray[i + 1] = P_s


# plot usingwith matplotlib package using pyplot
# setup figures
fig = plt.figure()
plt.plot(Zarray, Psarray,'-r', label='Ps_SSFM')
plt.legend(loc='best')
plt.suptitle("""Power of pump and signal laser in the Fiber """, fontweight ="bold")  
plt.xlabel("Positions in the Fiber (m)")
plt.ylabel("Optical powers (W)")
plt.grid()
plt.savefig('Ps_Pp.jpg')
plt.show()

fig = plt.figure()
plt.plot(tau*t0, pulse_s,'-r', label='Ps_SSFM')
plt.legend(loc='best')  
plt.xlabel("Time(s)")
plt.ylabel("Pulse")
plt.grid()
plt.show()