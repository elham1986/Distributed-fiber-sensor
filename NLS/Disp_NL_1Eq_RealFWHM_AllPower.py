# -*- coding: utf-8 -*-
"""
Created on Thu May  6 14:44:46 2021

@author: elham.rahmati
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May  6 10:53:14 2021

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
fiblen = 4           # fiber length in units of the dispersion length L_D
L_D = 1                 # in km
N = 1                   # soliton order
lambda_p = 1450e-12     # in km
lambda_s = 1550e-12     # in km
lambda0 = 1e-9          # in km
coeff1 = 0.5
alpha_dBs = 0.18
c = 3e-7                # Light speed in km/ps
vg_s_D = c/1.457        # Group velocity for laser pulse
ws = 2 * np.pi * c / lambda_s  # Laser frequency in rad/s
t0 = 5                   # in ps, in this case we can calculate FWHM = 2*np.sqrt(np.log(2))*t0 = 5e3 ps
vg_s = vg_s_D/L_D              # Group velocity for laser pulse Dimensionless
BetaS_D = 1/vg_s
BetaS2_D = 17                  # ps/nm * km
#omega0 = 0.02144
gamma = 2                      # 1/(W km)

# function that returns alpha (dB) to alpha (1/km)
def dB_To_mInv(dB):
    return (dB/20)*np.log(10)

# function that returns dBm to W
def dBm_To_W(dBm):
    return 10**(dBm/10)/1000

power_dBm_list = [30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20]
power_input_W = np.array([dBm_To_W(i) for i in power_dBm_list])

############### set simulation parameters ###############
nt = 1024                              # FFT points
Taumax = 150                           # window size
step_num = round(20 * fiblen * N**2)   # No. of Z steps
deltaξ = fiblen/step_num               # Step size in Z
dtau = (2*Taumax)/nt                   # Step size in tau

def pulse(t, width, tau0, amp):
    if t > width/2 :
        return amp*np.exp(-((t - width/2)/(np.sqrt(2)*tau0))**2)
    elif t < -width/2 :
        return amp*np.exp(-((t + width/2)/(np.sqrt(2)*tau0))**2)
    else:
        return amp

tau = np.arange(-nt/2,nt/2-1)*dtau                                   # Time array dimensionless
omega = scipy.fft.fftshift(np.arange(-nt/2,nt/2-1))*(np.pi/Taumax)   # Omega array

fig = plt.figure()
count = 0
P0_s = np.zeros((len(power_input_W)))
for p in power_input_W:
    Ps0 = p
    count = count + 1
    amp0 = np.sqrt(Ps0)
    print(p)
    pulse_s = np.array([pulse(t,0.0, t0, amp0) for t in tau])      # Tunable pulse shape (can be modified)

    uu_s = pulse_s / np.sqrt(Ps0)     # Field

    ############### plot input pulse shape and spectrum ###############
    temp_s = scipy.fft.fftshift(scipy.fft.ifft(uu_s))    # Fourier transform
    spect_s = abs(temp_s)**2                             # Input spectrum
    spect_s = spect_s/np.max(spect_s)                    # Normalize
    freq = scipy.fft.fftshift(omega)/(2*np.pi)           # Freq. Array
    
    ############### store dispersive phase shifts to speedup code ###############
    #dispersion_s = np.exp((-1j*omega**2*BetaS2/2- alpha_s/2) * deltaξ) # Phase factor
    dispersion_s = np.exp((1j*(omega)**2*BetaS2_D/2) * deltaξ) # Phase factor
    
    ############### begining of main loop ###############
    # Scheme : 1/2N -> D -> 1/2N    first half step nonlinear
    hhz_s = 1j * gamma * N**2* deltaξ    # Nonlinear phase factor
    temp_s = uu_s * np.exp(abs(uu_s)**2 * hhz_s/2)   # note hhz/2
    
    
    
    Zarray = np.zeros((step_num,1))
    Psarray = np.zeros((step_num,1))
    P0_s = max(abs(uu_s)**2)
    Zarray[0] = 0
    Psarray[0] = P0_s
    for i in np.arange(0,step_num-1):
        # print(i/(step_num - 2)*100)
        f_temp_s = scipy.fft.ifft(temp_s)*dispersion_s
        uu_s = scipy.fft.fft(f_temp_s)
        temp_s = uu_s * np.exp(hhz_s)
        P_s = max(abs(uu_s)**2)
        pulse_s = abs(uu_s) * np.sqrt(Ps0)
        Zarray[i + 1] = (i + 1) * deltaξ * L_D
        Psarray[i + 1] = P_s
        
    temp_s_f = temp_s * np.exp(-abs(uu_s)**2 * hhz_s/2)  # final field
    
    plt.plot(tau, pulse_s, label= "{0:.2f}".format(Ps0))
    plt.legend(loc='best')
    plt.xlabel("Time (s)")
    plt.ylabel("Field")
    plt.grid()

plt.show()
