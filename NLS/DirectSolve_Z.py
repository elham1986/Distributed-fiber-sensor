# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 20:56:12 2021

@author: Elham
"""


import numpy as np
import scipy
from mpmath import sech
import math
from matplotlib import pyplot as plt
import findiff as fdm
plt.close('all')

# SSFM code for solving the normalized NLS equation

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

# function that returns alpha (dB) to alpha (1/m)
def dB_To_mInv(dB):
    return (dB/10)*(10**-3)*np.log(10)

# Convert alpha (dB) to alpha (1/m) using function dB_To_mInv
alpha_s = dB_To_mInv(alpha_dBs)
alpha_p = dB_To_mInv(alpha_dBp)

def pulse(t, width, gamma, amp):
    if t > width/2 :
        return amp*np.exp(-(t - width/2)**2*gamma)
    elif t < -width/2 :
        return amp*np.exp(-(t + width/2)**2*gamma)
    else:
        return amp

def GetAmplitude(Ps0,Pp0, fiblen):
    nt = 1024                         # FFT points and window size
    Tmax = 32                         # FFT points and window size
    step_num = round(20*fiblen*N**2)  # No. of Z steps
    deltaZ = fiblen/step_num          # Step size in Z
    dtau = (2*Tmax)/nt                # Step size in tau
    tau = np.arange(-nt/2,nt/2 + 1)*dtau                 # Time array
    omega = np.fft.fftshift(np.arange(-nt/2,nt/2 + 1))*(np.pi/Tmax) # Omega array
    uu_s = np.array([pulse(t,25,0.1,np.sqrt(Ps0)) for t in tau])     # Tunable pulse shape (can be modified)
    uu_p = np.array([pulse(t,25,0.1,np.sqrt(Pp0)) for t in tau])   # Tunable pulse shape (can be modified)
    uu_s_0 = uu_s
    uu_p_0 = uu_p
    fdm_acc = 6
    d_dt = fdm.FinDiff(0, dtau, 1, acc=fdm_acc)
    Ds = -d_dt.matrix(tau.shape).toarray()/vg_s - alpha_s/2.0*np.identity(len(tau))
    Dp = -d_dt.matrix(tau.shape).toarray()/vg_p - alpha_p/2.0*np.identity(len(tau))
    uu_s_z = uu_s_0
    uu_p_z = uu_p_0
    Ns = g_s/2.0*np.diag(abs(uu_s_z)**2)
    Np = -g_p/2.0*np.diag(abs(uu_s_z)**2)
    
    for i in np.arange(0,100):
        uu_s_z_temp = uu_s_z
        uu_p_z_temp = uu_p_z
        uu_s_z = np.matmul( scipy.linalg.expm((Ds + Ns)*fiblen), uu_s_0)
        uu_p_z = np.matmul( scipy.linalg.expm((Dp + Np)*fiblen), uu_p_0)
        Ns = g_s / 2.0 * np.diag(abs(uu_p_z) ** 2)
        Np = -g_p / 2.0 * np.diag(abs(uu_s_z) ** 2)
        err_s = abs(np.amax(uu_s_z) - np.amax(uu_s_z_temp))
        err_p = abs(np.amax(uu_p_z) - np.amax(uu_p_z_temp))
        print([err_s, err_p])
        if max(err_s, err_p) < 0.0001:
            break
    Us = abs(uu_s_z)**2
    Up = abs(uu_p_z)**2
    print("")
    return [max(Us) , max(Up)]

ZArray = np.linspace(0.1,30000,30)
OUT = np.array([GetAmplitude(0.3,1, z) for z in ZArray])
print('')
print(OUT)
plt.figure()
plt.plot(ZArray, OUT[:,0] , '--b')
plt.plot(ZArray, OUT[:,1] , '--r')
plt.show()
