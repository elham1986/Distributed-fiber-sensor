# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:32:37 2021

@author: elham.rahmati
"""

import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

# Parameters
alpha_dB1 = 0.2
alpha_dB2 = 0.3
zmin = 0.0
zmax = 30000
n = 10000
Aeff = 100e-12
lambda_p = 1450e-9
lambda_s = 1550e-9
lambda0 = 1e-6
gR = (1e-13)*lambda0/(Aeff * lambda_p)
c = 3e8

# function that returns alpha_dB to alpha
def dB_To_mInv(dB):
    return (dB/20)*10**-3*np.log(10)

# Convert alpha
alpha_s = dB_To_mInv(alpha_dB1)
alpha_p = dB_To_mInv(alpha_dB2)
print(alpha_s)
print(alpha_p)
wp = 2 * np.pi * c / lambda_p
ws = 2 * np.pi * c / lambda_s

# function that returns dIs/dz  and  dIp/dz

def CoupleEq(z,Y):
    dY = np.zeros_like(Y)
    dY[0] = gR * Y[1] * Y[0] - alpha_s * Y[0]                # Y[0] means Ps
    dY[1] = -(wp / ws) * gR * Y[1] * Y[0]  - alpha_p * Y[1]  # Y[1] means Pp
    return dY

# initial condition
Y0 = np.array([10e-3,1])
B0 = 1.0

# scipy.integrate.solve_ivp(fun, t_span, y0, method='RK45', t_eval=None,
# dense_output=False, events=None, vectorized=False, args=None, **options)[source]
# It requires 3 inputs: Function, position or time points, and initial condition.
sol = solve_ivp(CoupleEq, [zmin,zmax], Y0, method='Radau', dense_output=True)

z = np.linspace(zmin, zmax, n+1) 
Y = sol.sol(z) 

Ps=Y[0,:]
Pp=Y[1,:]


# plot usingwith matplotlib package using pyplot
# setup figures
fig = plt.figure()
plt.plot(z, Y[0], 'b', label="Ps")
plt.plot(z, Y[1], '--r', label="Pp")
plt.legend(loc='best')
plt.suptitle("""Intensity of pump and signal laser """, fontweight ="bold")  
plt.xlabel("z (Propagation direction) (m)")
plt.ylabel("Power (W)")
plt.savefig('Is_Ip.jpg')
plt.show()