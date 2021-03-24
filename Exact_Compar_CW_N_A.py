# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 11:02:26 2021
General Numerical Solver for the comparison between Continuous Wave with numerical
calculation and the result of analytical equation.

The two first terms of Eq 2.3.41 in Book: Nonlinear Fiber Optics, GOVIND P. AGRAWAL.

author: Elham RAHMATI
email: elham.rahmati@febus-optics.com

Please feel free to use and modify this, but keep the above information. Thanks!
"""

import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

# Parameters
alpha_dB=0.2
zmin=0.0
zmax=800000
n=100

# function that returns alpha_dB to alpha
def dB_To_minv(dB):
    return (dB/20)*10**-3*np.log(10)

# Convert alpha
alpha = dB_To_minv(alpha_dB)

# function that returns dA/dt
def derivativeEq(z,A):
    dAdz=-(alpha/2)*A
    return dAdz

# initial condition
A0 = [1.0]
B0 = 1.0

# z points array
z=[zmin, zmax]
zList=np.linspace(zmin,zmax,num=n)
zArray = np.array(zList)

# scipy.integrate.solve_ivp(fun, t_span, y0, method='RK45', t_eval=None,
# dense_output=False, events=None, vectorized=False, args=None, **options)[source]
# It requires 3 inputs: Function, position or time points, and initial condition.
sol = solve_ivp(derivativeEq,z, A0,dense_output=True)
A = sol.sol(zArray)

# Analytical equation of A
B = B0*np.exp((-alpha/2)*zArray)

# plot usingwith matplotlib package using pyplot
# setup figures
fig = plt.figure()
plt.plot(zArray/10**5, B, 'b', label="Analytical") # plotting z, B separately 
plt.plot(zArray/10**5,A[0,:], '--r', label="Numerical") # plotting z, A separately 
plt.legend()
plt.suptitle("""Attenuation of Amplitude during propagation """, fontweight ="bold")  
plt.xlabel("z (Propagation direction) (km)")
plt.ylabel("Amplitude")
plt.savefig('Amplitude_attenuation.jpg')
plt.show()



