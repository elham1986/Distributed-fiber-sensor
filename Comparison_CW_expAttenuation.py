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

# Parameters
alpha_dB=0.2
zmin=0.0
zmax=800000

# function that returns alpha_dB to alpha
def alphaDBTOalpha(alpha_dB):
    alpha1 = (alpha_dB/10)*10**-3
    return alpha1

# function that returns variable alpha
def func_alpha(alpha_dB):
  alpha = alphaDBTOalpha(alpha_dB)
  return alpha

# function that returns dA/dt
def derivativeEq(z,A):
    alpha = func_alpha(alpha_dB)
    dAdz=-(alpha/2)*A
    return dAdz

# initial condition
A0 = [1.0]
B0 = 1.0
#A0=np.linspace(0.0,1.0)

# z points array
z=[zmin, zmax]
zList=np.linspace(zmin,zmax,num=1000)
zArray = np.array(zList)

# scipy.integrate.solve_ivp(fun, t_span, y0, method='RK45', t_eval=None,
# dense_output=False, events=None, vectorized=False, args=None, **options)[source]
# It requires 3 inputs: Function, position or time points, and initial condition.
A = solve_ivp(derivativeEq,z, A0)

# call function func_alpha for using variable alpha
alpha = func_alpha(alpha_dB)

# Analytical equation of A
B = B0*np.exp((-alpha/2)*zArray)


# plot usingwith matplotlib package using pyplot
plt.plot(A.t,A.y[0,:], 'r') # plotting z, A separately 
plt.plot(zArray, B, 'b') # plotting z, B separately 
plt.show()
