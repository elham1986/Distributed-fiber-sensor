# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 13:06:46 2021

General Numerical Solver for the Raman response function.Equation 2.3.36 Book: Nonlinear Fiber Optics, GOVIND P. AGRAWAL.

author: Elham RAHMATI
email: elham.rahmati@febus-optics.com

Please feel free to use and modify this, but keep the above information. Thanks!
"""

import numpy as np
from matplotlib import pyplot as plt

# parameters
tau1=12.2e-3
tau2=32e-3
tmin=0.0
tmax=1.4

# time points array
tlist=np.linspace(tmin,tmax,num=1000)
tarray=np.array(tlist)

# Equation or formula with given parameters
hR=((tau1**2 + tau2**2)/(tau1*tau2**2))*np.exp(-tarray/tau2)*np.sin(tarray/tau1)

# plot usingwith matplotlib package using pyplot
plt.plot(tarray,hR)
plt.show()