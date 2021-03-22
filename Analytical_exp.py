# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 10:29:30 2021
General Numerical Solver for the comparison between Continuous Wave with numerical
calculation and the result of analytical equation.

Using the two first terms of Eq 2.3.41 in Book: Nonlinear Fiber Optics, GOVIND P. AGRAWAL.


author: Elham RAHMATI
email: elham.rahmati@febus-optics.com

Please feel free to use and modify this, but keep the above information. Thanks!
"""
import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

# Parameters
alpha=0.05
A0=1.0
zmin=0.0
zmax=1000
z=[zmin, zmax]
zList=np.linspace(zmin,zmax,num=1000)

#zfloat = float('z')
zArray = np.array(zList)
# Analytical equation of A
A = A0*np.exp((-alpha/2)*zArray)

# plot usingwith matplotlib package using pyplot
plt.plot(zArray,A)
plt.show()
