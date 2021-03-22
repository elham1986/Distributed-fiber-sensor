# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 08:31:55 2021
General Numerical Solver for the Raman response function.Equation 2.3.36 Book: Nonlinear Fiber Optics, GOVIND P. AGRAWAL.

author: Elham RAHMATI
email: elham.rahmati@febus-optics.com

Please feel free to use and modify this, but keep the above information. Thanks!
"""

import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt

# Parameters
alpha=0.05

# function that returns dA/dt
def derivativeEq(A,z):
    dAdz=-(alpha/2)*A
    return dAdz

# initial condition
A0 = 1.0

# z points array
z = np.linspace(0.0,100,num=100)

# solve odeint: differential equations are solved in Python with Scipy.integrate package using ODEINT,
# It requires 3 inputs: Function, initial condition, and position or time points.
A = odeint(derivativeEq,A0,z)

# plot usingwith matplotlib package using pyplot
plt.plot(z,A)
plt.show()