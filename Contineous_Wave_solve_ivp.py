# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 09:11:25 2021
General Numerical Solver for the Raman response function.Equation 2.3.36 Book: Nonlinear Fiber Optics, GOVIND P. AGRAWAL.

author: Elham RAHMATI
email: elham.rahmati@febus-optics.com

Please feel free to use and modify this, but keep the above information. Thanks!
"""

import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

# Parameters
alpha=0.05
zmin=0.0
zmax=500
#z=[zmin, zmax]
# function that returns dA/dt
def derivativeEq(z,A):
    dAdz=-(alpha/2)*A
    return dAdz

# initial condition
A0 = [1.0]
#A0=np.linspace(0.0,1.0)

# z points array
#z = np.linspace(zmin,zmax,num=1000)
z=[zmin, zmax]

# scipy.integrate.solve_ivp(fun, t_span, y0, method='RK45', t_eval=None,
# dense_output=False, events=None, vectorized=False, args=None, **options)[source]
# It requires 3 inputs: Function, position or time points, and initial condition.
A = solve_ivp(derivativeEq,z, A0)

# plot usingwith matplotlib package using pyplot
plt.plot(A.t,A.y[0,:])
plt.show()
