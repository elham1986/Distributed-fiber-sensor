# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 22:37:04 2021

@author: Elham
"""


import numpy as np
import findiff as fdm
import matplotlib.pyplot as plt

fdm_acc = 2

z = np.linspace(0, 1, 20)
dz = z[1]-z[0]

t = np.linspace(0, 1, 30)
dt = t[1]-t[0]

Z, T = np.meshgrid(z, t, indexing='ij')
Fzt = Z - T

L = fdm.FinDiff(0, dz, 1, acc=fdm_acc) * fdm.FinDiff(1, dt, 1, acc=fdm_acc)

bc = fdm.BoundaryConditions(Fzt.shape)
bc[0, :] = 0
bc[:, 0] = 0

pde = fdm.PDE(L, Fzt.copy(), bc)
u = pde.solve()


plt.plot(t, u[-1,:])
plt.show()