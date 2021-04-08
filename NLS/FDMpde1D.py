# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 22:36:03 2021

@author: Elham
"""


import numpy as np
import findiff as fdm
import matplotlib.pyplot as plt

fdm_acc = 2

nt = 100
t = np.linspace(0, 10, nt)
dt = t[1]-t[0]

L = fdm.FinDiff(0, dt, 1, acc=fdm_acc)
f = np.cos(t)

bc = fdm.BoundaryConditions(f.shape)
bc[0] = 1

pde = fdm.PDE(L, f.copy(), bc)
u = pde.solve()


plt.plot(t, f)
plt.plot(t, u)
plt.show()