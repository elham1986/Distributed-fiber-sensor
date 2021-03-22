# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 20:27:09 2021

@author: Elham
"""

import numpy as np
from matplotlib import pyplot as plt


d = 0.4       # in micro meter
n0 = 1.0      # Refractive index
n1 = 2.0      # Refractive index
n2 = 1.5      # Refractive index
N = 50       
lamMin = 0.2
lamMax = 2
nlam = 10000


def Kz(lam,i):
    return n(i)*np.pi/lam

def n(i):
    if (i == 0):
        return n0
    res = i%2
    if (res == 0):
        nOut = n2
    elif (res == 1):
        nOut = n1
    return nOut

def t(i,j):
    return (n(i)*2)/(n(i)+n(j))

def r(i,j):
    return (n(i)-n(j))/(n(i)+n(j))
    
def Tj(lam,j):
    return np.array([[np.exp(-1j*Kz(lam,j)*d), 0],[0, np.exp(1j*Kz(lam,j)*d)]])

def Tij(i,j):
    return np.array([[1, r(i,j)],[r(i,j), 1]])/t(i,j)

def Tau(lam):
    Tout = Tij(0,1)
    for i in range(1,N+1):
        S = np.matmul(Tj(lam,i),Tij(i,i+1))
        Tout = np.matmul(Tout, S)
    return Tout

def ref(lam):
    tau = Tau(lam)
    return abs(tau[1,0]/tau[0,0])

lamArray = np.linspace(lamMin, lamMax, nlam)
rArray = np.zeros_like(lamArray)

for i in range(0,nlam):
    lam = lamArray[i]
    rArray[i] = ref(lam)

plt.plot(lamArray, rArray)