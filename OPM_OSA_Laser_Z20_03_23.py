# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 11:59:39 2021

@author: elham.rahmati
"""

import numpy as np
from matplotlib import pyplot as plt

opmArray = np.array([1.26, 1.58, 2.0, 2.52, 3.15, 3.82, 4.33, 4.60])  #  We put P2_OPM column, this connectioin is E2000, just before the FC connector
osaArray = np.array([0.675, 0.845, 1.07, 1.35, 1.67, 2.04, 2.46, 2.60]) #  We put P3_OSA column, just zafter the FC connector

fit = np.polyfit(opmArray[0:5], osaArray[0:5], 1)
f = np.poly1d(fit)
fit_osaArray = f(opmArray[0:5])

# plot usingwith matplotlib package using pyplot
# setup figures
fig = plt.figure(1)
plt.plot(opmArray[0:5],osaArray[0:5], '.b', label="Ps")
plt.plot(opmArray[0:5],fit_osaArray[0:5], 'r', label='({:.1E}) + ({:.1E}) x'.format(fit[1], fit[0]))
plt.legend(loc='best')  
plt.ylabel("OSA (mW)")
plt.xlabel("OPM (mW)")
plt.savefig('OSA_As_OPMLaser_Z20.jpg')
plt.show()