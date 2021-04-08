# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 14:48:05 2021

@author: elham.rahmati
"""

import numpy as np
from matplotlib import pyplot as plt

opmArray = np.array([32.03, 40.35, 50.83])
osaArray = np.array([13.27, 16.51, 20.61])

fit = np.polyfit(opmArray, osaArray, 1)
f = np.poly1d(fit)
fit_osaArray = f(opmArray)

# plot usingwith matplotlib package using pyplot
# setup figures
fig = plt.figure(1)
plt.plot(opmArray,osaArray, '.b', label="Ps")
plt.plot(opmArray,fit_osaArray, 'r', label='({:.1E}) + ({:.1E}) x'.format(fit[1], fit[0]))
plt.legend(loc='best')  
plt.ylabel("OSA (mW)")
plt.xlabel("OPM (mW)")
plt.savefig('OSA_As_OPMPump_Z20.jpg')
plt.show()