# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 23:51:41 2021

@author: Elham
"""


# libraries
import matplotlib.pyplot as plt
import numpy as np
 
# create data
x = np.linspace(-10,10,1000)
y = np.linspace(-3,3,300)
 
# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
nbins=300
X, Y = np.meshgrid(x,y)
Z =np.sin(X+Y)
 
# Make the plot
plt.pcolormesh(X, Y, Z, shading='auto')
plt.show()
