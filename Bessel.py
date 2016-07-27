# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 13:49:05 2016

@author: Ryan-Rhys
"""

import numpy as np
import scipy.special as sp
from scipy.integrate import quad
import matplotlib.pyplot as plt

a = 2

Bessel = lambda(x): sp.jv(0,x)

x_vals = range(0,100000)
Bessel_vals = []

for x in range(0,100000):
    Bessel_vals.append(sp.jv(0,x))
    
plt.figure(1)
plt.plot(x_vals, Bessel_vals)


def Integrand_Bessel1(y):
    return (sp.jv(0,y))
    
Bessel1, err = quad(Integrand_Bessel1, 0, np.inf)

