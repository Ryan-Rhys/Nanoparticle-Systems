# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 21:59:22 2016

@author: Ryan-Rhys
"""

import numpy as np
import matplotlib.pyplot as plt

k = 1.38*(10**-23) 
T = 273.0
hbar = 1.0545718*(10**-34)
pi = np.pi

# Matsubara coefficient

Matsubara = (2*pi*k*T)/(hbar)
rad_s = 1.05*(10**11)
eV_conversion = 6.55*(10**-5)

# Conversion from rad s^-1 to electronvolts required.

Matsubara_coeff = lambda n: ((Matsubara/rad_s)*eV_conversion)*n

# fj is dimensionless - Pendry parameters

fj_1 = 0.2975
fj_2 = 0.6061
fj_3 = 1.4782
fj_4 = 0.9284
fj_5 = 0.8600
fj_6 = 0.4280
fj_7 = 0.5323
fj_8 = 0.3951
fj_9 = 0.1167
fj_10 = 0.3728

# wj in electronvolts

wj_0 = 8.7663
wj_1 = 2.7260
wj_2 = 3.1091
wj_3 = 3.9415
wj_4 = 5.2526
wj_5 = 7.9749
wj_6 = 10.4650
wj_7 = 14.1798
wj_8 = 21.2070
wj_9 = 30.0898
wj_10 = 42.9033

# gj in electronvolts

gj_0 = 0.039299
gj_1 = 0.3274
gj_2 = 0.6066
gj_3 = 1.3318
gj_4 = 2.3931
gj_5 = 3.7913
gj_6 = 4.9667
gj_7 = 7.8967
gj_8 = 6.2030
gj_9 = 6.9297
gj_10 = 43.5640



Epsilon_Pendry = lambda n: 1 + (wj_0**2)/((Matsubara_coeff(n)**2) + (Matsubara_coeff(n)*gj_0)) + (fj_1*(wj_1**2))/((-Matsubara_coeff(n)**2) - (wj_1**2) - (gj_1*Matsubara_coeff(n))) + (fj_2*(wj_2**2))/((-Matsubara_coeff(n)**2) - (wj_2**2) - (gj_2*Matsubara_coeff(n)))    + (fj_3*(wj_3**2))/((-Matsubara_coeff(n)**2) - (wj_3**2) - (gj_3*Matsubara_coeff(n)))    + (fj_4*(wj_4**2))/((-Matsubara_coeff(n)**2) - (wj_4**2) - (gj_4*Matsubara_coeff(n)))    + (fj_5*(wj_5**2))/((-Matsubara_coeff(n)**2) - (wj_5**2) - (gj_5*Matsubara_coeff(n)))    + (fj_6*(wj_6**2))/((-Matsubara_coeff(n)**2) - (wj_6**2) - (gj_6*Matsubara_coeff(n)))    + (fj_7*(wj_7**2))/((-Matsubara_coeff(n)**2) - (wj_7**2) - (gj_7*Matsubara_coeff(n)))    + (fj_8*(wj_8**2))/((-Matsubara_coeff(n)**2) - (wj_8**2) - (gj_8*Matsubara_coeff(n)))    + (fj_9*(wj_9**2))/((-Matsubara_coeff(n)**2) - (wj_9**2) - (gj_9*Matsubara_coeff(n)))    + (fj_10*(wj_10**2))/((-Matsubara_coeff(n)**2) - (wj_10**2) - (gj_10*Matsubara_coeff(n)))

print Epsilon_Pendry(1)   

                