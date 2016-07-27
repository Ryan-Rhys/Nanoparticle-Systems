# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 16:56:54 2016

@author: Ryan-Rhys
"""

import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt

k = 1.38*(10**-23) 
T = 0.01
hbar = 1.0545718*(10**-34)
pi = np.pi

omega_star = (k*T)/(hbar)

wj2 = [0.0, 3.87, 8.37, 23.46]
fj2 = [40.11, 59.61, 122.55, 1031.19]
gj2 = [0, 2.62, 6.41, 27.57]
    
wj2_Hz = [0.0, 935766000000000.1, 2023866000000000.0, 5672628000000000.0]
fj2_Hz = [2.3451209964e+30, 3.4852321764e+30, 7.165160262e+30, 6.02908332156e+31]
gj2_Hz = [0.0, 633516000000000.1, 1549938000000000.2, 6666426000000000.0]

Reduced_wj2 = []
Reduced_fj2 = []
Reduced_gj2 = []

for element in wj2_Hz:
    Reduced_wj2.append((element)/(omega_star))
    
    
for element in fj2_Hz:
    Reduced_fj2.append((element)/(omega_star**2))
    
    
for element in gj2_Hz:
    Reduced_gj2.append((element)/(omega_star))
    
Wwj = [5005260000000.0,
 16684200000000.002,
 22245600000000.0,
 48360000000000.0,
 101556000000000.03,
 1994850000000000.2,
 2418000000000000.0,
 2756520000000000.0,
 3143400000000000.0,
 3602820000000000.5,
 4473300000000000.5]
 
Wfj = [3.6542025e+25,
 2.0463534e+26,
 7.48380672e+25,
 3.1806178560000008e+25,
 7.893077400000002e+26,
 1.566922032e+29,
 3.315092508e+29,
 7.0160688e+29,
 1.537688412e+30,
 1.9761927119999998e+30,
 5.425759872e+30]
 
Wgj = [3627000000000.0005,
 9188400000000.0,
 6770400000000.0,
 6045000000000.0,
 13540800000000.0,
 123318000000000.02,
 212784000000000.0,
 372372000000000.0,
 495690000000000.0,
 715728000000000.0,
 1513668000000000.0]
 
Reduced_Wwj = []
Reduced_Wfj = []
Reduced_Wgj = []

for element in Wwj:
    Reduced_Wwj.append((element)/(omega_star))
    
for element in Wfj:
    Reduced_Wfj.append((element)/(omega_star**2))
    
for element in Wgj:
    Reduced_Wgj.append((element)/(omega_star))
