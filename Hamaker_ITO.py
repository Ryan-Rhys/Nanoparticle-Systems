# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 14:41:48 2016

@author: Ryan-Rhys
"""

# Script to calcluate the Hamaker coefficient using Drude-like models for ITO and Gold.
# 1 February - One extra problem that should be noted is that in using the Drude model
# the zero-frequency term becomes difficult to introduce.

import numpy as np

k = 1.38*(10**-23) 
T = 298.15
hbar = 1.0545718*(10**-34)
pi = np.pi

# Matsubara coefficient

Matsubara = (2*pi*k*T)/(hbar)
rad_s = 1.05*(10**11)
eV_conversion = 6.55*(10**-5)

# Conversion from rad s^-1 to electronvolts required.

Matsubara_coeff = (Matsubara/rad_s)*eV_conversion

n = 1.0
Hamaker = 0.0

while n < 10000.0:
    if n < 1:
        Hamaker += 0.5
        n += 1
    #elif n < 100:
        #eps_gold = 4.9752 + (78.618)/((Matsubara_coeff*n)**2 + (0.0379*Matsubara_coeff*n))
        #eps_water = 1.78
        #eps_gold_2 = 4.9752 + (78.618)/((Matsubara_coeff*n)**2 + (0.0379*Matsubara_coeff*n))
        #eps_ITO = 3.8375 + (2.4983)/((Matsubara_coeff*n)**2 + (0.1764*Matsubara_coeff*n))
        #Hamaker += ((eps_gold - eps_water)/(eps_gold + eps_water))*((eps_gold_2 - eps_water)/(eps_gold_2 + eps_water))
        #n += 1
    else:
        eps_gold = 1.0 + (78.618)/((Matsubara_coeff*n)**2 + (0.0379*Matsubara_coeff*n))
        eps_water = 1.0
        #eps_gold_2 = 1.0 + (78.618)/((Matsubara_coeff*n)**2 + (0.0379*Matsubara_coeff*n))
        eps_ITO = 1.0 + (2.4983)/((Matsubara_coeff*n)**2 + (0.1764*Matsubara_coeff*n))
        Hamaker += ((eps_gold - eps_water)/(eps_gold + eps_water))*((eps_ITO - eps_water)/(eps_ITO + eps_water))
        n += 1

# 2 February - Addition made to the above formula. Calculating dielectric permittivity
# using Drude models for a certain range of frequencies and using Drude models
# with the first term equal to 1 after a certain range of frequencies.  
# Can test for the case of gold and gold to see what values of n
# to truncate are reasonable.   

Hamaker_coeff = (3.0/2.0)*Hamaker

print Hamaker_coeff

# 7 February - Hamaker coefficient compared for two different types of approximation
# in the Formula Combination Comparison Script, a value of ca. 20 was obtained for the
# Hamaker coefficient.