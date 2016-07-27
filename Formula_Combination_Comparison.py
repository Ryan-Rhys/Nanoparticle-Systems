# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 12:37:15 2016

@author: Ryan-Rhys
"""

import numpy as np

k = 1.38*(10**-23) 
T = 298.15
hbar = 1.0545718*(10**-34)
pi = np.pi

# Matsubara coefficient

Matsubara_constant = (2*pi*k*T)/(hbar)
rad_s = 1.05*(10**11)
eV_conversion = 6.55*(10**-5)

# Conversion from rad s^-1 to electronvolts required.

Matsubara_coeff = (Matsubara_constant/rad_s)*eV_conversion

Matsubara = lambda n: Matsubara_coeff*n

# Parameters for water found on page 266 of Parsegian's book.

Wwj = [2.07*(10**-2), 6.9*(10**-2), 9.2*(10**-2), 2.0*(10**-1), 4.2*(10**-1), 8.25, 
       10.0, 11.4, 13.0, 14.9, 18.5]
       
Wfj = [6.25*(10**-4), 3.5*(10**-3), 1.28*(10**-3), 5.44*(10**-4), 1.35*(10**-2), 
       2.68, 5.67, 12.0, 26.3, 33.8, 92.8]
       
Wgj = [1.5*(10**-2), 3.8*(10**-2), 2.8*(10**-2), 2.5*(10**-2), 5.6*(10**-2),
       0.51, 0.88, 1.54, 2.05, 2.96, 6.26]
       
# Johnson and Christy parameters for gold

wj3 = [0,3.0,4.8]
fj3 = [53.0,5.0,104.0]
gj3 = [1.8,0.8,4.4]
       
# February 5th - Formula for Epsilon_Water cross-referenced with the Hamaker script and
# appears to be fine.
       
# Again, it should be noted that the Debye term has been omitted as it is not compatible
# with the parameters for gold (there are no debye parameters for gold.)
       
# 8 February - 2.624 is the value of the square of the plasma frequency of ITO.
# 1.506 is the value of the plasma frequency of ITO.
# When the plasma frequency is decreased by an order of magnitude, the Hamaker coefficient
# decreases by 5 kT.

Epsilon_Water = lambda n: 1 + (Wfj[0])/((Wwj[0])**2 + Wgj[0]*Matsubara(n) + Matsubara(n)**2) + (Wfj[1])/((Wwj[1])**2 + Wgj[1]*Matsubara(n) + Matsubara(n)**2) + (Wfj[2])/((Wwj[2])**2 + Wgj[2]*Matsubara(n) + Matsubara(n)**2) + (Wfj[3])/((Wwj[3])**2 + Wgj[3]*Matsubara(n) + Matsubara(n)**2) + (Wfj[4])/((Wwj[4])**2 + Wgj[4]*Matsubara(n) + Matsubara(n)**2) + (Wfj[5])/((Wwj[5])**2 + Wgj[5]*Matsubara(n) + Matsubara(n)**2) + (Wfj[6])/((Wwj[6])**2 + Wgj[6]*Matsubara(n) + Matsubara(n)**2) + (Wfj[7])/((Wwj[7])**2 + Wgj[7]*Matsubara(n) + Matsubara(n)**2) + (Wfj[8])/((Wwj[8])**2 + Wgj[8]*Matsubara(n) + Matsubara(n)**2) + (Wfj[9])/((Wwj[9])**2 + Wgj[9]*Matsubara(n) + Matsubara(n)**2) + (Wfj[10])/((Wwj[10])**2 + Wgj[10]*Matsubara(n) + Matsubara(n)**2)  
Epsilon_Gold = lambda n: 1 + (fj3[0])/((wj3[0])**2 + gj3[0]*Matsubara(n) + Matsubara(n)**2) + (fj3[1])/((wj3[1])**2 + gj3[1]*Matsubara(n) + Matsubara(n)**2) + (fj3[2])/((wj3[2])**2 + gj3[2]*Matsubara(n) + Matsubara(n)**2)
Epsilon_ITO = lambda n: 1 + (2.624)/(Matsubara(n)**2 + 0.05*Matsubara(n)) + (10.919)/(23.232 + Matsubara(n)**2 - 0.82*Matsubara(n)) + (115.943)/(50.41 + Matsubara(n)**2 - 0.82*Matsubara(n))

n = 0
Hamaker = 0

while n < 10000.0:
    if n < 1.0:
        Hamaker += 0.5
        n += 1.0
    else:    
        Hamaker += ((Epsilon_Gold(n) - Epsilon_Water(n))/(Epsilon_Gold(n) + Epsilon_Water(n)))*((Epsilon_ITO(n) - Epsilon_Water(n))/(Epsilon_ITO(n) + Epsilon_Water(n)))
        n += 1.0
    
Hamaker_coeff = 1.5*Hamaker

print Hamaker_coeff

# 7 February - The value obtained for the Hamaker coefficient using the damped harmonic oscillator
# models for water and gold (omitting the second term) yields a value that is only
# ca. 7 units of KT higher than the Drude approximations.





