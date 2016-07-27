# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 11:26:59 2015

@author: Ryan-Rhys
"""

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

# Irani parameters

wj = [0.0,2.9,4.0,8.9]
fj = [9.7,4.95,41.55,207.76]
gj = [3.21,0.67,2.22,8.50]

# DESY parameters

wj2 = [0,3.87,8.37,23.46]
fj2 = [40.11,59.61,122.55,1031.19]
gj2 = [0,2.62,6.41,27.57]

# Johnson and Christy parameters

wj3 = [0,3.0,4.8]
fj3 = [53.0,5.0,104.0]
gj3 = [1.8,0.8,4.4]

# Parameters for water found on page 266 of Parsegian's book.

Wwj = [2.07*(10**-2), 6.9*(10**-2), 9.2*(10**-2), 2.0*(10**-1), 4.2*(10**-1), 8.25, 
       10.0, 11.4, 13.0, 14.9, 18.5]
       
Wfj = [6.25*(10**-4), 3.5*(10**-3), 1.28*(10**-3), 5.44*(10**-4), 1.35*(10**-2), 
       2.68, 5.67, 12.0, 26.3, 33.8, 92.8]
       
Wgj = [1.5*(10**-2), 3.8*(10**-2), 2.8*(10**-2), 2.5*(10**-2), 5.6*(10**-2),
       0.51, 0.88, 1.54, 2.05, 2.96, 6.26]


n = 1.0
Hamaker = 0.0

while n < 10000.0:
    Gi = 0
    Wi = 0
    Gold_eps = 0.0
    Water_eps = 0.0
    numerator = 0.0
    denominator = 0.0
    
    # Drude-Lorentz formula given below as well.
    # Conditional statement. Skip the parameters marked with a dashed line
    # in Parsegian's table if n = 0, because then the term goes to infinity.
    
    #if n == 0.0:
        #Gi += 1
        #while Gi < 3:
        
           #Gold_eps += (fj3[Gi])/((wj3[Gi])**2 + (Matsubara_coeff*n)**2 + (n*(Matsubara_coeff)*gj3[Gi]) )
        
           #Gi += 1
        
        #while Wi < 11:
        
           #Water_eps += (Wfj[Wi])/((Wwj[Wi])**2 + (Matsubara_coeff*n)**2 + (n*(Matsubara_coeff)*Wgj[Wi]) )
        
           #Wi += 1
    #else:
    
    #while Gi < 1:
    
        
    #Gold_eps += 4.9752 + ((78.618)/((Matsubara_coeff*n)**2 + (0.0379*Matsubara_coeff*n))) + ((22.81)/(12.96 + (Matsubara_coeff*n)**2  + (1.3*Matsubara_coeff*n) )) + ((7.464)/(7.84 + (Matsubara_coeff*n)**2  + (0.737*Matsubara_coeff*n) ))
        
    while Gi < 3:
        
           Gold_eps += (fj3[Gi])/((wj3[Gi])**2 + (Matsubara_coeff*n)**2 + (n*(Matsubara_coeff)*gj3[Gi]) )
        
           Gi += 1
        
    while Wi < 11:
        
        Water_eps += (Wfj[Wi])/((Wwj[Wi])**2 + (Matsubara_coeff*n)**2 + (n*(Matsubara_coeff)*Wgj[Wi]) )
        
        Wi += 1
        
    numerator = Gold_eps - Water_eps
    denominator = 2 + (Gold_eps + Water_eps)
    
    # Conditional Statement for the case of the zeroth order term.
    
    if n == 0.0:
        Hamaker += 0.5*((numerator/denominator)**2)
    else:
    
        Hamaker += (numerator/denominator)**2 
    
    #print n

    n += 1.0
    
Hamaker_coefficient_J = ((3.0*k*T)/(2.0))*Hamaker
Hamaker_coefficient_ergs = (Hamaker_coefficient_J/(1*(10**-7)))
Hamaker_kT = (3.0/2.0)*Hamaker

print Hamaker_kT

# 1 erg = 1*10^-7 joules
# For comparison Hamaker coefficient values see
# Parsegian note 24 paper in MSci folder.

# 8th November - currently able to reproduce the value given in
# the Parsegian reference using the parameters of the authors, when
# n is taken to 100,000.

# 13th November - The Drude-Loretnz calculation doesn't converge as I think
# that the two expressions are incompatible.

# 23rd November - Conditional case for the zeroth order (n = 0) term added.

# 26th January - The limiting factor in calculating a Hamaker coefficient
# is the fit to the absorption data. This dwarfs the influence of any contribution
# from zero-frequency screening.

# 5th February - I've just realised why I left out the first term in the damped oscillator
# model for water. I don't have similar terms available for gold
# in Parsegian's book.

