# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 15:39:16 2016

@author: Ryan-Rhys
"""

import numpy as np
import scipy.special as sp
from scipy.integrate import quad
import matplotlib.pyplot as plt

# 6th April 15:40 - Have copied the same constants over from the Image potential energy
# calculation script.
# Constants.

q = 1.602*(10)**-19 # Coulombs.
Epsilon_Zero = 8.854*(10)**-12 # Farads per metre.
kB = 1.38*(10)**-23 # joules per Kelvin.
NA = 6.022*(10)**23 # per mole.

# Quasi-constants

T = 298.15 # Kelvin.
Epsilon_R = 78.4 # Dimensionless.
Epsilon_Gold = 6.9
Epsilon_ITO = 3.3
Xi_Gold = (Epsilon_R/Epsilon_Gold)
Xi_ITO = (Epsilon_R/Epsilon_ITO)
Thomas_Fermi_ITO = 20.0 # Angstroms
Thomas_Fermi_Gold = 0.5 # Angstroms - for want of a better estimate at the moment.
K_ITO = 1.0/Thomas_Fermi_ITO # Angstroms^-1
K_Gold = 1.0/Thomas_Fermi_Gold # Angstroms^-1

# Ionic Strength will simply be the bulk concentration for a 1:1 electrolyte in mol/m^3.
# An example of a 1:1 electrolyte is NaCL.
# In the equation for the Debye length it should have units of mol/m^3
# in order to be dimensionally consistent with the units of the vacuum permittivity.

Bulk_Concentration = 63.0 # moles per metre cubed (or millimoles per dm^3) <- millimolar?.

L_B = (((q**2)/(np.pi*4*Epsilon_Zero*Epsilon_R*kB*T)))*10**10 # Angstroms.

# Debye_Length and Kappa in Angstroms and Angstroms^-1.
# Formula for Debye length taken from Wikipedia.

Debye_Length = (((1.0)/(np.sqrt((Bulk_Concentration*2.0*NA*(q**2))/(Epsilon_Zero*Epsilon_R*kB*T)))))*10**10
Kappa = (1.0/Debye_Length)

# 6th April - prefactor formula checked by calculator. I'm always wary of the order
# of operations with so many brackets!
# 6th April - Forgot that the formula is in Gaussian units. Hence, because my Bjerrum
# length is given in SI units I had to include a conversion factor. The prefactor is now
# just half the Bjerrum length in Angstrom.
# This should work out to be dimensionally consistent as the quantities within the
# integrand are in units of Angstroms^-1

# The prefactor is just the Bjerrum length in the case of the pair interaction.

prefactor = L_B

# Defining Geometric parameters a and z.

a = 2.0 # Angstroms
z = 0.0 # Angstroms

# 6th April 16:10 - First term checked manually using the fact that a zeroth order
# Bessel function evaluated with an argument of zero is equal to 1.
# 16:15 - Second term also checked manually.

First_Term = lambda K,p: sp.jv(0, K*p)*((K)/(np.sqrt(K**2 + Kappa**2)))*np.exp(-abs(z-a)*np.sqrt(K**2 + Kappa**2))
Second_Term = lambda K,p: sp.jv(0, K*p)*((K)/(np.sqrt(K**2 + Kappa**2)))*((Xi_Gold*np.sqrt(K**2 + Kappa**2) - np.sqrt(K**2 + K_Gold**2))/(Xi_Gold*np.sqrt(K**2 + Kappa**2) + np.sqrt(K**2 + K_Gold**2)))*np.exp(-(z+a)*np.sqrt(K**2 + Kappa**2))

rho_vals = []
Pair_Energy_vals = []

holder = 0.0

while holder < 10:
    rho_vals.append(holder)
    holder += 0.1
    
holder_2 = 0.0

while holder_2 < 10:
    result1, err = quad(First_Term, 0, np.inf, args = (holder_2,))
    result2, err = quad(Second_Term, 0, np.inf, args = (holder_2,))
    Pair_Energy_vals.append(prefactor*(result1 + result2))
    holder_2 += 0.1
    
plt.figure(1)
plt.plot(rho_vals, Pair_Energy_vals, label = 'Pair potential as a function of p')
plt.title('Pair Potential Energy as a function of cylindrical coordinate')
plt.xlabel('Distance of cylindrical coordinate p (Angstroms)')
plt.ylabel('Pair Potential Energy (kT)')












