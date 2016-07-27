# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 13:50:20 2015

@author: Ryan-Rhys
"""

import numpy as np
from matplotlib import pyplot as plt

# For most of the epxressions I have multiplied by factors in order to achieve
# the right dimensionality.

# Constants.

q = 1.602*(10)**-19 # Coulombs.
Epsilon_Zero = 8.854*(10)**-12 # Farads per metre.
kB = 1.38*(10)**-23 # joules per Kelvin.
NA = 6.022*(10)**23 # per mole.

# Quasi-constants

T = 298.15 # Kelvin.
Epsilon_R = 78.4 # Dimensionless.

# Ionic Strength will simply be the bulk concentration for a 1:1 electrolyte in mol/m^3.
# In the equation for the Debye length it should have units of mol/m^3
# in order to be dimensionally consistent with the units of the vacuum permittivity.
# Radius in nanometres.

Bulk_Concentration = 10.0 # moles per metre cubed (or millimoles per dm^3).
Radius = 20.0 # nanometres.

# surface charge density in Coulombs per square nanometre and Bjerrum length in nanometres.

no_charges = 200.0 # Dimensionless
sigma = ((no_charges*q)/(4*np.pi*Radius**2)) # Coulombs per square nanometre.
L_B = (((q**2)/(np.pi*4*Epsilon_Zero*Epsilon_R*kB*T)))*10**9 # nanometres.

# Debye_Length and Kappa in nanometres and nanometres^-1.
# Formula for Debye length taken from Wikipedia.

Debye_Length = (((1.0)/(np.sqrt((Bulk_Concentration*2.0*NA*(q**2))/(Epsilon_Zero*Epsilon_R*kB*T)))))*10**9
Kappa = (1.0/Debye_Length)

# sigma_star in Coulombs per square nanometre.
# 16th December - but will have dimensionality of moles as well! Source of error?

sigma_star = ((np.sqrt((8.0*Bulk_Concentration*Epsilon_Zero*Epsilon_R*NA*kB*T)/(q**2)))*10**-18)*q

# Zed in units of kT per nanometre.
# 16th December - I have a suspicion that something is wrong with Zed.
# I have introduced NA, Avogadro's number in the numerator in order to 
# make the units match.
# My order of magnitude seems to be completely off though.
# 12.30 pm update, problem seems to lie with the treatment of sigma and sigma star
# as either Coulombs per square nanometre or as area per charge. 

Zed = ((((64.0*(Bulk_Concentration*10**-27)*NA)/((L_B**2)*np.pi*(sigma**2)))*(1.0 - np.sqrt(1.0 + (sigma/sigma_star)**2) + (0.5*(sigma/sigma_star)**2))))*q**2

# Constants required for the VdW potential calculation.
# The Hamaker constant 'Hamaker' is taken from Johnson and Christie's data.
# This may be found in either Parsegian's book ca. pg. 66
# or from the his 1972 paper where he computes 3 different Hamaker constants
# from spectral data.

# Hamaker = 65.0 is that of bulk gold given in the 2015 paper from Iowa University.
# The VdW potential graph is very similar to that obtained by them.

Hamaker = 74.0
Hamaker_vals = []
PP_Vals = []


# Put xmin and xmax in nanometres.
# Initialising list of x values for the plots.
# L_vals - shorthand for values of the surface to surface separation in nanometres.

xmin = (0.25*Radius)
xmax = (2.0*Radius)
L_vals = np.arange(xmin,xmax,0.001)

# 26th February - Introucing extra line to compare against Silbey results.
#L_vals = np.arange(0.000,12.500,0.001)


# Initialising empty lists for storing the y values for the plots.
# These lists will have energy values appended to them in each iteration of the loop below.

Coulomb_Vals_1 = []
Coulomb_Vals_11 = []
Coulomb_Vals_111 = []
Coulomb_Vals_1111 = []

Coulomb_Vals_2 = []

Combined_Vals1 = []
Combined_Vals11 = []
Combined_Vals111 = []
Combined_Vals1111 = []
Combined_Vals2 = []
Combined_ValsPP = []


# Computing energy values for surface to surface separations in the range
# xmin - xmax.

for L in np.arange(xmin,xmax,0.001):
    
    U_1 = ((Radius)/(2))*(Zed)*(np.exp(-Kappa*L))
    Coulomb_Vals_1.append(U_1)
    
    Z = L + (2*Radius)
    U_h_sphere = (-Hamaker/3)*((Radius**2)/(Z**2 - 4*Radius**2)+(Radius**2)/(Z**2)+(0.5*np.log(1-(4*Radius**2/Z**2))))
    Combined_Vals1.append(U_1 + U_h_sphere)
    
no_charges = 500.0
Bulk_Concentration = 10.0
    
sigma = ((no_charges*q)/(4*np.pi*Radius**2))
Debye_Length = (((1.0)/(np.sqrt((Bulk_Concentration*2.0*NA*(q**2))/(Epsilon_Zero*Epsilon_R*kB*T)))))*10**9
Kappa = (1.0/Debye_Length)
sigma_star = ((np.sqrt((8.0*Bulk_Concentration*Epsilon_Zero*Epsilon_R*NA*kB*T)/(q**2)))*10**-18)*q
Zed = ((((64.0*(Bulk_Concentration*10**-27)*NA)/((L_B**2)*np.pi*(sigma**2)))*(1.0 - np.sqrt(1.0 + (sigma/sigma_star)**2) + (0.5*(sigma/sigma_star)**2))))*q**2
    
for L in np.arange(xmin,xmax,0.001):
    
    U_11 = ((Radius)/(2))*(Zed)*(np.exp(-Kappa*L))
    Coulomb_Vals_11.append(U_11)
    
    Z = L + (2*Radius)
    U_h_sphere = (-Hamaker/3)*((Radius**2)/(Z**2 - 4*Radius**2)+(Radius**2)/(Z**2)+(0.5*np.log(1-(4*Radius**2/Z**2))))
    Combined_Vals11.append(U_11 + U_h_sphere)
    
no_charges = 1000.0
Bulk_Concentration = 10.0
    
sigma = ((no_charges*q)/(4*np.pi*Radius**2))
Debye_Length = (((1.0)/(np.sqrt((Bulk_Concentration*2.0*NA*(q**2))/(Epsilon_Zero*Epsilon_R*kB*T)))))*10**9
Kappa = (1.0/Debye_Length)
sigma_star = ((np.sqrt((8.0*Bulk_Concentration*Epsilon_Zero*Epsilon_R*NA*kB*T)/(q**2)))*10**-18)*q
Zed = ((((64.0*(Bulk_Concentration*10**-27)*NA)/((L_B**2)*np.pi*(sigma**2)))*(1.0 - np.sqrt(1.0 + (sigma/sigma_star)**2) + (0.5*(sigma/sigma_star)**2))))*q**2

for L in np.arange(xmin,xmax,0.001):
    
    U_111 = ((Radius)/(2))*(Zed)*(np.exp(-Kappa*L))
    Coulomb_Vals_111.append(U_111)
    
    Z = L + (2*Radius)
    U_h_sphere = (-Hamaker/3)*((Radius**2)/(Z**2 - 4*Radius**2)+(Radius**2)/(Z**2)+(0.5*np.log(1-(4*Radius**2/Z**2))))
    Combined_Vals111.append(U_111 + U_h_sphere)
    
no_charges = 2000
Bulk_Concentration = 10.0
Radius = 20.0
    
sigma = ((no_charges*q)/(4*np.pi*Radius**2))
Debye_Length = (((1.0)/(np.sqrt((Bulk_Concentration*2.0*NA*(q**2))/(Epsilon_Zero*Epsilon_R*kB*T)))))*10**9
Kappa = (1.0/Debye_Length)
sigma_star = ((np.sqrt((8.0*Bulk_Concentration*Epsilon_Zero*Epsilon_R*NA*kB*T)/(q**2)))*10**-18)*q
Zed = ((((64.0*(Bulk_Concentration*10**-27)*NA)/((L_B**2)*np.pi*(sigma**2)))*(1.0 - np.sqrt(1.0 + (sigma/sigma_star)**2) + (0.5*(sigma/sigma_star)**2))))*q**2

for L in np.arange(xmin,xmax,0.001):
    
    U_1111 = ((Radius)/(2))*(Zed)*(np.exp(-Kappa*L))
    Coulomb_Vals_1111.append(U_1111)
     
    U_2 = ((sigma**2)/(q**2))*((16*((np.pi)**2)*Radius**2)/(Kappa**2))*(L_B/L)*((np.sinh(Kappa*Radius))**2)*(np.exp(-Kappa*L))
    Coulomb_Vals_2.append(U_2)
    
    # Setting Z, the centre-centre distance to be L + 2R, where L is the surface-surface
    # separation.
    
    Z = L + (2*Radius)

    U_PP = ((100**2)*L_B*np.exp(-Kappa*L))/(Z*(1 + Kappa*Radius))
    
    # U_PP is presumably the electrostatic potential of a point particle. - 26th Jan
    
    U_h_sphere = (-Hamaker/3)*((Radius**2)/(Z**2 - 4*Radius**2)+(Radius**2)/(Z**2)+(0.5*np.log(1-(4*Radius**2/Z**2))))
    
    U_h_planar = (-Hamaker/6)*(Radius/L + Radius/(2*Radius + L) + np.log(L/(2*Radius + L)))
    
    # 12th January - Formula for U_h_sphere checked manually and is fine.
    # 26th January - Formula for U_h_planar checked manually and is fine.
    
    Combined_Vals1111.append(U_1111 + U_h_sphere)
    Combined_Vals2.append(U_2 + U_h_sphere)
    Combined_ValsPP.append(U_PP + U_h_sphere)
    
    # 26th January - for Hamaker vals can either append sphere-sphere or sphere-plane
    # value.
    
    Hamaker_vals.append(U_h_sphere)
    PP_Vals.append(U_PP)

# Code for plotting the Coulomb Potential.

plt.figure(1)
plt.plot(L_vals, Coulomb_Vals_1, label = '200 Charges')
plt.plot(L_vals, Coulomb_Vals_11, label = '500 Charges')
plt.plot(L_vals, Coulomb_Vals_111, label = '1000 Charges')
plt.plot(L_vals, Coulomb_Vals_1111, label = '2000 Charges')
plt.xlabel('Surface to Surface Separation (nm)')
plt.ylabel('Coulomb Potential Energy (kT)')
plt.title('Coulomb Potential Energy as a Function of Separation', loc = 'right')
plt.xlim(xmin, xmax)

plt.legend(loc=9, bbox_to_anchor=(1.2, 0.65))
#plt.text(3,25, 'Radius = 6.5 nm')
#plt.text(3,20, 'NaCl Concentration = 63 mM')
#plt.text(3,15, 'Surface Charge Density = 400')

#plt.figure(2)
#plt.plot(L_vals, Hamaker_vals, label = 'Damped Oscillator and Drude-Lorentz')
#plt.xlabel('Surface to Surface Separation (nm)')
#plt.ylabel('VdW Potential (kT)')
#plt.title('VdW Potential as a Function of Separation', loc = 'right')
#plt.xlim(xmin, xmax)
#plt.text(1.7,-12, 'Radius = 6.5 nm')
#plt.text(1.7,-14, 'ITO planar surface and Gold NP')
#plt.text(1.7,-10, 'Different approximation for Dielectric Functions')


#plt.figure(2)
#plt.plot(L_vals, Coulomb_Vals_1)
#plt.xlabel('Surface to Surface Separation (nm)')
#plt.ylabel('Coulomb Potential (kT)')
#plt.title('Coulomb Potential as a Function of Separation', loc = 'right')
#plt.xlim(xmin, xmax)

# Code for plotting the combined interaction potential.

plt.figure(3)
plt.title('Total Interaction Potential as a Function of Separation (DESY)', y = 1.05)
#plt.plot(L_vals, Combined_Vals1, label = '200 Charges')
#plt.plot(L_vals, Combined_Vals11, label = '500 Charges')
plt.plot(L_vals, Combined_Vals111, label = '1000 Charges')
#plt.plot(L_vals, Combined_Vals1111, label = '2000 Charges')
#plt.xlim(xmin, xmax)
plt.ylim(-0.5, 0.5)
plt.xlabel('Surface to Surface Separation (nm)')
plt.ylabel('Coulomb Potential Energy (kT)')
plt.xlim(12.0, 32.0)
plt.legend(loc=9, bbox_to_anchor=(1.2, 0.65))
#plt.plot(L_vals, Combined_ValsPP)
#plt.xlabel('L (nm)')
#plt.ylabel('Interaction Potential (kT)')
#plt.title('Interaction Potential as a Function of Separation')
#plt.xlim(xmin, xmax)
#plt.text(10,-0.12, 'Radius = 6.5 nm')
#plt.text(10,-0.14, 'NaCl Concentration = 63 mM')
#plt.text(10,-0.16, 'Surface Charge Density = 400')
#plt.text(97,0.006,'Radius = ' + str(Radius) + ' nm')
#plt.text(97,0.004, 'charge number = ' + str(no_charges))
#plt.text(97,0.002, 'NaCl Concentration = ' + str(Bulk_Concentration/1000) + ' M')

#plt.figure(4)
#plt.plot(L_vals, Combined_Vals1)
#plt.xlabel('Surface to Surface Separation (nm)')
#plt.ylabel('Interaction Potential (kT)')
#plt.title('Interaction Potential as a Function of Separation', loc = 'right')
#plt.xlim(xmin, xmax)
#plt.text(30,0.008,'R = ' + str(Radius) + ' nm')
#plt.text(30,0.006, 'charge number = ' + str(no_charges))
#plt.text(30,0.004, 'NaCl Concentration = ' + str(Bulk_Concentration/1000) + ' M')

# 16th December - There is a bug in the code which I'm currently trying to fix.
# The order of magnitude of Israelachvili's Coulomb potential is wrong somehow.

#plt.figure(4)
#plt.plot(L_vals, PP_Vals)  

# 17th December - The point charge limit, now incorporated into this script drops
# off incredibly quickly.

# 13th January - I've now tried both Coulomb expressions. It would seem that the 
# one corresponding to that given by Israelachvili produces the correct behaviour
# when compared to the 2015 paper. 

# 22nd April - This is an important script. It would be good to plot the expressions in MathCad
# as well as here in order to troubleshoot the error with one of the expressions. Will have to
# dig out the derivation for this expression as well and write it up.


   
    

