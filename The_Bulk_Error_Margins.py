# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 18:25:56 2016

@author: Ryan-Rhys
"""

import numpy as np
import matplotlib.pyplot as plt


# Hamaker coefficient values taken from Parsegian and Weiss 1981.
# The values given in table IId in this paper were converted to units of kT.
# For reference: 1 erg = 1*10^-7 joules
# 4.04*10^-21 joules = 1 kT at room temperature.

# Hamaker coefficients for two gold nanospheres across water, no electrolyte.

DESY_Hamaker = 74.26
JC_Hamaker = 24.75

# Nanoparticle Radius in nanometres.

NP_Radius1 = 10.0
NP_Radius2 = 20.0
NP_Radius3 = 40.0


sep_vals = []

# Collecting the x-axis list (separation between nanoparticles in nm).
# Using largest surface to surface distance as 6 times the NP radius.
    
energy_vals_DESY1 = []
energy_vals_DESY2 = []
energy_vals_DESY3 = []

#energy_vals_JC = []

# Employing the centre-centre distance (holder_2) in the energy equation.

holder1 = (2*NP_Radius1) + 5.0
holder2 = (2*NP_Radius2) + 5.0
holder3 = (2*NP_Radius3) + 5.0

while holder1 < (6*NP_Radius1):
    
    sep_vals.append(holder1 - (2*NP_Radius1))
    # Formula below taken from Parsegian's 2006 book, page 155.
    energy_vals_DESY1.append((-DESY_Hamaker/3.0)*((NP_Radius1**2)/(holder1**2 - 4*NP_Radius1**2)+(NP_Radius1**2)/(holder1**2)+(0.5*np.log(1-(4*NP_Radius1**2/holder1**2)))))
    energy_vals_DESY2.append((-DESY_Hamaker/3.0)*((NP_Radius2**2)/(holder2**2 - 4*NP_Radius2**2)+(NP_Radius2**2)/(holder2**2)+(0.5*np.log(1-(4*NP_Radius2**2/holder2**2)))))
    energy_vals_DESY3.append((-DESY_Hamaker/3.0)*((NP_Radius3**2)/(holder3**2 - 4*NP_Radius3**2)+(NP_Radius3**2)/(holder3**2)+(0.5*np.log(1-(4*NP_Radius3**2/holder3**2)))))
    
    #energy_vals_JC.append((-JC_Hamaker/3.0)*((NP_Radius**2)/(holder**2 - 4*NP_Radius**2)+(NP_Radius**2)/(holder**2)+(0.5*np.log(1-(4*NP_Radius**2/holder**2)))))
    
    holder1 += 0.1
    holder2 += 0.1
    holder3 += 0.1
    
#energy_vals_DESY2 = []
#sep_vals = []
#holder7 = (2*NP_Radius2) + 0.1

#while holder7 < (6*NP_Radius1):
    
    #sep_vals.append(holder7 - (2*NP_Radius2))
    #energy_vals_DESY2.append((-DESY_Hamaker/3.0)*((NP_Radius2**2)/(holder7**2 - 4*NP_Radius2**2)+(NP_Radius2**2)/(holder7**2)+(0.5*np.log(1-(4*NP_Radius2**2/holder7**2)))))
    
    #holder7 += 0.1
    
plt.figure(1)
plt.title('Identical Gold NPs - Size Comparison')
plt.plot(sep_vals, energy_vals_DESY1, label = 'Radius of 10 nm')
plt.plot(sep_vals, energy_vals_DESY2, label = 'Radius of 20 nm')
plt.plot(sep_vals, energy_vals_DESY3, label = 'Radius of 40 nm')
#plt.plot(sep_vals, energy_vals_JC, label = 'Johnson and Christy')
plt.xlabel('Surface to Surface Separation (nm)')
plt.ylabel('Van der Waals Potential (kT)')
plt.legend(loc=9, bbox_to_anchor=(0.5, -0.2))
plt.xlim(5.0, 40.0)

