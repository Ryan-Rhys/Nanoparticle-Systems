# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 13:56:11 2015

@author: Ryan-Rhys
"""
import numpy
from numpy import *
import matplotlib
from matplotlib import pyplot as plt

# This script works. I have tested it using the Drude-Lorentz model proposed in the memo
# and have reproduced the same graphs as in the Mathcad file there.
# It is important in Python to avoid rounding by specifying numbers as floats!
# I noticed this when I seemed to be getting very square graphs!

lamda = 300.0

# Setting lamda = 300 nm, the beginning of the visible range of the spectrum.

# Omega has been replaced with (1240/lamda) at each instance of omega.
# Perhaps this would be a good idea to employ in my failed script where I tried to 
# plot Parsegian's data on gold on pg.266 of his book.

# Initialsing empty lists to store the data points as lamda is varied.

eps_vals_real = []
eps_vals_imag = []

# Values to plot on the x-axis, wavelengths in nm.

omega_vals = range(300,701)

# Compute values to plot on the y-axis for epsilon_of_omega evaluated at specific values of lamda.

while lamda < 701.0:
    
    first_term = (78.618)/((1240/lamda)*((1240/lamda) + 0.0379j))
    second_term = (22.81)/(12.96 - ((1240/lamda)**2) - (1.3j*(1240/lamda)))
    third_term = (7.464)/(7.84 - ((1240/lamda)**2) - (0.737j*(1240/lamda)))
    
    epsilon_of_omega = 4.9752 - first_term + second_term + third_term
    
    # Store real and imaginary components separately.
    
    eps_vals_real.append(epsilon_of_omega.real)
    eps_vals_imag.append(epsilon_of_omega.imag)
    
    lamda += 1.0
    
# Plot Re(epsilon_of_omega) in red
# Plot Im(epsilon_of_omega) in green

plt.title('Drude-Lorentz Model for Gold')
plt.xlabel('Omega (nm)')
plt.ylabel('Epsilon(omega)')
plt.plot(omega_vals, eps_vals_real, 'r', label = 'Re')
plt.plot(omega_vals, eps_vals_imag, 'g', label = 'Im')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
