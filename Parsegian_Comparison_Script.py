# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 11:08:46 2015

@author: Ryan-Rhys
"""
import numpy
import matplotlib
from matplotlib import pyplot as plt

wj = [0,2.9,4.0,8.9]
fj = [9.7,4.95,41.55,207.76]
gj = [3.21,0.67,2.22,8.50]

wj2 = [0,3.87,8.37,23.46]
fj2 = [40.11,59.61,122.55,1031.19]
gj2 = [0,2.62,6.41,27.57]

wj3 = [0,3.0,4.8]
fj3 = [53.0,5.0,104.0]
gj3 = [1.8,0.8,4.4]

omega_ranges = range(300,701)

epsilon_vals_real1 = []
epsilon_vals_imag1 = []
epsilon_vals_real2 = []
epsilon_vals_imag2 = []
epsilon_vals_real3 = []
epsilon_vals_imag3 = []

omega = 300.0

j = 1j

while omega < 701.0:
    i = 0
    x1, x2 = 0.0, 0.0
    y1, y2 = 0.0, 0.0
    while i < 4:
        x1 += (1.0 + (fj[i])/((wj[i])**2 - (1240.0/omega)**2 - (j*(1240.0/omega)*gj[i]) )).real
        y1 += (1.0 + (fj[i])/((wj[i])**2 - (1240.0/omega)**2 - (j*(1240.0/omega)*gj[i]) )).imag
        x2 += (1.0 + (fj2[i])/((wj2[i])**2 - (1240.0/omega)**2 - (j*(1240.0/omega)*gj2[i]) )).real
        y2 += (1.0 + (fj2[i])/((wj2[i])**2 - (1240.0/omega)**2 - (j*(1240.0/omega)*gj2[i]) )).imag
        i += 1
    epsilon_vals_real1.append(x1)
    epsilon_vals_imag1.append(y1)
    epsilon_vals_real2.append(x2)
    epsilon_vals_imag2.append(y2)
    omega += 1.0
    
omega2 = 300
    
while omega2 < 701.0:
    i = 0
    x3 = 0.0
    y3 = 0.0
    while i < 3:
        x3 += (1.0 + (fj3[i])/((wj3[i])**2 - (1240.0/omega2)**2 - (j*(1240.0/omega2)*gj3[i]) )).real
        y3 += (1.0 + (fj3[i])/((wj3[i])**2 - (1240.0/omega2)**2 - (j*(1240.0/omega2)*gj3[i]) )).imag
        i += 1
    epsilon_vals_real3.append(x3)
    epsilon_vals_imag3.append(y3)
    omega2 += 1.0   
    



plt.figure(1)
plt.title('Irani Fit')
plt.xlabel('Omega (nm)')
plt.ylabel('Epsilon(omega)')    
plt.plot(omega_ranges, epsilon_vals_real1, 'r', label = 'Re')
plt.plot(omega_ranges, epsilon_vals_imag1, 'g', label = 'Im')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
           
plt.figure(2)
plt.title('DESY Fit')
plt.xlabel('Omega')
plt.ylabel('Epsilon(omega)')    
plt.plot(omega_ranges, epsilon_vals_real2, 'r', label = 'Re')
plt.plot(omega_ranges, epsilon_vals_imag2, 'g', label = 'Im')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
           
plt.figure(3)
plt.title('Johnson and Christy Fit')
plt.xlabel('Omega')
plt.ylabel('Epsilon(omega)')    
plt.plot(omega_ranges, epsilon_vals_real3, 'r', label = 'Re')
plt.plot(omega_ranges, epsilon_vals_imag3, 'g', label = 'Im')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
