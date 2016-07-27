# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 10:35:52 2016

@author: Ryan-Rhys
"""

import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt


# Constants

k = 1.38*(10**-23) 
T = 298.15
hbar = 1.0545718*(10**-34)
pi = np.pi

# 18 February - omega_star value of 3.9*10^13 Hz.

omega_star = (k*T)/(hbar)

# Radius given in nanometres. CGS units of polarisability are cm^3.
# In this case radius would have to be given in units of centimetres for consistency.

# 8 February - Having looked at the formula for the energy given in the Silbey paper
# it would seem that as long as the dimensions for the interparticle distance and the
# particle radius are consistent, the energy shoudl come out in units of Joules.

Radius = 6.25

# Geometric Variables

interparticle_distance = 6.25
Image_Distance = 0.01
Surface_Distance = Image_Distance/2.0
psi = pi/2.0
xi = pi/2.0 + np.arccos(interparticle_distance/Image_Distance)

# Matsubara coefficient

Matsubara_constant = (2*pi*k*T)/(hbar)
rad_s = 1.05*(10**11)
eV_conversion = 6.55*(10**-5)

# Conversion from rad s^-1 to electronvolts required.
# Matsubara variables for use in conjunction with '#2' - the sum as opposed
# to the integral.

Matsubara_coeff = (Matsubara_constant/rad_s)*eV_conversion

Matsubara = lambda n: Matsubara_coeff*n

# 25 February - Alternative Parameters for gold

wj2 = [0.0, 3.87, 8.37, 23.46]
fj2 = [40.11, 59.61, 122.55, 1031.19]
gj2 = [0, 2.62, 6.41, 27.57]

# Johnson and Christy parameters for gold

# 19th February - Parameters below in eV.

# wj3 = [0,3.0,4.8]
# fj3 = [53.0,5.0,104.0]
# gj3 = [1.8,0.8,4.4]

# Parameters below in Hz.

wj3 = [0.0, 725400000000000.0, 1160640000000000.0]
fj3 = [3.09876372e+30, 2.923362e+29, 6.08059296e+30]
gj3 = [435240000000000.0, 193440000000000.0, 1063920000000000.1]


# 19th February - Code for converting values given above in electronvolts to Hz ('#5').

#5 wj3_Hz, fj3_Hz, gj3_Hz = [], [], []

#5 for element in wj3:
    #5 wj3_Hz.append(element*2.418*10**14)
    
#5 for element in fj3:
    #5 fj3_Hz.append(element*(2.418*10**14)*(2.418*10**14))
    
#5 for element in gj3:
    #5 gj3_Hz.append(element*2.418*10**14)

# 18 February - Computed dimensionless parameters by dividing by omega_star.
# Code for achieving this given as '#3' below. I have just extracted the resultant
# lists for ease of use.

#3 Reduced_wj3 = []
#3 Reduced_fj3 = []
#3 Reduced_gj3 = []

#3 for element in wj3:
    #3 Reduced_wj3.append((element)/(omega_star))
    
#3 for element in fj3:
    #3 Reduced_fj3.append((element)/(omega_star**2))
    
#3 for element in gj3:
    #3 Reduced_gj3.append((element)/(omega_star))
    
# 25th February - Alternative reduced parameters calculated in Junk_Script.
    
Reduced_wj2 = [0.0, 715095.9673904348, 1546602.9062165217, 4334922.84108]
Reduced_fj2 = [1369491653800.59, 2035287895364.0781, 4184273302748.998, 35208329555787.34]
Reduced_gj2 = [0.0, 484121.81771652185, 1184435.4395278264, 5094365.845207826]
    
Reduced_wj3 = [0.0, 18.59258625582396, 29.748138009318335]
Reduced_fj3 = [2035.6962194503851, 192.0468131556967, 3994.5737136384914]
Reduced_gj3 = [11.155551753494377, 4.958023001553056, 27.269126508541813]

# Parameters for water found on page 266 of Parsegian's book. Given in eV below.

# Wwj = [2.07*(10**-2), 6.9*(10**-2), 9.2*(10**-2), 2.0*(10**-1), 4.2*(10**-1), 8.25, 
       # 10.0, 11.4, 13.0, 14.9, 18.5]
       
# Wfj = [6.25*(10**-4), 3.5*(10**-3), 1.28*(10**-3), 5.44*(10**-4), 1.35*(10**-2), 
       # 2.68, 5.67, 12.0, 26.3, 33.8, 92.8]
       
# Wgj = [1.5*(10**-2), 3.8*(10**-2), 2.8*(10**-2), 2.5*(10**-2), 5.6*(10**-2),
       # 0.51, 0.88, 1.54, 2.05, 2.96, 6.26]
       
# 19th February - Parameters given in Hz.
       
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
 
# 19th February - Code for converting parameters in eV to Hz given below ('#6').
       
#6 Wwj_Hz, Wfj_Hz, Wgj_Hz = [], [], []

#6 for element in Wwj:
    #6 Wwj_Hz.append(element*2.418*10**14)
    
#6 for element in Wfj:
    #6 Wfj_Hz.append(element*(2.418*10**14)*(2.418*10**14))
    
#6 for element in Wgj:
    #6 Wgj_Hz.append(element*2.418*10**14)
       
# 18 February - Computed dimensionless parameters by dividing by omega_star.
# Code for achieving this given as '#4' below. I have just extracted the resultant
# lists for ease of use.
#4 Reduced_Wwj = []
#4 Reduced_Wfj = []
#4 Reduced_Wgj = []

#4 for element in Wwj:
    #4 Reduced_Wwj.append((element)/(omega_star))
    
#4 for element in Wfj:
    #4 Reduced_Wfj.append((element)/(omega_star**2))
    
#4 for element in Wgj:
    #4 Reduced_Wgj.append((element)/(omega_star))
    
# 19 February - Extracted Lists.
    
# 25th February - Trying omega star at 0 Kelvin (given second).

Reduced_Wwj = [0.12828884516518532,
 0.42762948388395117,
 0.5701726451786014,
 1.239505750388264,
 2.602962075815355,
 51.1296122035159,
 61.9752875194132,
 70.65182777213106,
 80.56787377523716,
 92.34317840392569,
 114.65428191091443]
 
#Reduced_Wwj = [3824.9319186,
 #12749.773062000002,
 #16999.697416,
 #36955.86394782609,
 #77607.3142904348,
 #1524429.3878478264,
 #1847793.1973913044,
 #2106484.245026087,
 #2402131.1566086956,
 #2753211.864113044,
 #3418417.4151739134]

Reduced_Wfj = [0.024005851644462086,
 0.1344327692089877,
 0.04916398416785835,
 0.020894693271339807,
 0.5185263955203812,
 102.93709185145343,
 217.78108611856007,
 460.9123515736721,
 1010.1662371989646,
 1298.2364569325096,
 3564.388852169731]
 
#Reduced_Wfj = [21339623.127034873,
 #119501889.5113953,
 #43703548.16416742,
 #18574007.969771158,
 #460935859.54395336,
 #91504303968.72554,
 #193593061008.4604,
 #409720764039.0696,
 #897971341185.6274,
 #1154046818710.046,
 #3168507241902.138]

Reduced_Wgj = [0.09296293127911981,
 0.23550609257377017,
 0.17353080505435697,
 0.154938218798533,
 0.34706161010871395,
 3.1607396634900735,
 5.4538253017083616,
 9.544194277989634,
 12.704933941479707,
 18.344685105746308,
 38.79652998715267]  
 
#Reduced_Wgj = [2771.689796086957,
 #7021.614150086956,
 #5173.8209526956525,
 #4619.482993478261,
 #10347.641905391305,
 #94237.45306695653,
 #162605.80137043478,
 #284560.1523982609,
 #378797.6054652174,
 #546946.786427826,
 #1156718.5415669566]
       
# 7th February - Not sure if damped oscillator will be compatible for continuous frequency
# input and not the discrete Matsubara frequency input. Can try anyway.
# probably best to use Thomas-Fermi like dielectric functions.
# Note that n is a substitute for kappa.
       
# 17th February - Trying a sum now instead.
# '#2' indicates code that that will compute the sum over the Matsubara frequencies
# when uncommented.
       
# 22nd February - Setting out parameters for ITO in eV, denoted by #8  
# 23rd February - Updating parameter lists with new ITO values from fit to
# refractive index.info data.
       
#8 omega_plasma = 1.80
#8 omega_tao = 0.05
#8 omega_1 = 4.82
#8 omega_2 = 7.1
#8 f_1 = 0.47
#8 f_2 = 2.3
#8 gamma_1 = 0.2
#8 gamma_2 = 0.1

#8 params_list_ITO_eV = [omega_plasma, omega_tao,omega_1, omega_2, f_1, f_2, gamma_1, gamma_2]
#8 params_list_ITO_Hz = []

#8 for element in params_list_ITO_eV:
    #8 params_list_ITO_Hz.append(element*2.418*10**14)

# 22nd February - Converting ITO parameters to reduced units (#9).
    
#9 params_list_ITO_Hz = [435240000000000.0,
 #9 12090000000000.0,
 #9 1165476000000000.2,
 #9 1716780000000000.0,
 #9 113646000000000.0,
 #9 556140000000000.0,
 #9 48360000000000.0,
 #9 24180000000000.0]
 
#9 Reduced_ITO = []

#9 for element in params_list_ITO_Hz:
    #9 Reduced_ITO.append((element)/(omega_star))
    
Reduced_ITO = [11.155551753494377,
 0.309876437597066,
 29.87208858435717,
 44.00245413878337,
 2.9128385134124204,
 14.254316129465037,
 1.239505750388264,
 0.619752875194132]
       
# 19th February - Defining K, the dimensionless variable:
       
K = lambda n: n/omega_star
 
#2 Epsilon_Water = lambda n: 1 + (Wfj[0])/((Wwj[0])**2 + Wgj[0]*Matsubara(n) + Matsubara(n)**2) + (Wfj[1])/((Wwj[1])**2 + Wgj[1]*Matsubara(n) + Matsubara(n)**2) + (Wfj[2])/((Wwj[2])**2 + Wgj[2]*Matsubara(n) + Matsubara(n)**2) + (Wfj[3])/((Wwj[3])**2 + Wgj[3]*Matsubara(n) + Matsubara(n)**2) + (Wfj[4])/((Wwj[4])**2 + Wgj[4]*Matsubara(n) + Matsubara(n)**2) + (Wfj[5])/((Wwj[5])**2 + Wgj[5]*Matsubara(n) + Matsubara(n)**2) + (Wfj[6])/((Wwj[6])**2 + Wgj[6]*Matsubara(n) + Matsubara(n)**2) + (Wfj[7])/((Wwj[7])**2 + Wgj[7]*Matsubara(n) + Matsubara(n)**2) + (Wfj[8])/((Wwj[8])**2 + Wgj[8]*Matsubara(n) + Matsubara(n)**2) + (Wfj[9])/((Wwj[9])**2 + Wgj[9]*Matsubara(n) + Matsubara(n)**2) + (Wfj[10])/((Wwj[10])**2 + Wgj[10]*Matsubara(n) + Matsubara(n)**2)  
#2 Epsilon_Gold = lambda n: 1 + (fj3[0])/((wj3[0])**2 + gj3[0]*Matsubara(n) + Matsubara(n)**2) + (fj3[1])/((wj3[1])**2 + gj3[1]*Matsubara(n) + Matsubara(n)**2) + (fj3[2])/((wj3[2])**2 + gj3[2]*Matsubara(n) + Matsubara(n)**2)      
#2 Epsilon_ITO = lambda n: 1 + (2.624)/(Matsubara(n)**2 + 0.05*Matsubara(n)) + (10.919)/(23.232 + Matsubara(n)**2 - 0.82*Matsubara(n)) 
 
# 19th February - '#1' below is the commented out code that uses the dimensional variable
# n.

#1 Epsilon_Water = lambda n: 1 + (Wfj[0])/((Wwj[0])**2 + Wgj[0]*n + n**2) + (Wfj[1])/((Wwj[1])**2 + Wgj[1]*n + n**2) + (Wfj[2])/((Wwj[2])**2 + Wgj[2]*n + n**2) + (Wfj[3])/((Wwj[3])**2 + Wgj[3]*n + n**2) + (Wfj[4])/((Wwj[4])**2 + Wgj[4]*n + n**2) + (Wfj[5])/((Wwj[5])**2 + Wgj[5]*n + n**2) + (Wfj[6])/((Wwj[6])**2 + Wgj[6]*n + n**2) + (Wfj[7])/((Wwj[7])**2 + Wgj[7]*n + n**2) + (Wfj[8])/((Wwj[8])**2 + Wgj[8]*n + n**2) + (Wfj[9])/((Wwj[9])**2 + Wgj[9]*n + n**2) + (Wfj[10])/((Wwj[10])**2 + Wgj[10]*n + n**2)  
#1 Epsilon_Gold = lambda n: 1 + (fj3[0])/((wj3[0])**2 + gj3[0]*n + n**2) + (fj3[1])/((wj3[1])**2 + gj3[1]*n + n**2) + (fj3[2])/((wj3[2])**2 + gj3[2]*n + n**2)      
#1 Epsilon_ITO = lambda n: 1 + (2.624)/(n**2 + 0.05*n) + (10.919)/(23.232 + n**2 - 0.82*n) 

Epsilon_Water = lambda K: 1 + (Reduced_Wfj[0])/((Reduced_Wwj[0])**2 + Reduced_Wgj[0]*K + K**2) + (Reduced_Wfj[1])/((Reduced_Wwj[1])**2 + Reduced_Wgj[1]*K + K**2) + (Reduced_Wfj[2])/((Reduced_Wwj[2])**2 + Reduced_Wgj[2]*K + K**2) + (Reduced_Wfj[3])/((Reduced_Wwj[3])**2 + Reduced_Wgj[3]*K + K**2) + (Reduced_Wfj[4])/((Reduced_Wwj[4])**2 + Reduced_Wgj[4]*K + K**2) + (Reduced_Wfj[5])/((Reduced_Wwj[5])**2 + Reduced_Wgj[5]*K + K**2) + (Reduced_Wfj[6])/((Reduced_Wwj[6])**2 + Reduced_Wgj[6]*K + K**2) + (Reduced_Wfj[7])/((Reduced_Wwj[7])**2 + Reduced_Wgj[7]*K + K**2) + (Reduced_Wfj[8])/((Reduced_Wwj[8])**2 + Reduced_Wgj[8]*K + K**2) + (Reduced_Wfj[9])/((Reduced_Wwj[9])**2 + Reduced_Wgj[9]*K + K**2) + (Reduced_Wfj[10])/((Reduced_Wwj[10])**2 + Reduced_Wgj[10]*K + K**2)
Epsilon_Gold = lambda K: 1 + (Reduced_fj3[0])/((Reduced_wj3[0])**2 + Reduced_gj3[0]*K + K**2) + (Reduced_fj3[1])/((Reduced_wj3[1])**2 + Reduced_gj3[1]*K + K**2) + (Reduced_fj3[2])/((Reduced_wj3[2])**2 + Reduced_gj3[2]*K + K**2)
#Epsilon_Gold = lambda K: 1 + (Reduced_fj2[0])/((Reduced_wj2[0])**2 + Reduced_gj2[0]*K + K**2) + (Reduced_fj2[1])/((Reduced_wj2[1])**2 + Reduced_gj2[1]*K + K**2) + (Reduced_fj2[2])/((Reduced_wj2[2])**2 + Reduced_gj2[2]*K + K**2) + (Reduced_fj2[3])/((Reduced_wj2[3])**2 + Reduced_gj2[3]*K + K**2)
Epsilon_ITO = lambda K: 1 + (Reduced_ITO[0]**2)/(K**2 + Reduced_ITO[1]*K) + (Reduced_ITO[2]*Reduced_ITO[4])/((Reduced_ITO[2])**2 + K**2 - Reduced_ITO[6]*K) + (Reduced_ITO[3]*Reduced_ITO[5])/((Reduced_ITO[3])**2 + K**2 - Reduced_ITO[7]*K)

# 22nd February - Taking out the Radius dependence to introduce dimensionlessness (#10).

Polarisability = lambda K: Epsilon_Water(K) * ((Epsilon_Gold(K) - Epsilon_Water(K))/(Epsilon_Gold(K) + 2*Epsilon_Water(K))) #10 * (Radius)**3

# 17 February - using discrete sum over Matsubara frequencies.
# Result for interparticle distance = 1 nm and image distance = 2 nm is
# -9.0487 *10^-8 in units of KT. roughly same as the integral (-1.556*10^-8)

#2 I_0 = 0
#2 I_1 = 0
#2 I_2 = 0

#2 for n in range(1,100000):
    #2 I_0 += Polarisability(n)**2
    #2 I_1 += ((Epsilon_ITO(n) - Epsilon_Water(n))/(Epsilon_ITO(n) + Epsilon_Water(n)))*Polarisability(n)**2
    #2 I_2 += (((Epsilon_ITO(n) - Epsilon_Water(n))/(Epsilon_ITO(n) + Epsilon_Water(n)))**2)*Polarisability(n)**2

# 8 February - using numerical quadrature for integration.

# 22nd February - Introducing A, dimensionless coefficient. Previous calcs given as (#11).

#11 def Integrand_I_0(n):
    #11 return Polarisability(n)**2
    
#11 I_0, err = quad(Integrand_I_0, 0, np.inf)

#11 I_0 = I_0*2.418*10**14

def Integrand_A1(K):
    return Polarisability(K)**2
    
A1, err = quad(Integrand_A1, 0, np.inf)

#11 def Integrand_I_1(n):
    #11 return ((Epsilon_ITO(n) - Epsilon_Water(n))/(Epsilon_ITO(n) + Epsilon_Water(n)))*Polarisability(n)**2
    
#11 I_1, err = quad(Integrand_I_1, 0, np.inf)

#11 I_1 = I_1*2.418*10**14

def Integrand_A2(K):
    return ((Epsilon_ITO(K) - Epsilon_Water(K))/(Epsilon_ITO(K) + Epsilon_Water(K)))*Polarisability(K)**2
    
A2, err = quad(Integrand_A2, 0, np.inf)

#11 def Integrand_I_2(n):
    #11 return (((Epsilon_ITO(n) - Epsilon_Water(n))/(Epsilon_ITO(n) + Epsilon_Water(n)))**2)*Polarisability(n)**2
    
#11 I_2, err = quad(Integrand_I_2, 0, np.inf)

#11 I_2 = I_2*2.418*10**14

def Integrand_A3(K):
    return (((Epsilon_ITO(K) - Epsilon_Water(K))/(Epsilon_ITO(K) + Epsilon_Water(K)))**2)*Polarisability(K)**2
    
A3, err = quad(Integrand_A3, 0, np.inf)

# Energy expression with angular dependence given below.

# 22nd February - Omitting old energy (#12)

#12 Energy = - ((3*hbar)/(pi*(interparticle_distance)**6))*I_0 + ((hbar*(2 - 3*np.cos(2*psi) - 3*np.cos(2*xi)))/(2*((Image_Distance)**3)*pi*((interparticle_distance)**3)))*I_1 - ((3*hbar)/(pi*(Image_Distance)**6))*I_2

# Energy expression with system orientated in an isosceles triangle given below.

# Energy = - ((3*hbar)/(pi*(interparticle_distance**6)))*I_0 - (((hbar*2)/(Image_Distance)**3)*2*pi*((interparticle_distance)**3))*I_1 - ((3*hbar)/(pi*(Image_Distance**6)))*I_2

# Energy = - ((3*hbar)/(pi*(interparticle_distance**6)))*I_0 - ((3*hbar)/(pi*(Image_Distance**6)))*I_2

# 9 February - First term in the above expression checked manually and is fine.

# 11 February - Second energy expression is for the formula without the term
# containing the angular dependence.

#12 kT_Energy = Energy/(k*T)  

#12 print kT_Energy

# 23rd February - Energy as a function of interparticle distance and Image_Distance.

Energy = ((Radius/interparticle_distance)**6)*A1 - (((5.0/3.0) - ((Image_Distance**2 - interparticle_distance**2)/(Image_Distance**2 + interparticle_distance**2)))/(2*((interparticle_distance/Radius)**3)*(((np.sqrt(interparticle_distance**2 + Image_Distance**2))/Radius)**3)))*A2 + (((Radius)/(np.sqrt(interparticle_distance**2 + Image_Distance**2)))**6)*A3
Energy = Energy*-1


Energy_bulk = (3.0/pi)*((Radius/interparticle_distance)**6)*A1
Energy_bulk = Energy_bulk*-1

Energy_Vals = []
Energy_bulk_vals = []
distance_vals = []

# 24th February - Deciding to set x, the surface-to-surface separation to 0 and will
# vary it up to 1 diameter. The diagram in Silbey is not so clear as to whether
# the interparticle separation refers to the centre-centre distance or the surface-
# surface distance.
# Remember factor of pi/3 in energy expression.

x = 0.001

while x < (2*Radius):
    distance_vals.append(x)
    x += 0.001

interparticle_distance = 2*Radius

Silbey_Energy = ((((Radius/interparticle_distance)**6)*A1)*-1)

while interparticle_distance < (4*Radius):
    Energy_Vals.append((3.0/pi)*(((Radius/interparticle_distance)**6)*A1 - (((5.0/3.0) - ((Image_Distance**2 - interparticle_distance**2)/(Image_Distance**2 + interparticle_distance**2)))/(2*((interparticle_distance/Radius)**3)*(((np.sqrt(interparticle_distance**2 + Image_Distance**2))/Radius)**3)))*A2 + (((Radius)/(np.sqrt(interparticle_distance**2 + Image_Distance**2)))**6)*A3)*-1)
    interparticle_distance += 0.001
    
interparticle_distance = 2*Radius
    
while interparticle_distance < (4*Radius):
    Energy_bulk_vals.append((3.0/pi)*(((Radius/interparticle_distance)**6)*A1)*-1)
    interparticle_distance += 0.001
    
    
Hamaker = 29.3

interparticle_distance = 2*Radius

Parsegian_Energy = -((Radius/interparticle_distance)**6)*(16.0/9.0)*Hamaker
Parsegian_Vals = []



while interparticle_distance < (4*Radius):
    Parsegian_Vals.append(-((Radius/interparticle_distance)**6)*(16.0/9.0)*Hamaker)
    interparticle_distance += 0.001

plt.figure(1)
#ITO = plt.plot(distance_vals, Energy_Vals, label = 'ITO')
Silbey = plt.plot(distance_vals, Energy_bulk_vals, label = 'Silbey')
Parsegian = plt.plot(distance_vals, Parsegian_Vals, label = 'Parsegian')
plt.xlabel('Surface to Surface Separation (nm)')
plt.ylabel('Energy (kT)')
plt.title('Interaction Potential (Large Separation Limit)')
#plt.text(4,-0.20,'ITO Distance = 0.5 Angstroms')
plt.legend()

