# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

CB2_V, CB2_C, CB1_V, CB1_C, CB_R, CB_uR = np.loadtxt('C:/Users/Brendon/Desktop/PHY224 Computational physics/Cell Battery Data.txt', skiprows=1, unpack=True)
# PS(-)_V(-) -- the first number represents the circuit layout option, the second represents the resistor it is plugged into

PS2_V1, PS2_C1, PS1_V1, PS1_C1, PS2_V2, PS2_C2, PS1_V2, PS1_C2, PS2_V3, PS2_C3, PS1_V3, PS1_C3, PS2_V4, PS2_C4, PS1_V4, PS1_C4 = np.loadtxt('C:/Users/Brendon/Desktop/PHY224 Computational physics/Power supply data.txt', skiprows=1, unpack=True)

#cell battery internal resistance is -9.52 mohms using multimeter
plt.figure()
plt.title ("V vs. I Plot: Cell Battery" )
plt.plot(CB2_C, CB2_V)
plt.plot(CB1_C, CB1_V)
plt.xlabel('Current')
plt.ylabel('Voltage')

# output resistance of cell battery using V_opencircuit= 6.41 V

#cell battery Circuit arrangement option 1
V_load1= CB1_V
R_load1= CB_R
I= V_load1/ R_load1 #this is current in the circuit for each resistor
V_int= 6.41 - CB1_V #this is voltage across the battery or terminal voltage when each resistor is connected
R_b_option1= V_int/ I # this is the output resistance of the battery for each resistor 
print(R_b_option1)

plt.figure()
plt.plot(I, V_load1)

#cell battery Circuit arrangement option 2
V_load2= CB2_V
R_load2= CB_R
I= V_load2/ R_load2 #this is current in the circuit for each resistor
V_int= 6.41 - CB2_V #this is voltage across the battery or terminal voltage when each resistor is connected
R_b_option2= V_int/ I # this is the output resistance of the battery for each resistor 
print(R_b_option2)

plt.figure()
plt.plot(I, V_load2)



plt.figure(figsize=(10,5))
plt.title("V vs. I Plot: Power Supply")
plt.plot(PS1_C1, PS1_V1, '-o', color='blue', label='Option1, Measurement1')
plt.plot(PS2_C1, PS2_V1, '-o', color='red', label='Option2, Measurement1')
plt.plot(PS1_C2, PS1_V2, '-o', color='green', label='Option1, Measurement2')
plt.plot(PS2_C2, PS2_V2, '-o', color='blueviolet', label='Option2, Measurement2')
plt.plot(PS1_C3, PS1_V3, '-o', color='gold', label='Option1, Measurement3')
plt.plot(PS2_C3, PS2_V3, '-o', color='darkslategrey', label='Option2, Measurement3')
plt.plot(PS1_C4, PS1_V4, '-o', color='deepskyblue', label='Option1, Measurement4')
plt.plot(PS2_C4, PS2_V4, '-o', color='m', label='Option2, Measurement4')
plt.legend()
plt.xlabel('Current')
plt.ylabel('Voltage')

# output resistance of power supply using V_opencircuit= 6.5 V
#10 mohm option 2 output resistance 
V_load= PS2_V1
R_load= 10
I= V_load/ R_load #this is current in the circuit for each resistor
V_int= 6.5 - PS2_V1 #this is voltage across the battery or terminal voltage when each resistor is connected
R_b= V_int/ I # this is the output resistance of the battery for each resistor 
print(R_b)

#10 mohm option 1 output resistance
V_load= PS1_V1
R_load= 10
I= V_load/ R_load #this is current in the circuit for each resistor
V_int= 6.5 - PS1_V1 #this is voltage across the battery or terminal voltage when each resistor is connected
R_b= V_int/ I # this is the output resistance of the battery for each resistor 
print(R_b)

