# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:58:46 2020

@author: Brendon
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#exercise 1:measure the stopping voltage (Vstop) for each of the 8 wavelengths provided.
f=np.array([1/390,1/455,1/505,1/535,1/590,1/615,1/640,1/935])*10**9*3*10**8
V_s1=[1.32,1.07,0.92,0.76,0.581,0.493,0.483,0.0]
"""sigma=np.array[40,40,30,30,10,10,10,10]"""

plt.plot(f,V_s1,'o')
plt.title("Stopping Voltage vs. Frequency of Diode")
plt.xlabel("Frequency(m)")
plt.ylabel("Voltage (V)")
e=1.602*10**-19
#1.measuring Plank's constant
def fx(freq,h,f0):
    return (h/e)*(freq)-f0

popt,pcov=curve_fit(fx,f,V_s1)
plt.plot(f,fx(f,*popt))
real_f0=popt[1]/(popt[0]/e)
pvar= np.diag(pcov)
print (pvar)

#2. measuring work function E_0

#3. measuring f0


#exercise 2: using variable LED, Graph Vstop and the photocurrent vs. the intensity setting.
V_s2=[0.992,0.994,0.996,1.001]
Intensity=[1,2,3,4]
I=np.array([0.199,0.210,0.223,0.251])/(100*10**3)
print(I)

plt.figure()
plt.scatter(Intensity,V_s2)
plt.title("Stopping Voltage vs. Intensity Setting")
plt.xlabel("Intensity setting")
plt.ylabel("Voltage (V)")

plt.figure()
plt.plot(Intensity,I,marker='o',color='red')
plt.title("Photocurrent vs. Intensity Setting")
plt.ylabel("Photocurrent(A)")
plt.xlabel("Intensity Setting")

#exercise 3: Using the work function that you determined in Exercise 1, estimate the time that
#an electron will need to absorb enough energy to escape the photocathode. 
#Compare this with the measured time delay between incident light and a
#photocurrent response
