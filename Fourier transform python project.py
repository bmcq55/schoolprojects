# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 09:31:50 2020

@author: Brendon
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#exercise 1: adding two sine waves together
N=2000 #number of sample points
t=np.arange(N)
print(t)
#Defined wave parameters
T1 =400.2 ; #period, seconds
amp1 = 3.04; # amplitude
T2 = 200; # period, seconds
amp2 = 1.13; # amplitude
T3=143
amp3=0.29

#Make signals
y1 = amp1 * np.sin(2.*(np.pi)*t/T1)
y2 = amp2 * np.sin(2.*np.pi*t/T2)
y3= amp3 * np.sin(2.*np.pi*t/T3)
signal = y1 + y2+y3

plt.figure(1)
plt.plot(t,signal)
plt.xlabel("time(s)")
plt.ylabel("amplitude")
plt.title("Sum of Two Sine Waves")
#plt.plot(t,y1)
#plt.plot(t,y2)

# fast fourier transformer of sigal
from scipy import fft
Y=fft(signal)

plt.figure(2)
plt.plot(np.imag(Y))
plt.title("FFT of the Sum of Two Sine Waves")
plt.xlabel("Absolute Value of the FFT")

#exercise 4: time dependent frequency function

amp=3
N=2000
t=np.arange(N)
T=0.04*t+9
y = amp * np.sin(2*np.pi*(t/T))

plt.figure(3)
plt.plot(t,y)
plt.title("frequency dependent sine wave function")
plt.xlabel("time")
plt.ylabel("amplitude")

Y_=fft(y)

plt.figure(4)
plt.plot(np.abs(np.imag(Y_)))
plt.title("FFT of the frequency dependent Sine Wave")
plt.xlabel("Absolute Value of the FFT")
