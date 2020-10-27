# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import cmath as cm

# #Q1 Newman,Exercise 7.2 

#7.2a
#load in sunspot data (months vs. sunspot count) and plot
data=np.loadtxt("C:/Users/brend/OneDrive/Desktop/PHY407 Comp Physics/sunspots.txt", float)
#print(data)

plt.figure()
plt.plotfile('sunspots.txt', delimiter=' ', cols=(0, 1), names=('months', 'sunspot count'))
plt.title("Sunspot Count Over the Course of many Months")
plt.xticks(np.arange(0,3200,400))
plt.show()
plt.savefig("sunspot graph.png")

#7.2b
#calculate Fourier transform of sunspot data, make a graph of magnitude squared of Fourier coefficients,|c_k|, 
#as a function of k-- AKA the power spectrum of sunpot signal


#discrete Fourier transform function
def dft(y):
    '''INPUT:function you want the Fourier transform of
    OUTPUT:coefficients of FT'''
    N = len(y)
    c = np.zeros(N//2+1, complex)
    for k in range(N//2+1):
        for n in range(N):
            c[k] += y[n]*np.exp(-2j*np.pi*k*n/N)
    return c

#take DFT of sunspot data
sunspots = data[:, 1]
c=dft(sunspots) # compute the Fourier coefficients

def mag_squared(a):
    return np.abs(a)**2

# plt.figure()
# plt.plot(mag_squared(c))
# plt.xlabel('k')
# plt.ylabel('|c_k|²')
# plt.show()


#Q1b: Newman 7.4:Fourier filtering and smoothing

#7.4a: read in data and plot
data2=np.loadtxt("C:/Users/brend/OneDrive/Desktop/PHY407 Comp Physics/dow.txt",unpack=True)
#print(data2)
#print(len(data2)) #find length of data

plt.figure()
plt.plot(data2, label='dow.txt data')
plt.title("Dow Jones Industrial Average from late 2006-2010")
plt.ylabel("Daily Closing Value")
plt.xlabel("Days")
plt.xticks(np.arange(0,1050,100))
plt.show()


#7.4b:use rfft to find coefficients of data2

c2 = np.fft.rfft(data2)
#print("the coefficients of the FT of the Dow data is" , c2)

#7.4c: set last 90% of c2 coefficients array to zero

#print(len(c2))
#513*0.1=51.3 coefficients

c2[51:] = 0
#print(c2)

#take inverse FT of 10% of c2
c2_10=np.fft.irfft(c2)

plt.plot(c2_10, label='10% of coefficients')
plt.show()


#7.4d:modify program so that it sets all but 2% of c2 to zero

#0.2*513=10 coefficients
c2[10:]=0

c2_2=np.fft.irfft(c2)

plt.plot(c2_2,label='2% of coefficients')
plt.show()
plt.legend()
plt.savefig("Dow Average.png")



