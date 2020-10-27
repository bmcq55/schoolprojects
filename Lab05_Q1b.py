# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 11:12:14 2020

@author: Brendon
"""
import numpy as np
import matplotlib.pyplot as plt

#Q1b: Newman 7.4:Fourier filtering and smoothing

#7.4a: read in data and plot
data2=np.loadtxt("C:/Users/brend/OneDrive/Desktop/PHY407 Comp Physics/dow.txt",unpack=True)
#print(data2)
#print(len(data2)) #find length of data

#plot dow.txt
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

#take inverse FT of 10% of c2 and plot
c2_10=np.fft.irfft(c2)
plt.plot(c2_10, label='10% of coefficients')
plt.show()


#7.4d:modify program so that it sets all but 2% of c2 to zero

#0.2*513=10 coefficients
c2[10:]=0
c2_2=np.fft.irfft(c2)

#plot inverse FT of 2% of coefficients
plt.plot(c2_2,label='2% of coefficients')
plt.show()
plt.legend()
plt.savefig("Dow Average.png")
