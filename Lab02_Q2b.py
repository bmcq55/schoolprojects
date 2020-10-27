# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 15:12:35 2020

@author: Brendon
"""

#Q2bi
#diffraction limit of telescope, ex. 5.4a
import numpy as np
import matplotlib.pyplot as plt
import scipy.special

def f(m,x,t): #integrand of Bessel function
    return np.cos(m*t - x*np.sin(t))

def J(m,x): #Bessel function
    if m<0:
        raise ValueError("m must be a nonnegative integer")
        
    a=0.0 #lower bound
    b=np.pi #upper bound
    N=1000 #number of steps
    if N % 2 == 1:
        raise ValueError("N must be an even integer.")
    h=(b-a)/N #step size
    
    #iteration of for loop
    s=f(m,x,a)+f(m,x,b)
    for n in range(1,N):
        t=a+n*h
        if n%2==1:
            s+=4.0*f(m,x,t)
        else:
            s+=2.0*f(m,x,t)
    return s*h/(3.0*np.pi)
    
#for x in range(21):   
   # print(x,J(0,x),J(1,x),J(2,x))
 
#plot of Bessel functions for m=0,1,2 (order)
v = np.linspace(0,20,100)     
plt.figure()
plt.plot(v,J(0,v),label='m=0',)
plt.plot(v,J(1,v),label='m=1')
plt.plot(v,J(2,v),label='m=2')
plt.legend()
plt.title("Bessel Functions for m=0,1,2")
plt.xlabel("x-value")
plt.ylabel("J(m,x)")
plt.savefig("Bessel Function #1.png")

#computing Bessel function values using scipy.special.jv 
x=np.linspace(0,20,1000)
j0=scipy.special.jv(0,x)
j1=scipy.special.jv(1,x)
j2=scipy.special.jv(2,x)
#print(x,j0,j1,j2)

#plotting the scipy.special.jv Bessel function for m=0,1,2 (order)
plt.figure()
plt.plot(x,j0,color='red',label='m=0')
plt.plot(x,j1,color='black',label='m=1')
plt.plot(x,j2,color='green',label='m=2')
plt.legend()
plt.title("Bessel Function Plot Using the Scipy Built-in Function")
plt.xlabel("x")
plt.ylabel("J(m,x)")
plt.savefig("Bessel scipy function.png")


#Lab 2, Q2b exercise 5.4b

from pylab import imshow,gray,show

wl=500*(10**-9) #wavelength in meters
k=(2*np.pi)/wl #k value

X = np.linspace(-1e-6, 1e-6, 500) #x-axis, 0.002 mm
Y = np.linspace(-1e-6, 1e-6, 500) #y-axis, 0.002 mm

#intensity of diffraction pattern 
I_2D = np.empty([len(X), len(Y)],float) #
for i in range(len(X)):
    x = X[i]
    for j in range(len(Y)):
        y = Y[j]
        r = np.sqrt(x**2 + y**2)
        I_2D[i,j] = ((scipy.special.jv(1,k*r))/(k*r))**2

#plot of diffraction pattern
plt.figure()
plt.pcolormesh(I_2D, cmap='hot', shading='auto', vmax = 0.01)
plt.xlabel("x Coordinate (0.002 micrometer scale)")
plt.ylabel("y Coordinate(0.002 micrometer scale)")
plt.show()
plt.savefig("diffraction pattern.png")
    
    


        