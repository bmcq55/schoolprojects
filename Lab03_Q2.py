# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 20:27:37 2020

@author: Brendon
"""

#quantum mechanical observables pseudocode
#create function called H(n,x)
#only permit integer n values greater than 0
#in the function, create if/elif/else statements for the first two Hermite polynomials
#for n>=2, write program that calculates hermite polynomial using the definition and returns value
#make another function for the wavefunction in the potential well
#plot the wavefunctions on a graph between x=[-4,4]

import numpy as np
import matplotlib.pyplot as plt
import math

#Q2a:
def h(n,x): #hermite polynomial function
    """INPUT: n=0,..∞ (nth hermite polynomial) as a function of x
    
        OUTPUT: value of nth hermite polynomial"""
    
    if n<0 or n!=int(n):
        raise ValueError("n must be greater than zero and an integer value")
    if n==0:
        return 1
    elif n==1:
        return 2*x
    else:
        return(2*x*(h(n-1,x))-2*(n-1)*(h(n-2,x)))

x=np.linspace(-4,4,100)
def psi(n,x):
    return 1 /np.sqrt(2**n*math.factorial(n)*np.sqrt(np.pi))*np.exp(-x**2/2)*h(n, x)


#plot of wavefunctions(n=0,1,2,3)
plt.figure()
plt.title("Quantum Harmonic Oscillator Wavefunctions(n=0,1,2,3)")
plt.xlabel("Position(x)")
plt.ylabel("Wavefunction(Ψ)")
plt.plot(x,psi(0,x),label='n=0')
plt.plot(x,psi(1,x),label='n=1')
plt.plot(x,psi(2,x),label='n=2')
plt.plot(x,psi(3,x),label='n=3')
plt.legend()
plt.show()
plt.savefig("quantum harmonic oscillator wavefunction(n=0,1,2,3).png" )


#Q2b: plot of wavefunction for n=30 for x=[-10,10]
x=np.linspace(-10,10,1000)
plt.figure()
plt.title("Quantum Harmonic Oscillator Wavefunctions(n=30)")
plt.xlabel("Position(x)")
plt.ylabel("Wavefunction(Ψ)")
plt.plot(x,psi(30,x),label="n=30")
plt.legend()
plt.show()
plt.savefig("quantum harmonic oscillator wavefunction(n=30).png" )

#Q2c: evaluate <x^2>, <p^2> and energy E using Gaussian quadrature on 100 points
from gaussxw import gaussxw
from gaussxw import gaussxwab

#root-mean-square(RMS) position
def rms_x(n):
    """ INPUT: nth energy level
    
        OUTPUT: value of the root-mean-square position(x) """
        
    def rms_xintegrand(z):
        def x(z):
            return z/(1-z) #change of variables: equation 5.67 

        return x(z)**2*np.abs(psi(n, x(z)))**2*(1/(1-z)**2) #integrand of <x>

    #integral of uncertainty of position 
    I=0.0
    N =100
    a=0
    b=1
    x,w = gaussxwab(N,a,b)
    for k in range(N):
        I+=w[k]*rms_xintegrand(x[k])
    return np.sqrt(2*I)

#prints momentum uncertainties for n=0,1,2...15    
for i in range(0,16,1): #print out of uncertainty in position for n=0,1,2,...15
    print("the uncertainty in position is:",rms_x(i), "for n=",i)


#RMS momentum
def rms_p(n):
    """ INPUT: nth energy level
    
    OUTPUT: value of the root-mean-square momentum(p) """
    def rms_pintegrand(z):
        def x(z):
            return z/(1-z)
        if n==0:
            return (np.abs((1/np.sqrt(2**n*math.factorial(n)*np.sqrt(np.pi))) * np.e**(-x(z)**2/2)*(-x(z)*h(0,x(z))+2*n*h(0,x(z)))))**2 #integrand of <p> for n=0
        elif n==1:
            return (np.abs((1/np.sqrt(2**n*math.factorial(n)*np.sqrt(np.pi))) * np.e**(-x(z)**2/2)*(-x(z)*h(1,x(z))+2*n*h(1,x(z)))))**2 #integrand of <p> for n=1
        else:
            return (np.abs((1/np.sqrt(2**n*math.factorial(n)*np.sqrt(np.pi))) * np.e**(-x(z)**2/2)*(-x(z)*h(n,x(z))+2*n*h((n-1),x(z)))))**2 #integrand of <p> for n>2
        
    #integral of uncertainty of momentum   
    I=0.0
    N =100
    a=0
    b=1
    x,w = gaussxwab(N,a,b)
    for k in range(N):
        I+=w[k]*rms_pintegrand(x[k])
    return np.sqrt(2*I)

#prints momentum uncertainties for n=0,1,2...15
for i in range(0,16,1):  
    print('the uncertainty in momentum is:',rms_p(i),'for n=',i)

#compute the total energy of the oscillator    
def energy(n):
    """ INPUT: nth energy level
    
    OUTPUT: total energy of quantum harmonic oscillator"""
    
    print("the energy is:" ,1/2*(rms_x(n)**2+rms_p(n)**2), "for n=",n)
    
##prints total enregy of oscillator for n=0,1,2...15
for k in range(0,16,1):
    print(energy(k))
    
print(rms_x(24))    
print(rms_p(24))
