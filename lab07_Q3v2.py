# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 15:49:41 2020

@author: Bmcq
"""

#lab 7, Q3: Hydrogen atom

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as pc

#Q3b. energy eigenvalue n=1,l=0 

m_e=pc.m_e #mass of electron
h_bar=pc.hbar #planks constant over 2pi
e=pc.e #electron charge
a_0=pc.physical_constants['Bohr radius'][0] #Bohr radius
eps_0=pc.epsilon_0
h=0.0002*a_0
r_inf=20.0*a_0

#energy eigenvalue at l=1,n=2
#potential function
def V(r):
    return -e**2/(4*np.pi*eps_0*r)

#change in R and S functions
def dRdS(u,r,E,l):
    R=u[0]
    S=u[1]
    fR=S/r**2
    fS=(2*m_e*r**2/h_bar**2)*(V(r)-E)*R+ (l*(l+1)*R)
    return np.array([fR,fS],float)

#calculate the wavefunction for a particular energy

def solve(f, E):
    r = np.arange(h, r_inf, h)
    R = np.zeros(len(r))
    S = np.zeros(len(r))
    R[0] = 0.0
    S[0] = 1.0
    u=np.array([R[0],S[0]],float)
    
    for i in range(1,len(r)):
        r0 = r[i-1]
        
        k1=h*f(u,r0,E)
        k2=h*f(u+0.5*k1,r0+0.5*h,E)
        k3=h*f(u+0.5*k2,r0+0.5*h,E)
        k4=h*f(u+k3,r0+h,E)
        u+=(k1+2*k2+2*k3+k4)/6
        
        R[i] = u[0]
        S[i] = u[1]
    
    return R

nvals=[1,2,2] #principal quantum number
lvals=[0,0,1] #azimuthal quantum number
E = np.zeros(len(nvals)) #estimated energies

for i in range(len(nvals)):
    n = nvals[i]
    l = lvals[i]
    F = lambda u, r, E: dRdS(u, r, E, l)
    
    E1=-15*e/n**2
    E2=-13*e/n**2
    R=solve(F, E1)
    S2 = R[-1]
    
    target=e/10000
    while abs(E1-E2)>target:
        R = solve(F, E2)
        S1,S2=S2,R[-1]
        E1,E2=E2,E2-S2*(E2-E1)/(S2-S1)
        
    E[i] = E2
    
    print("E at n=",n,",l=",l," is",E[i]/e, "eV")


    
#Q3c/d. calculate/plot normalized eigenfunction R for theoretical/simulated value

# These are the theoretical functions R(r)
def R10(r):
    coeff = 2/a_0**(3/2)
    return coeff*np.exp(-r/a_0)

def R20(r):
    coeff = 1/(2*a_0)**(3/2)
    return coeff*(2-r/a_0)*np.exp(-r/(2*a_0))

def R21(r):
    coeff = 1/(2*np.sqrt(6)*a_0**(3/2))
    return coeff*(r/a_0)*np.exp(-r/(2*a_0))

Rnl = [R10, R20, R21]

def simpson_samples_weights(a, b, N):
    
    if N % 2 == 1:
        raise ValueError("N must be an even integer.")
    w = np.ones(N+1)/3
    w[2:-1:2] = 2*w[2:-1:2]
    w[1:-1:2] = 4*w[1:-1:2]
    w = w*(b-a)/N
    x = np.linspace(a, b, N+1)
    return x, w

# Normalize wavefunctions and create plots
for i in range(len(E)):
    En = E[i]
    n = nvals[i]
    l = lvals[i]
    Rn = Rnl[i]
    
    # Simulated R(r)
    R = solve(F, En)
    
    r, weights = simpson_samples_weights(h, r_inf, len(R)-1)
    
    Rintegral = np.sum(R**2*weights)
    
    R = R/np.sqrt(Rintegral)
    
    # Theoretical R(r)
    Rtrue = Rn(r)
    Rintegral = np.sum(Rtrue**2*weights)
    
    Rtrue = Rtrue/np.sqrt(Rintegral)
    
    plt.figure()
    plt.title(r"Wavefunction for n = %d, l = %d"%(n,l))
    plt.xlabel("r Coordinate (Bohr Radii)")
    plt.ylabel("Amplitude")
    plt.plot(r/a_0, R*np.sqrt(a_0), 'b.', label='theoretical value') 
    plt.plot(r/a_0, Rtrue*np.sqrt(a_0), 'r',label='simulated value') 
    plt.legend()
    plt.show()    
    plt.savefig(r"Wavefunction_n%d_l%d.png"%(n,l))

