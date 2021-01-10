# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 15:21:39 2020

@author: BMCQ
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#pseudocode
#(step 1):assign constants in the problem
#(step 2):discretize TDSE in space/time arrays and create Hamiltonian array
#(step 3): implement CN scheme and solve the linear SOE
#(step 4): compute energy, normalization and expected position and plot results


def expval(op, psi): #expectation value
    '''iINPUT: operator,wavefunction
        OUTPUT: dot product of wavefunction and matrix multiplication
            of (op,psi)'''
    return np.real(np.vdot(psi, np.matmul(op, psi)))

def normalized(psi): #normalization
    '''INPUT:wavefunction
        OUTPUT:normalized wavefunction'''
    n = np.size(psi)
    psi_sq = expval(np.eye(n), psi)
    return psi/np.sqrt(psi_sq)

def animate(T,X,Y,Yname,figname): #animation
    '''INPUT:total time, X axis, Y axis, Y axis label, name of figure
        OUTPUT: animation of wavefunction'''
    plt.style.use('seaborn-pastel')
    
    xmin = np.min(X)
    xmax = np.max(X)
    xlim = (xmin,xmax)
    
    ymin = np.min(Y)
    ymax = np.max(Y)
    ylim = (ymin,ymax)
    
    fig = plt.figure()
    ax = plt.axes(xlim=xlim, ylim=ylim)
    plt.xlabel('Position (m)')
    plt.ylabel(Yname)
    
    line, = ax.plot([], [], lw=3)
    
    def init():
        line.set_data([], [])
        return line,
    def animate(i):
        line.set_data(X, Y[i,:])
        plt.title(r'%s at t = %.2E'%(Yname, T[i]))
        return line,
    
    anim = FuncAnimation(fig, animate, init_func=init,
                                   frames=len(T), interval=20, blit=True)
    
    
    anim.save(r'%s.gif'%(figname), writer='pillow')

#Time-dependent Schrodinger equation with Crank-Nicholson

#Q1a. square well potential

L=1*10**-8 #width of well in meters
m=9.109*10**-31 #mass of electron in kg
h_bar=6.626*10**-34 #Planck's constant

τ0=1*10**-18 #seconds
N0 = 3000
T0 = N0*τ0
N = 300
τ=T0/N

P = 1024 #grid spacing
a=L/P #step size
#print(a)

#initial condition
σ=L/25
κ=500/L
x_0=L/5
Ψ_0=1

def psi_0(x):
    return Ψ_0*np.exp(-(x-x_0)**2/(4*σ**2) +1j*κ*x)

#TDSE spatial discretization
x_p=np.linspace(-L/2,L/2,P+1)
#print(np.size(x_p))
psi=np.array(list(map(psi_0,x_p)),complex)
psi[0]=psi[P]=0

Psi = normalized(psi[1:P])
X = x_p[1:P]

def V_s(X):
    return 0*X

def H(V, X):
    #Hamiltonian matrix
    A = -h_bar**2/(2*m*a**2)
    B_p = V(X)-2*A #V(x)=0 in square well
    
    return np.diag(B_p) + A*np.eye(P-1,k=1) + A*np.eye(P-1,k=-1)

H_D = H(V_s, X)
#print(H_D)

#main loop for CN
I=np.eye(P-1) #rank-n Identity matrix

T=np.linspace(0,T0,N)
PsiVsT = np.zeros([N, P-1], dtype=complex)
EnergyVsT = np.zeros([N])
NormVsT = np.zeros([N])
PsiVsT[0] = Psi

R_mat=(I - 1j*τ/(2*h_bar)*H_D)
L_mat=(I + 1j*τ/(2*h_bar)*H_D)

for n in range(N-1):
    v = np.matmul(R_mat,PsiVsT[n]) #v=RΨ^n
    PsiNew = np.linalg.solve(L_mat,v) #LΨ^(n+1)=v             
    PsiVsT[n+1] = PsiNew

for n in range(N):
    EnergyVsT[n] = expval(H_D, PsiVsT[n])
    NormVsT[n] = expval(I, PsiVsT[n])

#plots/animations
plt.figure()
plt.plot(T, EnergyVsT, 'k')
plt.xlabel('Time (s)')
plt.ylabel('Energy')
plt.title('Energy vs Time(Square Well Potential')
plt.savefig('EnergyVsTimeQ1a')

plt.figure()
plt.plot(T, NormVsT, 'k')
plt.xlabel('Time (s)')
plt.ylabel('Normalization')
plt.title('Normalization vs Time(Square well potential)')
plt.savefig('NormVsTimeQ1a')

#q1.b. 
animate(T, X, np.real(np.conj(PsiVsT)*PsiVsT), r'Probability Density $|\psi|^2$', 'ProbDensQ1b')

#Q1.c using harmonic oscillator potential 

ω=3*10**15
N0 = 4000
T0 = N0*τ0
N = 300
τ=T0/N

def V_h(X):
    return 0.5*m*ω**2*X**2

def H(V, X):
    #Hamiltonian matrix
    A = -h_bar**2/(2*m*a**2)
    B_p = V(X)-2*A #V(x)=0 in square well
    
    return np.diag(B_p) + A*np.eye(P-1,k=1) + A*np.eye(P-1,k=-1)

H_D = H(V_h, X)

#main loop for CN
I=np.eye(P-1) #rank-n Identity matrix

T=np.linspace(0,T0,N)
PsiVsT = np.zeros([N, P-1], dtype=complex)
EnergyVsT = np.zeros([N])
NormVsT = np.zeros([N])
PsiVsT[0] = Psi

R_mat=(I - 1j*τ/(2*h_bar)*H_D)
L_mat=(I + 1j*τ/(2*h_bar)*H_D)

for n in range(N-1):
    v = np.matmul(R_mat,PsiVsT[n]) #v=RΨ^n
    PsiNew = np.linalg.solve(L_mat,v) #LΨ^(n+1)=v             
    PsiVsT[n+1] = PsiNew

for n in range(N):
    EnergyVsT[n] = expval(H_D, PsiVsT[n])
    NormVsT[n] = expval(I, PsiVsT[n])

#plots/animations
plt.figure()
plt.plot(T, EnergyVsT, 'k')
plt.xlabel('Time (s)')
plt.ylabel('Energy')
plt.title('Energy vs Time(harmonic oscillator)')
plt.savefig('EnergyVsTimeQ1c')

plt.figure()
plt.plot(T, NormVsT, 'k')
plt.xlabel('Time (s)')
plt.ylabel('Normalization')
plt.title('Normalization vs Time(harmonic oscillator)')
plt.savefig('NormVsTimeQ1c')

animate(T, X, np.real(np.conj(PsiVsT)*PsiVsT), r'Probability Density $|\psi|^2$', 'ProbDensQ1c')

#Q1d. using Double well potential

N0 = 6000
T0 = N0*τ0
N = 300
τ=T0/N
V_0=6*10**-17
x_1=L/4

#initial condition
σ=L/25
κ=500/L
x_0=L/3
Ψ_0=1

def psi_0(x):
    return Ψ_0*np.exp(-(x-x_0)**2/(4*σ**2) +1j*κ*x)

#TDSE spatial discretization
x_p=np.linspace(-L/2,L/2,P+1)
#print(np.size(x_p))
psi=np.array(list(map(psi_0,x_p)),complex)
psi[0]=psi[P]=0

Psi = normalized(psi[1:P])



def V_d(X): #double well potential
    return V_0*(X**2/x_1**2 -1)**2

def H(V, X):
    #Hamiltonian matrix
    A = -h_bar**2/(2*m*a**2)
    B_p = V(X)-2*A #V(x)=0 in square well
    
    return np.diag(B_p) + A*np.eye(P-1,k=1) + A*np.eye(P-1,k=-1)

H_D = H(V_d, X)

#main loop for CN
I=np.eye(P-1) #rank-n Identity matrix

T=np.linspace(0,T0,N)
PsiVsT = np.zeros([N, P-1], dtype=complex)
EnergyVsT = np.zeros([N])
NormVsT = np.zeros([N])
PsiVsT[0] = Psi

R_mat=(I - 1j*τ/(2*h_bar)*H_D)
L_mat=(I + 1j*τ/(2*h_bar)*H_D)

for n in range(N-1):
    v = np.matmul(R_mat,PsiVsT[n]) #v=RΨ^n
    PsiNew = np.linalg.solve(L_mat,v) #LΨ^(n+1)=v             
    PsiVsT[n+1] = PsiNew

for n in range(N):
    EnergyVsT[n] = expval(H_D, PsiVsT[n])
    NormVsT[n] = expval(I, PsiVsT[n])

#plots/animations
plt.figure()
plt.plot(T, EnergyVsT, 'k')
plt.xlabel('Time (s)')
plt.ylabel('Energy')
plt.title('Energy vs Time(double-well potential')
plt.savefig('EnergyVsTimeQ1d')

plt.figure()
plt.plot(T, NormVsT, 'k')
plt.xlabel('Time (s)')
plt.ylabel('Normalization')
plt.title('Normalization vs Time(double-well potential)')
plt.savefig('NormVsTimeQd')

animate(T, X, np.real(np.conj(PsiVsT)*PsiVsT), r'Probability Density $|\psi|^2$', 'ProbDensQ1d')