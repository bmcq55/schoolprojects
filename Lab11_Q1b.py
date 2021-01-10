# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 13:17:31 2020

@author: BMCQ
"""

import numpy as np
from math import sqrt,exp
from numpy import empty
from numpy.random import random,normal
import matplotlib.pyplot as plt

#Lab11_Q1b: Global minimum of a function
#part i)
#pseudocode:
#   (step 1): make function we are trying to find global minimum of, and if any, the constraints
#   (step 2): assign initial conditions/constants for the simulated annealling and make empty array to store x/y values
#   (step 3): make while loop that goes to the ground state of f/g, ie. the minimum 
#   (step 4): plot results for x and y separately
def f(x,y):
    '''function we are trying to find global minimum of'''
    return x**2-np.cos(4*np.pi*x)+(y-1)**2

#Gaussian distribution random number param.
mu=0 
sigma=1

#intial conditions/constants
Tmax = 1 #max temperature
Tmin = 1e-3 #minimum temperature
tau = 10e4 #cooling constant
fxy = f(2,2) #initial value of f(x,y)
t = 0
T = Tmax
x = 2
y = 2
xplot=[] #stores arrays of the "state" of the system
yplot=[]

#main loop
while T>Tmin:
 	
    t+=1
    T = Tmax*exp(-t/tau)
 	
    oldx = x
    oldy = y
    oldfxy = fxy
    
    r=np.random.normal(mu,sigma) #random number generator from gaussian distribution
    x+= r
    r=np.random.normal(mu,sigma)
    y+= r
    fxy = f(x,y)
    delta_fxy = fxy - oldfxy
 	
    if random() > exp(-delta_fxy/T): #cooling monte carlo simulation
        x=oldx
        y=oldy
        fxy= oldfxy

    xplot.append(x)
    yplot.append(y)
    
#plots
plt.figure()
plt.plot(xplot,'.')
plt.title('Finding Global Mimimum of f(x,y) Using Simulated Annealing(x value) ')
plt.xlabel('time')
plt.ylabel('x-value')
plt.show()
plt.savefig('global min f x value.png')

plt.figure()
plt.plot(yplot,'.')
plt.title('Finding Global Mimimum of f(x,y) Using Simulated Annealing(y value) ')
plt.xlabel('time')
plt.ylabel('y-value')
plt.show()
plt.savefig('global min f y value.png')

print('x = {}, y={} and f(x,y) = {}'.format(x,y,fxy))
    
#part ii)

def g(x,y):
    '''function we want to find global minimum of subject to x and y constraints,
    otherwise return large number outside of range of optimization'''
    if x>0 and x<100 and y>-20 and y<20:
        return np.cos(x)+np.cos(np.sqrt(2)*x)+np.cos(np.sqrt(3)*x)+(y-1)**2
    else:
        return 1e10
    
#print(g(2,2))

gxy = g(2,2)
t = 0
T = Tmax
x = 2
y = 2
x2plot=[]
y2plot=[]

#main loop
while T>Tmin:
 	
    t+=1
    T = Tmax*exp(-t/tau)
 	
    oldx = x
    oldy = y
    oldgxy = gxy
    
    r=np.random.normal(mu,sigma)
    x+= r
    r=np.random.normal(mu,sigma)
    y+= r
    gxy = g(x,y)
    delta_gxy = gxy - oldgxy
 	
    if random() > exp(-delta_gxy/T): #cooling
        x=oldx
        y=oldy
        gxy= oldgxy

    x2plot.append(x)
    y2plot.append(y)

#plots    
plt.figure()
plt.plot(x2plot,'.')
plt.title('Finding Global  of g(x,y) Using Simulated Annealing(x value) ')
plt.xlabel('time')
plt.ylabel('x-value')
plt.show()
plt.savefig('global min g x value.png')

plt.figure()
plt.plot(y2plot,'.')
plt.title('Finding Global Mimimum of g(x,y) Using Simulated Annealing(y value) ')
plt.xlabel('time')
plt.ylabel('y-value')
plt.show()
plt.savefig('global min g y value.png')

print('x = {}, y={} and g(x,y) = {}'.format(x,y,gxy))
    