# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 19:58:48 2020

@author: Brendon
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import math

#Lab 2,Q2ai
#Evaluating Dawson's function using trapezoidal and Simpson method for N=8 slices
# evaluate using scipy.special.dawsn 
#compare the three values

def trap(N): #trapezoidal integration method

    a=0.0
    b=4.0
    h=(b-a)/N
    x=np.linspace(a,b, N+1) #N+1 gives N intervals b/w a and b
    
    y=np.e**(x**2)
    y_r=y[1:] #right endpoint
    y_l=y[:-1]#left endpoint
    s=h/2 * np.sum(y_r+y_l)
    #return s
        
    print("the trapezoidal method gives the value", s*np.e**(-16), "for N=",N)

#evaluating with trapezoidal method N=8
trap(8)

def simps(N):#Simpson integration method
    
    if N % 2 == 1:
        raise ValueError("N must be an even integer.")
    a=0.0
    b=4.0
    h=(b-a)/N
    
    x=np.linspace(0,4,N+1)
    
    y = np.e**(x**2)
    S = h/3 * np.sum(y[0:-1:2] + 4*y[1::2] + y[2::2])
    value_2=np.e**(-16)*S

    print("the Simpson method for gives the value", value_2, "for N=",N)

#evaluating with simpson method N=8
simps(8)

#dawson using scipy   
x2=np.linspace(0,4,8, endpoint=True)   
value=special.dawsn(x2)
print("the scipy method for gives the value", value)

#Lab 2, Q2aii
#increase N for each method until value of integral is within 9 decimal point error of scipy value
#time each of the three methods

from time import time

start1=time()
trap(2**16)
end1=time()
diff1=end1-start1

print("the computational time of the trapezoidal method is",diff1,"seconds", "for N=65536")   

#simps(2**9)#does not give accuracy to 9th decimal point

start2=time() #timing the Simpson method
simps(2**10)
end2=time()

diff2=end2 - start2
print("the computational time of Simpsons method is", diff2, "seconds", "for N=1024")



x2=np.linspace(0,4,8, endpoint=True)   
start3=time() ##timing the scipy method

value=special.dawsn(x2)
print("the scipy method for gives the value", value)
end3=time()

diff3=end3-start3
print("the computational time of the scipy method is",diff3, "seconds")

#Lab 2,Q2aiii: errors on trapezoidal and Simpson's integration technique

h=(4.0-0.0)/(0.5*32)

I_2N=trap(64) #error in trapezoidal method for N=32
I_N=trap(32)
err_trap=(0.13194038496790617-0.13958009092677318)*1/(3) #err_trap=Ch^2=1/3(I_2N-I_N)
print("the error in the trapezoidal method is", abs(err_trap))

I_2K=simps(64) #error in Simpson's method for K=32
I_K=simps(32)
err_simps=(0.12939381631495048-0.13001118158375183)*1/(15) #err_simps=ch^4=1/15(I_2K-I_K)
print("the error in Simpson's method is",abs(err_simps))

av=(0.001001119613647461+0.0010135173797607422+0)/3
print(av)
