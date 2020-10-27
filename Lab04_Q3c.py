# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:26:45 2020

@author: Brendon
"""

import numpy as np
import matplotlib.pyplot as plt
from time import time
import scipy

#textbook problem 6.13(b)
#pseudocode:
#(step 1) given initial pair of points x1, x2, check that f(x1) and f(x2) have opposite signs and give a target accurary.
#(step 2) calculate midpoint x_prime=0.5*(x1+x2)
#(step 3) if f(x_prime) has same sign as f(x1), the x1=x_prime. Otherwise, x2=x_prime
#(step 4) if |x1-x2| > Ɛ, repeat step 2. Else, calculate 0.5*(x1+x2) and print out 

Ɛ=10**-6 #accurary
def f(x):
    return 5*np.e**(-x)+ x - 5 #nonlinear equation we want to find roots of

#plot of function we are trying to find the roots of
x=np.linspace(-1,8)
plt.figure()
plt.plot(x,f(x))
plt.title("Non-Linear Equation for Max Emitted Radition")
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.savefig("Wien's Displacement.png")
plt.show()

def binary_search(f,x1,x2,Ɛ):
    '''INPUT: non-linear equation, pair of points, accurary
    OUTPUT: roots of the equation'''
    
    #step 1
    if (f(x1)<0 and f(x2)<0) or (f(x1)>0 and f(x2)>0):
        return ValueError("f(x1) and f(x2) must have opposite signs")
    else:
        pass
    
    #step 2
    def same_sign(x,y):
        if x < 0 and y < 0 or x > 0 and y > 0:
            return True
        else:
            return False
    #step 3   
    def midpoint(x,y):
        return 0.5*(x+y) #value of x`
    
    count = 0
    while np.abs(x1-x2) > Ɛ:
        print("iteration Count:{}".format(count))
        count += 1 
        
        x_prime= midpoint(x1,x2)
        
        if same_sign(f(x1),f(x_prime)):
            x1=x_prime
        elif same_sign(f(x_prime),f(x2)):
            x2=x_prime
        elif np.abs(x_prime) < Ɛ:
            return x_prime
         
    return midpoint(x1,x2)
 
print("the binary search method gives the value of:" ,binary_search(f,2,7,Ɛ))


#exercise 6.13b using the relaxation method
#pseudocode
#(step 1): make an initial guess...an educated guess in this case
#(step 2) iterate until x converges to the root we found using the binary search method

x=2.0 #step 1
for k in range(8): #step 2
    x=5-5*np.e**-x
    print("the relaxation method gave a value of", x , "for iteration:",k)

#exercise 6.13b using Newton's method
#pseudocode
#(Step 1) Find f'(x) and define equation
#(step 2) Choose an initial value and use the formula x_prime=x-delta_x repeatedly to get better and better estimates of the root
#(step 3) if |x_prime-x| < Ɛ, stop iteration and output the roots/number of iterations

def Newton(x,Ɛ):
    '''INPUT: x is the initial value to start iteration, and Ɛ is the accuracy
    OUTPUT: the roots of the nonlinear equation and number of iterations'''
    for i in range(20):
        x_prime=x-(5*np.e**-x + x - 5)/(5*np.e**-x + 1) #step 1/2
        if abs(x_prime-x) < Ɛ: #step 3
            break
        x=x_prime
    print("the root is", x_prime, "after", i ,"iterations")
  
Newton(2,10**-6 )


#6.13c Estimate of surface temperature of the sun by optical pyrometry

λ=502*10**-9 #wavelength of suns emitted radiation in meters
c=299792458 #speed of light m/s
K_b=1.380649*10**-23 #boltzman constant Joules/Kelvin
h=6.62607015*10**-34 #planks constant J*s
x=4.9651141

b=h*c / (K_b*x)
print("the value of wien's displacement constant is", b)

#temperature of the sun
T=b/λ
print("The temperature of the sun is:", T,"Kelvin")
  

    

