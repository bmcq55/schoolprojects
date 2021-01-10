# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 08:04:26 2020

@author: BMCQ
"""

#Lab10 
#Q2: Volume of hypersphere (Newman 10.7)

#pseudocode
#(step 1) assign constants 
#(step 2) define hypersphere equations
#(step 3) do the mean value method in higher dimensional monte carlo integration to find volume
#(step 4) calculate the error in the integral via variance/std. deviation (eq. 10.31 and 10.32)

import numpy as np
import matplotlib.pyplot as plt

#step 1
N=1000000 #number of sample points
dim=10 #dimension of hypersphere
a=-1
b=1
k=0
k2=0


#step 2
def f(x):
    '''INPUT: ith dimensional vector, x
    OUTPUT: hypersphere in i dimensions: R=sqrt(sum(x_i^2))'''
    R_sq = np.zeros(x.shape[1],float)
 	
    for ri in x:
            R_sq+=ri**2
 	
    return R_sq<=1 #returns value on or inside sphere
 
x=(b-a)*np.random.random((dim,N)) -1
fx=f(x)
k=sum(fx) #average
k2=sum(fx**2) #variance
    
#step 3 
I=(b-a)**dim/N*np.sum(fx) #eq 10.33 in 10 dim
print('integal=',I)

#step 4
var = k2/N - (k/N)**2 # variance <f**2> - <f>**2
sigma_MV = (b-a)**dim*np.sqrt(var/N) # MV stands for Mean Value
print('error=',sigma_MV)



