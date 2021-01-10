# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 09:38:34 2020

@author: Nicholas Sullivan

Lab Partner: Brendon McHugh

Here, we compare integral estimates using mean value method and importance sampling
"""

import numpy as np
import matplotlib.pyplot as plt
import random

"""
Part a

Estimate of function x**(-1/2)/(1+np.exp(x)) from 0 to 1 using mean value 
method and importance sampling
"""

def mean_value_method(f, a, b, N):
    """
    INPUT: integrand f, integral bounds a and b, sample number N
    OUTPUT: integral estimate and error estimate
    """
    sum_f = 0
    sum_f2 = 0
    
    # For each random sample, update the sum of f and f**2
    for i in range(N):
        x = a + (b-a)*random.random()
        fx = f(x)
        sum_f = sum_f + fx
        sum_f2 = sum_f2 + fx**2
    
    # Calculate integral estimate, error estimate
    mean_f = sum_f/N
    int_f = mean_f*(b-a)
    
    mean_f2 = sum_f2/N
    var_f = mean_f2 - mean_f**2
    sigma_f = (b-a)*np.sqrt(var_f/N)
    
    return int_f, sigma_f

def importance_sampling_method(f, w, x, N):
    """
    INPUT: integrand f, normalized weight function w, random number generator x, sample number N
    OUTPUT: integral estimate and error estimate
    """
    sum_g = 0
    sum_g2 = 0
    
    # For each random sample, update the sum of g and g**2, where g = f/w
    for i in range(N):
        x0 = x()
        fx = f(x0)
        wx = w(x0)
        gx = fx/wx
        sum_g = sum_g + gx
        sum_g2 = sum_g2 + gx**2
    
    # Calculate integral estimate, error estimate of f
    mean_g = sum_g/N
    int_f = mean_g
    
    mean_g2 = sum_g2/N
    var_g = mean_g2 - mean_g**2
    sigma_f = np.sqrt(var_g/N)
    
    return int_f, sigma_f

# Integrand
def F(x):
    return x**(-1/2)/(1+np.exp(x))

# Normalized weight function
def W(x):
    return x**(-1/2)/2

# Random number generator with distribution from w
def X():
    z = random.random()
    return z**2

# We want to compare the integral estimates over 100 trials
trials = 100

I_mean_array = np.zeros([trials])
I_imp_array = np.zeros([trials])

for i in range(trials):
    # Mean value method
    I_mean, s_mean = mean_value_method(F, 0, 1, 10000)
    
    # Importance sampling
    I_imp, s_imp = importance_sampling_method(F, W, X, 10000)
    
    I_mean_array[i] = I_mean
    I_imp_array[i] = I_imp

# Histogram bins
bins = 10
bin_min = 0.80
bin_max = 0.88 

# Mean value histogram
plt.figure()
plt.hist(I_mean_array, bins, range=[bin_min, bin_max])
plt.xlabel('Integral Estimate')
plt.ylabel('Frequency')
plt.title('Histogram of Integral Estimates using Mean Value Method')
plt.savefig('Q3a_mean')

# Importance sampling histogram
plt.figure()
plt.hist(I_imp_array, bins, range=[bin_min, bin_max])
plt.xlabel('Integral Estimate')
plt.ylabel('Frequency')
plt.title('Histogram of Integral Estimates using Importance Sampling')
plt.savefig('Q3a_imp')


"""
Part b

Estimate of function np.exp(-2*np.abs(x-5)) from 0 to 10 using mean value 
method and importance sampling
"""

# Integrand
def F(x):
    return np.exp(-2*np.abs(x-5))

# Normalized weight function
def W(x):
    return np.exp(-(x-5)**2/2)/np.sqrt(2*np.pi)

# Random number from gaussian distribution 
def X():
    mu = 5
    sigma = 1
    return np.random.normal(mu, sigma)

# We want to compare the integral estimates over 100 trials
trials = 100

I_mean_array = np.zeros([trials])
I_imp_array = np.zeros([trials])

for i in range(trials):
    # Mean value 
    I_mean, s_mean = mean_value_method(F, 0, 10, 10000)
    
    # Importance sampling
    I_imp, s_imp = importance_sampling_method(F, W, X, 10000)
    
    I_mean_array[i] = I_mean
    I_imp_array[i] = I_imp

# Histogram bins
bins = 10
bin_min = 0.9
bin_max = 1.1

# Mean value histogram
plt.figure()
plt.hist(I_mean_array, bins, range=[bin_min, bin_max])
plt.xlabel('Integral Estimate')
plt.ylabel('Frequency')
plt.title('Histogram of Integral Estimates using Mean Value Method')
plt.savefig('Q3b_mean')

# Importance Sampling histogram
plt.figure()
plt.hist(I_imp_array, bins, range=[bin_min, bin_max])
plt.xlabel('Integral Estimate')
plt.ylabel('Frequency')
plt.title('Histogram of Integral Estimates using Importance Sampling')
plt.savefig('Q3b_imp')

