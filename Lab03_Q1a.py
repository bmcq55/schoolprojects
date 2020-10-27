# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:25:28 2020

@author: Brendon McHugh

Question 1 of Lab03

Here, we compare the performance of the Trapezoidal method,
Simpson's method, and Gaussian quadrature for calculating 
Dawson's function at x=4.
"""

"""
Question 1, Part (a) i

Here, we use the trapezoidal, Simpson and Gaussian methods
for calculating Dawson's function for x=4.
"""
# Import numpy, matplotlib, scipy, gaussxw
import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from gaussxw import gaussxwab


def trap_samples_weights(a, b, N):
    """
    Function to calculate the sample points and weights
    for trapezoidal integration over interval (a,b).

    Parameters
    ----------
    a : float
        Minimum of integration interval.
    b : float
        Maximum of integration interval.
    N : int
        Number of sample points.

    Returns
    -------
    x : float array
        Array of sample points.
    w : float array
        Array of weights.

    """
    w = np.ones(N+1)
    w[0] = 0.5
    w[-1] = 0.5
    w = w*(b-a)/N
    x = np.linspace(a, b, N+1)
    return x, w

def simpson_samples_weights(a, b, N):
    """
    Function to calculate the sample points and weights
    for Simpson integration over interval (a,b).

    Parameters
    ----------
    a : float
        Minimum of integration interval.
    b : float
        Maximum of integration interval.
    N : int
        Number of sample points.

    Returns
    -------
    x : float array
        Array of sample points.
    w : float array
        Array of weights.

    """
    if N % 2 == 1:
        raise ValueError("N must be an even integer.")
    w = np.ones(N+1)/3
    w[2:-1:2] = 2*w[2:-1:2]
    w[1:-1:2] = 4*w[1:-1:2]
    w = w*(b-a)/N
    x = np.linspace(a, b, N+1)
    return x, w


def dawson_integral(samples, weights):
    """
    Function to calculate integral of Dawson function integral
    based on the input sample points and weights

    Parameters
    ----------
    samples : float array
        Array of sample points.
    weights : float array
        Array of weights.

    Returns
    -------
    integral : float
        Integral needed for calculating Dawson's function.

    """
    
    dawson_int = lambda t: np.exp(t**2)
    integral = np.sum(weights*dawson_int(samples))
    return integral

# Calculates Dawson function using trapezoidal method at x
def dawson_trap(x, N):
    a=0.0
    b=x
    s, w = trap_samples_weights(a, b, N)
    
    coeff = np.exp(-x**2)
    integral = dawson_integral(s,w)
    return coeff*integral

# Calculates Dawson function using Simpson method at x
def dawson_simpson(x, N):
    
    if N % 2 == 1:
        raise ValueError("N must be an even integer.")
    a=0.0
    b=x
    s, w = simpson_samples_weights(a, b, N)
    
    coeff = np.exp(-x**2)
    integral = dawson_integral(s,w)
    return coeff*integral

# Calculates Dawson function using Gaussian Quadrature at x
def dawson_gauss(x, N):
    
    a=0.0
    b=x
    s, w = gaussxwab(N+1, a, b)
    
    coeff = np.exp(-x**2)
    integral = dawson_integral(s,w)
    return coeff*integral

# Calculate estimate for x=4, N between 8 and 2048, and each 
# integration method
x = 4
Nvals = 2**np.arange(3,12)
L = len(Nvals)

trap_est = np.zeros(L)
simpson_est = np.zeros(L)
gauss_est = np.zeros(L)

# Calculate Dawson function using each method, print table
value=special.dawsn(4)
print("True value (according to scipy): %.10f"%value)
print("N \tTrapezoidal \tSimpson \tGauss")
for k in range(L):
    N = Nvals[k]
    trap_est[k] = dawson_trap(x, N)
    simpson_est[k] = dawson_simpson(x, N)
    gauss_est[k] = dawson_gauss(x, N)
    print("%d \t%.10f \t%.10f \t%.10f"%(N, trap_est[k], simpson_est[k], gauss_est[k]))


"""
Question 1, Part (a) ii

Here, we compare the error in the trapezoidal, Simpson 
and Gaussian methods for calculating Dawson's function 
for x=4.
"""

# Calculate true error for x=4, N between 8 and 2048, and each 
# integration method
x = 4
Nvals = 2**np.arange(3,12)
L = len(Nvals)

trap_error = np.zeros(L)
simpson_error = np.zeros(L)
gauss_error = np.zeros(L)

# Perform calculation using each method, print table
value=special.dawsn(4)
print()
print("Error")
print("N \tTrapezoidal \tSimpson \tGauss")

for k in range(L):
    N = Nvals[k]
    trap_error[k] = abs(dawson_trap(x, N) - value)
    simpson_error[k] = abs(dawson_simpson(x, N) - value)
    gauss_error[k] = abs(dawson_gauss(x, N) - value)
    print("%d \t%.4E \t%.4E \t%.4E"%(N, trap_error[k], simpson_error[k], gauss_error[k]))


# Calculate error estimate for x=4, N between 8 and 2048, and each 
# integration method
trap_error_est = np.zeros(L)
simpson_error_est = np.zeros(L)
gauss_error_est = np.zeros(L)

# Perform calculation using each method, print table
value=special.dawsn(4)
print()
print("Error Estimate")
print("N \tTrapezoidal \tSimpson \tGauss")

for k in range(L):
    N = Nvals[k]
    trap_error_est[k] = abs(dawson_trap(x, N) - dawson_trap(x, N//2))/3
    simpson_error_est[k] = abs(dawson_simpson(x, N) - dawson_simpson(x, N//2))/15
    gauss_error_est[k] = abs(dawson_gauss(x, N) - dawson_gauss(x, N//2))
    print("%d \t%.4E \t%.4E \t%.4E"%(N, trap_error_est[k], simpson_error_est[k], gauss_error_est[k]))

# Create a figure with gridlines
plt.figure()
plt.grid(True)

# Create log-log plot of the error and error estimate
# for each method
plt.loglog(Nvals, trap_error, 'b-', label='Trapezoidal Error')
plt.loglog(Nvals, trap_error_est, 'b--', label='Trapezoidal Error Estimate')

plt.loglog(Nvals, simpson_error, 'r-', label='Simpson Error')
plt.loglog(Nvals, simpson_error_est, 'r--', label='Simpson Error Estimate')

plt.loglog(Nvals, gauss_error, 'g-', label='Gaussian Error')
plt.loglog(Nvals, gauss_error_est, 'g--', label='Gaussian Error Estimate')

# Label the x and y axes
plt.xlabel("Step Number N")
plt.ylabel("Absolute Error")

plt.legend()

# Set title of figure
plt.suptitle("Absolute Error of Integration Methods vs. Step Number")

# Save figure to file
plt.savefig('AbsError_N_1aii.png')


