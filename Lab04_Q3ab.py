# -*- coding: utf-8 -*-
"""
Created on Wed Oct 7 10:43:46 2020

@author: Brendon McHugh

Question 3a,b of Lab04

Here, we investigate and compare the performance of both the relaxation 
method and overrelaxation method. 
"""

"""
Question 3, Part (a)

Here, we investigate the performance of the relaxation method for the 
equation x = 1-e^(-cx), as a function of the parameter c. 
"""

# Exercise 6.10 (a)
# Import libraries
import numpy as np
import matplotlib.pyplot as plt


def relaxation(F, Fd, x0, tolerance):
    """
    Calculates solution to x = F(x) using relaxation method to a given 
    tolerance, starting at estimate x0.

    Parameters
    ----------
    F : function (float): float
        Function to which relaxation method is applied.
    Fd : function (float): float
        Derivative of F.
    x0 : float
        Initial estimate of x.
    tolerance : float
        Desired tolerance of final value.

    Returns
    -------
    x1 : float
        Estimate of solution to x = F(x).
    error : float
        Final (estimated) error of x.
    count : int
        Number of iterations required to produce estimate.

    """
    count = 0
    x1 = x0
    error = tolerance + 1
    while (np.abs(error)>tolerance):
        x0, x1 = x1, F(x1)
        if (np.abs(Fd(x1)) > 1e-6):
            error = (x0-x1)/(1-1/Fd(x1))
        else:
            error = 0
        count+=1
    return x1, error, count

# Define our function and its derivative for any given c. 
f = lambda x,c : 1 - np.exp(-c*x)
fd = lambda x,c : c*np.exp(-c*x)

# For c = 2
f1 = lambda x: f(x,2)
fd1 = lambda x: fd(x,2)

# Calculating estimate for x when c = 2, to a tolerance of 10^-6
tolerance = 1e-6
x = 1
xf, error, count = relaxation(f1, fd1, x, tolerance) 

print('For c = 2, Relaxation method calculates the estimate: ', xf, ' with error: ', error)

# Exercise 6.10 (b)

# Create array of c values and x values
c_array = np.linspace(0, 3, 301)
x_array = np.zeros(301)

# Calculate x values using relaxation method
for i in range(len(c_array)):
    c = c_array[i]
    f1 = lambda x: f(x,c)
    fd1 = lambda x: fd(x,c)
    
    tolerance = 1e-6
    x = 1
    
    xf, error, count = relaxation(f1, fd1, x, tolerance) 
    x_array[i] = xf
    #print("%.2f\t%.5f\t%d"%(c, xf, count))

# Plotting x as a function of c
plt.figure()
plt.plot(c_array, x_array, 'k-')
plt.grid(True)
plt.xlabel('c')
plt.ylabel('x estimate')
plt.title('Relaxation Estimate of x as a Function of Parameter c')
plt.savefig("Relaxation_c_vs_x")


"""
Question 3, Part (b)

Here, we investigate the performance of the overrelaxation method for
the same function. 
"""
# Exercise 6.11 (b)
# For c = 2
f1 = lambda x: f(x,2)
fd1 = lambda x: fd(x,2)

# Calculating estimate for x when c = 2, to a tolerance of 10^-6
tolerance = 1e-6
x = 1
xf, error, count = relaxation(f1, fd1, x, tolerance) 

print('For c = 2, Relaxation method converges in ', count, 'iterations.')

# Exercise 6.11 (c)
def overrelaxation(F, Fd, omega, x0, tolerance):
    """
    Calculates solution to x = F(x) using overrelaxation method to a given 
    tolerance, starting at estimate x0, with overrelaxation parameter omega.

    Parameters
    ----------
    F : function (float): float
        Function to which relaxation method is applied.
    Fd : function (float): float
        Derivative of F.
    omega : float
        Overrelaxation parameter.
    x0 : float
        Initial estimate of x.
    tolerance : float
        Desired tolerance of final value.

    Returns
    -------
    x1 : float
        Estimate of solution to x = F(x).
    error : float
        Final (estimated) error of x.
    count : int
        Number of iterations required to produce estimate.

    """
    count = 0
    x1 = x0
    error = tolerance + 1
    while (np.abs(error)>tolerance):
        x0, x1 = x1, x1+(1+omega)*(F(x1)-x1)
        error = (x0-x1)/(1-1/((1+omega)*Fd(x1)-omega))
        count+=1
    return x1, error, count


# Create array of omega values and count values
omega_array = np.linspace(0, 1, 101)
count_array = np.zeros(101)

# Calculate number of iterations required for each value of omega
c = 2
f1 = lambda x: f(x,c)
fd1 = lambda x: fd(x,c)
for i in range(len(omega_array)):
    omega = omega_array[i]
    tolerance = 1e-6
    x = 1
    xf, error, count = overrelaxation(f1, fd1, omega, x, tolerance) 
    count_array[i] = count

# Plotting count as a function of omega
plt.figure()
plt.plot(omega_array, count_array, 'k-')
plt.grid(True)
plt.xlabel(r'$\omega$')
plt.ylabel('Number of Iterations')
plt.title(r'Number of Iterations vs $\omega$ Parameter')
plt.savefig("Overrelaxation_count_vs_omega")
