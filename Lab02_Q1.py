# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:31:27 2020

@author: Nicholas Sullivan

Lab Partner: Brendon McHugh

Question 1 of Lab02

Here, we investigate the performance of the forward difference method 
and central difference method as approximations for the derivative
"""

"""
Question 1, Part (a)

Create function to calculate the forward difference method for a given step
size, and summarize the results when applied to the function f(x) = exp(-x^2)
"""

# Import numpy and matplotlib
import numpy as np
import matplotlib.pyplot as plt

# Define forward difference function
def forward_diff(f, x, h):
    """ 
    Function forward_diff calculates the forward difference of a given function
    
    INPUT:
        f - function to be differentiated
        x - point to be differentiated at
        h - forward increment
        
    OUTPUT:
        Forward difference of f between x and x+h.
    """
    
    return (f(x+h) - f(x))/h

# Define function f(x) = exp(-x^2)
f = lambda x: np.exp(-x**2)

# Define x = 0.5 as the point to differentiate at
x = 0.5

# Create array of h values from 10^-16 to 10^0
h_array = np.power(10.0, np.arange(-16, 1))

# Create empty array for storing forward difference values
forward_diff_array = np.empty(17)

# For each value in array of h values, calculate forward difference of f at x
for i in range(len(h_array)):
    h = h_array[i]
    forward_diff_array[i] = forward_diff(f, x, h)


"""
Question 1, Part (b)

Calculate the absolute error of the forward difference method as a function 
of step size, and print results
"""

# Define function g(x) = -2x*exp(-x^2) as derivative of f
g = lambda x: -2*x*np.exp(-x**2)

# Calculate true derivative of F at x
true_derivative = g(x)

# Calculate absolute error of forward difference values
forward_diff_error = np.abs(forward_diff_array - true_derivative)

# Print Results
print("Forward Difference")
print("Step Size \tDifference \tError")
for i in range(len(h_array)):
    print("%.2E \t%.6f \t%.2E"%(h_array[i],forward_diff_array[i], forward_diff_error[i]))


"""
Question 1, Part (c)

Plot the absolute error of the forward difference method as a function of
step size
"""

# Create a figure with gridlines
plt.figure()
plt.grid(True)

# Create log-log plot of the forward difference error vs step size
plt.loglog(h_array, forward_diff_error, 'k-', label='$X_0$ = 0.1')

# Label the x and y axes
plt.xlabel("Step Size")
plt.ylabel("Absolute Error")

# Set title of figure
plt.suptitle("Absolute Error of Forward Difference vs. Step Size")

# Save figure to file
plt.savefig('AbsError_StepSize_1c.png')

"""
Question 1, Part (d)

Create function to calculate the central difference method for a given step
size, calculate the absolute error when applied to the same function f(x),
and plot the comparison with the forward difference method
"""

# Define central difference function
def central_diff(f, x, h):
    """ 
    Function central_diff calculates the central difference of a given function
    
    INPUT:
        f - function to be differentiated
        x - point to be differentiated at
        h - central increment
        
    OUTPUT:
        Central difference of f between x-h/2 and x+h/2.
    """
    
    return (f(x+h/2) - f(x-h/2))/h

# Create empty array for storing central difference values
central_diff_array = np.empty(17)

# For each value in array of h values, calculate central difference of f at x
for i in range(len(h_array)):
    h = h_array[i]
    central_diff_array[i] = central_diff(f, x, h)

# Calculate absolute error of central difference values
central_diff_error = np.abs(central_diff_array - true_derivative)

# Print Results
print("Central Difference")
print("Step Size \tDifference \tError")
for i in range(len(h_array)):
    print("%.2E \t%.6f \t%.2E"%(h_array[i],central_diff_array[i], central_diff_error[i]))

# Create a figure with gridlines
plt.figure()
plt.grid(True)

# Create log-log plot of the forward difference error vs step size
plt.loglog(h_array, forward_diff_error, 'k-', label='Forward Difference')
# Create log-log plot of the central difference error vs step size
plt.loglog(h_array, central_diff_error, 'r-', label='Central Difference')

# Create legend
plt.legend()

# Label the x and y axes
plt.xlabel("Step Size")
plt.ylabel("Absolute Error")

# Set title of figure
plt.suptitle("Absolute Error of Forward and Central Difference vs. Step Size")

# Save figure to file
plt.savefig('AbsError_StepSize_1d.png')
