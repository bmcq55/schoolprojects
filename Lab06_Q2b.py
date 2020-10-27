# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 21:13:32 2020

@author: Brendon McHugh

Question 2 of Lab06

Here, we have pseudocode for two particles interacting with Lennard-Jones
potential.
"""

"""
Question 2, Part (b)

Pseudocode for simulating interaction of two particles with Lennard-Jones
potential, using the verlet method
"""

# Import numpy and matplotlib libraries

# Define function (acceleration) that takes the x and y coordinates of 
# both particles as input, and returns the x and y components of 
# acceleration for both particles.
    # rx = x1 - x2
    # ry = y1 - y2
    # r = sqrt(rx^2 + ry^2)
    # coeff = 4/r2*(12/r^12 - 6/r^6)
    # ax = coeff*rx
    # ay = coeff*ry
    # return ax, ay, -ax, -ay

# Define function (verlet) that takes as input the acceleration function,
# the initial time, position and velocity, as well as the time step and number
# of iterations. Returns array of time, position and velocity of each particle
# at each iteration step
    # T = t0:t0+h*N
    # v = v0 + 0.5*h*acceleration(r0,t0)
    # R = zero array length N+1
    # V = zero array length N+1
    # R[0] = r0
    
    # for i from 1:N
        # R[i+1] = R[i] + h*v
        # k = h*acceleration(R[i+1], T[i+1])
        # V[i+1] = v + 0.5*k
        # v = v + k
    # return T, R, V

# Call verlet for (x1, y1, x2, y2) = (4, 4, 5.2, 4), h = 0.01, N = 100
# and plot trajectory

# Call verlet for (x1, y1, x2, y2) = (4.5, 4, 5.2, 4), h = 0.01, N = 100
# and plot trajectory

# Call verlet for (x1, y1, x2, y2) = (2, 3, 3.5, 4.4), h = 0.01, N = 100
# and plot trajectory
