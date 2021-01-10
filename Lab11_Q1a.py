# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 17:15:54 2020

@author: BMCQ
"""

#Lab_11, Q1: Simulated Annealing 
import numpy as np
from math import sqrt,exp
from numpy import empty
from random import random,randrange,seed
import matplotlib.pyplot as plt

N = 25
R = 0.02
Tmax = 10.0
Tmin = 1e-3
tau = 1e5

# Function to calculate the magnitude of a vector
def mag(dr):
    return sqrt((dr[0])**2+(dr[1])**2)

# Function to calculate the total length of the tour
def distance(path, r):
    s = 0.0
    for i in range(N):
        s += mag(r[path[i+1]] - r[path[i]])
    return s

# Seed value for initial city locations
c_s=10
seed(c_s)

# Choose N city locations
r = empty([N, 2], float)
for i in range(N):
    r[i,0]=random()
    r[i,1]= random()

# Choose initial path, including first point repeated at the end
path = np.arange(N+1, dtype=int)
path[N] = 0

# Calculate initial distance
D = distance(path, r)
print(r'The initial distance is %.4f'%(D))

# Plot intiial path
plt.figure() 
plt.scatter(r[path, 0], r[path, 1],color='red')
plt.plot(r[path, 0], r[path, 1])
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Initial Travelling Salesman Path with Initial City seed{}'.format(c_s))
plt.savefig('TSMP initial path seed{}'.format(c_s))

# Seed for simulation - change to get different results
s_s=25
seed(s_s)
# Main loop
t = 0
T = Tmax
while T>Tmin:

    # Cooling
    t += 1
    T = Tmax*exp(-t/tau)

    # Choose two cities to swap (excluding the first) and make sure they are distinct
    i,j = randrange(1,N),randrange(1,N)
    while i==j:
        i,j = randrange(1,N),randrange(1,N)

    # Swap them and calculate the change in distance
    oldD = D
    path[i],path[j] = path[j],path[i]
    D = distance(path, r)
    deltaD = D - oldD

    # If the move is rejected, swap them back again
    if random()>exp(-deltaD/T):
        path[i],path[j] = path[j],path[i]
        D = oldD

# Calculate final distance
D = distance(path, r)
print(r'The final distance is %.4f'%(D))

# Plot final path
plt.figure() 
plt.scatter(r[path, 0], r[path, 1],color='red')
plt.plot(r[path, 0], r[path, 1])
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Final Travelling Salesman Path with Initial City Seed{} and Simulation Seed{}'.format(c_s,s_s))
plt.savefig('TSMP{} and{}.png'.format(c_s,s_s))
