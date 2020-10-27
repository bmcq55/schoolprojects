# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 09:28:43 2020

@author: Brendon McHugh

Question 2 of Lab01

Investigation of the evolution of a population of insects that grows 
according to the logistic map x_p -> r(1-x_p)x_p. Further exploration of 
bifurcation diagrams and chaotic systems.

"""

"""
Question 2, Part (a)
"""

# Import numpy for array manipulations
# Import matplotlib for plotting functions

# From keyboard, read number of years to simulate and save as variable pmax
# From keyboard, read initial number of insects and save as variable n0
# From keyboard, read maxmimum number of insects and save as variable nmax
# From keyboard, read maximum reproduction rate and save as variable r


# Define x0 = n0/nmax as the population as a fraction of nmax
# Initialize xp, array of length pmax+1 filled with zeros, to store population data
# Set xp[0] = x0 as the initial population


# Define for loop to update population for years p from 1 to pmax
    # Update population as xp[p] = r*(1-xp[p-1])*xp[p-1]

# Plot array xp with respect to array 0:pmax

"""
Question 2, Part (b)
"""

# Import numpy for array manipulations
import numpy as np
# Import matplotlib for plotting functions
import matplotlib.pyplot as plt

def update_population(x0, r, pmax):
    """ Define function update_population accepting as arguments the 
    initial normalized population x0, the maximum reproduction rate r,
    and the number of years to simulate pmax. It will return an array of 
    length pmax+1 that contains the population at each year from 0 to pmax."""
    
    # Initialize xp, array of length p+1 filled with zeros,
    # to store population data
    xp = np.zeros(pmax+1)
    
    # Set xp[0] = x0 as the initial population
    xp[0] = x0

    # Define for loop to update population for years p in range 1:pmax
    for p in range(1,pmax+1):
        # Update population as xp[p] = r*(1-xp[p-1])*xp[p-1]
        xp[p] = r*(1-xp[p-1])*xp[p-1]
    
    # Return array of population data
    return xp

"""
Question 2, Part (c)
"""

X0 = 0.1
R_array = np.arange(8)*(2/8) + 2
PMAX = 50

fig, axs = plt.subplots(2, 4)


for n in range(len(R_array)):
    R = R_array[n]
    pop_array = update_population(X0, R, PMAX)
    time_array = np.arange(PMAX+1)
    print(n//4," ", n%4)
    axs[n//4, n%4].plot(time_array, pop_array, 'k-')
    axs[n//4, n%4].set_title('r = %3.2f'%R)
    axs[n//4, n%4].grid(True)
    
for ax in axs.flat:
    ax.set(xlabel="Time (Years)", ylabel="Normalized Population")

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

fig.suptitle("Normalized Insect Population with Respect to Time in Years")

fig.savefig('PopVsTime.png')


"""
Question 2, Part (d)
"""

X0 = 0.1
R_array = np.arange(2,4,0.005)
PMAX = 2000

plt.figure()
plt.grid(True)

for n in range(len(R_array)):
    R = R_array[n]
    pop_array = update_population(X0, R, PMAX)
    time_array = np.arange(PMAX+1)
    if (R<3):
        plt.scatter((R*np.ones(PMAX+1))[-100:-1], pop_array[-100:-1], s = 0.01, color=(0,0,0))
    else:
        plt.scatter((R*np.ones(PMAX+1))[-1000:-1], pop_array[-1000:-1], s = 0.01, color=(0,0,0))
    
plt.xlabel("Population Growth Parameter")
plt.ylabel("Final Normalized Population")

plt.suptitle("Bifurcation Diagram")

plt.savefig('Bifurc_d.png')


"""
Question 2, Part (e)
"""

X0 = 0.1
R_array = np.arange(3.738,3.745,1e-5)
PMAX = 2000

plt.figure()
plt.grid(True)

for n in range(len(R_array)):
    R = R_array[n]
    pop_array = update_population(X0, R, PMAX)
    time_array = np.arange(PMAX+1)
    if (R<3):
        plt.scatter((R*np.ones(PMAX+1))[-100:-1], pop_array[-100:-1], s = 0.01, color=(0,0,0))
    else:
        plt.scatter((R*np.ones(PMAX+1))[-1000:-1], pop_array[-1000:-1], s = 0.01, color=(0,0,0))
    
plt.xlabel("Population Growth Parameter")
plt.ylabel("Final Normalized Population")

plt.suptitle("Bifurcation Diagram")

plt.savefig('Bifurc_e.png')


"""
Question 2, Part (f)
"""

X0_1 = 0.1
eps = 0.0001
X0_2 = X0_1 + eps
R = 3.71
PMAX = 50

time_array = np.arange(PMAX+1)

plt.figure()
plt.grid(True)

pop_array_1 = update_population(X0_1, R, PMAX)
plt.plot(time_array, pop_array_1, 'k-', label='$X_0$ = 0.1')

pop_array_2 = update_population(X0_2, R, PMAX)
plt.plot(time_array, pop_array_2, 'r-', label='$X_0$ = 0.1001')

    
plt.xlabel("Time (Years)")
plt.ylabel("Normalized Population")
plt.legend()

fig.suptitle("Population Evolution Between Similar Initial Conditions")

fig.savefig('PopVsTime_f.png')


"""
Question 2, Part (g)
"""

pop_diff = np.abs(pop_array_1 - pop_array_2)
a = eps
lbda = 0.4
spread_fit = a*np.exp(lbda*time_array)

plt.figure()
plt.grid(True)

plt.semilogy(time_array, pop_diff, 'k-', label='Population Divergence')
plt.semilogy(time_array, spread_fit, 'r--', label='Exponential Fit')
    
plt.xlabel("Time (Years)")
plt.ylabel("Normalized Population")
plt.legend()

fig.suptitle("Difference in Population Evolution Between Similar Initial Conditions")

fig.savefig('PopVsTime_g.png')

