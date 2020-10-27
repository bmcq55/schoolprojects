# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 21:13:32 2020

@author: Brendon McHugh

Question 3 of Lab06

Here, we use the verlet method to simulate the Lennard-Jones interactions of
N particles, both with and without periodic boundary conditions.

3a --> Line 23
3b --> Line 200
3c --> Line 333
Verlet Method --> Line 86
Acceleration Function --> Line 147

"""

"""
Question 3, Part (a)

Pseudocode for simulation of N particles without periodic boundary conditions
"""

# Import numpy and matplotlib libraries

# Define function (acceleration) that takes the x and y coordinates of 
# all N particles as input, and returns the x and y components of 
# acceleration for all particles.
    # ax = zero array of length N
    # ay = zero array of length N
    
    # For each pair (i,j):
        # Find distance between i,j
        # rx = x[i] - x[j]
        # ry = y[i] - y[j]
        # r = sqrt(rx^2 + ry^2)
        
        # Calculate accelerations
        # coeff = 4/r2*(12/r^12 - 6/r^6)
        # ax[i] += coeff*rx
        # ay[i] += coeff*ry
        # ax[j] += -coeff*rx
        # ay[j] += -coeff*ry
        
    # return ax, ay

# Define function (verlet) that takes as input the acceleration function,
# the initial time, position and velocity, as well as the time step (h) and number
# of iterations (N). Returns array of time, position and velocity of each particle
# at each iteration step
    # Initialize time array
    # T = t0:t0+h*N
    # v = v0 + 0.5*h*acceleration(r0,t0)
    # R = zero array length N+1
    # V = zero array length N+1
    # R[0] = r0
    
    # Perform verlet method for each time step
    # for i from 1:N
        # R[i+1] = R[i] + h*v
        # k = h*acceleration(R[i+1], T[i+1])
        # V[i+1] = v + 0.5*k
        # v = v + k
        
    # return T, R, V

# Define x0, y0 as 4x4 grid of positions
# Define v0 as zero array

# Call verlet for (x0, y0), v0, h = 0.01, N = 1000
# and plot trajectory

"""
Here, we implement the functions as described in the above pseudocode.

"""

# Import libraries
import numpy as np
import matplotlib.pyplot as plt

def verlet(f, t0, r0, v0, h, N):
    """
    Implements verlet method for a given acceleration function

    Parameters
    ----------
    f : function 
        Computes accleration given coordinates.
    t0 : float
        Initial time.
    r0 : float array
        Initial value of general coordinate array.
    v0 : float array
        Initial value of general velocity array.
    h : float
        Time step.
    N : int
        Number of iterations.

    Returns
    -------
    T : float array 
        Array of time values.
    R : 2d float array
        2D array of general coordinates at each time.
    V : 2d float array
        2D array of general velocities at each time.

    """
    
    # Intialize T, R and V arrays to hold values
    T = np.linspace(t0, t0+h*N, N+1)
    dof = len(r0)
    R = np.zeros([N+1, dof])
    R[0] = r0
    V = np.zeros([N+1, dof])
    V[0] = v0
    
    # Implement verlet's method for N time steps
    vc = v0 + 0.5*h*f(r0, T[0])
    for i in range(N):
        R[i+1] = R[i] + h*vc
        k = h*f(R[i+1], T[i+1])
        V[i+1] = vc + 0.5*k
        vc = vc + k
    # Return values
    return T, R, V

# Unpacks coordinate array into x and y coordinates
def unpack(R):
    N = len(R[0])//2
    x = R[:,:N]
    y = R[:,N:]
    return x, y

# Defines separation between two particle based on their coordinates
# Trivial in this case
def separation_nonperiodic(x1, x2):
    return x1 - x2

 
def LJ_many(sep, r, t, sigma, eps, m): 
    """
    Calculates acceleration of N particles with Lennard-Jones potential,
    given their coordinate arrays

    Parameters
    ----------
    sep : function
        Calculates separation between two particles based on their coordinates.
    r : float array
        General coordinate array.
    t : float
        Time.
    sigma : float
        Length scale of LJ potential.
    eps : float
        Energy scale of LJ potential.
    m : float
        Particle mass.

    Returns
    -------
    a : float array 
        General acceleration array.

    """
    # Unpacking coordinates into x and y values
    N = len(r)//2
    x = r[:N]
    y = r[N:]
    
    # Intializing acceleration arrays
    ax = np.zeros(N)
    ay = np.zeros(N)
    # Iterating over each pair (i,j)
    for i in range(N):
        for j in range(i+1,N):
            # Calculating distance between i, j
            rx = sep(x[i], x[j])
            ry = sep(y[i], y[j])
            r2 = rx**2 + ry**2
            
            # Updating acceleration of each particle
            factor = sigma**2/r2
            coeff = (4*eps)/(m*r2)*(12*factor**6 - 6*factor**3)
            
            ax[i] = ax[i] + coeff*rx
            ax[j] = ax[j] - coeff*rx
            ay[i] = ay[i] + coeff*ry
            ay[j] = ay[j] - coeff*ry
    
    # Return accleration array
    a = np.reshape([ax, ay], [2*N])
    return a

# Sets initial conditions for N particles arranged in a rectangle with 
# dimensions Lx, Ly. Code taken from assignment document.
def particle_setup(N, Lx, Ly):
    dx = Lx/np.sqrt(N)
    dy = Ly/np.sqrt(N)
    
    x_grid = np.arange(dx/2, Lx, dx)
    y_grid = np.arange(dy/2, Ly, dy)
    xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
    x_initial = xx_grid.flatten()
    y_initial = yy_grid.flatten()
    
    r_initial = np.reshape([x_initial, y_initial], [2*N])
    v_initial = np.zeros(2*N)
    return r_initial, v_initial

# Kinetic energy of single particle
def KE_particle(v, m): 
    return 0.5*m*v**2

# Potential energy of particle pair, given the squared distance between them
def PE_particle_pair(r2, sigma, eps):
    factor = sigma**2/r2
    return 4*eps*(factor**6 - factor**3)

# Kinetic and potential energy, using m = eps = sigma = 1.
KE1 = lambda v: KE_particle(v, 1)
PE1 = lambda r2: PE_particle_pair(r2, 1, 1)

# Kinetic energy of N particles, given array of velocities length 2N
def kinetic(KE_single, V):
    KE_array = KE_single(V)
    return np.sum(KE_array, axis=1)

# Potential energy of N particles, given array of coordinates length 2N,
# as well as separation function.
def potential(PE_pair, sep, R):
    N = len(R[0])//2
    x, y = unpack(R)
    
    PE_array = np.zeros(len(R))
    for i in range(N):
        for j in range(i+1,N):
            rx = sep(x[:,i], x[:,j])
            ry = sep(y[:,i], y[:,j])
            r2 = rx**2 + ry**2
            
            PE_array = PE_array + PE_pair(r2)
    return PE_array


# Plots trajectories of particles, as well as the kinetic, potential and total
# energy as a function of time. 
def plot_LJ(mapping, KE, PE, T, R, V, particles, suffix):
    # Unpack the x,y positions and x,y velocities
    x, y = unpack(R)
    Vx, Vy = unpack(V)
    
    # Map positions into where they will be plotted
    xmap = mapping(x)
    ymap = mapping(y)
    
    # Plot particle trajectories
    plt.figure()
    for i in particles:
        plt.plot(xmap[:,i], ymap[:,i], ',', markersize=1, label='Particle '+str(i+1))
    plt.grid(True)
    plt.xlabel(r"x Position ($\sigma$)")
    plt.ylabel(r"y Position ($\sigma$)")
    plt.title("Particle Trajectories")
    plt.savefig("Trajectories_LJ_"+suffix)
    
    # Calculate the kinetic and potential energies as function of time
    total_kinetic = KE(V)
    total_potential = PE(R)
    
    # Plot kinetic, potential and total energy
    plt.figure()
    plt.plot(T, total_kinetic, label='Kinetic')
    plt.plot(T, total_potential, label='Potential')
    plt.plot(T, total_kinetic + total_potential, label='Total Energy')
    plt.grid(True)
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.title("Kinetic, Potential and Total Energy")
    plt.savefig("Energy_LJ_"+suffix)
    
    # Plot deviation of total energy from initial
    plt.figure()
    plt.plot(T, (total_kinetic + total_potential)/(total_kinetic + total_potential)[0] - 1, label='Total Energy')
    plt.grid(True)
    plt.xlabel("Time")
    plt.ylabel("Fluctuation of Energy")
    plt.title("Fluctuation of Total Energy")
    plt.savefig("Energy_Fluctuation_LJ_"+suffix)

"""
Question 3, Part (b)

Here, we perform simulation of 16 particles interacting with Lennard Jones 
potential, without periodic boundary conditions
"""
# 16 particles, in a 4.0 by 4.0 square
N = 16
Lx = 4.0
Ly = 4.0

# Set initial positions and velocites of particles
r0, v0 = particle_setup(N, Lx, Ly)

# Define acceleration function by using a normal (non-periodic) separation function
LJ_acceleration_nonperiodic = lambda r, t: LJ_many(separation_nonperiodic, r, t, 1, 1, 1)
# Identity mapping for purposes of plotting
identity = lambda x: x

# These are the particles we want to plot
plot_particles = np.arange(N)

# There are the functions for total kinetic and potential energy, using 
# non-periodic separation function
totalKE = lambda V: kinetic(KE1, V)
totalPE = lambda R: potential(PE1, separation_nonperiodic, R)

# Perform numerical simulation using Verlet's method, using our 
# acceleration function, and store the array of times, positions and velocities
T, R, V = verlet(LJ_acceleration_nonperiodic, 0, r0, v0, 0.01, 1000)
# PLot our results
plot_LJ(identity, totalKE, totalPE, T, R, V, plot_particles, "3i")

"""
Question 3, Part (c)

Here, we perform simulation of 16 particles interacting with Lennard Jones 
potential, with periodic boundary conditions
"""
# 16 particles, in a 4.0 by 4.0 square
N = 16
Lx = 4.0
Ly = 4.0

# Set initial positions and velocites of particles
r0, v0 = particle_setup(N, Lx, Ly)

# Define separation function that takes into account periodic boundary conditions
def separation_periodic(x1, x2, L):
    return np.mod(x1-x2+L/2, L) - L/2

sep = lambda x1, x2: separation_periodic(x1, x2, Lx)

# Define acceleration function by using a periodic separation function
LJ_acceleration_periodic = lambda r, t: LJ_many(sep, r, t, 1, 1, 1)

# Periodic mapping for purposes of plotting
periodic = lambda x: np.mod(x,Lx)

# These are the particles we want to plot
plot_particles = np.arange(N)

# There are the functions for total kinetic and potential energy, using 
# periodic separation function
totalKE = lambda V: kinetic(KE1, V)
totalPE = lambda R: potential(PE1, sep, R)

# Perform numerical simulation using Verlet's method, using our 
# acceleration function, and store the array of times, positions and velocities
T, R, V = verlet(LJ_acceleration_periodic, 0, r0, v0, 0.01, 1000)
# PLot our results
plot_LJ(periodic, totalKE, totalPE, T, R, V, plot_particles, "3ii")


