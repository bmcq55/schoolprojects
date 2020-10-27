# -*- coding: utf-8 -*-
"""
Created on Wed Oct 7 10:43:46 2020

@author: Brendon McHugh

Question 2 of Lab04

Here, we calculate the bound states of a quantum particle in a linear potential 
inside an infinite potential well. 
"""

"""
Question 2, Part (c)

Here, we calculate the eigenvalues of the 10x10 approximation of the 
Hamiltonian matrix, and define the functions needed.
"""

# Import libraries
import scipy.constants as const
import numpy as np
import matplotlib.pyplot as plt


def Linear_Potential_Hamiltonian(m, n, L, M, a):
    """
    Hamiltonian matrix element H_mn in units of electron volts for a linear 
    potential inside infinite potential well.

    Parameters
    ----------
    m : int
        First matrix index of Hamiltonian.
    n : int
        Second matrix index of Hamiltonian.
    L : float
        Well width in units of meters.
    M : float
        Mass in units of units of kilograms*(eV/Joule).
    a : float
        Energy constant in units of eV.

    Returns
    -------
    entry : float
        Hamiltonian element H_mn in units of eV.

    """
    
    # Planck's constant in enits of eV*s
    hbar = const.hbar/const.eV
    entry = 0
    # Calculates entry depending on m, n
    if (m==n): 
        # If m = n
        entry = (hbar**2*np.pi**2*n**2)/(2*M*L**2) + a/2
    elif ((m+n)%2 == 1):
        # If one of m, n is odd and one is even
        entry = -8*a*m*n/(np.pi**2*(m**2-n**2)**2)
    else:
        # Otherwise, matrix element is 0
        entry = 0
    return entry


def Hamiltonian_Matrix(H_Element, Nmax):
    """
    Calulculates the Hamiltonian matrix of size (Nmax x Nmax) using
    function H_element.

    Parameters
    ----------
    H_Element : function (int, int) : (float)
        Function calculating Hamiltonian element at (m,n).
    Nmax : int
        Size of matrix.

    Returns
    -------
    H_mat : TYPE
        Hamiltonian matrix of size (Nmax x Nmax).

    """
    H_mat = np.zeros([Nmax,Nmax])
    for m in range(1, Nmax+1):
        for n in range(1, Nmax+1):
            H_mat[m-1,n-1] = H_Element(m,n)
    return H_mat


# Length in units of meters, time in units of seconds, energy in units of eV 
# Mass is in units of eV*s^2/m^2, or eV*(kg/J)
L = 5*const.angstrom
a = 10.0
M = const.electron_mass/const.eV

# Our specific potential
H1 = lambda m, n: Linear_Potential_Hamiltonian(m, n, L, M, a)

# Calculate 10x10 Hamiltonian matrix
Nmax = 10
H_mat = Hamiltonian_Matrix(H1, Nmax)

# Calculate eigenvalues and eigenvectors
D,V = np.linalg.eigh(H_mat)

# Print results
print("10x10 Approximation")
print("Energy Level\tEnergy (eV)")
for i in range(10):
    print("%d\t%.8f"%(i, D[i]))
print()


"""
Question 2, Part (d)

Here, we calculate the eigenvalues of the 100x100 approximation of the 
Hamiltonian matrix.
"""

# Calculate 100x100 Hamiltonian matrix
Nmax = 100
H_mat = Hamiltonian_Matrix(H1, Nmax)

# Calculate eigenvalues and eigenvectors
D,V = np.linalg.eigh(H_mat)

# Print results
print("100x100 Approximation")
print("Energy Level\tEnergy (eV)")
for i in range(10):
    print("%d\t%.8f"%(i, D[i]))
print()


"""
Question 2, Part (e)

Here, we plot the wavefunctions associated with the first three energy levels
in the 100x100 approximation to the Hamiltonian matrix. Wavefunction 
normalized using Simpson method.
"""


def Wavefunction(V, L, N_points):
    """
    Calculating wavefunction at N_points points in the interval [0,L], using 
    sinusoidal coefficients contained in array V. Returns array of points
    that were sampled, and array of wavefunction amplitudes.

    Parameters
    ----------
    V : array of type float
        Coefficients of sinusoidal terms.
    L : float
        Length of interval to compute wavefunction.
    N_points : int
        Number of points for wavefunction to be calculated.

    Returns
    -------
    x : numpy array of type float
        Array of points that were sampled.
    psi : numpy array of type float
        Array of wavefunction amplitudes at these sample points.

    """
    # Convert length into units of Angstroms
    L = L/const.angstrom 
    x = np.linspace(0,L,N_points)
    psi = np.zeros(N_points)
    for m in range(1, len(V)+1):
        psi += V[m-1]*np.sin(np.pi*m*x/L)
    
    # Normalize using trapezoidal integration
    w = np.ones(N_points)
    w[0] = 0.5
    w[-1] = 0.5
    w = w*L/(N_points-1)
    A = np.sum(w*np.abs(psi)**2)
    psi = psi/np.sqrt(A)
    return x, psi


# Plotting Wavefunction Amplitude
plt.figure()
x, psi = Wavefunction(V[:,0], L, 1000)
plt.plot(x, psi, 'r-', label='Energy Level 0')
x, psi = Wavefunction(V[:,1], L, 1000)
plt.plot(x, psi, 'g-', label='Energy Level 1')
x, psi = Wavefunction(V[:,2], L, 1000)
plt.plot(x, psi, 'b-', label='Energy Level 2')

plt.grid(True)
plt.legend()
plt.xlabel(r'Position ($\AA$)')
plt.ylabel('Amplitude')
plt.title('Wavefunction Ampitude vs Position')
plt.savefig("Wavefunction_vs_Position")

# Plotting Probability Density
plt.figure()

x, psi = Wavefunction(V[:,0], L, 1000)
plt.plot(x, np.abs(psi)**2, 'r-', label='Energy Level 0')
x, psi = Wavefunction(V[:,1], L, 1000)
plt.plot(x, np.abs(psi)**2, 'g-', label='Energy Level 1')
x, psi = Wavefunction(V[:,2], L, 1000)
plt.plot(x, np.abs(psi)**2, 'b-', label='Energy Level 2')

plt.grid(True)
plt.legend()
plt.xlabel(r'Position ($\AA$)')
plt.ylabel('Probability Density')
plt.title('Probability Density vs Position')
plt.savefig("Probability_vs_Position")







