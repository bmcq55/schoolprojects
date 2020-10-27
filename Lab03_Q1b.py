# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:25:28 2020

@author: Brendon McHugh

Question 1(b) of Lab03

Here, we use gaussian quadrature to calculate the 
probability of blowing snow based on the average 
hourly windspeed at 10m, the snow surface age
and the average hourly temperature.
"""

# Import numpy, matplotlib, gaussxw
import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxwab


def prob_blowing_snow(u10, Ta, th, N):
    """
    Function that uses gaussian quadrature to calculate 
    the probability of blowing snow

    Parameters
    ----------
    u10 : float
        Average hourly windspeed at 10m.
    Ta : float
        Average hourly temperature
    th : float
        Snow surface age.
    N : int
        Number of sample points.

    Returns
    -------
    float
        Probability of blowing snow in these conditions.

    """
    # Mean and standard deviation wind speed
    ubar = 11.2 + 0.365*Ta + 0.00706*Ta**2 + 0.9*np.log(th)
    delta = 4.3 + 0.145*Ta + 0.00196*Ta**2
    
    # Integrand as a function u, wind speed
    integrand = lambda u : np.exp(-(ubar-u)**2/(2*delta**2))
    
    # Calculation of integral using gaussian quadrature
    s, w = gaussxwab(N, 0, u10)
    integral = np.sum(w*integrand(s))
    
    # Normalize integral to calculate probability
    return integral/(np.sqrt(2*np.pi)*delta)
    
# Number of sample points
N = 100

# Values of u10 and th to investigate
u10_vals = [6,8,10]
th_vals = [24,48,72]

# Values of Ta to plot
Ta = np.linspace(-30, 10, 100)

# Array to store probability values
P_bs = np.zeros(len(Ta))

# Colours and line types for plotting
colours = ('r', 'g', 'b')
lines = ('.', '-', ':')

# Create plot of proability of blowing snow vs temperature
plt.figure()
for (u10, colour) in zip(u10_vals, colours):
    for (th, line) in zip(th_vals, lines):
        P_bs = prob_blowing_snow(u10, Ta, th, N)
        plot_str = colour + line
        plt.plot(Ta, P_bs, plot_str, label=r'$u_{10}$ = %.4f, $t_h$ = %.4f'%(u10, th))
        
plt.xlabel('Average Hourly Temperature (Degrees C)')
plt.ylabel('Probability of Blowing Snow')
plt.legend()
plt.grid(True)

plt.title('Probability of Blowing Snow vs Temperature')
plt.savefig('ProbSnow_Temp_1b')


# Create log plot of proability of blowing snow vs temperature
plt.figure()
for (u10, colour) in zip(u10_vals, colours):
    for (th, line) in zip(th_vals, lines):
        P_bs = prob_blowing_snow(u10, Ta, th, N)
        plot_str = colour + line
        plt.semilogy(Ta, P_bs, plot_str, label=r'$u_{10}$ = %.4f, $t_h$ = %.4f'%(u10, th))
        
plt.xlabel('Average Hourly Temperature (Degrees C)')
plt.ylabel('Probability of Blowing Snow')
plt.legend()
plt.grid(True)

plt.title('Probability of Blowing Snow vs Temperature')
plt.savefig('ProbSnow_Temp_Semilogy_1b')
