# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 11:25:28 2020

@author: Brendon McHugh

Question 3 of Lab03

Here, we plot height as a function of position for the
Lake Geneva region, calculate the gradient, and make 
a relief plot of the region.
"""

"""
Question 3, Part (a)

Pseudocode for loading the height file, calculating gradient,
and producing relief plot.
"""
# Import relevent libraries (numpy, matplotlib, struct)
# Load file containing height values for Lake Geneva region

# for i,j in range(1201, 1201)
    # Transfer height value into i,j index of w (height array)

# Create gx and gy, gradient arrays

# Calculate gradients for most values using central
# difference method

# Calculate x gradient on west edge using forward difference method
# Calculate x gradient on east edge using backward difference method

# Calculate y gradient on north edge using forward difference method
# Calculate y gradient on south edge using backward difference method

# Calculate 'intensity' using equation on page 212

# Plot height array, using cutoff to hide bad values
# Plot relief map using 'intensity' array

"""
Question 3, Part (b)

Here, we use the pseudocode to produce height plot and relief
plots.
"""

# Import struct, matplotlib, numpy
import struct
import matplotlib.pyplot as plt
import numpy as np

# Load file containing Lake Geneva data
filename = 'N46E006.hgt'
f = open(filename, 'rb')

# Load file into height array
L = 1201
w = np.zeros([L,L])

for i in range(L):
    for j in range(L):
        # Read height value, write into array
        buf = f.read(2) 
        w[i,j] = struct.unpack('>h', buf)[0] 

# Latitude and Longitude ranges
long_min = 6
long_max = 7
lat_min = 47
lat_max = 46

# Calculate gradient from altitude array

# Distance between grid points in meters
h = 83

# Create arrays for x and y gradients
gx = np.zeros([L,L])
gy = np.zeros([L,L])

# Calculate gradient for center values
gx[1:(L-1), :] = (w[2:L,:] - w[0:(L-2),:])/(2*h)
gy[:,1:(L-1)] = (w[:,2:L] - w[:,0:(L-2)])/(2*h)

# Calculate x gradient on east and west edges
gx[0,:] = (w[1,:] - w[0,:])/(h)
gx[L-1,:] = (w[L-1,:] - w[L-2,:])/(h)

# Calculate y gradient on north and south edges
gy[:,0] = (w[:,1] - w[:,0])/(h)
gy[:,L-1] = (w[:,L-1] - w[:,L-2])/(h)

# The sun is shining from the southwest
phi = np.pi/6

# Calculate 'intensity' of light from the sun for relief plot
I = np.zeros([L,L])
I = -(np.cos(phi)*gx + np.sin(phi)*gy)/(np.sqrt(gx**2 + gy**2 + 1))

# Function to plot height as function of position
def plot_height(X, Y, title, figtitle='Height_Position'):
    xmin = X[0]
    xmax = X[-1]
    ymin = Y[0]
    ymax = Y[-1]
    
    # Range of longitudes and latitudes to plot
    long_range = long_min + (long_max-long_min)*X/(L-1)
    lat_range = lat_min + (lat_max-lat_min)*Y/(L-1)
    
    # Create figure of height as a function of position
    plt.figure()
    plt.pcolormesh(long_range, lat_range, w[xmin:(xmax+1),ymin:(ymax+1)], vmin = 0, cmap = 'jet', shading='auto')
    plt.xlabel(r'Longitude ($\degree$ E)')
    plt.ylabel(r'Latitude ($\degree$ N)')
    plt.gca().set_aspect('equal')
    plt.title(title)
    plt.colorbar()
    plt.savefig(figtitle)

# Function to make relief plot as function of position
def plot_relief(X, Y, title, figtitle='Relief'):
    xmin = X[0]
    xmax = X[-1]
    ymin = Y[0]
    ymax = Y[-1]
    
    # Range of longitudes and latitudes to plot
    long_range = long_min + (long_max-long_min)*X/(L-1)
    lat_range = lat_min + (lat_max-lat_min)*Y/(L-1)
    
    # Create relief plot based on 'Intensity' values
    plt.figure()
    plt.pcolormesh(long_range, lat_range, I[xmin:(xmax+1),ymin:(ymax+1)], cmap = 'gray', shading='auto')
    plt.xlabel(r'Longitude ($\degree$ E)')
    plt.ylabel(r'Latitude ($\degree$ N)')
    plt.gca().set_aspect('equal')
    plt.title(title)
    plt.savefig(figtitle)

# Region to zoom in on
X1 = np.arange(0,L)
Y1 = np.arange(0,L)

# Plotting height, relief map
plot_height(X1, Y1, 'Plot of Height (m) vs Position', figtitle='Height_Position_3b')
plot_relief(X1, Y1, 'Relief plot of Lake Geneva', figtitle='ReliefPlot_3b')

# Region to zoom in on
X2 = np.arange(300,400)
Y2 = np.arange(1000,1100)

# Plotting relief map for this region
plot_relief(X2, Y2, 'Relief plot of Lake Geneva, detail', figtitle='ReliefPlot2_3b')

