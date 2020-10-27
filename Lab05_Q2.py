# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 11:05:28 2020

@author: Brendon McHugh

Question 2 of Lab05

[Description of File]
"""

"""
Question 2, Part (a)

[Description]
"""

import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
import struct

# Load file containing blurred image
filename = 'blur.txt'
f = open(filename, "r")
lines = f.readlines()
f.close()

# Turn string array into numpy float array
data = np.array([np.fromstring(line, sep=' ') for line in lines])

# Show figure
plt.figure()
plt.imshow(data, cmap='gray')


"""
Question 2, Part (b)

[Description]
"""
#  
L = len(data)

# x and y coordinates, with negative values in the second half
x = [np.arange(0,L/2), np.arange(-L/2,0)]
y = [np.arange(0,L/2), np.arange(-L/2,0)]

# 2-D arrays of x and y coordinates, and r^2
X,Y = np.meshgrid(x, y)
R2 = X**2 + Y**2

# Gaussian point spread function
sigma = 25
gaussian = np.exp(-R2/(2*sigma**2))

# Show point spread function
plt.figure()
plt.imshow(gaussian, cmap='gray')


"""
Question 2, Part (c)

[Description]
"""
# Take fourier transform of data, gaussian
data_ft = fft.rfft(data)
gaussian_ft = fft.rfft(gaussian)

# If value is less that 1e-3, set equal to 1/L^2 (so that it won't have any affect)
gaussian_ft[gaussian_ft<1e-3] = 1/L**2

# Divide fourier transforms to yield unblurred fourier transformed data
unblur_ft = data_ft/(L**2*gaussian_ft)

# Take inverse fourier transform to yield the unblurred image
unblur = fft.irfft(unblur_ft)

# Show unblurred image
plt.figure()
plt.imshow(unblur, cmap='gray')


"""
Question X, Part (d)

[Description]
"""


"""
Question X, Part (e)

[Description]
"""