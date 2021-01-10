# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 06:40:32 2020

@author: BMCQ
"""

#Lab 10: 

"""
Part a

Here, we plot the trajectory of particle undergoing simulated brownian 
motion in a domain of size Lp = 101.
"""

#Q1 Brownian Motion (Newman 10.3)

#pseudocode
#(step1): create a function that let's a particle move in random directions 1 unit at a time
#(step 2): assign constants and write a main program that has a particle starting at the centre of 
#           grid, and moves randomly over the course of Nt steps
#(step 3): plot the randm trajectory or walk of the particle

import numpy as np
import random
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation


def nextmove(x, y):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = left, 3 = right"""
    
    direction=random.randrange(4)
  
  
    if direction == 0:  # move up
        y += 1
    elif direction == 1:  # move down
        y -= 1
    elif direction == 2:  # move right
        x += 1
    elif direction == 3:  # move left
        x -= 1
    else:
        print("error: direction isn't 0-3")

    return x, y

#main program

Lp = 101  # size of domain
Nt = 5000  # number of time steps
position_x=np.empty([Nt]) # arrays to record the trajectory of the particle
position_y=np.empty([Nt])

centre_point = (Lp-1)//2  # middle point of domain
xp = centre_point
yp = centre_point

position_x[0] = xp #initial starting point at centre of grid
position_y[0] = yp

for i in range(Nt-1):
    xpp, ypp = nextmove(position_x[i], position_y[i])
    position_x[i+1] = np.mod(xpp, Lp)
    position_y[i+1] = np.mod(ypp, Lp)

#plot of trajectory
plt.figure()
plt.title('Brownian Motion of a Particle')
plt.xlabel('x-axis')
plt.ylabel('y-axis')
plt.plot(position_x,position_y)
plt.show()
plt.savefig('Q1a_brownianmotion.png')
    

#Q1.b Diffusion Limited Aggregation(DLA)

"""
Part b

Here, we plot the diffusion limited aggregation of particles in a 101x101 domain
"""

#pseudocode
#(step 1): modify program in part a so that the particle performs a random walk until it reaches
#           a point on the edge and anchors to it. Have another particle start at the centre and repeat this
#         process until it either sticks to another particle or the edge of the system
#           Store the anchored particles in an array. At each step of the random walk, check the particle
#           position using array that stores which points neighbour an anchor point
#           In anchor_neighbour array, if element = 2, square has an anchored point,
#           if element = 1, square neighbours an edge or anchored point and element = 0 otherwise.

def animate_scatter(T,X,Y,Xname,Yname,figname): #animation
    '''INPUT:total time, X axis, Y axis, X axis label, Y axis label, name of figure
        OUTPUT: scatter-plot animation of DLA process with frames of each particle after anchoring'''
    plt.style.use('seaborn-pastel')
    
    xmin = np.min(X)
    xmax = np.max(X)
    xlim = (xmin,xmax)
    
    ymin = np.min(Y)
    ymax = np.max(Y)
    ylim = (ymin,ymax)
    
    fig = plt.figure()
    ax = plt.axes(xlim=xlim, ylim=ylim)
    plt.xlabel(Xname)
    plt.ylabel(Yname)
    
    Nt = len(T)
    colors = np.arange(Nt)
    scat = ax.scatter(X[:i], Y[:i], c = colors, marker='+', cmap='viridis')
    
    def init():
        return scat,
    def animate(i):
        scat.set_offsets(np.swapaxes([X[:i], Y[:i]],0,1))
        scat.set_array(colors[:i])
        plt.title(r'DLA run for %d particles'%(T[i]))
        return scat,
    
    anim = FuncAnimation(fig, animate, init_func=init,
                                   frames=len(T), interval=20, blit=True)
    
    
    anim.save(r'%s.gif'%(figname), writer='pillow')


def nextmove_DLA(x, y):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = left, 3 = right"""
    
    direction=random.randrange(4)
    
    if direction == 0:  # move up
        y += 1
    elif direction == 1:  # move down
        y -= 1
    elif direction == 2:  # move right
        x += 1
    elif direction == 3:  # move left
        x -= 1
    else:
        print("error: direction isn't 0-3")

    return x, y

def anchor_particle(anchor_neighbour, x,y):
    """ INPUT: array containing data of which squares neighbour anchored particles, 
    x and y position of new anchored particle
    OUTPUT: None, array is modified based on position of new anchored particle"""
    xp = np.mod(x+1,Lp)
    yp = np.mod(y+1,Lp) 
    xm = np.mod(x-1,Lp)
    ym = np.mod(y-1,Lp) 
    # Neighbouring points now neighbour an anchor point
    anchor_neighbour[xp, y] = 1
    anchor_neighbour[xm, y] = 1
    anchor_neighbour[x, yp] = 1
    anchor_neighbour[x, ym] = 1
    # This square now has an anchor point
    anchor_neighbour[x, y] = 2

def DLA_anchor(anchor_neighbour, x0, y0):
    """ INPUT: array containing data of which squares neighbour anchored particles, 
    initial x and y position of particle (placed in center)
    OUTPUT: final x,y anchor point, AND array is modified based on position of new anchored particle"""
    xpp = x0
    ypp = y0
    #run if particle does not attach to anchored particle
    while anchor_neighbour[xpp, ypp] == 0: 
        xpp, ypp = nextmove_DLA(xpp, ypp)
    
    # Update anchor particle array
    anchor_particle(anchor_neighbour, xpp, ypp)
    return xpp, ypp

def DLA_animate(Lp, N, figname):
    """ INPUT: Size of domain Lp, number of particles to simulate N, name of figure
    OUTPUT: None, animation and figure are produced based on results"""
    # anchor neighbour array 
    anchor_neighbour = np.zeros([Lp,Lp])
    anchor_neighbour[0,:] = 1
    anchor_neighbour[:,0] = 1
    anchor_neighbour[Lp-1,:] = 1
    anchor_neighbour[:,Lp-1] = 1
    
    # list to represent x and y positions of anchored points
    anchored_point_x = np.empty([N])
    anchored_point_y = np.empty([N])
    
    centre_point = (Lp-1)//2  # middle point of domain
    
    # for each particle, release it from center; once it has anchored, store in array
    for j in range(N):
        xp = centre_point
        yp = centre_point
        
        xpp, ypp = DLA_anchor(anchor_neighbour, xp, yp)
        anchored_point_x[j] = xpp
        anchored_point_y[j] = ypp
    
    # Produce scatter plot animation using anchor array 
    animate_scatter(np.arange(1,N+1), anchored_point_x, anchored_point_y, '$x$', '$y$', figname+'_anim')
    
    # Produce still figure with color gradient from newer to older particles
    plt.figure()
    plt.scatter(anchored_point_x, anchored_point_y, c = np.arange(N), marker='+', cmap='viridis')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'DLA Run for %d particles'%N)
    plt.savefig(figname+'_fig')

# Simulate 100 particles in a domain of size 101
Lp = 101  # size of domain
N = 100  # number of particles

DLA_animate(Lp, N, 'Q1b_DLA_aggregateL101_N100')

"""
Part c

Here, we simulate DLA run for domain of size Lp = 151 until 
the center is anchored. 
"""

def DLA_animate_stopcenter(Lp, figname):
    """ INPUT: Size of domain Lp, name of figure
    OUTPUT: None, animation and figure are produced based on results when 
    anchored points reach center"""
    # anchor neighbour array 
    anchor_neighbour = np.zeros([Lp,Lp])
    anchor_neighbour[0,:] = 1
    anchor_neighbour[:,0] = 1
    anchor_neighbour[Lp-1,:] = 1
    anchor_neighbour[:,Lp-1] = 1
    
    # list to represent x and y positions of anchored points (length Lp^2 to accomodate 
    # worst case scenario)
    anchored_point_x = np.empty([Lp**2])
    anchored_point_y = np.empty([Lp**2])
    
    centre_point = (Lp-1)//2  # middle point of domain
    xp = centre_point
    yp = centre_point
    
    
    # for each particle, release it from center; once it has anchored, store in array
    N = 0
    while(anchor_neighbour[centre_point, centre_point] != 2):
        xp = centre_point
        yp = centre_point
        
        xpp, ypp = DLA_anchor(anchor_neighbour, xp, yp)
        anchored_point_x[N] = xpp
        anchored_point_y[N] = ypp
        N += 1
    
    # Produce scatter plot animation using anchor array (only up to N)
    animate_scatter(np.arange(1,N+1), anchored_point_x[:N], anchored_point_y[:N], '$x$', '$y$', figname+'_anim')
    
    # Produce still figure with color gradient from newer to older particles (only up to N)
    plt.figure()
    plt.scatter(anchored_point_x[:N], anchored_point_y[:N], c = np.arange(N), marker='+', cmap='viridis')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'DLA Run for %d particles'%N)
    plt.savefig(figname+'_fig')


# Simulate in a domain of size 101 until anchored particles reach center
Lp = 151  # size of domain

DLA_animate_stopcenter(Lp, 'Q1c_DLA_aggregateL151')






