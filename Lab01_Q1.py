# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:39:44 2020

@author: Brendon
"""
import numpy as np #import numpy
import matplotlib.pyplot as plt #import graphing visualization functions


#Lab 1, Question 1.b:

#pseudocode of the Euler-Cromer method:
#Assign a value to the constants, G, M_s, alpha
#create a function that takes in parameters u_x0,u_y0,v_x0,v_y0 (u=position,v=velocity)
#Assign a value to the step size, delta_t, create a timevalue
#create a for loop that updates the position variables, (x,y), when iterating
#in this for loop, have the equations for Newtonian planetary orbits for both x and y
#return u_x,u_y,v_x,v_y
#plot (v_x,t),(v_y,t) and (x,y) 

#Lab 1, Question 1.c:

# Constants in SI units
#G=6.67*10**(-11)
#M_s=2.0*10**30

# Constants in AU, M_s, yr
G=39.5
M_s=1


def eulercromer(u_x0,u_y0,v_x0,v_y0): #u_x0/u_y0 is initial position, v_x0/v_y0 is initial velocity, delta_t is step size
    
    r=np.sqrt(u_x0**2+u_y0**2) #r^2=x^2+y^2 is the distance between mercury and the sun
    t0 = 0
    tf = 1
    Nt = 10001
    delta_t=(tf-t0)/(Nt-1)
    t=np.linspace(t0,tf,Nt)
    N=len(t)-1
    
    ux=np.empty(N+1) #u is position
    uy=np.empty(N+1)
    vx=np.empty(N+1) #v is velocity
    vy=np.empty(N+1)
    
    #initial conditions
    ux[0]=u_x0
    uy[0]=u_y0
    vx[0]=v_x0
    vy[0]=v_y0
    
    #integration
    for n in range(N):
        r=np.sqrt(ux[n]**2+uy[n]**2) #r^2=x^2+y^2 is the distance between mercury and the sun
        vx[n+1]=vx[n]-((G*M_s)/r**3)*delta_t*ux[n]
        vy[n+1]=vy[n]-((G*M_s)/r**3)*delta_t*uy[n]
        ux[n+1]=ux[n]+ vx[n+1]*delta_t
        uy[n+1]=uy[n]+ vy[n+1]*delta_t
        
    #print ("the x position is:",ux, "the y position is:",uy,vx,vy,t)
    
    #plots
    plt.figure()
    plt.plot(t,vx)
    plt.title("Velocity vs. time(x-coordinate)")
    plt.xlabel("time(s)")
    plt.ylabel("velocity(m/s)")
    plt.savefig("velvstime(x).png")
    
    plt.figure()
    plt.plot(t,vy)
    plt.title("Velocity vs. time(y-coordinate)")
    plt.xlabel("time(s)")
    plt.ylabel("velocity(m/s)")
    plt.savefig("velvstime(y).png")
    
    plt.figure()
    plt.plot(ux,uy)
    plt.title("x-position vs. y-position")
    plt.xlabel("x-position")
    plt.ylabel("y-position")
    plt.savefig("xvsy.png")
    
    spec_ang_mom = ux*vy - uy*vx
    
    plt.figure()
    plt.plot(t,spec_ang_mom)
    plt.title("Specific Angular Momentum vs. time")
    plt.xlabel("Time (yr)")
    plt.ylabel("Specific Angular Momentum ($AU^2$/yr)")
    plt.ylim(3.5, 4.0)
    plt.savefig("angmomvstime.png")
    

#using our example of mercury around the sun

eulercromer(0.47,0.0,0.0,8.17)


#Lab 1,Question 1.d:

alpha=0.01

def eulercromer_GR(u_x0,u_y0,v_x0,v_y0): #eulercromer function using the general relvativity form of the gravitational force
    
    r=np.sqrt(u_x0**2+u_y0**2) #r^2=x^2+y^2 is the distance between mercury and the sun
    t0 = 0
    tf = 10
    Nt = 10001
    delta_t=(tf-t0)/(Nt-1)
    t=np.linspace(t0,tf,Nt)
    N=len(t)-1
    
    ux=np.empty(N+1) #u is position
    uy=np.empty(N+1)
    vx=np.empty(N+1) #v is velocity
    vy=np.empty(N+1)
    
    #initial conditions
    ux[0]=u_x0
    uy[0]=u_y0
    vx[0]=v_x0
    vy[0]=v_y0
    
    #integration
    for n in range(N):
        r=np.sqrt(ux[n]**2+uy[n]**2) #r^2=x^2+y^2 is the distance between mercury and the sun
        vx[n+1]=vx[n]-((G*M_s)/r**3)*delta_t*(1-alpha/r**2)*ux[n]
        vy[n+1]=vy[n]-((G*M_s)/r**3)*delta_t*(1-alpha/r**2)*uy[n]
        ux[n+1]=ux[n]+ vx[n+1]*delta_t
        uy[n+1]=uy[n]+ vy[n+1]*delta_t
        
    #print ("the x position is:",ux, "the y position is:",uy,vx,vy,t)       
    
    #plot
    plt.figure()
    plt.plot(ux,uy)
    plt.title("position x vs. position y in GR")
    plt.show()

#showing the precession of mercury using 
eulercromer_GR(0.47,0.0,0.0,8.17)



