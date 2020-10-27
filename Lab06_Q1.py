  # -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

#Newman Excercise 8.8b: Space garbage

#equations of motion for ball bearing in xy-plane

#psdeudocode: (step 1) assign values to constants/initial conditions. 
            #(step 2) convert two 2nd order equations of motion into 4 equations
            # (step 3) write function that uses 4th order RungeKutta (R4) method to solve system
            #(step 4) plot the orbit of the ball bearing
#step 1:
G=1 #constants
M=10
L=2
x_i=1.0 #initial conditions
y_i=0.0
v_xi=0.0
v_yi=1.0



#(step 2):convert two 2nd order equations of motion into 4 equations
#dx/dt=v_x, dv_x/dt=f(x,v,t)=d²x/dt² and similar for dy/dt=v_y
#(step 3): function
def f(r,t): #f(r,t)=f(x,v_x,y,v_y,t) split vector in components

    ''' INPUT:vector componens r, and time t
    OUTPUT: returns the array of vector component values'''
     
    x=r[0] #x/y=position
    v_x=r[1] #v_x/v_y=velocity
    y=r[2]
    v_y=r[3]
    
    dis=np.sqrt(x**2+y**2) #r=sqrt(x^2+y^2)
    a_x=(-G*M*x)/(dis**2*np.sqrt(dis**2+L**2/4)) #a_x/a_y=acceleration
    a_y=(-G*M*y)/(dis**2*np.sqrt(dis**2+L**2/4))
    return np.array([v_x,a_x,v_y,a_y],float) #return array

a=0.0 #starting value t=0
b=10.0 #ending value t=10
N=1000 #number of steps
h=(b-a)/N #size of single step

tpoints=np.arange(a,b,h) #list of values of t at each step of the RK4 method
xpoints=[] #empty arrays to append values of x/y position's
ypoints=[]
r=np.array([x_i,v_xi,y_i,v_yi],float) #r vector initial conditions

#Runge-Kutta fourth order method
for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[2])
    k1=h*f(r,t)
    k2=h*f(r+0.5*k1,t+0.5*h)
    k3=h*f(r+0.5*k2,t+0.5*h)
    k4=h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6


# (step 4) plot x vs y position
plt.figure()
plt.title("Space Garbage Orbit")
plt.xlabel("position x (m)")
plt.ylabel("position y (m)")
plt.plot(xpoints,ypoints)
plt.show()    
plt.savefig("space garbage plot.png")
    
    
    
    
    
    
    
    
    