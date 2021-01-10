# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 09:52:14 2020

@author: Bmcq
"""
import numpy as np
import matplotlib.pyplot as plt
from time import time

#Lab08
#Q1a:(Capacitor) Adapted from Newman 9.3


'''Gauss-Seidel method for solving PDES pseudocode:
    Step 1: assign constants, grid spacing M, and target accuracy. Create grid arrays.
    Step 2:make a while loop that iterates over the grid of phi values(or whichever function is 
    satisfying the PDE), and include the boundary conditions in the loop using an if statement and an else
    statement to update the function with delta and continue onto the next iteration if delta is not
    greater than the target accuracy.
    Step 3: End iterating through the while loop when the target accuracy is reached. Plot the results.'''
    
#constants
M=100   #grid squares on a side in mm
V=1   #voltage on plates
target=1e-6 #target accuracy in V

#create arrays to hold potential values
phi=np.zeros([M+1,M+1],float)
dphi=np.zeros([M+1,M+1],float)

#main loop
delta=1.0
start=time()
while delta> target:
    
    #delta=0.0
    for i in range(1,M):
        for j in range(1,M):
            
            if j==20 and i>20 and i<80: #capacitor plate +V,position in mm
                phi[i,j]=V
            elif j==80 and i>20 and i<80: #capacictor plate -V, position in mm
                phi[i,j]= -V 
            elif i==0 or i==M or j==0 or j==M: #boundary of system
                phi[i,j]=phi[i,j]
            else:
                dphi[i,j] = (phi[i+1,j]+phi[i-1,j]+phi[i,j+1]+phi[i,j-1])/4 - phi[i,j]
                phi[i,j] = phi[i,j] + dphi[i,j]
    
    delta=np.max(np.abs(dphi))
end=time()
diff1=end-start
print('the time for the GS method is',diff1)

x=np.linspace(0,10,M+1)
y=np.linspace(0,10,M+1)
X,Y=np.meshgrid(x,y)

#contour plot
plt.figure()           
plt.contourf(X, Y, phi,cmap='plasma')
plt.title('Capacitors Gradient of Potential(Gauss-Seidel Method)')
plt.xlabel('position-x (cm)')
plt.ylabel('position-y(cm)')
plt.show()
plt.savefig('capacitor contour 1.png')

#streamline plot
Ey, Ex = np.gradient(-phi, y, x) # careful about order
fig = plt.figure()
strm = plt.streamplot(X, Y, Ex, Ey, color=phi, linewidth=2, cmap='plasma')
cbar = fig.colorbar(strm.lines)
cbar.set_label('Potential $V$')
plt.title('Electric field lines(Gauss-Seidel Method)')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.axis('equal')
plt.tight_layout()
plt.show()
plt.savefig('capacitor streamline 1.png')

#Q2b: (Capacitor) Overelaxation ω=0.1

#main loop
delta=1.0 
ω=0.1 #overrelaxation parameter
start=time()
while delta>target: #limit of accuracy
    
    delta=0.0
    for i in range(1,M):
        for j in range(1,M):
            
            if j==20 and i>20 and i<80: #capacitor plate +V,position in mm
                phi[i,j]=V
            elif j==80 and i>20 and i<80: #capacictor plate -V, position in mm
                phi[i,j]= -V 
            elif i==0 or i==M or j==0 or j==M: #boundary of system
                phi[i,j]=0
            else:
                Δphi= (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])/4 - phi[i,j] #change in phi on iteration
                phi[i,j] = phi[i,j] + (1+ω)*Δphi
                    
            if Δphi>delta: 
                delta=np.max(np.abs(Δphi))
                
end=time()
diff2=end-start
print('the time for the GS method with overrelaxation (omega=0.1) is',diff2)
#make x,y points for contour/streamline plot
x=np.linspace(0,10,M+1) 
y=np.linspace(0,10,M+1)
X,Y=np.meshgrid(x,y)            
  
#contour plot
plt.figure()           
plt.contourf(X,Y,phi, cmap='autumn')
plt.title('Capacitors Gradient of Potential (Gauss-Seidel w/ overrelaxation)')
plt.xlabel('position-x (cm)')
plt.ylabel('position-y(cm)')
plt.show()
plt.savefig('capacitor contour 2.png')

#streamline plot
Ey, Ex = np.gradient(-phi, y, x) 
fig = plt.figure()
strm = plt.streamplot(X, Y, Ex, Ey, color=phi, linewidth=2, cmap='autumn')
cbar = fig.colorbar(strm.lines)
cbar.set_label('Potential $V$')
plt.title('Electric field lines (Gauss-Seidel w/ overrelaxation)')
plt.xlabel('$x$ (cm)')
plt.ylabel('$y$ (cm)')
plt.axis('equal')
plt.tight_layout()
plt.show()
plt.savefig('capacitor streamline 2.png')
           

#Q2b: (Capacitor) Overelaxation omega=0.5    

'''Gauss-Seidel method for solving PDES pseudocode:
    Step 1: assign constants, grid spacing M, and target accuracy. Create grid arrays.
    Step 2:make a while loop that iterates over the grid of phi values(or whichever function is 
    satisfying the PDE), and include the boundary conditions in the loop using an if statement and an else
    statement to update the function with Δphi, which "overshoots" the value of phi by an amount controlled by 
    the parameter ω. Continue onto the next iteration if delta is not greater than the target accuracy.
    Step 3: End iterating through the while loop when the target accuracy is reached. Plot the results.'''
#main loop    
delta=1.0 
ω=0.5 #overrelaxation parameter
start=time()
while delta>target: #limit of accuracy
    
    delta=0.0
    for i in range(1,M):
        for j in range(1,M):
            
            if j==20 and i>20 and i<80: #capacitor plate +V,position in mm
                phi[i,j]=V
            elif j==80 and i>20 and i<80: #capacictor plate -V, position in mm
                phi[i,j]= -V 
            elif i==0 or i==M or j==0 or j==M: #boundary of system
                phi[i,j]=0
            else:
                Δphi= (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])/4 - phi[i,j] #change in phi on iteration
                phi[i,j] = phi[i,j] + (1+ω)*Δphi
                    
            if Δphi>delta: 
                delta=np.max(np.abs(Δphi))
end=time()
diff3=end-start
print('the time for the GS method with overrelaxation (omega=0.5) is',diff3)              

#make x,y points for contour/streamline plot
x=np.linspace(0,10,M+1) 
y=np.linspace(0,10,M+1)
X,Y=np.meshgrid(x,y)            
  
#contour plot
plt.figure()           
plt.contourf(X,Y,phi, cmap='spring')
plt.title('Capacitors Gradient of Potential (Gauss-Seidel w/ overrelaxation)')
plt.xlabel('position-x (cm)')
plt.ylabel('position-y(cm)')
plt.show()
plt.savefig('capacitor contour 3.png')

#streamline plot
Ey, Ex = np.gradient(-phi, y, x) 
fig = plt.figure()
strm = plt.streamplot(X, Y, Ex, Ey, color=phi, linewidth=2, cmap='spring')
cbar = fig.colorbar(strm.lines)
cbar.set_label('Potential $V$')
plt.title('Electric field lines (Gauss-Seidel w/ overrelaxation )')
plt.xlabel('$x$ (cm)')
plt.ylabel('$y$(cm)')
plt.axis('equal')
plt.tight_layout()
plt.show()
plt.savefig('capacitor streamline 3.png')



