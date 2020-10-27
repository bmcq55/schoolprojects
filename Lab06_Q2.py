# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 21:13:32 2020

@author: Nicholas Sullivan

Lab Partner: Brendon McHugh

Question X of LabXX

[Description of File]
"""

"""
Question X, Part (a)

[Description]
"""

import numpy as np
import matplotlib.pyplot as plt


def verlet(f, t0, r0, v0, h, N):
    T = np.linspace(t0, t0+h*N, N+1)
    vc = v0 + 0.5*h*f(r0, T[0])
    dof = len(r0)
    
    R = np.zeros([N+1, dof])
    R[0] = r0
    V = np.zeros([N+1, dof])
    V[0] = v0
    
    for i in range(N):
        R[i+1] = R[i] + h*vc
        k = h*f(R[i+1], T[i+1])
        V[i+1] = vc + 0.5*k
        vc = vc + k
    
    return T, R, V

def plot_LJ(T, R, V, suffix):
    N = len(R[0])//2
    x = R[:,:N]
    y = R[:,N:]
    
    # Vx = V[:,:N]
    # Vy = V[:,N:]
    
    # plt.figure()
    # for i in range(N):
    #     plt.plot(T, x[:,i], label='Particle '+str(i+1)+' x position')
    #     plt.plot(T, y[:,i], label='Particle '+str(i+1)+' y position')
    # plt.grid(True)
    # plt.legend()
    # plt.xlabel("Time")
    # plt.ylabel("Position")
    # plt.title("Position vs Time in LJ Potential")
    # plt.savefig("Pos_vs_Time_LJ_"+suffix)
    
    # plt.figure()
    # for i in range(N):
    #     plt.plot(T, Vx[:,i], label='Particle '+str(i+1)+' x velocity')
    #     plt.plot(T, Vy[:,i], label='Particle '+str(i+1)+' y velocity')
    # plt.grid(True)
    # plt.legend()
    # plt.xlabel("Time")
    # plt.ylabel("Velocity")
    # plt.title("Velocity vs Time in LJ Potential")
    # plt.savefig("Vel_vs_Time_LJ_"+suffix)
    
    plt.figure()
    for i in range(N):
        plt.plot(x[:,i], y[:,i], 'o', markersize=1, label='Particle '+str(i+1))
    plt.grid(True)
    plt.legend()
    plt.xlabel(r"x Position ($\sigma$)")
    plt.ylabel(r"y Position ($\sigma$)")
    plt.title("Particle Trajectories in LJ Potential")
    plt.savefig("Trajectories_LJ_"+suffix)

def LJ_two(r, t, sigma, eps, m): 
    rx = r[0] - r[1]
    ry = r[2] - r[3]
    r2 = rx**2 + ry**2
    factor = sigma**2/r2
    return (4*eps)/(m*r2)*(12*factor**6 - 6*factor**3)*np.array([rx, -rx, ry, -ry])

LJ_acceleration = lambda r, t: LJ_two(r, t, 1, 1, 1)


x1 = 4
y1 = 4
x2 = 5.2
y2 = 4 
T, X1, V1 = verlet(LJ_acceleration, 0, [x1, x2, y1, y2], [0, 0, 0, 0], 0.01, 100)
plot_LJ(T, X1, V1, "2bi")

x1 = 4.5
y1 = 4
x2 = 5.2
y2 = 4 
T, X2, V2 = verlet(LJ_acceleration, 0, [x1, x2, y1, y2], [0, 0, 0, 0], 0.01, 100)
plot_LJ(T, X2, V2, "2bii")

x1 = 2
y1 = 3
x2 = 3.5
y2 = 4.4
T, X3, V3 = verlet(LJ_acceleration, 0, [x1, x2, y1, y2], [0, 0, 0, 0], 0.01, 100)
plot_LJ(T, X3, V3, "2biii")



