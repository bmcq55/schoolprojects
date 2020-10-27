# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 13:43:00 2020

@author: Brendon
"""
import numpy as np
from numpy.linalg import solve

#Q1a. modify GaussElim function to make a function that incorporates partial pivoting
def GaussElim(A_in, v_in):
    """Implement Gaussian Elimination. This should be non-destructive for input
    arrays, so we will copy A and v to
    temporary variables
    IN:
    A_in, the matrix to pivot and triangularize
    v_in, the RHS vector
    OUT:
    x, the vector solution of A_in x = v_in """
    # copy A and v to temporary variables using copy command
    A = np.copy(A_in)
    v = np.copy(v_in)
    N = len(v)

    for m in range(N):
        # Divide by the diagonal element
        div = A[m, m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m+1, N):
            mult = A[i, m]
            A[i, :] -= mult*A[m, :]
            v[i] -= mult*v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = np.empty(N, dtype=v.dtype)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i]*x[i]
    return x
    print("the value that the GaussElim function gives is:", x)
    

def partial_pivot(A_in, v_in):
    """ IN: 
    A_in, the matrix to pivot and triangularize
    v_in, the RHS vector
    OUT:
    x, the vector solution of A_in x = v_in """
    # copy A and v to temporary variables using copy command
    A = np.copy(A_in)
    v = np.copy(v_in)
    N = len(v)

    for m in range(N):
        
        #pivoting: check if A[m,m] is larger than elements below
        #if not, swap rows for both the row vector and the column vector v
        for i in range(m+1,N):
            if A[m,m] < A[i,m]:
                A[[m,i],:] = A[[i,m],:]	
                v[[m,i]] = v[[i,m]]
            
        # Divide by the diagonal element
        div = A[m,m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m+1, N):
            mult = A[i, m]
            A[i, :] -= mult*A[m, :]
            v[i] -= mult*v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = np.empty(N, float)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i]*x[i]
    return x
    print("the value that the partial_pivot function gives is:",x)
    
#check if functions gives answer from text
#GaussElim(np.array([[2,1,4,1],
#           [3,4,-1,-1],
#           [1,-4,1,5],
#           [2,-2,1,3]],float),
#        np.array([-4,3,9,7],float))
#
#partial_pivot(np.array([[2,1,4,1],
#           [3,4,-1,-1],
#           [1,-4,1,5],
#           [2,-2,1,3]],float),
#        np.array([-4,3,9,7],float))
    
#Q2b:     
def LU_dec(A_in,v_in):
    x=np.linalg.solve(A_in, v_in)
    return x
    print("the value that the LU_dec function gives is:",x)
    
    
#test LU_dec to ensure it gives the same answer as the above
#LU_dec(np.array([[2,1,4,1],
#           [3,4,-1,-1],
#           [1,-4,1,5],
#           [2,-2,1,3]],float),
#        np.array([-4,3,9,7],float))   

from time import time
import matplotlib.pyplot as plt


def accuracy_timing(N):
    v_in=np.random.rand(N)
    A_in=np.random.rand(N,N)
    

    #Gaussian Elimintation timing
    start1=time()
    x1=GaussElim(A_in,v_in)
    end1=time()
    time1=end1-start1
    #return time1
    #print("the timing of the GaussElim function is:", time1)
    v_sol1=np.dot(A_in,x1)
    err1=np.mean(abs(v_in-v_sol1))
    #print("the error for GaussElim function is:",err1)
    #return err1
    
    
    #partial pivoting timing
    start2=time()
    x2=partial_pivot(A_in, v_in)
    end2=time()
    time2=end2-start2
    #return time2
    #print("the timing of the partial_pivot function is:", time2)
    v_sol2=np.dot(A_in,x2)
    err2=np.mean(abs(v_in-v_sol2))
    #return err2
    #print("the error for the partial_pivot function is:",err2)
    
    #LU_dec timing
    start3=time()
    x3=LU_dec(A_in,v_in)
    end3=time()
    time3=end3-start3
    #return time3
    #print("the timing of the LU_dec function is:", time3)
    v_sol3=np.dot(A_in,x3)
    err3=np.mean(abs(v_in-v_sol3))
    #return err3
    #print("the error for the partial_pivot function is:",err3)
    
    # Return all times and errors
    return time1, time2, time3, err1, err2, err3

N_array = np.arange(5,305,20)
N_length = len(N_array)

# Arrays of storing timing and error data
gauss_timing = np.zeros(N_length)
gauss_error = np.zeros(N_length)
pivot_timing = np.zeros(N_length)
pivot_error = np.zeros(N_length)
LU_timing = np.zeros(N_length)
LU_error = np.zeros(N_length)

for i in range(N_length):
    N = N_array[i]
    print(N)
    
    time1, time2, time3, err1, err2, err3 = accuracy_timing(N)
    
    gauss_timing[i] = time1
    gauss_error[i] = err1
    pivot_timing[i] = time2
    pivot_error[i] = err2
    LU_timing[i] = time3
    LU_error[i] = err3

# Plot Timing Data
plt.figure()
plt.plot(N_array, gauss_timing, 'r-', label='Gaussian Elimination')
plt.plot(N_array, pivot_timing, 'g-', label='Partial Pivot')
plt.plot(N_array, LU_timing, 'b-', label='LU Decomposition')
plt.xlabel('Matrix Dimension N')
plt.ylabel('Time')
plt.legend()
plt.grid(True)
plt.title('Timing of Various Approaches vs Matrix Dimension')
plt.savefig('Timing_vs_N')

# Plot Error Data
plt.figure()
plt.plot(N_array, gauss_error, 'r-', label='Gaussian Elimination')
plt.plot(N_array, pivot_error, 'g-', label='Partial Pivot')
plt.plot(N_array, LU_error, 'b-', label='LU Decomposition')
plt.xlabel('Matrix Dimension N')
plt.ylabel('Error')
plt.legend()
plt.grid(True)
plt.title('Error of Various Approaches vs Matrix Dimension')
plt.savefig('Error_vs_N')

#Q1c. Newman Exercise 6.2
    
from cmath import polar, phase
import math

##Complex circuit
R1 = R3 = R5 = 1000 # resistance in ohms
R2 = R4 = R6 = 2000 
C1 = 10**-6 # capacitance in farads
C2 = 0.5*10**-6 # in farads
x_p= 3 # in volts
ω = 1000 # in hertz

#array of circuit elements
#each row is found using Kirchoff's loop law, and apply Ohm's/Capacitance law with the constants x1,x2,x3
A_in = np.array([[1/R1 + 1/R4 + 1j*ω*C1,-1j*ω*C1,0 ],
            [-1j*C1, 1/R2 + 1/R5 + 1j*ω*C1 + 1j*ω*C2,-1j*ω*C2],
            [0,-1j*ω*C2,1/R3 + 1/R6 + 1j*ω*C2]], complex)
#the RHS array 
v_in = np.array([x_p/R1, x_p/R2, x_p/R3 ], complex)

#solve for x1,x2,x3 and their respective ω
x = LU_dec(A_in,v_in)
for i in range(len(x)):
    r,theta = polar(x[i]) #finds the amplitude and phase of voltages in radians of x1,x2, and x3
    print(r, 'V', theta, 'rads') #r=voltage(x1,x2,x3), theta is ω  (the phase of the x1,x2,x3)
    print(r, 'V', theta*(180/np.pi), 'degrees') #converts phase in radians to degrees

def V(t, V0, omega):
    return V0*np.exp(1j*omega*t)

#calculates V_i=x_i*e^(iωt) for i=1,2,3
t=np.linspace(0,0.05,1000) #time axis
# def V1(t): 
#     return 1.8594351717288344*np.e**((-0.32028664093499326j)*t)
# def V2(t):
#     return 0.9167737698470294*np.e**((-0.25666984698199286j)*t)
# def V3(t):
#     return 1.9908106484802868*np.e**((-0.18042737305109224j)*t)



#plot volage vs. time for each of the voltages and their respective phases
plt.figure()
plt.title("V1,V2 and V3 vs. time")
plt.xlabel("time(s)")
plt.ylabel("amplitude(V)")
#plt.plot(t,V1(t), label='V_1')
#plt.plot(t,V2(t),label='V_2')
#plt.plot(t,V3(t),label='V_3')
plt.plot(t,V(t,x[0],ω), label='V_1')
plt.plot(t,V(t,x[1],ω),label='V_2')
plt.plot(t,V(t,x[2],ω),label='V_3')
plt.legend()
plt.show()
plt.savefig("V1,V2 and V3 vs. time.png")


#replacing resistor R6 with an inductor L=R6/ω
R1 = R3 = R5 = 1000 # resistance in ohms
R2 = R4 = R6 = 2000 
C1 = 10**-6 # capacitance in farads
C2 = 0.5*10**-6 # in farads
x_p= 3 # in volts
ω = 1000 # in hertz
L=R6/ω

#array of circuit elements
#each row is found using Kirchoff's loop law, and apply Ohm's/Capacitance law with the constants x1,x2,x3
A_in = np.array([[1/R1 + 1/R4 + 1j*ω*C1,-1j*ω*C1,0 ],
            [-1j*C1, 1/R2 + 1/R5 + 1j*ω*C1 + 1j*ω*C2,-1j*ω*C2],
            [0,-1j*ω*C2,1/R3 + 1/L + 1j*ω*C2]], complex)
##the RHS array 
v_in = np.array([x_p/R1, x_p/R2, x_p/R3 ], complex)

#solve for x1,x2,x3 and their respective ω
x = LU_dec(A_in,v_in)
for i in range(len(x)):
    r,theta = polar(x[i])
    print(r, 'V', theta, 'rads')
    print(r, 'V', theta*(180/np.pi), 'degrees')

#calculates V_i=x_i*e^(iωt) for i=1,2,3    
# t=np.linspace(0,50,1000) #time axis
# def V1_L(t):
#     return 1.9605003337204925*np.e**((-0.44552291170495756j)*t)
# def V2_L(t):
#     return 0.707392457519951*np.e**((-0.7820580716915961j)*t)
# def V3_L(t):
#     return 0.00650486597476126*np.e**((0.07607725045005106j)*t)

#plot of the circuit with inductor
plt.figure()
plt.title("V1,V2 and V3 with Inductor vs. time")
plt.xlabel("time(s)")
plt.ylabel("amplitude(V)")
#plt.plot(t,V1(t), label='V_1')
#plt.plot(t,V2(t),label='V_2')
#plt.plot(t,V3(t),label='V_3')
plt.plot(t,V(t,x[0],ω), label='V_1')
plt.plot(t,V(t,x[1],ω),label='V_2')
plt.plot(t,V(t,x[2],ω),label='V_3')
plt.legend()
plt.show()
plt.savefig("V1,V2 and V3 with inductor vs. time.png")
 
    
    
    
    
    
    


    