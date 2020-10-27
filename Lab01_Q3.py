# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 09:33:25 2020

@author: Brendon
"""

# Import numpy for array functions
# Import matplotlib for plotting functions
# Import time for timing
import numpy as np
from matplotlib import pyplot as plt
from time import time


#q3 timing maxtrix multiplication


def time_mm(N): #timing matrix multiplication with argument N (matrix size)
    start=time()

    
    A=np.ones([N,N],float)*2 #constant matrix A (all 2 entries)
    B=np.ones([N,N],float)*5 #constant matrix B (all 5 entries)
    C=np.zeros([N,N],float)

    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i,j]+= A[i,k]*B[k,j]
                
    
    end=time()
    diff=end-start
    return diff

Nmax = 200     
            
N_array=np.arange(2,Nmax+1,1)
t_array=np.zeros(Nmax-1)

for N in N_array:
    t_array[N-2]=time_mm(N)

print(N_array)
print(t_array)

#plots
plt.figure()
plt.plot(N_array, t_array, 'k-')
plt.xlabel("Matrix Size, $N$")
plt.ylabel("Computation Time (s)")
plt.grid(True)
plt.suptitle("Computation Time vs Matrix Size")
plt.savefig('TimevsN.png')


plt.figure()
plt.plot(N_array**3, t_array, 'k-')
plt.xlabel("Cube of Matrix Size, $N^3$")
plt.ylabel("Computation Time (s)")
plt.grid(True)
plt.suptitle("Computation Time vs Cube of Matrix Size")
plt.savefig('TimevsN3.png')


plt.figure()
plt.loglog(N_array, t_array, 'k-')
plt.xlabel("Matrix Size, $N$")
plt.ylabel("Computation Time (s)")
plt.grid(True)
plt.suptitle("Computation Time vs Matrix Size")
plt.savefig('TimevsNLogLog.png')


plt.figure()
plt.plot(np.log(N_array), np.log(t_array), 'k-')
plt.xlabel("Log Matrix Size, $\log(N)$")
plt.ylabel("Log Computation Time")
plt.grid(True)
plt.suptitle("Log Computation Time vs Log Matrix Size")
plt.savefig('LogTimevsLogN.png')


def time_mm_np(): #timing the np.dot module to do the same matrix multiplication A dot B
   
    start=time()
    
    A=np.ones([N,N],float)*2 #constant matrix A (all 2 entries)
    B=np.ones([N,N],float)*5 #constant matrix B (all 5 entries)
    
    np.dot(A,B)
    
    end=time()

    diff_2=end-start
    return diff_2
    
    

N_max = 500     
            
Narray=np.arange(2,N_max+1,1)
tarray=np.zeros(N_max-1)

for N in Narray:
    tarray[N-2]=time_mm_np()

print(Narray)
print(tarray)

#plots 
plt.figure()
plt.plot(Narray, tarray, 'k-')
plt.xlabel("Matrix Size, $N$")
plt.ylabel("Computation Time (s)")
plt.grid(True)
plt.suptitle("Computation Time vs Matrix Size for Numpy Dot Product")
plt.savefig('TimevsN (Numpy Dot Product).png')


plt.figure()
plt.plot(Narray**3, tarray, 'k-')
plt.xlabel("Cube of Matrix Size, $N^3$")
plt.ylabel("Computation Time (s)")
plt.grid(True)
plt.suptitle("Computation Time vs Cube of Matrix Size for Numpy Dot Product")
plt.savefig('TimevsN3(Numpy Dot Product).png')


