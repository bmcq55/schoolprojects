# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 09:28:25 2020

@author: Brendon and Lin
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#1.thermistor ,R(ohms),T(C)
#T,R=[(1, 6.45M),(11,4.58M),(18,2.70M),(22,2.95M),(32,1.60M),(37,1.38M),(47,0.8M),(52,0.60M),(59,0.40M),(70,0.23M),(80,120k),(91,43k)]

#1/T=a1 + a2*x +a3*x*x +a4*x*x*x
T_data=np.array([1,11,18,22.0,32.0,37,47,52,59,70,80,91])+273.15
R_data=np.array([6.45,4.58,2.95,2.70,1.60,1.38,0.8,0.6,0.4,0.23,0.120,0.043])*10**6
x0=[0.0,0.0,0.0,0.0]
T_data2=1/T_data
unc=[]
for i in T_data:
    unc.append((0.5/i)*(1/i))
    
R_data2=np.log(R_data)
def f(x,a1,a2,a3,a4):#model fuction
    return a1 + a2*x +a3*(x**2) +a4*(x**3)

plt.scatter(R_data2,T_data2)
plt.ylim(0.0025,0.004)
plt.title("T and R curve of the thermistor")
plt.xlabel("Ln(Resistance),Ohm")
plt.ylabel("1/Temperature,K")

popt,pcov= curve_fit(f, R_data2, T_data2,x0,sigma=unc)
plt.plot(R_data2,f(R_data2,*popt))
print("the fit parameters are (a1,a2,a3,a4,)",popt)
plt.show()

#user input function
R= float(input("input resistance (ohms):"))
print("the temperature at the given resistance is:",1/f(np.log(R),*popt)-273.15,"K")

#chisquared 
resi=T_data2-f(R_data2,*popt)
chisq = sum((resi / unc) ** 2)
print("the value of Chi squared is",chisq)

#2.Diode direction 1 I(mA),V(V)
#I,V=[(0.16,0.717),(0.21,0.730),(0.42,0.761),(0.52,0.767),(0.80,0.784),(0.91,0.787),(0.96,0.788)]
I_data=np.array([0.16,0.21,0.42,0.52,0.80,0.91,0.96])*.001
V_data=np.array([0.717,0.730,0.761,0.767,0.784,0.787,0.788])
x0=[0]
q=1.6*10**-19
T=293.15
def g(V,I0):
    return I0*(np.e**((q*V)/(1.3*10**-23*T))-1)

plt.scatter(V_data,I_data)
plt.ylim(0,0.001)
popt1,pcov1=curve_fit(g,V_data,I_data,x0)
vdata=np.linspace(.7,.8,20)
plt.plot(vdata,g(vdata,*popt1))
plt.title("I and V curve across the diode")
plt.xlabel("V,V")
plt.ylabel("I,A")
print("the value of I0 is",popt1,"A")
#The fitting program of curve_fit was unable to find a fit for k, thus it has to be manually provided.

#chisquared of diode
resid=I_data-g(V_data,*popt)
chi2=sum((resid/unc))**2


#diode direction 2 I(mA),V(V)
#I,V=[(0.001,4.83),(0.002,15.17)]
#we could not gather sufficient data for analysis due to the limitation of the multimeter at low current, however it is clear that this is the reverse direction from the low current