# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import statistics as sts

l0=282.8
#1.calibration
l=[587.6,656.0,667.8,706.5,504.8,501.6,492.2,471.3,447.1,438.8]
y= [8.61, 7.65,7.20,6.77,11.05,11.22,11.65,12.70,14.25,14.90]
l.sort()
y.sort(reverse=True)
print(y)

#2.Hartman Relation fit and user input
def f(x,m,b):
    return (m/(x-282.8))+b
xdata = np.linspace(430,720,10)
plot=plt.plot(l,y)
popt,pcov = curve_fit(f,xdata,y)
plt.plot(xdata,f(xdata,*popt))
print(popt)

#user input function
def g(y):
    return (1882.62945896/(y-3.36434391)) + 282.8
    """return (282.8*3.36434391 -282.8*y -1882.62945896 )/(3.36434391-y)"""

#3.rydberg constant
y2=[17.7, 15.4,12.0,7.4]#hydrogen spectral line measurements (from violet to red)
y21=[]
for i in y2:
    y21.append(g(i))
plt.plot(y21,y2)



print ((g(17.7), g(15.4),g(12.0),g(7.4)))

R1=1/(y21[3]*((1/4)-(1/3**2)))
R2=1/(y21[2]*((1/4)-(1/4**2)))
R3=1/(y21[1]*((1/4)-(1/5**2)))
R4=1/(y21[0]*((1/4)-(1/6**2)))
r=print ([R1,R2,R3,R4])

print((0.009608981421209509 + 0.010649488303161057 + 0.0108417053690311 + 0.010866285390655475)/4) #mean rydberg

#verifying balmer series 

print(0.010491615121014285*(1/4-1/9)) #n=3
print(0.010491615121014285*(1/4-1/16)) #n=4
print(0.010491615121014285*(1/4-1/25))#n=5
print(0.010491615121014285*(1/4-1/36)) #n=6

print(1/414.12495974657554)
print(1/439.2210081180543)
print(1/500.80653469052174)
print(1/749.2989823154134)

#4. gas identification (tube #11)
y3= [7.65,7.70,7.85,8.10,8.20,8.75,9.90,10.0,12.05]#unknown spectral line measurements(from red to blue)
y31=[]
for i in y3:
    y31.append(g(i))
plt.plot(y31,y3)

print(g(7.65), g(7.7),g(7.85),g(8.10), g(8.2),g(8.75), g(9.9), g(10.0),g(12.05))

print((6.62607015*10**-34)*(299792458)*(0.010491615121014285))

