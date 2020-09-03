# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 17:27:21 2019

@author: Brendon
"""
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

#RC circuit, DC voltage, q1, ch1- applied voltage
t_ch1_RC,V_ch1_RC= np.loadtxt('C:/Users/Brendon/Desktop/RC_DC_exp_1_ch_1.txt', skiprows=2, unpack=True)
plt.figure()
plt.plot(t_ch1_RC, V_ch1_RC)
plt.title('RC Circuit Channel 1 (DC): Applied Voltage')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')

#RC circuit, DC voltage, q1, ch2- voltage across capacitor
t_ch1_RC,V_ch1_RC= np.loadtxt('C:/Users/Brendon/Desktop/RC_DC_exp_1_ch_2.txt', skiprows=2, unpack=True)
plt.figure()
plt.plot(t_ch1_RC, V_ch1_RC)
plt.title('RC Circuit Channel 2(DC): Voltage Across Capacitor')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')

def test_func(t_ch1_RC, a, b):
    return a*(1-np.exp(-t_ch1_RC *b))

popt,pcov= curve_fit(test_func,t_ch1_RC, V_ch1_RC, p0=[2,2])

print(popt)

#RC circuit, exp 1, q2, ch1-applied voltage
t_ch1_RC,V_ch1_RC= np.loadtxt('C:/Users/Brendon/Desktop/RC_ch_1_ex1.txt', skiprows=2, unpack=True)
plt.figure()
plt.plot(t_ch1_RC, V_ch1_RC)
plt.title('RC Circuit Channel 1: Applied Voltage')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')


#RC circuit, exp 1, q2, ch2-resistor voltage
t_ch2_RC,V_ch2_RC= np.loadtxt('C:/Users/Brendon/Desktop/RC_ch_2_ex1.txt', skiprows=2, unpack=True)
plt.figure()
plt.plot(t_ch2_RC, V_ch2_RC)
plt.title('RC Circuit Channel 2: Resistor Voltage')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')

#RC circuit,exp.1 q2, resistor voltage fit
(time, measured_voltage)= np.loadtxt('C:/Users/Brendon/Desktop/sample from RC circuit, time, resistor voltage.txt', skiprows=2, unpack=True)
(measured_V_applied)= np.loadtxt('C:/Users/Brendon/Desktop/sample from RC circuit, applied voltage, q2,exp1.txt', skiprows=2, unpack=True)

"""def voltage_resistor(t,a,b)
    #RC=0.01034 #470 kOhm x 22 nF
    return a * np.exp()


g=[0.0000001,0.001]
n=len(time)
y=np.empty(n)
for i in range(n):
    y[i]=voltage_resistor(time[i],g[0])
    
t=time.values
V_m=measured_voltage.values
popt,pcov= curve_fit(voltage_resistor, t, V_m,g)
print(popt)


plt.figure()
plt.plot(time,measured_voltage)
#ime,y, 'ro')

popt,pcov= curve_fit(voltage_resistor,time, measured_voltage, measured_V_applied)
plt.figure()
plt.plot(voltage_resistor)"""



#LR circuit, exp 1, q3, ch1
t_ch1_LR,V_ch1_LR,sample= np.loadtxt('C:/Users/Brendon/Desktop/LR_Ch_1_Ex1.txt', skiprows=2, unpack=True)
print(t_ch1_LR,V_ch1_LR,sample)
plt.figure()
plt.plot(t_ch1_LR,V_ch1_LR)
plt.title('LR Circuit Channel 1: Resistor Voltage')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')

#LR ciruit, exp.1 q3, Ch2
t_ch2_LR,V_ch2_LR= np.loadtxt('C:/Users/Brendon/Desktop/LR_ch_2_ex1.txt', skiprows=2, unpack=True)

plt.figure()
plt.plot(t_ch2_LR, V_ch2_LR)
plt.title('LR circuit Channel 2: Applied Voltage')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')


#LC circuit, exp.1, q4, ch1 voltage across resistor
t_ch1_LC,V_ch1_LC= np.loadtxt('C:/Users/Brendon/Desktop/LC_ch_1_ex1.txt', skiprows=2, unpack=True)
plt.figure()
plt.plot(t_ch1_LC,V_ch1_LC)
plt.title('LC Circuit Channel 1: Inductor Voltage')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')

#LC circuit, exp.1, q4, ch2, applied voltage
t_ch2_LC,V_ch2_LC= np.loadtxt('C:/Users/Brendon/Desktop/LC_ch_2_ex1.txt', skiprows=2, unpack=True)
plt.figure()
plt.plot(t_ch2_LC,V_ch2_LC)
plt.title('LC Circuit Channel 2: Applied Voltage')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')

#LC circuit, exp.1, q4. Ch1, applied voltage for capacitor
t_ch1_LC,V_ch1_LC= np.loadtxt('C:/Users/Brendon/Desktop/LC _ch_1_ex1_V_app.txt', skiprows=2, unpack=True)
plt.figure()
plt.plot(t_ch1_LC,V_ch1_LC)
plt.title('LC Circuit Channel 1: Applied Voltage (Capacitor) ')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')

#LC circuit, exp.1, q4,ch2, voltage across capacitor
t_ch2_LC,V_ch2_LC= np.loadtxt('C:/Users/Brendon/Desktop/LC_ch_2_exp1_V_C.txt', skiprows=2, unpack=True)
plt.figure()
plt.plot(t_ch2_LC,V_ch2_LC)
plt.title('LC Circuit Channel 1: Voltage Across Capacitor ')
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')

#1. LR coding

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def f(x, A, B):
    return B*np.exp(A * x)

def g(x, y):
    return y
    
decay_data_CH1 = np.loadtxt('C:/Users/Brendon/Desktop/LR_Ch_1_Ex1.txt', skiprows=2, unpack=True, delimiter = ' ')
decay_data_CH2 = np.loadtxt('C:/Users/Brendon/Desktop/LR_ch_2_ex1.txt', skiprows=2, unpack=True ,delimiter = ' ')

T1 = decay_data_CH1.T[0]
V1 = decay_data_CH1.T[1]
T2 = decay_data_CH2.T[0]
V2 = decay_data_CH2.T[1]

plt.plot(T1,V1,label="CH1")
plt.plot(T2,V2,label="CH2")
popt1, pcov1 = curve_fit(f,T1,V1)
popt2, pcov2 = curve_fit(g,T2,V2)

y1 = f(T1,*popt1)
y2 = f(T2,*popt2)

plt.plot(T1, y1, "red")
plt.plot(T2, y2, "green")
plt.title("The Voltage vs. Time Plot of Half Period")
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')

#2. RC coding

def f(x, A, B):
    return B*np.exp(A * x)

decay_data_CH1 = np.loadtxt('C:/Users/Brendon/Desktop/RC_ch_1_ex1.txt', skiprows=2, unpack=True, delimiter = ' ')
decay_data_CH2 = np.loadtxt('C:/Users/Brendon/Desktop/RC_ch_2_ex1.txt', skiprows=2, unpack=True, delimiter = ' ')

T1 = decay_data_CH1.T[0]
V1 = decay_data_CH1.T[1]
T2 = decay_data_CH2.T[0]
V2 = decay_data_CH2.T[1]

plt.plot(T1,V1,label="CH1: total voltage")
plt.plot(T2,V2,label="CH2: voltage of resistor")

popt1, pcov1 = curve_fit(f,T1,V1)
popt2, pcov2 = curve_fit(f,T2,V2)

y1 = f(T1,*popt1)
y2 = f(T2,*popt2)

plt.plot(T1, y1, "red")
plt.plot(T2, y2, "green")

print(popt2[0])

plt.title("The Voltage vs. Time Plot of Half Period")
plt.xlabel('Time(s)')
plt.ylabel('Voltage(V)')


















