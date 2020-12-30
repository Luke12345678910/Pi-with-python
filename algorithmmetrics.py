#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparison of:
Chudnovsky algorithm
Machin's formula
brent-salamin formula

To compute Pi

Using GMP arbitrary-precision arithmetic 
Plotted by matplotlib.pyplot

@author: luke
"""

import time
import gmpy2
from gmpy2 import mpz
from gmpy2 import mpfr
import matplotlib.pyplot as plt

#not needed if unless curve fitting
#import numpy as np
#from scipy.optimize import curve_fit


#######################################################################################################################    

def binarysplit(a, b):
    P = mpz(1)
    Q = mpz(1)
    T = mpz(1)
    a = mpz(a)
    b = mpz(b)
    if (b - a) == 1:
        if(a != 0):
            P = (6*a-5)*(2*a-1)*(6*a-1)
            Q = a**3*10939058860032000
        T = P * (13591409 + a*545140134)
        if (a%2) == 1:
            T *= -1
    else:
        mid = (a + b) // 2
        P1, Q1, T1 = binarysplit(a, mid)
        P2, Q2, T2 = binarysplit(mid, b)
        P = P1 * P2
        Q = Q1 * Q2
        T = Q2 * T1 + P1 * T2
    return P, Q, T

def calc_pi_chud(digits):
    digits_per_iteration = 14.1816474627254776555
    iterations = mpz(digits/digits_per_iteration + 1)
    P, Q, T = binarysplit(0, iterations)
    pi = mpz(Q*426880*gmpy2.isqrt(mpz(10)**(digits*2)*10005)) // T
    return pi
    
#########################################################################################################################
    
#Calculate arctan(1/x) to specified accuracy with binary splitting the series
def arctanbs(a, b, base):
    R = mpz(1)
    Q = mpz(1)
    P = mpz(1)
    a = mpz(a)
    b = mpz(b)
    q = mpz(base)
    if (b - a) == 1:
        R = (2*a+3)
        Q = (2*a+3)*(q**(2*(b-a)))
        P = (-1)**(a+1)*R*q**(2*(b-a-1))/(2*a+3)

    else:
        mid = (a + b) // 2
        R1, Q1, P1 = arctanbs(a, mid, q)
        R2, Q2, P2 = arctanbs(mid, b, q)
        R = R1 * R2
        Q = Q1 * Q2
        P = Q2 * P1 + R1 * P2
    return R, Q, P
#calculate pi using machin's formula
def calc_pi_mach(digits):
    #set accuracy of divisions
    gmpy2.get_context().precision = int(3.4*digits)
    iterations = mpfr(digits)/1.38
    R, Q, P = arctanbs(0, iterations, 5)
    R2, Q2, P2 = arctanbs(0, iterations, 239)
    pi = 16*mpfr(1)/5*(1+mpfr(P)/Q) - 4*mpfr(1)/239*(1+mpfr(P2)/Q2)
    return pi

#Calculate arctan(1/x) to specified accuracy
def arctan(x, digits):
    #1/x - 1/x^3 + 1/x^5 + ...
    x = mpz(x)
    xsquared = -1*x**2
    xterm = x
    #1 + 3 + 5 + 7 + ...
    constterm = mpz('1')
    term = mpz('1')
    total = mpz('0')
    #fixed point divisions
    thing = mpz(10)**(digits+3)
    #main loop
    while(abs(term) > 0):
        term = thing//(xterm*constterm)
        total += term
        xterm = xterm*xsquared
        constterm += 2
    return total

def calc_pi_mach2(digits):
    #pi = 16*arctan(1/5) - 4*arctan(1/239)
    pi = 16*arctan(5, digits) - 4*arctan(239, digits)
    return pi

#########################################################################################################################
    
def calc_pi_brent(iterations):
    digits = 1.38*2**iterations
    gmpy2.get_context().precision = int(digits*3.4)
    a = mpfr('1')
    b = mpfr('1')/(gmpy2.sqrt(2))
    d = mpfr('0')
    for i in range(iterations):
        wa = a
        wb = b
        a = (wa+wb)/2
        b = gmpy2.sqrt(wa*wb)
        d += (2**(i+1))*((a**2)-(b**2))
    pi = (4*(a**2))/(1-2*d)
    return pi
    
#########################################################################################################################
    
def collect_data(algorithm, stop, increment, timing, space):
    timedatax = []
    timedatay = []
    memorydatax = []
    memorydatay = []
    i = 10
    while (i<=stop):
        if timing:
            start_time = time.time()
        if space:
            tracemalloc.start()
        if algorithm == 0:
            calc_pi_chud(i)
        elif algorithm == 1:
            calc_pi_mach(i)
        elif algorithm == 2:
            calc_pi_mach2(i)
        if timing:
            duration = time.time() - start_time
            timedatax.append(i)
            timedatay.append(duration)
        if space:
            current, peak = tracemalloc.get_traced_memory()
            peak = peak/10**6
            memorydatax.append(i)
            memorydatay.append(peak)
            tracemalloc.stop()
        print("done with " + str(i))
        i+=increment
    if timing:
        return timedatax, timedatay
    if space:
        return memorydatax, memorydatay
    
def collect_data2(iterations, increment, timing, space):
    timedatax = []
    timedatay = []
    memorydatax = []
    memorydatay = []
    i = increment
    while (i<=iterations):
        if timing:
            start_time = time.time()
        if space:
            tracemalloc.start()
        calc_pi_brent(i)
        if timing:
            duration = time.time() - start_time
            timedatax.append(1.364*2**i)
            timedatay.append(duration)
        if space:
            current, peak = tracemalloc.get_traced_memory()
            peak = peak/10**6
            memorydatax.append(1.364*2**i)
            memorydatay.append(peak)
            tracemalloc.stop()
        
        print("done with " + str(i))
        i+=increment
    if timing:
        return timedatax, timedatay
    if space:
        return memorydatax, memorydatay
    

    
    
    

fig, ax = plt.subplots(1,2,figsize=(20,10))     
    
timedatax, timedatay = collect_data(1, 100010, 10000, True, False)
print("Machin")
print(timedatax)
print(timedatay)
print("--------------------------")
ax[0].plot(timedatax, timedatay)
timedatax, timedatay = collect_data2(22, 1, True, False)
print("Brent")
print(timedatax)
print(timedatay)
print("--------------------------")
ax[0].plot(timedatax, timedatay)
timedatax, timedatay = collect_data(0, 6000010, 600000, True, False)
print("Chud")
print(timedatax)
print(timedatay)
print("--------------------------")
ax[0].plot(timedatax, timedatay)














'''
timedatax, timedatay = collect_data(0, 6000010, 600000, False, True)
ax[1].plot(timedatax, timedatay)
timedatax, timedatay = collect_data(1, 10010, 1000, False, True)
ax[1].plot(timedatax, timedatay)
timedatax, timedatay = collect_data2(22, 1, False, True)
ax[1].plot(timedatax, timedatay)

'''

'''
timedatax, timedatay = collect_data(0, 1000000, 100000, True, False)
for i in range(0, len(timedatax)):
    timedatax[i] = timedatax[i]/10000
for i in range(0, len(timedatay)):
    timedatay[i] = timedatay[i]
x_data = np.array(timedatax) 
y_data = np.array(timedatay) 
print(x_data)
print(y_data)
def fit_func(x, a, b):
    return a*x*(np.log2(x)) + b

params = curve_fit(fit_func, x_data, y_data)

[a, b] = params[0]
print([a, b])
plt.scatter(x_data, y_data)
x = np.linspace(0,100)
plt.plot(x, a*x*np.log2(x)+b)
'''
'''
timedatax, timedatay = collect_data2(22, 1, True, False)
for i in range(0, len(timedatax)):
    timedatax[i] = timedatax[i]/100000
for i in range(0, len(timedatay)):
    timedatay[i] = timedatay[i]
x_data = np.array(timedatax) 
y_data = np.array(timedatay) 
print(x_data)
print(y_data)
def fit_func(x, a, b):
    return a*x*(np.log2(x)**2) + b

params = curve_fit(fit_func, x_data, y_data)

[a, b] = params[0]
print([a, b])
plt.scatter(x_data, y_data)
x = np.linspace(0,60)
plt.plot(x, a*x*np.log2(x)**2+b)
    '''











