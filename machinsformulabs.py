#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Machin Formula to Compute Pi
Using Arctan Taylor Polynomial and
pi/4 = 4arctan(1/5) - arctan(1/239)
Implemented Binary Splitting Optimization

Using GMP arbitrary-precision arithmetic 

@author: luke
"""

import gmpy2
from gmpy2 import mpz
from gmpy2 import mpfr

#if overflowing precise text file, set to false
preccheck = True

#Binary splitting optimization
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
    gmpy2.get_context().precision = int(4*digits)
    iterations = mpfr(digits)/1.38
    R, Q, P = arctanbs(0, iterations, 5)
    R2, Q2, P2 = arctanbs(0, iterations, 239)
    pi = 16*mpfr(1)/5*(1+mpfr(P)/Q) - 4*mpfr(1)/239*(1+mpfr(P2)/Q2)
    return pi

    #number of decimal digits to calculate
pi = calc_pi_mach(10000)

#check against select digits from verified computations
pi = str(pi)

if (not preccheck):
    with open("selectpidigits.txt") as f:
        content = f.readlines()
    content = [x.strip() for x in content] 
    content = [x.replace(" ", "") for x in content] 
    content = [x.replace(",", "") for x in content] 
    
    x = True
    for group in content:
        number_of_digits = int(group[50:(len(group))])
        if x:
            if (len(pi) < number_of_digits):
                print("Failed at " + str(number_of_digits))
                x = False
            if x:
                for i in range(0, 49):
                    if (group[i] != pi[number_of_digits-49+(i+1)]):
                        if x:
                            print("Failed at " + str(number_of_digits))
                            x = False
else:
    with open('pidigits.txt', 'r') as f:
        data = f.read()
    x = True
    i = 0
    while x:
        pi = str(pi)
        if (data[i] != pi[i]) or (i == len(pi)-1):
            x = False
        i = i + 1      
    print("Decimal digits of accuracy: " + str(i-2))
