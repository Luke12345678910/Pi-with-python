#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Brent-salamin formula to compute Pi
Using GMP arbitrary-precision arithmetic 
"""


import gmpy2
from gmpy2 import mpfr

#if overflowing precise text file, set to false
preccheck = True

#Calculation
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


#number of iteration to use
pi = calc_pi_brent(16)


pi = str(pi)
#check against a select digits from verified computations
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
                    
