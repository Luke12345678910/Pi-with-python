#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Using arctan taylor polynomial and
pi/4 = 4arctan(1/5) - arctan(1/239)

Using GMP arbitrary-precision arithmetic 

@author: luke
"""

from gmpy2 import mpz

#if overflowing precise text file, set to false
preccheck = True


def arctan(x, digits):
    x = mpz(x)
    xsquared = -1*x**2
    xterm = x
    constterm = mpz('1')
    term = mpz('1')
    total = mpz('0')
    thing = mpz(10)**(digits+3)
    while(abs(term) > 0):
        term = thing//(xterm*constterm)
        total += term
        xterm = xterm*xsquared
        constterm += 2
    return total

def calc_pi_mach(digits):
    pi = 16*arctan(5, digits) - 4*arctan(239, digits)
    return pi

pi = calc_pi_mach(10000)

#check against select digits from verified computations
pi = str(pi)
pi = pi[:1] + '.' + pi[1:]

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

