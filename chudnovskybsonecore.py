#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chudnovski Algorithm to compute Pi
Using GMP arbitrary-precision arithmetic 
@author: luke
"""
import gmpy2
from gmpy2 import mpz
from gmpy2 import mpfr

#if overflowing precise text file, set to false
preccheck = False


#binary splitting
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

def calc_pi(digits):
    #experimentally calculated digits per iteration
    digits_per_iteration = 14.1816474627254776555
    iterations = mpz(digits/digits_per_iteration + 1)
    P, Q, T = binarysplit(0, iterations)
    #final fixed point division
    pi = mpz(Q*426880*gmpy2.isqrt(mpz(10)**(digits*2)*10005)) // T
    return pi


pi = calc_pi(2000000)

#check against select digits from previous computations
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

