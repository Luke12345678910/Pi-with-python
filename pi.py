#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Estimating pi through randomly placing points inside a
square inscribed by a circle

Very slow but kind of cool

@author: luke
"""

import random
count = 0
iterations = 10000000
for i in range(iterations):
    a = random.random()
    b = random.random()
    if ((a*a + b*b) <= 1):
        count+=1
    
pi = count/(i+1) * 4
print(pi)