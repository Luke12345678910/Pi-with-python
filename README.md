# Pi-with-python

Various Scripts to Calculate Pi using python

Requires gmpy2 software(a python version of GMP Multi-precision reals and large integer arithmetic)

brentsalamin.py, chudnovskybsonecore.py, machinsformula.py, and machinsformulabs.py calculate pi using various algorithms.

chudnovskybsonecore.py is the fastest in most cases. 
An optimisation completed by storing the large integers used in the calculation in a fully or partially factored form has been shown to speed up the calculation, however it is not implemented due to the increased space requirements. 
Other optimisations such as multi-threading and decimal to binary conversions(used in software such as y-cruncher) are not implemented as well due to the fact they rely on the specific system.)

algorithmetrics.py compares BrentSalamin, Chudnovsky, and Machin algorithms.

Pi.py calculates pi by generating random number inside a square inscribed with a circle and calculating how many land in the circle. 
