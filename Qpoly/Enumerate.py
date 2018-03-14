#!/usr/bin/python


import numpy as np
import pickle
from fractions import Fraction
import math
import sympy
import time
import itertools

def is_prime(n):
  if n == 2 or n == 3: return True
  if n < 2 or n%2 == 0: return False
  if n < 9: return True
  if n%3 == 0: return False
  r = int(n**0.5)
  f = 5
  while f <= r:
    #print('\t',f)
    if n%f == 0: return False
    if n%(f+2) == 0: return False
    f +=6
  return True

maxvert = 70
#primes = [i for i in range(maxvert) if is_prime(i)]
primes = [i for i in range(6) if is_prime(i)]
#maxval = int(math.log(1000))/math.log(2))+1
maxval = 2
t0 = time.time()
Matrices = [np.array([[x[0], x[6], x[7]],[x[1], x[3], x[4]+x[6]-x[7]],[x[2], x[4], x[5]]]) for x in itertools.product(range(maxval),range(maxval),range(maxval),range(maxval),range(maxval),range(maxval),range(maxval),range(maxval)) if x[4]+x[6]-x[7]>=0]
t1 = time.time()
print('Step 1:%0.f seconds' % int(t1-t0))

Matrices = [mat for mat in Matrices if mat.max()>0]
t2 = time.time()
print('Step 2:%0.f seconds' % int(t2-t1))

dictionary = {0:np.array([[1,1,1],[1,1,1],[1,1,1]])}
i = 0
o = np.ones((3,1))
for prime in primes:
    print('Starting prime: %0.f' % prime)
    newmax = math.log(maxvert)/math.log(prime)
    newmats = [prime**mat for mat in Matrices if mat.max()<newmax]
    newdict = {}
    for newm in newmats:
        for indx in dictionary:
            pot = dictionary[indx]*newm
            if np.dot(pot,o).max()<maxvert/2:
                newdict[i+1] = pot
                i+=1
    for key in newdict:
        dictionary[key] = newdict[key]
    print('Dictionary size: %0.f' % len(dictionary))

t3 = time.time()
print('Step 3:%0.f seconds' % int(t3-t2))

for mat in range(len(dictionary)):
    sums = np.dot(dictionary[mat],o)
    if sums[0]+1==sums[1] and sums[1]==sums[2]:
        print(dictionary[mat])
        print('\n')
    else:
        dictionary.pop(mat)
print('Dictionary size: %0.f' % len(dictionary))