#!/usr/bin/python

import numpy as np
import math
from itertools import combinations

### Initialize OA
def Initialize_OA(n,k):
    O = np.zeros([n*n,2])
    for row in range(n*n):
        O[row,0] = row//n
        O[row,1] = (row % n)
    return O

def Prime_OA(n):
    ### Only works if n is prime
    O = np.zeros([n*n,n+1])
    for row in range(n*n):
        O[row,0] = row//n
        O[row,1] = row % n
    for alpha in range(1,n):
        for e in range(n*n):
            O[e,alpha+1] = (O[e,0]+O[e,1]*alpha) % n
    return O

def Build_adj(O):
    n = int(math.sqrt(O.shape[0]))
    k = O.shape[1]
    A = np.zeros([int(n*n),n*n])
    for col in range(k):
        for j in range(n*n):
            for h in range(j+1,n*n):
                if O[j,col] == O[h,col]:
                    A[j,h] = 1
                    A[h,j] = 1
    return A

def Build_next_col(O):
    n = int(math.sqrt(O.shape[0]))
    k = O.shape[1]
    A = Build_adj(O)

n=17
O = Prime_OA(n)

success = 0
combs = list(combinations([i for i in range(n-1)],9))
indx = 0
win = []
while indx < len(combs):
    if indx % 100 == 0:
        print(indx/len(combs))
    columns = [0,1]+[i+2 for i in combs[indx]]
    A = Build_adj(O[:,columns])
    C = 17-((A[n:,0:n] @ A[0:n,n:]) * A[n:,n:])
    u = n-np.unique(C)
    #print(u)
    if len(np.unique(u % 2)) == 1 and np.unique(u % 2)[0] == 0:
         win.append(columns)
         print(win)
    else:
        indx+=1
