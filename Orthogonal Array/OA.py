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

check_complement = 1

success = 0
combs = list(combinations([i for i in range(n-1)],7))
indx = 0
wineven = []
winodd = []

if check_complement == 0:
    print('checking main graph\n')
    while indx < len(combs):
        if indx % 100 == 0:
            print(indx/len(combs))
        columns = [0,1]+[i+2 for i in combs[indx]]
        A = Build_adj(O[:,columns])
        C = (A[n:,0:n] @ A[0:n,n:])[A[n:,n:]==1]
        u = np.unique(C)
        #print(u)
        #print(np.unique(u % 2))
        if len(np.unique(u % 2)) == 1:
            if np.unique(u % 2)[0] == 0:
                wineven.append(columns)
            elif np.unique(u % 2)[0] ==1:
                winodd.append(columns)
            print(columns)
            print(np.unique(u % 2))
        indx+=1
else:
    print('checking complement\n')
    while indx < len(combs):
        if indx % 100 == 0:
            print(indx/len(combs))
        columns = [0,1]+[i+2 for i in combs[indx]]
        A = Build_adj(O[:,columns])
        total_index = {i for i in range(n+1)}
        acceptable = total_index-set(columns)
        next_square = acceptable.pop()
        coclique = [i for i in range(n*n) if O[i,next_square] == 1]
        remaining = [i for i in range(n*n) if i not in coclique]
        sub_A = A[:,coclique]
        sub_A = sub_A[remaining,:]
        rem_A = A[remaining,:]
        rem_A = rem_A[:,remaining]-1*np.eye(len(remaining))
        C = (sub_A @ np.transpose(sub_A))[rem_A==0]
        u = np.unique(C)#[:-1]
        #if np.unique(C)[-1] != 9:
        #    print(np.unique(C))
        #if len(u)<3:
        print(u)
        #print(u)
        #print(np.unique(u % 2))
        if len(np.unique(u % 2)) == 1:
            if np.unique(u % 2)[0] == 0:
                wineven.append(columns)
            elif np.unique(u % 2)[0] ==1:
                winodd.append(columns)
            print(columns)
            print(np.unique(u % 2))
        indx+=1
    
print('even')
print(wineven)
print('odd')
print(winodd)
