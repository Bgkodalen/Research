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

def smalldivisors(n):
    if n == 0:
        return [0]
    else:
        return [i for i in range(1,int(math.sqrt(n))+1) if (n % i == 0)]

def commondivisors(n,m):
    x = min(m,n)
    if m*n == 0:
        return [0]
    else:
        return [i for i in range(1,int(x)+1) if (m % i == 0 ) and (n % i == 0)]

def isint(x):
    tol = 10**-8
    if (np.round(x)-x)**2<tol:
        return 1
    else:
        return 0

def Lsm(Q):
    ### Generate L* from Q
    ### Ls = [q^i_{j,k}]_{i,k}
    X = sum(Q[0,:])
    P = np.linalg.inv(Q)*X
    d = Q.shape[0]
    Ls = np.zeros((d,d,d))
    for i in range(d):
        for j in range(d):
            for k in range(d):
                Ls[j,i,k] = sum([P[0,l]*Q[l,i]*Q[l,j]*Q[l,k] for l in range(d)])/(Q[0,i]*sum(P[0,:]))
    return Ls

maxvert = 92
def check(v1,p111,p112,p113,p121,p122,p123,p131,p132,p133):
    valid = 1
    dummy =0
    if min(v1,p111,p112,p113,p121,p122,p123,p131,p132,p133)<0:
        valid = 0
    #v2 = p112*v1/p121
    #v3 = p113*v1/p131
    #p212 = p122*v2/v1
    #p213 = p123*v2/v1
    #if p112+p212+p213 != v2:
    #    valid = 0
    #p312 = p132*v3/v1
    #p313 = p133*v3/v1
    #if p312 != p213:
    #    valid = 0
    #if p133*v3/v1 != p313:
    #    valid = 0
    #if p113+p312+p313 != v3:
    #    valid = 0
    if v1+v2+v3+1>maxvert:
        valid = 0
    return valid

def checkeig(L1,X):
    valid = 1
    tol = 10**-8
    [p,Qt] = np.linalg.eig(L1)
    if np.abs(Qt[0,:]).min()>tol:
        mults = [sum(X)/np.dot(Qt[:,i]*Qt[:,i]/(Qt[0,i]**2),X) for i in range(4)]
    else:
        mults = [0]
        valid = 0
    for value in mults:
        if np.isreal(value):
            if isint(np.real(value))==0:
                valid = 0
        else:
            valid = 0
    if valid == 1:
        mults = np.round(mults)
        indx = [i for i in range(4) if mults[i]==1]
        if len(indx) == 1:
            Qt[:, 0], Qt[:, indx[0]] = Qt[:, indx[0]], Qt[:, 0].copy()
            mults[0],mults[indx] = mults[indx],mults[0]
            Q = np.array([[Qt[i,j]*mults[j]/Qt[0,j] for j in range(4)] for i in range(4)])
            if X[1] == 20 and X[2] == 30:
                print(Q)
            [e,t] = np.linalg.eig(Q)
            if np.abs(e).min()>tol:
                Ls = Lsm(Q)
                if Ls.min()<-tol:
                    valid = 0
            else:
                valid = 0
        else:
            valid = 0
    return valid

### Try 2

maxv1 = 21
numschemes = 0
for v1 in range(1,maxv1):
    for p111 in range(v1-1):
        for p112 in range(1,v1-p111-1):
            p113 = v1-p111-p112-1
            for p121 in smalldivisors(v1*p112):
                if p121>0:
                    v2 = int(np.round(v1*p112/p121))
                    denom = Fraction(v2/v1).limit_denominator(v1).denominator
                    for mult in range(int(np.round(v1/denom)+1)):
                        p122 = mult*denom
                        p123 = v1-p121-p122
                        if p123 > 0:
                            for v3 in commondivisors(p113*v1,p123*v2):
                                p131 = p113*v1/v3
                                p132 = p123*v2/v3
                                p133 = v1-p131-p132
                                if check(v1,p111,p112,p113,p121,p122,p123,p131,p132,p133):
                                    L1 = np.array([[0,v1,0,0],[1,p111,p112,p113],[0,p121,p122,p123],[0,p131,p132,p133]])
                                    if checkeig(L1,[1,v1,v2,v3]):
                                        numschemes+=1
                                        if v1 == 20 and v2 == 30:
                                            print(L1)
                                        #print(L1)
print(numschemes)