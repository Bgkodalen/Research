#!/usr/bin/python


import numpy as np
import pickle
from fractions import Fraction
import math
import sympy
import time
import itertools
import Qpoly


tol = 10**-5

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

def geg1(n):
    k1 = n**3*(n**2-1)*((n**7+2*n**6-3*n**4-17*n**3+45*n**2+14*n-76)/(-2*(n**4-13*n**3+15*n**2+12*n-32)*(n**2-1)) + math.sqrt(n**10+4*n**9+6*n**8+2*n**7-35*n**6+22*n**5+145*n**4-72*n**2+32*n+16)/(-2*(n**4-13*n**3+15*n**2+12*n-32)))
    k2 = ((3*n**2-5)*n**4)/2
    return max(k1,k2)

def oldk(n):
    return int(((1/2)*n**5+n**4-3*n**2*(1/2)-3*n*(1/2)+1/2+(1/2)*math.sqrt(n**10+4*n**9+4*n**8-6*n**7-18*n**6-10*n**5+13*n**4+14*n**3-n**2-2*n+1))*n**2)

def generatelist(n,k0,k1):
    s=-n**2
    #maxk = geg1(n)
    maxk = oldk(n)
    possible = {}
    newpossible = {}
    for k in range(maxk):
        if k%n == 0:
            possible[k] = {}
            newpossible[k] = {}
    print(maxk)
    #for kn in range(k0,k1):
    for kn in range(1,maxk//n):
        k = kn*n
        if k%1000 == 0:
            print(k)
        if k<(9*n**2*(1/4)-13/4+3*math.sqrt(9*n**4-26*n**2+17)*(1/4))*n**4:
            lowbnd = max(2*k/(3*n**2)-n**2/3,1)
        else:
            lowbnd = (-n**4+math.sqrt(n**8-4*k*n**6+6*k*n**4+k**2)+k)/(2*n**2)
        upbnd = (k-n**2)/(n**2+n)
        for r in range(int(lowbnd),int(upbnd)+1):
            mu = int(r*s+k)
            v = int((k-r)*(k-s)/(k+r*s))
            f = int((n**2-1)*k*(k-s)/(mu*(r-s)))
            m = int((k-r)*(-s)/mu)
            g = int((k+f*r)/(n**2))
            if ((k-r)*(k-s)) % (k+r*s) == 0: #Integral v
                if ((k-r)*(n**2)) % mu == 0: #Integral m
                    if ((n**2-1)*k*(k-s)) % (mu*(r-s))==0: #Integral f
                        if (k+f*r) % n**2 == 0: #Integral g
                            if (k-r*n**2) %2 == 0: #p^1_21 is integral
                                if (mu+n*(r-n))*(n+1) % (2*n) == 0: #p^1_11 is integral
                                    if 1+f<=m*(m+1)/2: #absolute bound
                                        if (1+f+g<g*(g+1)/2) or ((k+s)*r**2+2*r*(k-s**2)-s*(k+s)==0 and 1+f<g*(g+1)/2): #absolute bound
                                            if v<=m*f: #absolute bound
                                                if mu>=n*(r+n): #positive p^1_31
                                                    possible[k][r] = 1
                                                    if 15*n**4*(2*n**2-3)*r**2+(n**6-45*k*n**2+76*k)*n**2*r+k*(16*k+n**6)*(n**2-2)>=0:# New bound
                                                        newpossible[k][r]=1
    posslist = [[k,r] for k in possible for r in possible[k]]
    newposslist = [[k,r] for k in newpossible for r in newpossible[k]]
    return posslist,newposslist