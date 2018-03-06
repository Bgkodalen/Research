#!/usr/bin/python

### This file is intended to be a one-off method for examining the various Q-polynomial association schemes.
# This file no longer imports data from Williford's tables. This is not intended to be a workplace for the
# Schemes found in various locations.

from urllib.request import urlopen
import re
import pickle
import numpy as np
from fractions import Fraction
from math import sqrt
import sympy
import time

def Qm(P):
    X = sum(P[0,:])
    Q = np.linalg.inv(P)
    return Q*X

def Lsm(P):
    ### Generate L* from Q
    ### Ls = [q^i_{j,k}]_{i,k}
    Q = Qm(P)
    d = Q.shape[0]
    Ls = np.zeros((d,d,d))
    for i in range(d):
        for j in range(d):
            for k in range(d):
                Ls[j,i,k] = sum([P[0,l]*Q[l,i]*Q[l,j]*Q[l,k] for l in range(d)])/(Q[0,i]*sum(P[0,:]))
    return Ls

def Lm(P):
    ### Generate L from P
    ### L = [p^i_{j,k}]_{i,k}
    tol = 10**-14

    Q = Qm(P)
    d = P.shape[0]
    L = np.zeros((d,d,d))
    for i in range(d):
        for j in range(d):
            for k in range(d):
                L[j,i,k] = sum([Q[0,l]*P[l,i]*P[l,j]*P[l,k] for l in range(d)])/(P[0,i]*sum(P[0,:]))
    Lt = np.int64(np.rint(L))
    if np.absolute(L-Lt).max()>tol:
        return L
    return Lt

def Gegproj(Ls,k=10,verbose=0,binary=0):
### Here, we take in L* and compute our Gegenbauer projections for degrees 0 up to k.
    tol = 10**(-14)
    if len(Ls.shape) == 3: ### If you gave all of L*
        Ls = Ls[1,:,:]
    side = Ls.shape[0]
    n = Ls[0,1]
    Projections = np.zeros((side,k+1))
    Projections[0,0] = 1
    Projections[1,1] = 1/n
    for i in range(2,k+1):
        Projections[:,i] = (2*i+n-4)/(i+n-3)/n*np.dot(Ls,Projections[:,i-1])-(i-1)/(i+n-3)*Projections[:,i-2]
    if verbose ==1:
        print(Ls)
        print(n)
        print(Projections)
    if binary == 1:
        projections2 = np.chararray((side,k+1))
        for i in range(side):
            for j in range(k+1):
                if Projections[i,j]>tol:
                    projections2[i,j] = '+'
                elif Projections[i,j]<-tol:
                    projections2[i,j] = '-'
                elif Projections[i,j]>=-tol and Projections[i,j]<=tol:
                    projections2[i,j] = '0'
        return projections2
    return Projections

def Gnk(n,k,dim=1):
### This function builds the Gegenbauer polynomials recursively. Note this is the normalized two-term recurrence where G_1(t) = t.
    if dim == 1:
        t = sympy.Symbol('t')
    else:
        t = sympy.MatrixSymbol('t',dim,dim) ### allows for t to be a matrix input

    if k == 0:
        return 1
    elif k == 1:
        return t
    else:
        return sympy.simplify((2*k+n-4)/(k+n-3)*t*Gnk(n,k-1) + (k-1)/(k+n-3)*Gnk(n,k-2))