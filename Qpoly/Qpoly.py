#!/usr/bin/python

### This file is intended to be a one-off method for examining the various Q-polynomial association schemes.
# This file no longer imports data from Williford's tables. This is not intended to be a workplace for the
# Schemes found in various locations.

from urllib.request import urlopen
import re
import pickle
import numpy as np
from fractions import Fraction
import math
import sympy
import time
import itertools

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

def kreinarray(X,vals):
    d = int(len(X)/2)
    L1s = np.zeros((d+1,d+1))
    for i in range(d):
        L1s[i,i+1] = X[i]
        L1s[i+1,i] = X[i+d]
        if i>0:
            L1s[i,i] = L1s[0,1]-L1s[i,i+1]-L1s[i,i-1]
    L1s[d,d] = L1s[0,1]-L1s[d,d-1]
    [q,Pt] = np.linalg.eig(L1s)
    order = np.argsort(q)
    order = order[::-1]
    Pt[:,[i for i in range(d+1)]] = Pt[:,order]
    for i in range(d+1):
        Pt[:,i] = Pt[:,i]*vals[i]/Pt[0,i]
    return Pt

def Gegproj(Ls,k=10,verbose=0,binary=0):
### Here, we take in L* and compute our Gegenbauer projections for degrees 0 up to k.
    tol = 10**(-12)
    q = 1
    if len(Ls.shape) == 3: ### If you gave all of L*
        Ls = Ls[q,:,:]
    side = Ls.shape[0]
    n = Ls[0,q]
    Projections = np.zeros((side,k+1))
    Projections[0,0] = 1
    Projections[q,1] = 1/n
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

def Equiangular(Q):
    ### Look for Equiangular lines
    success = 1
    string = ''
    tol = 10**-8
    d = Q.shape[0]-1
    indx = Q[0,:].argmax()
    goodindx = [i for i in range(d+1)]
    goodindx.pop(indx)
    numlines = sum(Q[0,:])
    A = np.array(Q[1:,goodindx],dtype='float')
    B = np.linalg.inv(A)
    for innerprods in list(itertools.product([-1,1],repeat = d)):
        rank = sum(Q[0,goodindx])
        if sum(innerprods) == d:
            continue
        x = np.dot(B,innerprods)
        if sum(x<-tol) == 0:
            for i in range(len(x)):
                if x[i]<tol:
                    rank = rank-Q[0,goodindx[i]]
            ip = np.dot(Q[0,goodindx],x)
            if abs((ip-1)/2-np.rint((ip-1)/2))>tol:
                success = 0
            string = string+ "%0.f lines in %0.f with angle 1/%0.3f using coeffs: %s\n" % (numlines, rank, np.dot(Q[0,goodindx],x),innerprods)
    if len(string) == 0:
        return ["No possible Equiangular lines",success]
    else:
        return[string,success]


    
