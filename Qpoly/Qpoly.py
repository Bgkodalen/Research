#!/usr/bin/python

### This file is intended to be a one-off method for examining the various Q-polynomial association schemes. Ideally, there will be a paired 
# GUI file which will call the various functions used here. As such, we try to keep any code restricted to function calls, allowing the file
# to be imported without much difficulty.

from urllib.request import urlopen
import re
import pickle
import numpy as np
from fractions import Fraction
from math import sqrt
import sympy

class r:
### This class is used to define constants used for Irrational schemes. Each attribute allows for a new constant to be defined.
#These are recycled for each scheme, so only one is needed, though two extra are provided in case it is needed later.
    def __init__(self):
        self.a1 = ''
        self.a2 = ''
        self.a3 = ''


def GetSchemes():
### This function is used to generate the dictionary holding all the information from Williford's table.
# The first part grabs all possible paramater sets, while the second generates the desired dictionary.
    fivebipstr = urlopen('http://www.uwyo.edu/jwilliford/data/qbip5_table.html').read().decode('utf-8')
    fivebip = re.findall(r'<\d+,\d+>[a-z]?',fivebipstr)

    fourbipstr = urlopen('http://www.uwyo.edu/jwilliford/data/qbip4_table.html').read().decode('utf-8')
    fourbip = re.findall(r'<\d+,\d+>[a-z]?',fourbipstr)

    threeprimstr = urlopen('http://www.uwyo.edu/jwilliford/data/qprim3_table.html').read().decode('utf-8')
    threeprim = re.findall(r'<\d+,\d+>[a-z]?',threeprimstr)


    Associationschemes = {3:threeprim,4:fourbip,5:fivebip}

    for key in Associationschemes:
        print('starting '+str(key)+'-class schemes\n')
        schemes = Associationschemes[key]
        Associationschemes[key] = {}
        for scheme in schemes:
            if scheme[-1] == '>':
                params = scheme.strip('<>').split(',')
            else:
                scheme = scheme.replace('>','')
                params = scheme.strip('<').split(',')
            S = Schemeparams(params[0],params[1],key)
            Associationschemes[key][scheme] = S


    pickle.dump(Associationschemes,open("schemes.p","wb"))
    return

def Schemeparams(v,m,d):
    # This function searches through Willifords tables to find the parameters of the d-class Q-poly scheme
    # with v vertices and multiplicity m. 

    uptol = 10**8

    base = {3:'http://www.uwyo.edu/jwilliford/data/qprim3/qpd3',
    4:'http://www.uwyo.edu/jwilliford/data/qbip4/qbd4',
    5:'http://www.uwyo.edu/jwilliford/data/qbip5/qbd5'}

    schemetxt = urlopen(base[d]+'v'+str(v)+'m'+str(m)+'.txt').read().decode('utf-8')
    
    constants = re.findall(r'r\.\d = sqrt\(\d*\)',schemetxt)
    for item in constants:
        exec(item.replace('.','.a'))

    schemetxt = schemetxt.replace('.','.a')
    schemetxt = re.sub(r'\s*\+\s*','+',schemetxt)


### I think these search calls can be simplified to simply check for brackets, but it is working now so I am not going to mess with it.
    Ptxt = re.search(r'P[\s\w]*=[\s\w\.\*()\+/\\\-\[\],]*\n\n',schemetxt)
    Prows = re.findall(r'\[[\.\*(\+)\s\d\w\-\\/]*\]',Ptxt[0].replace('\n',''))
    P = np.array([[int(num) if Fraction(eval(num)).denominator == 1 else Fraction(num) if Fraction(eval(num)).denominator<uptol else float(eval(num)) for num in Prows[i].strip('\[\]').split()] for i in range(d+1)])
    
    Qtxt = re.search(r'Q[\s\w]*=[\s\w\.\*()\+/\\\-\[\],]*\n\n',schemetxt)
    Qrows = re.findall(r'\[[\.\*(\+)\s\d\w\-\\/]*\]',Qtxt[0].replace('\n',''))
    Q = np.array([[int(num) if Fraction(eval(num)).denominator == 1 else Fraction(num) if Fraction(eval(num)).denominator<uptol else float(eval(num)) for num in Qrows[i].strip('\[\]').split()] for i in range(d+1)])

    Ltxt = re.search(r'L[\s\w]*=[\s\w\[\]\-,]*\]\n\n',schemetxt)
    Lrows = re.findall(r'\d[\s\d\-]*\d',Ltxt[0].replace('\n',''))
    L = np.array([[[int(num) for num in Lrows[i].split()] for i in range((d+1)*j,(d+1)*(j+1))] for j in range(d+1)])

    Lstxt = re.search(r'L\*[\.\s\w]*=[\s\w\.\*()\+/\\\-\[\],]*\](\n\n|$)',schemetxt)
    Lsrows = re.findall(r'\d[\.\s\d\-\\/]*\d',Lstxt[0].replace('\n',''))
    Ls = np.array([[[int(num) if Fraction(eval(num)).denominator == 1 else Fraction(num) if Fraction(eval(num)).denominator<uptol else float(eval(num)) for num in Lsrows[i].strip('\[\]').split()] for i in range((d+1)*j,(d+1)*(j+1))] for j in range(d+1)])
    
### Each scheme currently lists the P, Q, L, and L* matrices. More can be added later if needed.    
    scheme = {'P':P,'Q':Q,'L':L,'L*':Ls}

    return scheme

def Gegproj(Ls,k=10,verbose=0):
### Here, we take in L* and compute our Gegenbauer projections for degrees 0 up to k.
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