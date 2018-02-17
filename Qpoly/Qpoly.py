#!/usr/bin/python

from urllib.request import urlopen
import re
import pickle
import numpy as np
from fractions import Fraction
from math import sqrt

class r:
    def __init__(self):
        self.a1 = ''
        self.a2 = ''
        self.a3 = ''

def GetSchemes():

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
    
    
    scheme = {'P':P,'Q':Q,'L':L,'L*':Ls}

    return scheme
