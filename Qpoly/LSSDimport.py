#!/usr/bin/python

### This file imports data from Williford's tables.

import xlrd
import pickle
import numpy as np
from math import sqrt


workbook = xlrd.open_workbook("C:/Users/bgkodalen/Desktop/Research2/LinkedSystem/Examples_1035.xlsx")
worksheet = workbook.sheet_by_index(0)
### params = [[v,k,lambda,w]]
params = [[worksheet.cell(i,j).value for j in [1,2,3,14,15,16]] for i in range(1,38)]

schemes = pickle.load(open("augschemesi.p",'rb'))
schemes[0] = {}

for paramset in params:
    v = paramset[0]
    k = paramset[1]
    l = paramset[2]
    w = paramset[3]
    comm = paramset[4]
    tpe = paramset[5]
    s = sqrt(k-l)
    P = np.array([[1, k*(w-1), v-1, (v-k)*(w-1)],[1, s*(w-1), -1, -s*(w-1)],[1, -s, -1, s],[1,-k,v-1,k-v]])
    schemes[0]['<%0.f,%0.f,%0.f;%0.f>'%(v,k,l,w)] = {'P':P,'exists':tpe,'comments':comm,'irrational':0}

pickle.dump(schemes,open('Schemes.tot','wb'))

