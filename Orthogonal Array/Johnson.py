#!/usr/bin/python

import numpy as np
import math
from itertools import combinations

def Johnson_adj(n,k):
    combs = list(combinations([i for i in range(n)],k))
    bincombs = np.array([[1 if i in tup else 0 for i in range(n)] for tup in combs])
    A = np.dot(bincombs,np.transpose(bincombs))
    A2 = np.ones(A.shape)-A.copy()
    A[A!=1] = 0
    A2[A2!=1] = 0
    return A,A2