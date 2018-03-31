#!/usr/bin/python

import numpy as np
import sys
import pickle
import Qpoly
import math

tol = 10**-8

def isint(x):
    if (np.round(x)-x)**2<tol:
        return 1
    else:
        return 0


def absolute(Ls):
    f = [Ls[i,0,i] for i in range(5)]
    valid = 1
    for i in range(5):
        for j in range(5):
            lower = sum([f[k] for k in range(5) if Ls[i,k,j]>tol])
            if i==j:
                upper = 1/2*f[i]*(f[i]+1)
            else:
                upper = f[i]*f[j]
            if lower>upper+tol:
                valid = 0
    return valid



if __name__ == '__main__':

    if len(sys.argv) == 1:
        print('need max k parameter')
        sys.exit(-1)
    
    schemes = pickle.load(open("alldata2.p",'rb'))
    maxk = int(sys.argv[1])
    schemes[4]['bipartite'] = {}

    for k in range(maxk):
        for n in range(2,int(math.sqrt(k))):
            s = -n**2
            for r in range(1,int(k/(-s))):
                valid = 1
                mu = k+r*s
                f = (s+1)*k*(k-s)/(mu*(s-r))
                v = (k-r)*(k-s)/mu
                g = v-1-f
                if isint(v) and isint(f) and g>0:
                    P = np.array([[1,k,2*(v-1-k),k,1],[1,k/n,0,-k/n,-1],[1,r,-2*(1+r),r,1],[1,-n,0,n,-1],[1,s,-2*(s+1),s,1]])
                    try:
                        i = np.linalg.inv(P)
                    except:
                        valid = 0
                    if valid == 1:
                        L = Qpoly.Lm(P)
                        Ls = Qpoly.Lsm(P)
                        Q = Qpoly.Qm(P)
                        X = sum(P[0,:])
                        m = Q[0,1]
                        valid = 1
                        for ii in range(4):
                            for jj in range(4):
                                for kk in range(4):
                                    valid = valid*(Ls[ii,jj,kk] > -tol)*isint(L[ii,jj,kk])*(L[ii,jj,kk]>-tol)
                        if valid == 1 and absolute(Ls):
                            name = '<%0.f,%0.f>'%(-s,k)
                            if Ls[4,4,2]>tol:
                                print('<%0.f,%0.f>'%(X,k))
                                schemes[4]['bipartite'][name] = {'P':P,'comments':''}
                                if (Qpoly.Gegproj(Ls)).min()<-tol:
                                    schemes[4]['bipartite'][name]['exists'] = '-'
                                else:
                                    schemes[4]['bipartite'][name]['exists'] = '?'
                            else:
                                schemes[4]['B&A'][name] = {'P':P,'exists':'?','comments':''}
                    

    pickle.dump(schemes,open('newdata.p','wb'))