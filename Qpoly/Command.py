### This file is purely for repetitive commands.
import time
import Qpoly
import pickle
schemes = pickle.load(open('schemes.p','rb'))
obs = []
f = open('obstructions.txt','w')
for d in range(3,6):
    f.write('%d-class: \n' % d)
    dschemes = schemes[d]
    for key in dschemes:
        Proj = Qpoly.Gegproj(dschemes[key]['L*'])
        if Proj.min()<-10**(-8):
            print(key+': \n')
            print(Proj.min())
            print('\n\n')
            obs.append(key)
            f.write('\t%s\n' % key)