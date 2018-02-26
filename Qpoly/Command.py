### This file is purely for repetitive commands.
import pickle
import Qpoly
schemes = pickle.load(open('augschemes.p','rb'))
for cla in schemes:
    print('\n')
    for scheme in schemes[cla]:
        print(scheme)
        params = scheme.strip('<>').split(',')
        schemes[cla][scheme]['irrational'] = Qpoly.CheckRational(params[0],params[1],cla)
pickle.dump(schemes,open('augschemesi.p','wb'))