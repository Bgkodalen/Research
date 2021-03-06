import re
import pickle
import numpy as np
from fractions import Fraction
from tkinter import *
from tkinter.ttk import *
import math
import Qpoly

### Load in all the current information.
schemes = pickle.load(open('alldata2.p','rb'))
#schemes = pickle.load(open('alldata.p','rb'))
#schemes = pickle.load(open('bigfile','rb'))

### There is an annoying issue with the keys here where the lettered schemes dont have a trailing >.
# Below is a temporary patch to this so that I don't have to redo the data compilation.
for d in schemes:
    for prim in schemes[d]:
        for scheme in schemes[d][prim]:
            if scheme[-1]!= '>':
                schemes[d][prim][scheme+'>'] = schemes[d][prim].pop(scheme)
### The actual GUI
root = Tk()
root.title("Cometric Association schemes")
root.geometry('700x200')

lbl = Label(root, text="This project is intended to help analyze properties of Q-polynomial association schemes.",font = ("Arial Bold",10))
lbl.grid(row=0,column=0,columnspan=100,sticky=W+E+N+S)

### Close button
def quit(window):
    window.destroy()
btn = Button(root,text = "Close",command = lambda: quit(root))
btn.grid(column = 99, row = 99)

### Ability to enter class number --- Change this to radio buttons?
numclasses = IntVar()
numclasses.set(3)
Radiobutton(root, text = "3", variable = numclasses, value = 3).grid(column=1,row=1,sticky=W)
Radiobutton(root, text = "4", variable = numclasses, value = 4).grid(column=1,row=2,sticky=W)
Radiobutton(root, text = "5", variable = numclasses, value = 5).grid(column=1,row=3,sticky=W)

imprim = StringVar()
imprim.set('primitive')
Radiobutton(root, text = "primitive", variable = imprim, value = 'primitive').grid(column=2,row=1,sticky=W)
Radiobutton(root, text = "antipodal", variable = imprim, value = 'antipodal').grid(column=2,row=2,sticky=W)
Radiobutton(root, text = "bipartite", variable = imprim, value = 'bipartite').grid(column=2,row=3,sticky=W)

def Matrixfrmt(mat,name,window,r,c,string=0):
    dim = mat.shape
    if len(name)>0:
        Label(window,text = name+'=').grid(row=r+int(dim[0]/2),column=c)
    if len(dim)>2:
        for i in range(dim[0]):
            [t,c] = Matrixfrmt(mat[i],'',window,r,c,string)
        return [t,c]
    else:
        rowpos = r
        for i in range(dim[0]):
            colpos = c+2
            for j in range(dim[1]):
                if string == 1:
                    Label(window,text = mat[i,j]).grid(row=rowpos,column = colpos)
                elif (type(mat[i,j]) == bytes) or (type(mat[i,j]) == str):
                    Label(window,text = mat[i,j]).grid(row=rowpos,column=colpos)
                elif round(mat[i,j],2) == int(mat[i,j]):
                    Label(window,text = int(mat[i,j])).grid(row=rowpos,column=colpos)
                else:
                    Label(window,text = format(float(mat[i,j]),'.2f')).grid(row=rowpos,column=colpos)
                colpos+=1         
            rowpos+=1
        
        for b in range(dim[0]):
            Label(window,text = '|').grid(row = r+b,column = c+1)
            Label(window,text = '|').grid(row = r+b,column = colpos)
        Label(window,text = ' ').grid(row = rowpos,column = colpos+1)
        Label(window,text = ' ').grid(row = rowpos,column = colpos+2)
        return [rowpos,colpos+2]


### Various filters on the schemes
irrat = IntVar()
irr = Checkbutton(root,text="Irrational schemes", variable = irrat)
irr.grid(column = 1,row = 4,sticky=W)

geg = IntVar()
ge = Checkbutton(root,text="Spherical bound", variable = geg)
ge.grid(column = 1,row = 5,sticky=W)

SD = IntVar()
SDcheck = Checkbutton(root,text="Spherical design", variable = SD)
SDcheck.grid(column = 2,row = 4,sticky=W)
numdesign = Combobox(root,width = 2)
numdesign.grid(column = 3,row = 4,sticky = W)
numdesign['values'] = (3,4,5,6,7,8,9)
numdesign.current(0)

equi = IntVar()
equicheck = Checkbutton(root, text = "Equiangular lines", variable = equi)
equicheck.grid(column = 2, row = 5)

def sortedlist(lst):
    x = [[int(string.split(',')[0].strip('<')),int(string.split(',')[1].split('>')[0].strip('abcdefghijklmnopqrstuvwxyz')),i] for i,string in enumerate(lst) if (string.split(',')[0].strip('<')).isdigit()]
    x.sort()
    l = [lst[indx[2]] for indx in x]
    return l

def parameterlist():
    selectedschemes = schemes[numclasses.get()][imprim.get()]
    schemelist = [scheme for scheme in selectedschemes]
    tol=10**(-14)
    if irrat.get():
        schemelist = [scheme for scheme in schemelist if selectedschemes[scheme]['irrational'] == 1]
    if geg.get():
        schemelist = [scheme for scheme in schemelist if (Qpoly.Gegproj(Qpoly.Lsm(selectedschemes[scheme]['P']))).min()<-tol]
    if SD.get():
        schemelist = [scheme for scheme in schemelist if sum(np.absolute(Qpoly.Gegproj(Qpoly.Lsm(selectedschemes[scheme]['P']))[0,1:(int(numdesign.get())+1)]))<tol]
    if equi.get():
        newlist = []
        tableexamples = []
        optim = []
        nearop = []
        for scheme in schemelist:
            use = 0
            [a,success,examples] = Qpoly.Equiangular(Qpoly.Qm(selectedschemes[scheme]['P']))
            if success:
                for example in examples:
                    if example[1]<example[2]**2 and example[0] > round(example[1]*(1-1/(example[2]**2))/(1-example[1]/(example[2]**2)))-2:
                        use=1
                        tableexamples.append([scheme, example[0],example[1],example[2],round(example[1]*(1-1/(example[2]**2))/(1-example[1]/(example[2]**2))-example[0])])
                        if abs(round(example[1]*(1-1/(example[2]**2))/(1-example[1]/(example[2]**2))-example[0]))<1:
                            optim.append([scheme, example[0],example[1],example[2]])
                        else:
                            nearop.append([scheme, example[0],example[1],example[2]])
            if use:
                newlist.append(scheme)
        schemelist = [scheme for scheme in newlist]
        f = open('schemetableop.txt','w')
        g = open('schemetablenearop.txt','w')
        halfop = len(optim)//2
        for i in range(halfop):
            f.write('$\\left<%s\\right>$ & %0.f & %0.f & %0.f & \qquad &$\\left<%s\\right>$ & %0.f & %0.f & %0.f \\\\\n' % (optim[i][0][1:-1],optim[i][1],optim[i][2],optim[i][3],optim[i+halfop][0][1:-1],optim[i+halfop][1],optim[i+halfop][2],optim[i+halfop][3]))
        halfnear = len(nearop)//2
        for i in range(halfnear):
            g.write('$\\left<%s\\right>$ & %0.f & %0.f & %0.f & \qquad &$\\left<%s\\right>$ & %0.f & %0.f & %0.f \\\\\n' % (nearop[i][0][1:-1],nearop[i][1],nearop[i][2],nearop[i][3],nearop[i+halfnear][0][1:-1],nearop[i+halfnear][1],nearop[i+halfnear][2],nearop[i+halfnear][3]))
        if len(optim) > halfop * 2:
            f.write('$\\left<%s\\right>$ & %0.f & %0.f & %0.f \\\\\n' % (optim[-1][0][1:-1],optim[-1][1],optim[-1][2],optim[-1][3]))
        if len(nearop)> halfnear *2:
            g.write('$\\left<%s\\right>$ & %0.f & %0.f & %0.f \\\\\n' % (nearop[-1][0][1:-1],nearop[-1][1],nearop[-1][2],nearop[-1][3]))
        f.close()
        g.close()
    if numclasses.get() == 4:
        schemelist = [scheme+' '+str(-selectedschemes[scheme]['P'][4,1]) for scheme in schemelist]
    #    schemelist = [scheme for scheme in schemelist if sum(Qpoly.Qm(selectedschemes[scheme]['P'])[0,:])-Qpoly.Qm(selectedschemes[scheme]['P']).max()<140]
    if len(schemelist) == 0:
        schemelist = ['None']
    else:
        #+selectedschemes[scheme]['exists']
        if numclasses.get()==4:
            schemelist = sorted(schemelist, key=lambda x: int(re.search(r'\d+$',x).group()))
        else:
            schemelist = sortedlist([scheme for scheme in schemelist])
    
    #def sortThird(val):
    #    return val[2]
    #if numclasses.get()==4 and imprim.get() == 'bipartite':
    #    schemelist = [[selectedschemes[scheme]['P'][0,1],selectedschemes[scheme]['P'][2,1]] for scheme in schemelist if selectedschemes[scheme]['P'][3,1] == -3]
    #    schemelist.sort()
    return tuple(schemelist)


### Shows the available parameters based on the Radiobutton input.
params = Combobox(root)
params.grid(column = 90,row = 1,columnspan=9)
params['values'] = parameterlist()
params['height'] = 50

params.current(0)
temp = params['values'][0]
def update_comb(temp):
    params['values'] = parameterlist()
    if temp != params['values'][0]:
        params.current(0)
    temp = params['values'][0]
    if params['values'][0]!='None':
        num['text'] = str(len(params['values']))+' Schemes'
    else:
        num['text'] = '0 Schemes'
    root.after(1000,lambda: update_comb(temp))

### The Examine Scheme window.
def examine():
    Details = Toplevel(root)
    scheme = re.findall(r'<[,\d\w;]*>',params.get())[0]
    
    Data = schemes[numclasses.get()][imprim.get()][scheme]

    ### Load in the P matrix and calculate the rest
    P = Data['P']
    Q = Qpoly.Qm(P)
    L = Qpoly.Lm(P)
    Ls = Qpoly.Lsm(P)
    Geg = Qpoly.Gegproj(Ls,10,0,1)



    ### Display the various information
    fp = Frame(Details)
    fp.grid(row = 2, column = 2,sticky = W)
    [r,c] = Matrixfrmt(P,'P',fp,2,2)
    [r,c] = Matrixfrmt(Q,'Q',fp,2,c+1)
    [r,t] = Matrixfrmt(Geg,'G',fp,2,c+1,1)
    if imprim.get()=="bipartite" and numclasses.get() == 4:
        #Gegsrg = Qpoly.Gegproj(Qpoly.Lsm(Data['psrg']),10,0,1)
        #[r,t] = Matrixfrmt(Gegsrg,'Gsrg',fp,2,t+1,1)
        Label(Details, text = ("(k,r,s) = ("+str(P[0,1])+","+str(P[2,1])+","+str(P[3,3])+")")).grid(row = 100,column = 1)

    fl = Frame(Details)
    fl.grid(row = 4, column = 2,sticky = W)
    [r,t] = Matrixfrmt(L,'L',fl,r+1,2)

    fls = Frame(Details)
    fls.grid(row = 5, column = 2,sticky = W)
    [r,t] = Matrixfrmt(Ls,'L*',fls,r+1,2)

    ### Extra information
    if 'exists' in Data:
        exists = Data['exists']
    else:
        exists = '?'
    if 'Comments' in Data:
        comments = Data['Comments']
    else:
        comments = ''

    if exists == '-':
        color = 'red'
    elif exists == '?':
        color = 'yellow'
    elif exists == 'O':
        color = 'blue'
    elif exists == 'P':
        color = 'yellow'
    else:
        color = 'green'
    Label(Details, text = (scheme+' '+exists),background=color).grid(row = 1,column = 1)
    Label(Details, text = comments).grid(row = 1,column = 2)
    
    ### Temporary Labels
    if Ls[1,1,1]>0:
        m = Ls[1,0,1]
        a = Ls[1,1,1]**2
        b = Ls[1,2,2]/Ls[1,1,1]
        c = Ls[1,1,2]*Ls[1,2,1]
        d = 4*m*(2*m-3)/(m+6)
        Label(Details, text = '%0.2f + (2+%0.2f)*%0.2f > %0.2f' % (a,b,c,d)).grid(row = 99,column = 1,columnspan = 5)
        Label(Details, text = '%0.2f > %0.2f' % (a+(2+b)*c,d)).grid(row = 100,column = 1,columnspan = 5)
    
    Label(Details, text = Qpoly.Equiangular(Q)[0]).grid(row = 6, column = 2)

    l = Q[:,1]/Q[0,1]
    s = l[1]
    X = sum(Q[0,:])

    Label(Details, text= "Max degree "+str(((math.log((1+s**2)*X*(X-2)))/(-2*math.log(s))))).grid(row=7,column = 4)

    exclose = Button(Details,text = "Close",command = lambda: quit(Details))
    exclose.grid(column = 99, row = 99)

ex = Button(root, text = "Examine Scheme", command = examine)
ex.grid(column = 99, row = 1)

num = Label(root,text = '')
num.grid(row = 2,column = 99)



root.after(1000,lambda: update_comb(temp))
root.mainloop()


