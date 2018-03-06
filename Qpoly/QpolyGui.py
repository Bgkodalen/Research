import re
import pickle
import numpy as np
from fractions import Fraction
from tkinter import *
from tkinter.ttk import *
import math
import Qpoly

### Load in all the information from Williford's tables.
schemes = pickle.load(open('Schemes.tot','rb'))
### There is an annoying issue with the keys here where the lettered schemes dont have a trailing >.
# Below is a temporary patch to this so that I don't have to redo the data compilation.
for key in schemes:
    for scheme in schemes[key]:
        if scheme[-1]!= '>':
            schemes[key][scheme+'>'] = schemes[key].pop(scheme)
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
Radiobutton(root, text = "3-class primitive", variable = numclasses, value = 3).grid(column=1,row=1,sticky=W)
Radiobutton(root, text = "3-class antipodal", variable = numclasses, value = 0).grid(column=1,row=2,sticky=W)
Radiobutton(root, text = "4-class bipartite", variable = numclasses, value = 4).grid(column=1,row=3,sticky=W)
Radiobutton(root, text = "5-class bipartite", variable = numclasses, value = 5).grid(column=1,row=4,sticky=W)

### The Examine Scheme window.
def examine():
    Details = Toplevel(root)
    scheme = re.findall(r'<[,\d\w;]*>',params.get())[0]
    #p = scheme.strip('<>').split(',')
    
    P = schemes[numclasses.get()][scheme]['P']
    #Q = schemes[numclasses.get()][scheme]['Q']
    Q = Qpoly.Qm(P)
    #L = schemes[numclasses.get()][scheme]['L']
    L = Qpoly.Lm(P)
    #Ls = schemes[numclasses.get()][scheme]['L*']
    Ls = Qpoly.Lsm(P)
    Geg = Qpoly.Gegproj(Ls,10,0,1)

    fp = Frame(Details)
    fp.grid(row = 2, column = 2,sticky = W)
    [r,c] = Matrixfrmt(P,'P',fp,2,2)
    [r,c] = Matrixfrmt(Q,'Q',fp,2,c+1)
    [r,t] = Matrixfrmt(Geg,'G',fp,2,c+1,1)

    fl = Frame(Details)
    fl.grid(row = 4, column = 2,sticky = W)
    [r,t] = Matrixfrmt(L,'L',fl,r+1,2)

    fls = Frame(Details)
    fls.grid(row = 5, column = 2,sticky = W)
    [r,t] = Matrixfrmt(Ls,'L*',fls,r+1,2)

    ### Extra information
    if 'exists' in schemes[numclasses.get()][scheme]:
        exists = schemes[numclasses.get()][scheme]['exists']
    else:
        exists = '?'
    if 'Comments' in schemes[numclasses.get()][scheme]:
        comments = schemes[numclasses.get()][scheme]['Comments']
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
    Label(Details, text = comments).grid(row = 1,column = 3)
    
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

    exclose = Button(Details,text = "Close",command = lambda: quit(Details))
    exclose.grid(column = 99, row = 99)



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

irrat = IntVar()
irr = Checkbutton(root,text="Irrational schemes", variable = irrat)
irr.grid(column = 1,row = 10,sticky=W)

geg = IntVar()
ge = Checkbutton(root,text="Spherical bound", variable = geg)
ge.grid(column = 1,row = 11,sticky=W)


SD = IntVar()
SDcheck = Checkbutton(root,text="Spherical design", variable = SD)
SDcheck.grid(column = 2,row = 1,sticky=W)
numdesign = Combobox(root,width = 2)
numdesign.grid(column = 3,row = 1,sticky = W)
numdesign['values'] = (3,4,5,6,7,8,9)
numdesign.current(0)

equi = IntVar()
equicheck = Checkbutton(root, text = "Equiangular lines", variable = equi)
equicheck.grid(column = 2, row = 10, sticky = W)


ex = Button(root, text = "Examine Scheme", command = examine)
ex.grid(column = 99, row = 1)


def parameterlist():
    schemelist = [scheme for scheme in schemes[numclasses.get()]]
    tol=10**(-14)
    if irrat.get():
        schemelist = [scheme for scheme in schemelist if schemes[numclasses.get()][scheme]['irrational'] == 1]
    if geg.get():
        schemelist = [scheme for scheme in schemelist if (Qpoly.Gegproj(Qpoly.Lsm(schemes[numclasses.get()][scheme]['P']))).min()<-tol]
    if SD.get():
        schemelist = [scheme for scheme in schemelist if sum(np.absolute(Qpoly.Gegproj(Qpoly.Lsm(schemes[numclasses.get()][scheme]['P']))[0,1:(int(numdesign.get())+1)]))<tol]
    if equi.get():
        schemelist = [scheme for scheme in schemelist if sum(Qpoly.Qm(schemes[numclasses.get()][scheme]['P'])[0,:])-Qpoly.Qm(schemes[numclasses.get()][scheme]['P']).max()<140]
    if len(schemelist) == 0:
        schemelist = ['None']
    else:
        schemelist = [scheme+schemes[numclasses.get()][scheme]['exists'] for scheme in schemelist]
    return tuple(schemelist)


### Shows the available parameters based on the Radiobutton input.
params = Combobox(root)
params.grid(column = 96,row = 1,columnspan=3)
params['values'] = parameterlist()

params.current(0)
temp = params['values'][0]
def update_comb(temp):
    params['values'] = parameterlist()
    if temp != params['values'][0]:
        params.current(0)
    temp = params['values'][0]
    if params['values'][0]!='None':
        Label(root, text = str(len(params['values']))+' Schemes').grid(row = 2,column = 99)
    else:
        Label(root, text = '0 Schemes').grid(row = 2,column = 99)
    root.after(1000,lambda: update_comb(temp))

root.after(1000,lambda: update_comb(temp))
root.mainloop()


