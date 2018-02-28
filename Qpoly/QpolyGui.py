import re
import pickle
import numpy as np
from fractions import Fraction
from tkinter import *
from tkinter.ttk import *
import Qpoly

### Load in all the information from Williford's tables.
schemes = pickle.load(open('augschemesi.p','rb'))
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
Radiobutton(root, text = "3-class primitive", variable = numclasses, value = 3).grid(column=1,row=1)
Radiobutton(root, text = "4-class bipartite", variable = numclasses, value = 4).grid(column=1,row=2)
Radiobutton(root, text = "5-class bipartite", variable = numclasses, value = 5).grid(column=1,row=3)

### The Examine Scheme window.
def examine():
    Details = Toplevel(root)
    scheme = re.findall(r'<[,\d\w]*>',params.get())[0]

    
    P = schemes[numclasses.get()][scheme]['P']
    Q = schemes[numclasses.get()][scheme]['Q']
    L = schemes[numclasses.get()][scheme]['L']
    Ls = schemes[numclasses.get()][scheme]['L*']
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
    exists = schemes[numclasses.get()][scheme]['exists']
    comments = schemes[numclasses.get()][scheme]['Comments']

    if exists == '-':
        color = 'red'
    elif exists == '?':
        color = 'yellow'
    else:
        color = 'green'
    Label(Details, text = (scheme+' '+exists),background=color).grid(row = 1,column = 1)
    Label(Details, text = comments).grid(row = 1,column = 3)
    
    
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
irr.grid(column = 1,row = 10)



ex = Button(root, text = "Examine Scheme", command = examine)
ex.grid(column = 99, row = 1)




### Shows the available parameters based on the Radiobutton input.
params = Combobox(root)
params.grid(column = 97,row = 1,columnspan = 2)
if irrat.get():
    params['values'] = tuple([scheme+schemes[numclasses.get()][scheme]['exists'] for scheme in schemes[numclasses.get()] if schemes[numclasses.get()][scheme]['irrational'] == 1])
else:
    params['values'] = tuple([scheme+schemes[numclasses.get()][scheme]['exists'] for scheme in schemes[numclasses.get()]])
if len(params['values']) == 0:
    params['values'] = ['None']
params.current(0)
temp = numclasses.get()
def update_comb(temp):
    if irrat.get()==1:
        params['values'] = tuple([scheme+schemes[numclasses.get()][scheme]['exists'] for scheme in schemes[numclasses.get()] if schemes[numclasses.get()][scheme]['irrational'] == 1])
    else:
        params['values'] = tuple([scheme+schemes[numclasses.get()][scheme]['exists'] for scheme in schemes[numclasses.get()]])
    if len(params['values']) == 0:
        params['values'] = ['None']
    if temp != numclasses.get():
        params.current(0)
    temp = numclasses.get()
    root.after(1000,lambda: update_comb(temp))

root.after(1000,lambda: update_comb(temp))
root.mainloop()


