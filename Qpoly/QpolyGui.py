import re
import pickle
import numpy as np
from fractions import Fraction
from tkinter import *
from tkinter.ttk import *
import Qpoly

### Load in all the information from Williford's tables.
schemes = pickle.load(open('augschemes.p','rb'))
testmat = np.array([[1,2,3,4],[1,2,3,5], [0,1,1,2], [1,1,1,1]])
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
    Label(Details, text = params.get()).grid(row = 1,column = 1)
    
    P = schemes[numclasses.get()][params.get()]['P']
    Q = schemes[numclasses.get()][params.get()]['Q']
    L = schemes[numclasses.get()][params.get()]['L']
    Ls = schemes[numclasses.get()][params.get()]['L*']
    Geg = Qpoly.Gegproj(Ls,10,0,1)

    fp = Frame(Details)
    fp.grid(row = 2, column = 2,sticky = W)
    [r,c] = Matrixfrmt(P,'P',fp,2,2)
    [r,c] = Matrixfrmt(Q,'Q',fp,2,c+1)
    [r,t] = Matrixfrmt(Geg,'G',fp,2,c+1)

    fl = Frame(Details)
    fl.grid(row = 4, column = 2,sticky = W)
    [r,t] = Matrixfrmt(L,'L',fl,r+1,2)

    fls = Frame(Details)
    fls.grid(row = 5, column = 2,sticky = W)
    [r,t] = Matrixfrmt(Ls,'L*',fls,r+1,2)

    
    
    exclose = Button(Details,text = "Close",command = lambda: quit(Details))
    exclose.grid(column = 99, row = 99)


def Matrixfrmt(mat,name,window,r,c):
    dim = mat.shape
    if len(name)>0:
        Label(window,text = name+'=').grid(row=r+int(dim[0]/2),column=c)
    if len(dim)>2:
        for i in range(dim[0]):
            [t,c] = Matrixfrmt(mat[i],'',window,r,c)
        return [t,c]
    else:
        rowpos = r
        for i in range(dim[0]):
            colpos = c+2
            for j in range(dim[1]):
                if (type(mat[i,j]) == bytes) or (type(mat[i,j]) == str):
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



ex = Button(root, text = "Examine Scheme", command = examine)
ex.grid(column = 99, row = 1)




### Shows the available parameters based on the Radiobutton input.
params = Combobox(root)
params.grid(column = 97,row = 1,columnspan = 2)
params['values'] = tuple(list(schemes[numclasses.get()]))
params.current(0)
temp = numclasses.get()
def update_comb(temp):
    params['values'] = tuple(list(schemes[numclasses.get()]))
    if temp != numclasses.get():
        params.current(0)
    temp = numclasses.get()
    root.after(1000,lambda: update_comb(temp))

root.after(1000,lambda: update_comb(temp))
root.mainloop()


