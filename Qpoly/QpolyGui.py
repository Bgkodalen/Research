import re
import pickle
import numpy as np
from fractions import Fraction
from tkinter import *
from tkinter.ttk import *

### Load in all the information from Williford's tables.
schemes = pickle.load(open('schemes.p','rb'))

### The actual GUI
window = Tk()
window.title("Cometric Association schemes")
window.geometry('700x200')

lbl = Label(window, text="This project is intended to help analyze properties of Q-polynomial association schemes.",font = ("Arial Bold",10))
lbl.grid(row=0,column=0,columnspan=100,sticky=W+E+N+S)

### Close button
def quit():
    window.destroy()
btn = Button(window,text = "Close",command = quit)
btn.grid(column = 99, row = 99)

### Ability to enter class number
numclasses = DoubleVar()
Label(window,text = "Number of classes:").grid(row=1,sticky=W)
d = Entry(window,width = 10, textvariable = numclasses)
d.grid(column=1,row=1)
d.insert(0,3)

### Function to grab updates
update_in_progress = False
def update_d(*args):
    global update_in_progress
    if update_in_progress: return
    try:
        classes = numclasses.get()
    except ValueError:
        return
    update_in_progress = True
    numclasses.set(classes)
    update_in_progress = False

numclasses.trace("w",update_d)



params = Combobox(window)
if type(d.get())==int:
    params['values'] = tuple(list(schemes[numclasses.get()]))
else:
    params['values'] = tuple(list(schemes[3]))
params.current(0)
params.grid(column = 97,row = 1,columnspan = 2)

window.mainloop()


