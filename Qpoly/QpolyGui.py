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
Label(window,text = "Number of classes:").grid(row=1,sticky=W)
d = Entry(window,width = 10)
d.grid(column=1,row=1)
d.insert(0,3)

### Function to grab updates
def update_d(*args):
    global update_in_progress
    if update_in_progress: return
    try:
        temp_f_float = temp_f_number.get()
    except ValueError:
        return
    new_temp_c = round((temp_f_float - 32) * 5 / 9, 2)
    update_in_progress = True
    temp_c_number.set(new_temp_c)
    update_in_progress = False



params = Combobox(window)
if type(d.get())==int:
    params['values'] = tuple(list(schemes[d.get()]))
else:
    params['values'] = tuple(list(schemes[3]))
params.current(0)
params.grid(column = 97,row = 1,columnspan = 2)

window.mainloop()


