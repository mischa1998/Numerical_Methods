from tkinter import *
from sys import *
from PIL import Image, ImageTk

def get_value(entryWidget):
    value = entryWidget.get()
    return int(value)




def solve(var):
    f = open('data/pp.txt', 'w')
    N = get_value(nl)
    T = get_value(t)
    Moment = get_value(moment)
    f.write(str(N) + ' ')
    f.write(str(T) + ' ')
    ind = lwlistbox.curselection()
    f.write(str(ind[0]) + ' ')
    f.write(str(var.get()) + ' ')
    f.write(str(Moment))
    f.close()
    root.destroy()





root = Tk()
root.title("Численные методы ЛР2")
img = Image.open("data/LW1.jpg")
render = ImageTk.PhotoImage(img)
initil = Label(root, image=render)

nl = Entry(root)
nlab = Label(root, text="Число разбиений по Х")
t = Entry(root)
tlab = Label(root, text="Число разбиений по Т")
time = Label(root, text="Момент времени t")
moment = Entry(root)
var=IntVar()
rbutton1=Radiobutton(root,text='1 порядок точности 2 точки',variable=var,value=0)
rbutton2=Radiobutton(root,text='2 порядок точности 3 точки',variable=var,value=1)
rbutton3=Radiobutton(root,text='2 порядок точности 2 точки',variable=var,value=2)
q = Button(root, text="Solve", command= lambda: solve(var))
scheme = Label(root, text="Схема")

mode = ["Явная", "Неявная", "Кранк-Николсон"]
lwlistbox = Listbox(selectmode="SINGLE")
for modes in mode:
    lwlistbox.insert(END, modes)


initil.pack()
nl.pack()
nlab.pack()
t.pack()
tlab.pack()
moment.pack()
time.pack()
scheme.pack()
lwlistbox.pack()
rbutton1.pack()
rbutton2.pack()
rbutton3.pack()
q.pack()


root.mainloop()
