#
# Polynomial interpolation
# Vincent Legat - 2018
# Ecole Polytechnique de Louvain
#
 
from tkinter import *
from numpy   import polyfit,polyval
 
 
# =========================================================================
 
def mouse(event):
  global X,U
  print("New data : " + str(event.x) + "," + str(event.y))
  x1,y1 = (event.x - 5),(event.y - 5)
  x2,y2 = (event.x + 5),(event.y + 5)
  myCanvas.create_oval(x1,y1,x2,y2,fill='red')
  X.append(event.x)
  U.append(event.y)
 
# =========================================================================
 
def cleanCallBack():
  global X,U
  X,U=[],[]
  myCanvas.delete('all')
 
# =========================================================================
 
def computeCallBack():
  global X,U,myWindow
  myCanvas.delete('curve')
  if (len(X) >= 2):
    myPolynomial = polyfit(X,U,len(X)-1,rcond=2e-64)
    curve = []
    for x in range(0,myWindow.winfo_width(),1):
      uh = polyval(myPolynomial,x)
      curve.append((x,uh))
    myCanvas.create_line(curve,fill='blue',smooth=1,tag='curve')
 
# ============================= mainProgram ===============================
 
X,U=[],[]
 
myWindow = Tk()
myWindow.geometry('570x420')
myWindow.title('Interpolation')
myCanvas = Canvas(myWindow,width=570,height=420,bg='white')
myCanvas.bind('<Button-1>',mouse)
myCanvas.pack(fill='both',expand=True)
 
computeButton = Button(myWindow,text='Compute',command=computeCallBack)
computeButton.place(x=20,y=30,width=120,height=25)
cleanButton = Button(myWindow,text='Clean',command=cleanCallBack)
cleanButton.place(x=20+130,y=30,width=120,height=25)
exitButton = Button(myWindow,text='Exit',command=myWindow.quit)
exitButton.place(x=20+260,y=30,width=120,height=25)
 
myWindow.mainloop()
 
# =========================================================================