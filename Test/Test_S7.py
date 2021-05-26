
from scipy.interpolate import CubicSpline 
import numpy as np

def b(t,T,i,p):
  if p == 0:
    return (T[i] <= t)*(t < T[i+1])
  else:
    u  = 0.0 if T[i+p ]  == T[i]   else (t-T[i])/(T[i+p]- T[i]) * b(t,T,i,p-1)
    u += 0.0 if T[i+p+1] == T[i+1] else (T[i+p+1]-t)/(T[i+p+1]-T[i+1]) * b(t,T,i+1,p-1)
    return u

    
def periodicBspline(X,Y,t,p): 

  T = np.arange(-p, len(X) + 1 + p, 1)
  n = len(T) - 1 
  #print(len(T))
  X = [*X, *X[0:p]]
  Y = [*Y, *Y[0:p]]


  B_spl = np.zeros((n-p,len(t)))

  for i in range(n-p):
    B_spl[i,:] = b(t, T, i, p)

  x = X @ B_spl 
  y = Y @ B_spl

  return x, y




X = [0,1,3,5,6]
Y = [1,2,3,4,7]
t = np.linspace(0,len(X),len(X)*100 + 1)
      
x,y = periodicBspline(X,Y,t,2)


#
# -2- Un joli dessin :-)
#

import matplotlib.pyplot as plt
import matplotlib 
matplotlib.rcParams['toolbar'] = 'None'  


fig = plt.figure("Approximation avec des B-splines")
plt.plot(X,Y,'.r',markersize=10)
plt.plot([*X,X[0]],[*Y,Y[0]],'--r')
plt.plot(x,y,'-b')
plt.axis("equal")#; plt.axis("off")
plt.show()
