# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 20-21
# Problème 4
#
# Script de test
# Vincent Legat
#
# -------------------------------------------------------------------------
# 

import numpy as np
import time

start_time = time.time()

# ============================================================
# FONCTIONS A MODIFIER [begin]
#

def b(t,T,i,p):
  if p == 0:
    return (T[i] <= t)*(t < T[i+1])
  else:
    u  = 0.0 if T[i+p ]  == T[i]   else (t-T[i])/(T[i+p]- T[i]) * b(t,T,i,p-1)
    u += 0.0 if T[i+p+1] == T[i+1] else (T[i+p+1]-t)/(T[i+p+1]-T[i+1]) * b(t,T,i+1,p-1)
    return u

    
def bspline(X,Y,t):
  p = 3 

  T = np.arange(-3, len(X) + 1 + 3, 1)
  #T = [*(np.linspace(T[0], T[0],3)), *T, *(np.linspace(T[-1], T[-1],3))]
  n = len(T) - 1 
  #print(len(T))
  X = [*X, *X[0:3]]
  Y = [*Y, *Y[0:3]]


  B_spl = np.zeros((n-p,len(t)))

  for i in range(n-p):
    B_spl[i,:] = b(t, T, i, p)

  x = X @ B_spl 
  y = Y @ B_spl

  return x, y
   
#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
# -1- Approximation d'un rectangle :-)     
#

X = [0,1,3,3,3,5,6,3,3,3]
Y = [3,6,4,4,4,6,3,0,0,0]
t = np.linspace(0,len(X),len(X)*100 + 1)
      
x,y = bspline(X,Y,t)
print(round(time.time()-start_time, 3))

#
# -2- Un joli dessin :-)
#

import matplotlib.pyplot as plt
import matplotlib 
matplotlib.rcParams['toolbar'] = 'None'

p = 3
T = np.arange(-3, len(X) + 4, 1)
#T = [*(np.linspace(T[0], T[0],3)), *T, *(np.linspace(T[-1], T[-1],3))]
t_B = np.linspace(0, len(T)-1-p, len(T) * 100 + 1)

fig_Bspline = plt.figure("spline de base cubique")

for i in range(len(T)-p-1):
  B_spline = b(t_B, T, i, p)
  plt.plot(t_B, B_spline)
  

plt.plot(T, np.zeros(len(T)), '.r', markersize=10)
plt.axis(ymax=1) #, xmax=6, xmin=-2)
  


fig = plt.figure("Approximation avec des B-splines")
plt.plot(X,Y,'.r',markersize=10)
plt.plot([*X,X[0]],[*Y,Y[0]],'--r')
plt.plot(x,y,'-b')
plt.axis("equal")#; plt.axis("off")
plt.show()
