#
# PYTHON for DUMMIES 20-21
# Problème 3
# Fonction interpolation avec des exposants négatifs
#
# Vincent Legat
#
# -------------------------------------------------------------------------
# 
 
from matplotlib import pyplot as plt
plt.rcParams['toolbar'] = 'None'
import numpy as np
from numpy.linalg import solve
import time
 
start_time = time.time()
 
# ============================================================
# FONCTIONS A MODIFIER [begin]
#
# Ici, on renvoie toujours bêtement : uh(x) = x - pi :-)
#
 
def interpolation(X,U,x):

  A = np.zeros((len(X), len(X)))
  n = len(X) //2

  for i, val_X in enumerate(X):
    for exp in range(-n, n+1):
      A[i, exp+n] = val_X**exp 
  
  coef = solve(A, U)
  uh = np.zeros(len(x))

  for k, x_val in enumerate(x):
    x_vect = np.zeros(2*n+1)
    for l in range(-n,n+1):
      x_vect[l+n] = x_val ** l 

    uh[k] = np.dot(x_vect, coef)

  return uh 


#
# FONCTIONS A MODIFIER [end]
# ============================================================
 
 
 
 
#
# -1- Test de la fonction interpolation
#     On considère un jeu des 3 fonctions u(x)
#
 
n = 2; m = 100
x = np.linspace(1,6,m)
X = np.linspace(1,6,2*n+1)
 
functions = [lambda x : 1 + 6*(x-1)*(x-6)*np.exp(-x), 
             lambda x : -2.7 + 1/(x-0.82),
             lambda x : np.sign(x-3.5)]
  
for u in functions:
 
  plt.figure()
  plt.plot(x,u(x),'-b',label='Fonction u')
  U = u(X)
  uh = np.polyval(np.polyfit(X,U,len(X)-1),x)
  plt.plot(x,uh,'-g',label='Interpolation polynomiale')
  uh = interpolation(X,U,x)
  plt.plot(x,uh,'-r',label='Mon interpolation :-)')
  plt.plot(X,U,'ob')
  plt.xlim((0.8,6.2)); plt.ylim((-3,3))
  plt.title('Interpolation avec exposants négatifs : 2n+1 = %d ' % len(X))
  plt.legend(loc='upper right')

print(time.time()-start_time)  
plt.show()
print(time.time()-start_time)