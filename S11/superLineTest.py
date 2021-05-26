# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 20-21
# Problème 9
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

import numpy as np
import time

start_time = time.time()

#
# -------- PARTIE A MODIFIER -------------------------------------------
#
#
#
# -1- Droite au sens de moindres carrés usuels
# 
  
def superLineInitialGuest(X,U):

#
# A MODIFIER ..... [begin]
#
  A = np.array([X, np.ones(len(X))])
  a = A@A.T
  b = A@U

  alpha = np.linalg.solve(a,b)

  
#
# A MODIFIER ..... [end]
# 
    
  return alpha
  
#
# -------------------------------------------------------------------------
#
# -2- Iteration de Newton-Raphson pour obtenir
#     les coefficients d'une droite minimisant la "vraie distance" entre 
#     les données et la droite  
# 

def superLineIterate(X,U,alpha):

#
# A MODIFIER ..... [begin]
#
  a,b = alpha
  #F_ab = sum(U -a*X -b)**2 / (1+a**2)

  aa = (1+a**2)
  C = (U -a*X-b)
  C2 = C**2

  dfa = (-2*sum(C*X)*aa - 2*a*sum(C2))

  df = np.zeros(2)
  ddf = np.zeros((2,2))

  df[0] = dfa / (aa**2)
  df[1] = -2*sum(C)/aa

  ddf[1,0] = (2*sum(X)*aa + 4*a*sum(C)) / aa**2
  ddf[0,1] = ddf[1,0]
  ddf[0,0] = ( ((2*sum(X*X)*aa + 2*sum(C*X)*a*2 - 2*sum(C2) + 4*a*sum(C*X))*aa +  8*a*a*sum(C2))) / aa**3
  ddf[1,1] = 2*len(X) / aa

  #dfa = sum(-2*U*X + 2*a**2*U*X - 2*a*X**2 + 2*b*X - 2*a**2*b*X - 2*a**2*b*X - 2*a*U**2 - 2*a*b**2 + 4*a*b*U)/ aa**2
  #dfb = -sum(2*U - 2*a*X - 2*b)/aa
  #df = -np.array([dfa, dfb])
  #
  ##Mtn on doit fair la matrice hessienne
  #ddf = np.zeros((2,2))
  #
  #ddf[0,0] = sum(12*a*U*X - 4*a**3*U*X + 2*X**2 - 6*a**2*X**2 -12*a*b*X + 4*a**3*b*X-2*U**2+6*a**2*U**2 -2*b**2 + 6*a**2*b**2 + 4*b*U - 12*a**2*b*U)/ aa**3
  #ddf[0,1] = sum(2*X - 2*a**2*X - 4*a*b + 4*a*U)/aa**2
  #ddf[1,0] = sum(-2*X + 2*a**2*X - 4*a*b + 4*a*U)/aa**2
  #ddf[1,1] = 2/ aa
  #
  dalpha = -np.linalg.solve(ddf, df)
  
#
# A MODIFIER ..... [end]
# 

  return dalpha
  
#
# -------------------------------------------------------------------------
#
#  
# -1- Schéma de Newton-Raphson pour calculer les coefficients 
#     d'une droite minimisant la "vraie distance" entre 
#     les données et la droite
#  
#     superLineInitialGuest = fournit le candidat initial
#     superLineIterate = calcul l'incrément pour une itération
#

def superLine(X,U):
  alpha  = superLineInitialGuest(X,U)
  nmax   = 30
  iter   = 0
  tol    = 1e-6
  dalpha = 1
  error = np.zeros(nmax)
  while (iter < nmax) and (np.linalg.norm(dalpha) > tol):
    dalpha = superLineIterate(X,U,alpha)
    alpha = alpha + dalpha
    error[iter] = np.linalg.norm(dalpha)
    iter = iter + 1
    print('   Iteration %i : %14.7e (a =%14.7e b =%14.7e)' % (iter,np.linalg.norm(dalpha),alpha[0],alpha[1]))


#
#     Estimation numérique du taux de convergence de la méthode
#     en comparant la diminution successive des erreurs :-)
#     10-2,10-4,10-8 => erreur quadratique !
#

  if (iter > 1) :
    rate = np.mean(np.log(error[1:iter])/np.log(error[:iter-1]))
    print('   Observed rate of convergence : %.2f ' % rate)
    print('   Theoretical rate             : %.2f ' % 2.0)
    
  return alpha

#
# -2- Un ensemble de données
#

X = np.array([0.2,-0.7,-0.8, 0.4, 0  ,-0.5, 0.6])
U = np.array([1.0, 2.0, 2.0, 1.4, 1.0, 0.9, 1  ])
m = len(X)


#
# -3- Calcul des deux approximations
#     Et un joli dessin avec la représentation des résidus 
#     dans les deux cas :-)
#
#     Observer le calcul pour obtenir l'intersection (XSuper,USuper) d'une droite
#     et de la perpendiculaire à cette droite passant par le point (X,U)
#     Un joli exercice de géométrie pour l'examen d'entrée :-)
#


pu = superLineInitialGuest(X,U)
x = np.linspace(-1,1,100)
uLine = pu[0] * x + pu[1]
ULine = pu[0] * X + pu[1]
print(" === Usual least squares linear regression         : u(x) = %.3f x + %.3f " % (pu[0],pu[1]))

pu = superLine(X,U)
print("Time to run = {:e}s".format(time.time()-start_time) )
uSuper = pu[0] * x + pu[1]
t = -(pu[0] * X - U + pu[1])/(1 + pu[0]*pu[0])
XSuper = X + pu[0] * t
USuper = U - t
print(" === Minimizing the distance between line and data : u(x) = %.3f x + %.3f " % (pu[0],pu[1]))

from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['toolbar'] = 'None'
#plt.rcParams['figure.facecolor'] = 'lavender'
#plt.rcParams['axes.facecolor'] = 'lavender'
fig = plt.figure("Super line approximation")

plt.plot(x,uLine ,'-b',X,ULine,'.b',linewidth=0.5)
plt.plot(x,uSuper,'-g',XSuper,USuper,'.g',linewidth=0.5)
for i in range(m) :
  plt.plot([X[i],X[i]]     ,[U[i],ULine[i]] ,'-b',linewidth=0.5)
  plt.plot([X[i],XSuper[i]],[U[i],USuper[i]],'-g',linewidth=0.5)
plt.plot(X,U,'or')

plt.axis('equal')

plt.show()




