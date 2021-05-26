# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 20-21
# Problème 0
# Calcul du nombre d'or
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 


# ============================================================
# FONCTIONS A MODIFIER [begin]
#
# Ici, on renvoie toujours la meme liste quelque soit la valeur de phi :-(
#

def solveRatio(phi) :
  if phi == 0.5:
    return phi
  x1 = (phi)/((2*phi)-1)
  x2 = 1 - x1
  return [x1, x2]

#
# FONCTIONS A MODIFIER [end]
# ============================================================

# -------------------------------------------------------------------------
#
# -1- Test de la fonction solveRatio
#     Essayer comme argument phi = 0.5 :-)
#

import numpy as np

phi = (1 + np.sqrt(5.0)) / 2
x = solveRatio(phi)
print("Solution with phi = %.6f is given by : \n  " % phi,end='')
print(x)

#
# -2- Un petit design moderniste :-)
#     Typiquement, le genre de profil de fenêtre que Lecorbusier utilisait !
#

X = np.array([0,x[0],x[0],0,0]) + x[1]
Y = np.array([0,0,x[1],x[1],0])

import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['toolbar'] = 'None'
plt.rcParams['figure.facecolor'] = 'moccasin'
plt.figure('Golden ratio :-)')
plt.axis('equal')
plt.axis('off')
plt.text(0.5,-0.0625*phi,'$\phi = %.6f$' % phi) 
plt.plot([0,phi,phi,0,0,1,1],[0,0,1,1,0,0,1],'-b')
plt.plot(X,Y,'-r')
plt.plot(Y,X,'-r')
plt.show()
