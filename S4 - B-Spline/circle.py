#
# Splines cubiques
# Vincent Legat - 2018
# Ecole Polytechnique de Louvain
#

import matplotlib 
from matplotlib import pyplot as plt
matplotlib.rcParams['toolbar'] = 'None'

from numpy import *
from scipy.interpolate import CubicSpline as spline

#
# -1- Dessiner le cercle
#

T = pi * arange(0,13) / 6 
t = linspace(T[0],T[-1],100)
plt.plot(sin(T),cos(T),'ob')
plt.plot(sin(t),cos(t),':b')

#
# -1- Interpoler un quart de cercle
#

T = pi * arange(0,4) / 6
X = sin(T)
Y = cos(T)

x = linspace(X[0],X[-1],100)
yx = spline(X,Y)(x)

t = linspace(T[0],T[-1],100)
xt = spline(T,X)(t)
yt = spline(T,Y)(t)

plt.plot(x,yx,'-r')
plt.plot(xt,yt,'-b')
plt.plot(X,Y,'or')

plt.axis('off')
plt.axis('equal')


plt.show()