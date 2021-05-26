import numpy as np
from scipy.interpolate import CubicSpline as spline
from matplotlib import pyplot as plt

k=np.linspace(0,3,4)

X_k = np.sin(k * np.pi/6) 
Y_k = np.cos(k * np.pi/6) 

T_k = (np.pi * k/6)

#point 1
x_new = np.linspace(0,1,100)

y_h = spline(X_k, Y_k)(x_new)

#point 2
t_new = np.linspace(0,np.pi/2,100)

x_ht = spline(T_k, X_k)(t_new)
y_ht = spline(T_k, Y_k)(t_new)

#point 3
T_circle = np.array([k *np.pi/30 for k in range(0,61)])
X_circle = np.sin(T_circle)
Y_circle = np.cos(T_circle)

plt.figure().add_subplot(111).set_aspect('equal', adjustable='box')
plt.plot(x_new, y_h,'r--',label="interpolation en x")
plt.plot(X_circle, Y_circle,'g',label="vrai cercle")
plt.plot(x_ht,y_ht,'b--',label="interpolation en t")
plt.legend()
plt.show()