#Maxime Delacroix

import numpy as np

def cubicSolve(U, h):
    """
    Résoudre le système pour avoir les coef cubique pour chaque interval de pas h
    Donne un liste de longeur n+1 avec des listes des coef de l'equation cubique
    """
    A = []
    b = []
    for i in range(len(U)):
        x = np.linspace(0,0, len(U))
        if i == 0:
            x[0] = (1*(h**2))/6
            A.append(x)
            b.append(0)
        elif i == (len(U)-1):
            x[-1] = (1*(h**2))/6
            A.append(x)
            b.append(0)
        else:
            x[i-1] = (1*(h**2))/6
            x[i] = (4*(h**2))/6
            x[i+1] = (1*(h**2))/6
            A.append(x)
            b.append(U[i-1] - (2*U[i]) + U[i+1])
    return np.linalg.solve(A, b)
    


def cubicSpline(X, h, U, ddU, x):
    """
    Donne une liste Oord de longueur x qui est l'interpretation par spline cubique
    des ces points
    """

    Oord = []
    for x_elem in x:                
        if x_elem >= X[-1]:             
            i = len(U)-1                
            Oord.append(((ddU[i-1]/(6*h))*((X[i]-x_elem)**3))               
            +((ddU[i]/(6*h))*((x_elem-X[i-1])**3)               
            +(((U[i-1]/h)-((h*ddU[i-1])/6))*(X[i]-x_elem))              
            +(((U[i]/h)-((ddU[i]*h)/6))*(x_elem-X[i-1]))))                  
        elif x_elem <= X[0]:                
            i = 1               
            Oord.append(((ddU[i-1]/(6*h))*((X[i]-x_elem)**3))               
            +((ddU[i]/(6*h))*((x_elem-X[i-1])**3)               
            +(((U[i-1]/h)-((h*ddU[i-1])/6))*(X[i]-x_elem))              
            +(((U[i]/h)-((ddU[i]*h)/6))*(x_elem-X[i-1]))))                  
        else:
            for j in range(1,len(X)):
                if x_elem >= X[j-1] and x_elem < X[j]:
                    i = j
                    Oord.append(((ddU[i-1]/(6*h))*((X[i]-x_elem)**3)) 
                    +((ddU[i]/(6*h))*((x_elem-X[i-1])**3) 
                    +(((U[i-1]/h)-((h*ddU[i-1])/6))*(X[i]-x_elem)) 
                    +(((U[i]/h)-((ddU[i]*h)/6))*(x_elem-X[i-1]))))                       

    return Oord


#
# -1- Test du splines cubiques
#
u = lambda x : X**2 + 1
Xstart = -5.0
Xend = 5.0
n = 10
X,h = np.linspace(Xstart,Xend,n+1,retstep=True)
#
# -2- Calcul des splines cubiques sur un intervalle quelconque
#
Xstart = -4.2
Xend = 9.7
U = u(X)
ddU = cubicSolve(U,h)
x = np.linspace(-4.2,9.7,200)
uh = cubicSpline(X,h,U,ddU,x)
print("==== Computing the cubic splines curve")
print("ddU = ",end=''); print(ddU)


#
# -3- Et zou, un petit dessin :-)
#
from matplotlib import pyplot as plt

plt.rcParams['toolbar'] = 'None'
plt.figure("Cubic splines :-)")
plt.plot(X,U,'or',markerSize='5',label='Data points')
plt.plot(x,uh,'-b',label='Cubic splines interpolation (or extrapolation :-)')
plt.legend(loc='lower left')
plt.show()