# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 20-21
# Problème 6
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 

import numpy as np
import time 

start_time = time.time()
# ============================================================
# FONCTIONS A MODIFIER [begin]
#
def I_ij(i, j, I_ij1, I_i1j1):
  return ( ((2**(2*j)*(I_ij1))-I_i1j1) / (2**(2*j)-1) ) 


def romberg(f,a,b,n,nmax,tol):

  I = 0.0
  errorEst = 0.0 

  x_n, h = np.linspace(a, b, n+1, retstep=True)
  I_0 = np.trapz(u(x_n), dx=h)
  
  n *= 2
  I_i1 = np.zeros(2)

  i_n, h = np.linspace(a, b, n+1, retstep=True)
  I_i1[0] = np.trapz(u(i_n), dx=h)

  I_i1[1] = I_ij(1, 1, I_i1[0], I_0)
  errorEst = abs(I_i1[1] - I_i1[0])

  i=2

  while n <= nmax and errorEst > tol:
    I_i = np.zeros(i+1) #I+1 pour avoir la prochaine longeur
    
    n = n*2
    i_n, h = np.linspace(a, b, n+1, retstep=True)
    I_i[0] = np.trapz(u(i_n), dx=h)

    for j in range(1,i+1):
      I_i[j] = I_ij(i, j, I_i[j-1], I_i1[j-1])


    errorEst = abs(I_i[i] - I_i[i-1])
    
    I_i1 = I_i
    i+=1

  I = I_i1[i-1]

  return I,n,errorEst  
  
  
#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
 
def u(x):
  return (x*x+x+1)*np.cos(x)
  
a = 0
b = np.pi/2
n = 1


I,n,errorEst = romberg(u,a,b,n,100,1e-8)
print(time.time()-start_time)
errorExact = abs(2.03819742706724 - I)
print("\n ======= Integral of (x*x+x+1)*cos(x) between 0 and pi/2 = %21.14e " % I)
print("  True error = %14.7e" % errorExact)
print("  Est. error = %14.7e" % errorEst)
print("  Number of intervals = %d" % n)


from matplotlib import pyplot as plt
import matplotlib 
matplotlib.rcParams['toolbar'] = 'None'

plt.figure("Mastering the trapezoids integration rule :-)")
x = np.array([a,b])
plt.plot(x,u(x),'.k',markersize=5) 
x = np.linspace(a,b,200)
plt.fill(np.append(x,[0]),np.append(u(x),[0]),'xkcd:sky blue')
x = np.linspace(-0.5,2.0,300)
plt.plot(x,u(x),'-k')
plt.title('Integrating (x*x+x+1)*cos(x) between 0 and pi/2')
plt.gca().axhline(y=0,color='k',linewidth=1.0)
plt.text(0.1,0.1,"I = %8.6f" % I,fontsize=12)
plt.show()