# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 20-21
# Problème 5
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 


from matplotlib import pyplot as plt
import numpy as np
import time

# ============================================================
# FONCTIONS A MODIFIER [begin]
#

def trapezeEasy(f,a,b,n):
  I = 0.0
  i_n, h = np.linspace(a, b, n+1, retstep=True)
  U_n = f(i_n)
  U_n = U_n * 2
  U_n[0] = f(a)
  U_n[n] = f(b)
  
  I = np.sum(U_n) * (h/2)
  
  return I

  
def trapezeFun(f,a,b,n,nmax,tol):
  I = trapezeEasy(f, a, b, n//2)
  I_new = trapezeEasy(f, a, b, n)
  errorEst = abs(I_new - I)

  while errorEst > tol and n*2 < nmax:
    I = I_new
    n = n*2
    I_new = trapezeEasy(f, a, b, n)
  
    errorEst = abs(I_new - I)
  

  return I_new,n,errorEst 
  
#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
def u(x):
  return np.cos(x)
  
a = 0
b = np.pi/2
n = 10

start_time = time.time()

I = trapezeEasy(u,a,b,n)
errorExact = abs(1.0 - I)
print(" ======= Integral of sinus between 0 and pi/2 = %21.14e " % I)
print("  True error = %14.7e" % errorExact)
print("  Number of intervals = %d" % n)
print("\n")

I,n,errorEst = trapezeFun(u,a,b,n,200000,1e-12)
print(time.time()-start_time)
errorExact = abs(1.0 - I)
print(" ======= Integral of sinus between 0 and pi/2 = %21.14e " % I)
print("  True error = %14.7e" % errorExact)
print("  Est. error = %14.7e" % errorEst)
print("  Number of intervals = %d" % n)



plt.figure("Discovering the trapezoids integration rule :-)")
x = [a,b]
plt.plot(x,u(x),'.k',markersize=5) 
x = np.linspace(a,b,200)
plt.fill(np.append(x,[0]),np.append(u(x),[0]),'xkcd:sky blue')
x = np.linspace(-np.pi/2,np.pi,300)
plt.plot(x,u(x),'-k')
plt.title('Integrating sinus between 0 and pi/2')
plt.gca().axhline(y=0,color='k',linewidth=1.0)
plt.text(0.1,0.1,"I = %6.4f" % I,fontsize=12)
plt.show()