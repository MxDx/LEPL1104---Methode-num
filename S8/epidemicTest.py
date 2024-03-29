# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 20-21
# Problème 7
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

import numpy as np 
import time

start = time.time()

# -------------------------------------------------------------------------
#
# -1- Equations SIR 
#     offertes généreusement par l'équipe didactique !
#
#           dSdt(t) = - beta*S(i)*I(t)
#           dIdt(t) =   beta*S(i)*I(t) - gamma*I(t)
#           dRdt(t) =                  + gamma*I(t)
# 
#

def f(u,beta,gamma):
  dSdt = - beta * u[0] * u[1]
  dIdt =   beta * u[0] * u[1] - gamma * u[1]
  dRdt =                        gamma * u[1]  
  return np.array([dSdt,dIdt,dRdt])
 
 
# -------------------------------------------------------------------------    
#
# -2- Schema de d'Euler explicite d'ordre 1
# 
  
def epidemicEuler(Xstart,Xend,Ustart,n,beta,gamma):
  X,h = np.linspace(Xstart,Xend,n+1, retstep=True)  
  
#
# A MODIFIER ..... [begin]
#
  U = np.zeros((n+1,3))
  U[0] = Ustart
  for i in range(n):
    U[i+1] = U[i] + h*f(U[i,:], beta, gamma)

#
# A MODIFIER ..... [end]
# 

  return X,U
 

# -------------------------------------------------------------------------    
#
# -2- Schema de Taylor classique d'ordre 4
# 
  

def epidemicTaylor(Xstart,Xend,Ustart,n,beta,gamma):
  X,h = np.linspace(Xstart,Xend,n+1, retstep=True)
  
#
# A MODIFIER ..... [begin]
#
  U = np.zeros((n+1,3))
  U[0] = Ustart


  for i in range(n):
    du    = f(U[i], beta, gamma)

    d2u   = np.array([-beta*(du[0] * U[i, 1] + U[i,0] * du[ 1]), 
                      beta*(du[0] * U[i, 1] + U[i,0] * du[1]) - gamma*du[1], 
                      gamma*du[1]] )

    d3u   = np.array([-beta*(d2u[0] * U[i, 1] + 2*du[0]*du[1] + U[i,0] * d2u[1]), 
                      beta*(d2u[0] * U[i, 1] + 2*du[0]*du[1] + U[i,0] * d2u[1]) - gamma*d2u[1], 
                      gamma*d2u[1]] )

    d4u = np.array([-beta*(d3u[0] * U[i, 1] + 3*d2u[0]*du[1] + 3*du[0]*d2u[1] + U[i,0] * d3u[1]), 
                    beta*(d3u[0] * U[i, 1] + 3*d2u[0]*du[1] + 3*du[0]*d2u[1] + U[i,0] * d3u[1]) - gamma*d3u[1], 
                    gamma*d3u[1]] )

    U[i+1] = U[i] + h*du + (h**2)/2 * d2u + (h**3)/6 * d3u + (h**4)/24 * d4u
#
# A MODIFIER ..... [end]
#

  return X,U

# -------------------------------------------------------------------------    
#
# -3- Schema de Runge-Kutta d'ordre 4
# 
  
def epidemicRungeKutta(Xstart,Xend,Ustart,n,beta,gamma):
  X,h = np.linspace(Xstart,Xend,n+1,retstep=True)
  
#
# A MODIFIER ..... [begin]
#
  U = np.zeros((n+1,3))
  U[0] = Ustart

  for i in range(n):
    K1 = f(U[i,:], beta, gamma)
    K2 = f(U[i,:] + h/2 * K1, beta, gamma)
    K3 = f(U[i,:] + h/2 * K2, beta, gamma)
    K4 = f(U[i,:] + h * K3, beta, gamma)
    U[i+1] = U[i] + h/6 *(K1 + 2*K2 + 2*K3 + K4)
 
    
#
# A MODIFIER ..... [end]
# 
  return X,U

# -------------------------------------------------------------------------


# ============================= a so simple class ! =======================

class EpidemicExplMethods(object):
  "Class of explicit integrators for SIR equations"
#
# Un constructeur qui est exécutée automatiquement lorsqu'on instancie un
# nouvel objet de la classe : elle s'appelle obligatoirement __init__()
#
  def __init__(self,name,order,integrator):
    self.name = name
    self.order = order
    self.f = integrator

integrators = [EpidemicExplMethods("Explicit Euler",1,epidemicEuler),
               EpidemicExplMethods("Taylor",4,epidemicTaylor),
               EpidemicExplMethods("Runge-Kutta",4,epidemicRungeKutta)]


# =========================================================================


#
# -1- Paramètres de la simulation
#     Nombre initial de susceptibles = 10000
#     Nombre initial de contaminés = 2
#

Ustart = [10000,2,0]
beta  = 0.00001315
gamma = 0.01798476

#
# -2- Analyse de la convergence des trois méthodes
#

Xstart = 0; Xend = 50; Uref = 541.867329232525
print(" ============ Exact confirmed cases as reference  : Uref(%d) = %14.7e " % (Xend,Uref))

for integrator in integrators:
  print("\n ================")
  error = np.zeros(4)
  for j in range(4):
    n = 50*pow(2,j); h = (Xend - Xstart)/n  
    X,U = integrator.f(Xstart,Xend,Ustart,n,beta,gamma)
    error[j] = abs(U[-1,1]-Uref)
    print(" ==== %4s (order=%d) n=%6d h=%6.3f :    U(%d) = %14.7e : eh(Xend) = %8.2e " 
            % (integrator.name.rjust(15),integrator.order,n,h,Xend,U[-1,1],error[j]))
  order = np.mean(np.log(error[:-1]/error[1:])/np.log(2))
  print(" ================ Estimated order of %s : %.4f " % (integrator.name,order))
  
#
# -3- Comparons notre simulation avec des données réelles !
#

Xstart = 0; Xend = 150
m = 10; n = Xend * m

#
#  Pour comprendre pourquoi se confiner est utile :-)
#  Se confiner permet de réduire la valeur de beta
#
#  Tester les deux lignes qui précèdent !
#  Ensuite, on peut imaginer de faire varier beta dans le temps pour être plus réaliste !
#   
#     Xend = 300; beta = beta/2
#     Xend = 1200; beta = beta/4
#

X,U = epidemicRungeKutta(Xstart,Xend,Ustart,n,beta,gamma)

print("Time : {:e} [s]".format(time.time() - start))


#
# -3.1- Données du Japon à partir du 22 janvier 2020
#       obtenues ici : [HDX](https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases)
#

dataConfirmed = np.array([2,1,2,2,4,4,7,7,11,15,20,20,20,22,22,45,25,25,26,26,26,28,28,29,43,59,66,
                       74,84,94,105,122,147,159,170,189,214,228,241,256,274,293,331,360,420,461,
                       502,511,581,639,639,701,773,839,825,878,889,924,963])
dataRecovered = np.array([0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,4,9,9,9,9,12,12,12,13,18,18,22,22,
                       22,22,22,22,22,22,32,32,32,43,43,43,46,76,76,76,101,118,118,118,118,118,
                       144,144,144,150,191])
days = np.arange(dataConfirmed.shape[0])

#
# -3.2- Une jolie figure : et hop :-) 
#       Le choix du nombre de 10000 pour estimer le mélange homogène équivalent
#       est assez arbitraire et rend donc la prédiction assez aléatoire !
#
#       Toutefois : jouer sur les facteurs beta et gamma
#       permet d'avoir une compréhension intuitive de la dynamique d'une épidemie
#       Se confiner permet de réduire beta : refaites la simulation en diminuant beta
#       d'un facteur deux : qu'observez-vous ?
#
#       AVERTISSEMENT : ceci est un modèle jouet pour illustrer le cours LEPL1104
#       Vous n'êtes pas encore des épidémiologistes avertis même pas dilettantes....
#       Mais, cela devrait vous permettre de comprendre intuitivement la dynamique 
#       d'une épidémie !
#

from matplotlib import pyplot as plt

# 
#       Un petit truc pour éviter d'avoir la barre de commandes et un fond blanc
#       lorsqu'on fait une capture d'écran 
#

import matplotlib
matplotlib.rcParams['toolbar'] = 'None'
plt.rcParams['figure.facecolor'] = 'lavender'
plt.rcParams['axes.facecolor'] = 'lavender'

fig = plt.figure("SIR Equations")
X=X[::m]; U=U[::m,:]
plt.plot(X,U[:,0],'-b',X,U[:,1],'-r',X,U[:,2],'-g')
plt.plot(days[::5],dataConfirmed[::5],'or')
plt.plot(days[::5],dataRecovered[::5],'og')
plt.text(0,8800,"Susceptibles",color='blue',fontsize=12)
plt.text(38,3500,"Infectious",color='red',fontsize=12)
plt.text(120,7800,"Recovered",color='green',fontsize=12)
plt.show()



