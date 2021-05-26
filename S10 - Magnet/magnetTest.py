# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 20-21
# Problème 8
#
# Script de test
#  Vincent Legat
#
# Largement inspiré du programme de Nicolas Roisin :-)
# Ou les méthodes numériques pour obtenir la solution du projet P2 !
#
# -------------------------------------------------------------------------
#


import numpy as np
from scipy.spatial import Delaunay

import time

# ------------------------------------------------------------------------------------
#
# Intégration du flux magnétique
#
#    Xmagnet,Ymagnet : tableaux numpy des coordonnées x,z des sommets  de l'aimant 
#    Zmagnet : scalaire contenant la hauteur du l'aimant
#    Xcoil,Ycoil : tableaux numpy des coordonnées x,z des sommets de la bobine
#    triangles : tableau contenant les indices des 3 sommets de chaque élément
#    Xshift : tableau numpy contenant les translation de l'aimant sur une période
#    mu0 : perméabilité du vide
#    mu  : valeur absolue de la composante z du momemt magnétique du l'aimant = [0,0,-mu]
#   
#  La fonction renvoie un vecteur phi contenant le flux du champs magnétique intercepté
#  par une spire exprimé en [T cm2]
#


def g(X,Y,Coord):
  rep = np.zeros((len(Coord), 3))

  for i,coord in enumerate(Coord):

    X1,X2,X3 = X[Coord[i,:]]
    Y1,Y2,Y3 = Y[Coord[i,:]]
    rep[i,0] = (X1+X2+X3)/3
    rep[i,1] = (Y1+Y2+Y3)/3
    rep[i,2] = abs(np.linalg.det(np.array([[X2-X1, Y2-Y1], [X3-X1, Y3-Y1]]) )) /2

  return rep  


def magnetComputeInduction(Xmagnet,Ymagnet,Zmagnet,Xcoil,Ycoil,triangles,Xshift,mu0,mu) :
  
  m     = len(Xshift)
  phi = np.zeros(m)
  
#
# A MODIFIER ..... [begin]
#  
  gCoil = g(Xcoil, Ycoil, triangles)
  gMagnet = g(Xmagnet, Ymagnet, triangles)
  
  gXcoil, gYcoil, CoilArea = gCoil[:,0], gCoil[:,1], gCoil[:,2]
  gXmagnet, gYmagnet, MagnetArea = gMagnet[:,0], gMagnet[:,1], gMagnet[:,2]

  nNode = len(triangles)
  cst_1 = mu0/(4*np.pi)
  mu_prop = mu/np.sum(MagnetArea)

  for i in range(m):
    phi_i = 0.0
    shift = Xshift[i]

    for j in range(nNode):
      B_Xj = 0
      v_posj = np.array([gXcoil[j], gYcoil[j], 0])

      for f in range(nNode):        
        v_posf = np.array([gXmagnet[f]-shift, gYmagnet[f], Zmagnet])
        r_f = 1/(np.linalg.norm(v_posf-v_posj)**3)

        v_mu = np.array([gXmagnet[f]-shift, gYmagnet[f],-mu_prop*MagnetArea[f]])
        r = v_posf-v_posj
        norm_v_mu = r / np.linalg.norm(r)
        mu_norm = (3*np.dot(v_mu,norm_v_mu)*norm_v_mu-v_mu)[2]

        B_Xj += cst_1*r_f*mu_norm

      phi_i += B_Xj * CoilArea[j]

    phi[i] = phi_i

    print(" Iteration %2d  : shift = %6.3f [cm] : phi = %.8f" % (i,Xshift[i],phi[i]))
     
#
# A MODIFIER ..... [end]
#
   
  return phi

# ------------------------------------------------------------------------------------ 
#
# Script de test
#
#
# -0- Paramètres matériels
#
# ------------------------------------------------------------------------------------


mu0     = 4e-7*np.pi*1e-2     # permeabilité du vide en [H/cm] 
Rmagnet = 1.25             # rayon de l'aimant [cm]
Hmagnet = 0.6              # épaisseur de l'aimant [cm]
Zmagnet = 0.5              # position verticale de l'aimant en [cm]
Br      = 1.4              # magnetisation residuelle du Néodyme fer bore (NdFeB) en [T] ou [kg/(A s)]
mu      = Rmagnet**2*Hmagnet*np.pi*Br / mu0    
                           # moment magnétique de l'aimant [A cm2]
Rcoil   = 1                # rayon de la bobine [cm]
nSpires = 200



# ------------------------------------------------------------------------------------
#
# -1- Construction d'un maillage de triangles pour un cercle de rayon unitaire
#
# ------------------------------------------------------------------------------------


nR      = 6
nTheta  = 6
nNode   = 1 + sum(np.arange(1,nR))*nTheta
R     = np.zeros(nNode)
Theta = np.zeros(nNode)

index = 1; dR = 1.0/(nR-1)
for i in range(1,nR):
    dTheta = 2*np.pi/(i*nTheta)
    for j in range(0,i*nTheta):
        R[index]     = i*dR
        Theta[index] = j*dTheta; index += 1

X       = R*np.cos(Theta)
Y       = R*np.sin(Theta)

triangles = Delaunay(np.stack((X,Y),1)).simplices
nElem = len(triangles)

print(" Number of triangles : %d " % nElem)
print(" Number of nodes     : %d " % nNode)


# ------------------------------------------------------------------------------------
#
# -2- Calcul du flux et de la tension induite dans la bobine
#
# ------------------------------------------------------------------------------------

m       = 41
Xstart  = -5                        # [cm]
Xstop   =  5                        # [cm]
Xshift  = np.linspace(Xstart,Xstop,m)
Tstart  = 0                         # [s]
Tstop   = 0.5                       # [s]
T,delta = np.linspace(Tstart,Tstop,m,retstep=True)

Xmagnet = Rmagnet*R*np.cos(Theta)
Ymagnet = Rmagnet*R*np.sin(Theta) 
Xcoil   = Rcoil*R*np.cos(Theta)
Ycoil   = Rcoil*R*np.sin(Theta) 
    
start_time = time.time()

phi     = magnetComputeInduction(Xmagnet,Ymagnet,Zmagnet,Xcoil,Ycoil,triangles,
                                                                   Xshift,mu0,mu)  
print(time.time()-start_time)
phi     = phi * nSpires    
voltage = - np.diff(phi) / (delta*10)

# ------------------------------------------------------------------------------------
#
# -3- Quelques jolis plots et animation
#
# ------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.rcParams['toolbar'] = 'None'

def frame(i):
  plt.clf()

  n = 50
  X,Z =np. meshgrid(np.linspace(-2,2,n),np.linspace(-2,2,n))
  Y  = np.zeros_like(X)
  Bx = np.zeros(np.shape(X))
  Bz = np.zeros(np.shape(X))

  for iElem in range(nElem):
    Xp = X - Xdipole[iElem] - Xshift[i]
    Yp = Y - Ydipole[iElem]
    Zp = Z - Zmagnet
    r     = np.sqrt(Xp*Xp + Yp*Yp + Zp*Zp)
    coeff = -(mu0*mu) / (4*np.pi*r**5)
    Bx   += coeff * (3*Zp*Xp)
    Bz   += coeff * (3*Zp*Zp - r*r)
  plt.streamplot(X,Z,Bx,Bz, density=1.4, linewidth=None, color='blue')

  x = np.array([-Rmagnet,Rmagnet,Rmagnet,-Rmagnet,-Rmagnet]) + Xshift[i]
  y = np.array([0,0,Hmagnet,Hmagnet,0])+Zmagnet-Hmagnet/2.0
  plt.fill(x,y,facecolor='blue',alpha=1)

  x = [-Rcoil,Rcoil]
  y = [0,0]
  plt.plot(x,y,"-r",linewidth=4)
  
  plt.xlim((-2,2)); plt.ylim((-2,2))
  plt.title('Electromagnetic Field')   

# ------------------------------------------------------------------------------------

fig=plt.figure("Maillage de l'aimant")
plt.plot(Xmagnet,Ymagnet,'or')
plt.triplot(Xmagnet,Ymagnet,triangles,'-k')
Xdipole = np.mean(Xmagnet[triangles[:,:]],axis=1)  
Ydipole = np.mean(Ymagnet[triangles[:,:]],axis=1)  
plt.plot(Xdipole,Ydipole,'ob')  

plt.axis("equal")
plt.axis("off")

plt.figure("Flux et tension induite sur une période")
plt.plot(T,phi,'-r')
plt.plot(T[1:],voltage,'-b')
plt.text(0.01,-100,"$N\phi(t)$ [T cm$^2$]",color='red',fontsize=12)
plt.text(0.3,100,r"$-N\dfrac{\partial \phi}{\partial t}(t)$ [mV]",color='blue',fontsize=12)
plt.text(0.4,-210,"time [s]",color='black',fontsize=12)

plt.figure("Un joli plot pour le coordinateur :-)",figsize=(10, 10))
frame(20)

movie = animation.FuncAnimation(plt.figure("Claude's project",figsize=(10,10)),frame,41,interval=20,repeat=False)
plt.show()









