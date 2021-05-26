import numpy as np
from scipy.spatial import Delaunay
import time

def magnetComputeInduction(Xmagnet,Ymagnet,Zmagnet,Xcoil,Ycoil,triangles,Xshift,mu0,mu) :


    m     = len(Xshift)
    n= len(triangles)
    phi = np.zeros(m)
    X_magnet=Xmagnet[triangles]
    Y_magnet=Ymagnet[triangles]
    surfaces= 1/2* ((X_magnet[:,1]-X_magnet[:,0]) * (Y_magnet[:,2]-Y_magnet[:,0]) - (X_magnet[:,2]-X_magnet[:,0])*(Y_magnet[:,1]-Y_magnet[:,0]))
    SMagnettotal=sum(surfaces)
    mu_e=-mu*surfaces/SMagnettotal
    X_magnetc = 1/3* (X_magnet[:,1]+X_magnet[:,0]+ X_magnet[:,2])
    Y_magnetc = 1/3* (Y_magnet[:,1]+Y_magnet[:,0]+ Y_magnet[:,2])

    X_coil=Xcoil[triangles]
    Y_coil=Ycoil[triangles]
    surf_coil= 1/2* ((X_coil[:,1]-X_coil[:,0]) * (Y_coil[:,2]-Y_coil[:,0]) - (X_coil[:,2]-X_coil[:,0])*(Y_coil[:,1]-Y_coil[:,0]))
    X_coilc = 1/3* (X_coil[:,1]+X_coil[:,0]+ X_coil[:,2])
    Y_coilc = 1/3* (Y_coil[:,1]+Y_coil[:,0]+ Y_coil[:,2])
     # boucle pour l'evolution dans le temps du magnet
    coef0=mu0/(4*np.pi)

    for i in range(m):        
        X_magnetc_t=X_magnetc+Xshift[i]
      
        Xc , Xm =np.meshgrid(X_coilc,X_magnetc_t)
        Yc , Ym =np.meshgrid(Y_coilc,Y_magnetc)
        r1=Xc-Xm
        r2=Yc-Ym
        r3=-Zmagnet
        norm_r = (r1*r1+r2*r2+r3*r3)**(1/2)
        coef=coef0/(norm_r * norm_r * norm_r)
        rchapeau3=-Zmagnet/norm_r
        dotprod=mu_e*rchapeau3
        BoCoil=np.sum(coef*(3*dotprod*rchapeau3-mu_e), axis=0)
        phi[i]=np.vdot(surf_coil,BoCoil) 
        print(" Iteration %2d  : shift = %6.3f [cm] : phi = %.8f" % (i,Xshift[i],phi[i]))
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
