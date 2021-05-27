from numpy import *
from scipy.sparse.linalg import spsolve 
import scipy.sparse as sparse
from timeit import default_timer as timer

#
# - La solution du prof :-)
#   Y a sans doute moyen de faire mieux....
#   Faudrait aussi écrire des commentaires : but I am very lazy !
#
#   Globalement : les trucs :-)
#    -1- Ne pas résoudre les variables contraintes en introduisant une numérotation ad-hoc
#        dans le tableau "map"
#    -2- Ne résoudre qu'une partie du problème symétrique....
#        Attention : indexFull permet d'écrire les conditions de symétrie 
#    -3- Construire directement la matrice csr sans passer par une matrice dok
#    -4- A la fin, reconstruire la partie symétrique de la solution....
#
#   C'est en partie inspiré de la fonction MATLAB numgrid et delsq de Moler :-)
#   et d'un clone de ces fonctions en Python trouvé sur le net
#   https://stackoverflow.com/questions/21097657/numpy-method-to-do-ndarray-to-vector-mapping-as-in-matlabs-delsq-demo
#
#   Une vraie mine d'or, ce site au passage :-)
#


def poissonSolve(nCut) :

  n = 2*nCut + 1;  m = 0; h = 2/(n-1)   
  map = zeros((n,n),dtype=int)
  for i in range(nCut-1):
    map[i+1,1:-(1+i)] = arange(n-2-i)+1 + m
    m += (n-2-i)  
  map = map.flatten() 
  index = where(map)[0]
  indexFull = concatenate([index,(arange(nCut-1)+2)*(n-1)])
  
  i = j =  arange(m)
  coeff = 4*ones(m)

  for iNeigh in [-1,1,n,-n]:
     mapNeigh = map[indexFull+iNeigh]
     indexNeigh = where(mapNeigh)[0]
     mNeigh = len(indexNeigh)
     i = concatenate([i,map[indexFull[indexNeigh]]-1])
     j = concatenate([j,mapNeigh[indexNeigh]-1])
     coeff = concatenate([coeff,-ones(mNeigh)])

  A = sparse.csr_matrix((coeff,(i,j)),(m,m))
  B = ones(m)
  X = (spsolve(A,B))*(h*h) 

    
  U = zeros((n,n)); i,j = unravel_index(index,(n,n))
  U[i,j] = X; U[::-1,::-1][j,i] = X
  return U

#
# -2- Visualisation de la solution
#
 
n = 20
tic = timer()
U = poissonSolve(n) 
print("      Elapsed time : %f seconds" % (timer() - tic))
print(" ==== Maximum value of U : %10.8f " % amax(U))
print(" ==== Minimum value of U : %10.8f " % amin(U))
 
 
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['toolbar'] = 'None'
myColorMap = matplotlib.cm.jet
 
 
X = linspace(-1,1,2*n+1); U = abs(U)
plt.figure("Python as Matlab clone...")
plt.contourf(X,X,U,10,cmap=myColorMap)
plt.contour(X,X,U,10,colors='k',linewidths=1)
plt.hlines(X,X.min(),X.max(),color='white',linewidths=0.5)
plt.vlines(X,X.min(),X.max(),color='white',linewidths=0.5)
plt.axis("off"); plt.axis("equal")
plt.show()
        


        

