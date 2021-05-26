import numpy as np

def compute(n,t,U):
    #
    # -1- Calcul du noeud qui precede Xestim
    # Gestion des cas critiques (extrapolation `a droite et `a gauche)
    #
    k = int(np.floor(t))
    k = n-2 if (k+2 > n) else k
    k = -n+1 if (k-1 < -n) else k
    #
    # -2- Estimation avec des polynomes de Lagrange
    # xi = coordonnee locale avec Xlocal = [-1 0 1 2]
    #
    Ulocal = U[n+k-1:n+k+3]
    xi = t - k
    phi = np.array([ -xi*(xi-1)*(xi-2) ,3*(xi+1)*(xi-1)*(xi-2),
    -3*(xi+1)*xi*(xi-2),(xi+1)*xi*(xi-1)]) / 6
    print(Ulocal)
    print(phi)
    return Ulocal @ phi 


def myCompute(n, t, U):
    
    #On vérifie si t est une extropolation
    #ou une interpolation

    index = int(t)
    index = n-2 if(index+2 > n) else index
    index = -n+1 if (index-1 < -n) else index

    Ulocal = U[n+index-1:n+index+3]

    return Ulocal