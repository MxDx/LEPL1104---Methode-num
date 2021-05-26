def diophantine(a,b,c,nmax,normRequest) :
    best_dist = normRequest
    sol = [0,0,0]
    for i in range(1, nmax+1):
        for j in range(1, nmax+1):
            for k in range(1, nmax+1):
                if (i**a) + (j**b) == (k**c):
                    dist_norm = abs(normRequest-((i**2)+(j**2)+(k**2)))
                    if dist_norm <= best_dist and i > sol[0]:
                        best_dist = dist_norm
                        sol = [i, j, k]

    return sol
#
# -1- Test de la fonction diophantine : on calcule la solution de l’´equation x**a + y**b = z**c
# dont la somme des carr´es est la plus proche de 1800
# On effectue la recherche pour tous les entiers compris entre 0 et 30 :-)
#
a = 2; b = 3; c = 4
xmax = 30
norm = 800
x = diophantine(a,b,c,xmax,norm)
print("==== Computing the solution of the diophantine equation :")
print(" === [a,b,c] = [%d,%d,%d] xmax = %d norm = %d " %(a,b,c,xmax,norm))
print(" === Solution is ",x)
#
# -2- Le plus grand triangle rectangle de c^ot´es entiers de longueur inf´erieure ou ´egale `a 20 :-)
# dont la somme des carr´es est la plus proche de 200
#
x = diophantine(2,2,2,20,200)
from matplotlib import pyplot as plt

plt.rcParams['toolbar'] = 'None'
plt.rcParams['figure.facecolor'] = 'moccasin'
plt.figure("Diophantine equations :-)")
plt.axis('equal')
plt.axis('off')
plt.text(x[0]/2 ,-0.5 ,'x = %d'%x[0],weight='bold',color='red')
plt.text(x[0]+0.5,x[1]/2,'y = %d'%x[1],weight='bold',color='red')
plt.text(x[0]/4 ,x[1]/2,'z = %d'%x[2],weight='bold',color='red')
plt.text(-1 ,x[1] ,'%d + %d = %d*%d + %d*%d = %d*%d = %d'%
(x[0]*x[0],x[1]*x[1],x[0],x[0],x[1],x[1],x[2],x[2],x[2]*x[2]),weight='bold',color='red')
plt.fill([0,x[0],x[0],0],[0,0,x[1],0],'-b')
plt.show()
