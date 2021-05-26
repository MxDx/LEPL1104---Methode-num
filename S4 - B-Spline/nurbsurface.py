#
# Une surface NURBS
# Vincent Legat - 2018
# Ecole Polytechnique de Louvain
#

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt 
import mpl_toolkits.mplot3d 

# =========================================================================

def b(t,T,i,p):
  if p == 0:
    return (T[i] <= t)*(t < T[i+1])
  else:
    u  = 0.0 if T[i+p ]  == T[i]   else (t-T[i])/(T[i+p]- T[i]) * b(t,T,i,p-1)
    u += 0.0 if T[i+p+1] == T[i+1] else (T[i+p+1]-t)/(T[i+p+1]-T[i+1]) * b(t,T,i+1,p-1)
    return u

# ============================= mainProgram ===============================

X = [[1,2,3],[1,2,3],[1,2,3]]
Y = [[1,1,1],[2,2,2],[3,3,3]]
Z = [[1,2,3],[2,3,1],[3,1,4]]
T = [0.1,0,0,1,1,1.1]
S = [0.1,0,0,1,1,1.1]
p = 2; n = 5

t = np.arange(T[p],T[n-p]+0.05,0.1)
s = np.arange(S[p],S[n-p]+0.05,0.1)
Bt = np.zeros((n-p,len(t)))
Bs = np.zeros((n-p,len(s)))
for i in range(0,n-p):
  Bt[i,:] = b(t,T,i,p)
  Bs[i,:] = b(s,S,i,p)

x = Bs.T @ X @ Bt
y = Bs.T @ Y @ Bt
z = Bs.T @ Z @ Bt

matplotlib.rcParams['toolbar'] = 'None'
myColorMap = matplotlib.cm.jet

plt.figure("Surface NURBS")
ax = plt.axes(projection='3d',proj_type='ortho',azim=235)
ax.plot_surface(x,y,z,alpha=1,cmap=myColorMap,linewidth=0.5,edgecolors='k')
ax.set_xticks([1,1.5,2,2.5,3])
ax.set_yticks([1,1.5,2,2.5,3])
ax.set_zticks([1,1.5,2,2.5,3,3.5,4])
ax.margins(x=0,y=0,z=0)
plt.tight_layout()
plt.show()


# =========================================================================
