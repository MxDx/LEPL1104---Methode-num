#
# A nice "potje voor mij moeder" from NURBS
# Vincent Legat - 2018
# Ecole Polytechnique de Louvain
#

from numpy import *
import matplotlib 
import matplotlib.pyplot as plt 
import mpl_toolkits.mplot3d 


def b(t,T,i,p):
  if p == 0:
    return (T[i] <= t)*(t < T[i+1])
  else:
    u  = 0.0 if T[i+p ]  == T[i]   else (t-T[i])/(T[i+p]- T[i]) * b(t,T,i,p-1)
    u += 0.0 if T[i+p+1] == T[i+1] else (T[i+p+1]-t)/(T[i+p+1]-T[i+1]) * b(t,T,i+1,p-1)
    return u


# ============================= mainProgram ===============================


T = [0,0,0,1,2,3,4]
S = [0,0,0,1,1,2,2,3,3,4]
R = [0,10,1,15]
H = [0,0,10,15]


a = sqrt(3)
Xc = [1,0,1/2,1,3/2,2,1]
Yc = [0,0,a/2,a,a/2,0,0]
Xc = Xc - mean(Xc[0:6])
Yc = Yc - mean(Yc[0:6])
Zc = ones(shape(Xc))
Wc = [1,0.5,1,0.5,1,0.5,1]
X = outer(Xc,R)
Y = outer(Yc,R)
Z = outer(Zc,H)
W = outer(Wc,[1,1,1,1])

p = 2; h = 0.05

nt = size(T)-1
t = arange(T[p],T[nt-p]+h,h)
Bt = zeros((nt-p,len(t)))
for i in range(0,nt-p):
  Bt[i,:] = b(t,T,i,p)

ns = size(S)-1
s = arange(S[p],S[ns-p]+h,h)
Bs = zeros((ns-p,len(s)))
for i in range(0,ns-p):
  Bs[i,:] = b(s,S,i,p)
  
w = (Bs.T @ W @ Bt)
x = (Bs.T @ (W * X) @ Bt) / w
y = (Bs.T @ (W * Y) @ Bt) / w
z = (Bs.T @ (W * Z) @ Bt) / w

matplotlib.rcParams['toolbar'] = 'None'
myColorMap = matplotlib.cm.jet
plt.figure("Mooie potje voor moeder Legat")
ax = plt.axes(projection='3d')
ax.plot_surface(x,y,z,alpha=1,cmap=myColorMap,linewidth=0.5,edgecolors='k')
ax.set_axis_off()
ax.set_aspect('auto') 

# ax.view_init(90, 0)     # vue de haut
# ax.view_init(180, 0)    # vue de face
ax.view_init(45,0)        # vue normale

#
# Petit truc pour forcer les trois axes a être égaux :-)
# On fixe manuellement les limites du plot
#

max_range = array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
mid_x = (x.max()+x.min()) * 0.5
mid_y = (y.max()+y.min()) * 0.5
mid_z = (z.max()+z.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.show()



