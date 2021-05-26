from numpy import *
from scipy.interpolate import CubicSpline 



def periodicCubicSpline(X,Y,t,nLoop):
  x = linspace(1,0,len(t))
  y = linspace(1,1,len(t))  
  return x,y
