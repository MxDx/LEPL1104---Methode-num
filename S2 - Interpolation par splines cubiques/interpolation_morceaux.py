import numpy as np
from scipy.interpolate import CubicSpline as spline 
import matplotlib.pyplot as plt

X = np.arange(-55, 75, 10)
U =  [3.25,3.37,3.35,3.20,3.12,3.02,3.02,3.07,3.17,3.32,3.30,3.22,3.10] 

x = np.linspace(X[0],X[-1],100)
uh = spline(X,U)(x)


