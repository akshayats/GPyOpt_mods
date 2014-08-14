from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

'''
Benchmark of one dimensional functions interesting to optimize. 

List of avaiable functions so far:
- forrester

The classes are oriented to create a python function which contain.
- *.f : the funtion itself
- *.plot: a plot of the function if the dimension is <=2.
- *.sensitivity: The Sobol coefficient per dimension when these are available.
- *.min : value of the global minimum(s) for the default parameters.

NOTE: the imput of .f must be a nxD numpy array. The dimension is calculated within the function.

Javier Gonzalez August, 2014
'''

class function1d:
	def plot(self):
		X = np.arange(0.0, 1.0, 0.01)
		Y = self.f(X)
		fig = plt.figure()
		plt.plot(X, Y, lw=2)
		plt.show()

class forrester(function1d):
	def __init__(self,sd=None):
		self.D = 1		
		if sd==None: self.sd = 0
		else: self.sd=sd
		self.min = 0.78 ## approx
		self.fmin = -6 ## approx
                
	def f(self,X):
		X = X.reshape((len(X),1))
		n = X.shape[0]
		fval = ((6*X -2)**2)*np.sin(12*X-4)
		if self.sd ==0:
			noise = np.zeros(n).reshape(n,1)
		else:
			noise = np.random.normal(0,self.sd,n).reshape(n,1)
		return fval.reshape(n,1) + noise


