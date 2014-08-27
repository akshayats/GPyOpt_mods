import numpy as np
from .util.general import reshape, ProjNullSpace
from random import random
from numpy.linalg inport norm

## 
## functions used in sampling steps 
##
def samples_multimensional_uniform(bounds,num_data):
	'''
	Sample uniformly from a box constrain defined by bounds
	'''
        dim = len(bounds)
        Z_rand = np.zeros(shape=(num_data,dim))
        for k in range(0,dim): Z_rand[:,k] = np.random.uniform(low=bounds[k][0],high=bounds[k][1],size=num_data)
	return Z_rand

def  Slice_ShrinkRank(xx, P, dP, s0):
	'''
	Multivariate slice sampling with shrinking rank covariance adaptation.
 	Based on the implementation of Philipp Henning.

	Inputs:
	param xx: the last sample of the Markov Chain
	param P : function returning the probabilty at xx
	param dP: function returning the derivarive of the probability at xx	
	param s0: initial proposal width (free parameter). 
	param input_dim: input dimension
	'''
	input_dim = len(xx)
	x = xx.reshape(input_dim,1)
	f = P(x)
	logf = np.log(f)
	logy = np.log(rand()) + logf

	theta = 0.95
	k = 0
	s = s0
	c = np.zeros(input_dim).reshape(input_dim,1)
	J = np.array([])

	while True:		
		c = np.append(ProjNullSpace(J,v), x + s[k] * np.random.random((input_dim, 1))
		sx = 1/sum(1/s)
		mx = sx * sum( np.dot((1/s), c-x))
		xk = x + ProjNullSpace(J, mx + sx * np.random.random((input_dim, 1)))
     
		fk = P(xk)
		dfk = dP(xk)
		logfk = np.log(fk)	
		dlogfk = dfk / fk
		if logfk > logy:
			return  xxk		 
		else:
			g = ProjNullSpace(J, dlogfk)

			if (J.shape[1] < input_dim - 1)	 and (g*dlogfk > 0.5 * norm(g) * norm(dlogfk)):
				J = np.append(J, g/norm(g))
				s = np.append(s,s[k])  ####
			else:
				s = np.append(s, theta * s[k]) 
				if s[k+1] < np.spacing(1)
					print 'bug found: contracted down to zero step size, still not accepted.\n' 
					return None
		k=+1
	return xxk    #check  


