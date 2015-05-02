# Copyright (c) 2015, Javier Gonzalez
# Copyright (c) 2015, the GPy Authors (see GPy AUTHORS.txt)
# Licensed under the BSD 3-clause license (see LICENSE.txt)

import GPy

from ..core.acquisition import AcquisitionEI, AcquisitionMPI, AcquisitionLCB, AcquisitionEL1 
from ..core.bo import BO
from ..util.general import samples_multidimensional_uniform
import warnings
warnings.filterwarnings("ignore")


class BayesianOptimization(BO):
    def __init__(self, f, bounds=None, kernel=None, X=None, Y=None, model_data_init = None, model_optimize_interval=1, acquisition='EI', acquisition_par= 0.01, model_optimize_restarts=2, 
        sparseGP=False, num_inducing=None, normalize=False, true_gradients=True, exact_feval=False, verbosity=0):
        '''
        Bayesian Optimization using EI, MPI and LCB (or UCB) acquisition functions.
    
        This is a thin wrapper around the methods.BO class, with a set of sensible defaults
        :param *f* the function to optimize. Should get a nxp numpy array as imput and return a nx1 numpy array.
        :param bounds: Tuple containing the box constrains of the function to optimize. Example: for [0,1]x[0,1] insert [(0,1),(0,1)].
        :param kernel: a GPy kernel, defaults to rbf + bias.  
        :param X: input observations. If X=None, some  points are evaluated randomly.
        :param Y: output values. If Y=None, f(X) is used.
        :param model_data_init: number of initial random evaluations of f is X and Y are not provided (default, 3*input_dim).  
        :param model_optimize_interval: number of iterations after which the parameters of the model are optimized (1, Default). 
        :param acquisition: acquisition function ('EI': Expec. Improvement. 'MPI': Maximum Prob. Improvement. 'EL1': Expected Loss. LCB: Lower Confidence Bound). Default, EI.      
        :param acquisition_par: parameter of the acquisition function.    
        :param model_optimize_restarts: number of initial points for the GP parameters optimization (2, default)
        :param sparseGP: whether to use an sparse GP (False, default).
        :param num_inducing: number of inducing points for a Sparse GP (None, default)
        :param normalize: whether to normalize the Y's for optimization (False, default).
        :param true_gradients: whether the true gradients of the acquisition function are used for optimization (True, default). 
        :param exact_feval: set the noise variance of the GP if True (False, default).
        :param verbosity: whether to show (1) or not (0, default) the value of the log-likelihood of the model for the optimized parameters.
    
        '''
        # ------- Get default values 
        self.model_data_init = model_data_init  
        self.num_inducing = num_inducing
        self.sparseGP = sparseGP
        self.input_dim = len(bounds)
        self.normalize = normalize
        self.exact_feval = exact_feval
        self.model_optimize_interval = model_optimize_interval
        self.model_optimize_restarts = model_optimize_restarts
        self.verbosity = verbosity 
        if f==None: 
            print 'Function to optimize is required.'
        else:
            self.f = f  

        # ------- Initialize model 
        if bounds==None: 
            raise 'Box constraints are needed. Please insert box constrains.' 
        else:
            self.bounds = bounds
        if  model_data_init ==None:
            self.model_data_init = 3*self.input_dim
        else:
            self.model_data_init = model_data_init
        if X==None or Y == None:
            self.X = samples_multidimensional_uniform(self.bounds, self.model_data_init)
            self.Y = f(self.X)
        else:
            self.X = X
            self.Y = Y
        if kernel is None: 
            self.kernel = GPy.kern.RBF(self.input_dim, variance=.1, lengthscale=.1)  + GPy.kern.Bias(self.input_dim)
        else:
            self.kernel = kernel
        self._init_model()
        

        # ------- Initialize acquisition function
        self.acqu_name = acquisition
        if  acquisition_par == None:
            self.acquisition_par = 0
        else:
            self.acquisition_par = acquisition_par
        
        if acquisition==None or acquisition=='EI': 
            acq = AcquisitionEI(acquisition_par)
        elif acquisition=='MPI':
            acq = AcquisitionMPI(acquisition_par)
        elif acquisition=='LCB':
            acq = AcquisitionLCB(acquisition_par)
        elif acquisition=='EL1':
            acq = AcquisitionEL1(acquisition_par)
        else:   
            print 'The selected acquisition function is not valid. Please try again with EI, MPI, EL1 or LCB'
        if (acquisition=='EI' or acquisition=='MPI' or acquisition =='EL1' or acquisition =='LCB'):
            super(BayesianOptimization ,self).__init__(acquisition_func=acq)
    
    
    def _init_model(self):
        '''
        Initializes the Gaussian Process over *f*.
        :param X: input observations.
        :param Y: output values.

        ..Note : X and Y can be None. In this case model_data_init*input_dim data are uniformly generated to initialize the model.
        
        '''
        if self.sparseGP == True:
            if self.num_inducing ==None:
                raise 'Sparse model, please insert the number of inducing points'
            else:           
                self.model = GPy.models.SparseGPRegression(self.X, self.Y, kernel=self.kernel, num_inducing=self.num_inducing)
        else:       
            self.model = GPy.models.GPRegression(self.X,self.Y,kernel=self.kernel)
            
        if self.exact_feval == True:
            self.model.Gaussian_noise.constrain_fixed(1e-4) #to avoid numerical problems
        else:
            self.model.Gaussian_noise.constrain_bounded(1e-4,1e6) #to avoid numerical problems
            
            
            

