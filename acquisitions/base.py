# Copyright (c) 2016, the GPyOpt Authors
# Licensed under the BSD 3-clause license (see LICENSE.txt)

from ..core.task.cost import constant_cost_withGradients

class AcquisitionBase(object):
    """
    Base class for acquisition functions in Bayesian Optimization

    :param model: GPyOpt class of model
    :param space: GPyOpt class of domain
    :param optimizer: optimizer of the acquisition. Should be a GPyOpt optimizer

    """

    analytical_gradient_prediction = False

    def __init__(self, model, space, optimizer, cost_withGradients=None):
        self.model = model
        self.space = space
        self.optimizer = optimizer
        self.analytical_gradient_acq = self.analytical_gradient_prediction and self.model.analytical_gradient_prediction # flag from the model to test if gradients are available

        if cost_withGradients is  None:
            self.cost_withGradients = constant_cost_withGradients
        else:
            self.cost_withGradients = cost_withGradients

    @staticmethod
    def fromDict(model, space, optimizer, cost_withGradients, config):
        raise NotImplementedError()

    def acquisition_function(self,x):
        """
        Takes an acquisition and weights it so the domain and cost are taken into account.
        """
        import copy
        import numpy as np
        import matplotlib.pyplot as plt
        # Vanilla acquisition function
        f_acqu       = self._compute_acq(x)
        # Subject that to barriers
        f_acqu_new   = copy.deepcopy(f_acqu)
        #####TSA:: Check if there are barriers and then force values to low
        if len(self.barriers):
            print('Using barriers to barricade')
            # Visit the barriers one by one
            invalidities    = (f_acqu!=np.inf)
            for b in self.barriers:
                # Within all the barriers.
                # Higher than down lim and lower than upper lim
                invalidities   = invalidities &(x>b[0]) &(x<b[1])

            # Force to low value
            f_acqu_new[invalidities]   = 0
            if len(f_acqu_new)==1:
                print("Invalidated an optimum pick")


        # Plotting for sanity check
        if len(x)>1:
            # plt.figure()
            Idxs    = x.flatten().argsort()
            print(Idxs.shape)
            print(f_acqu_new.flatten()[Idxs].shape)
            print(x.flatten()[Idxs].shape)
        # plt.plot(x.flatten()[plotorder], f_acqu_new.flatten()[plotorder] - 1)

        # print(np.shape(x))
        # print(np.shape(f_acqu_new))
        #
        # plotorder    = x.flatten().argsort()
        # plt.figure()
        # plt.plot(x.flatten()[plotorder], f_acqu_new.flatten()[plotorder])
        # plt.show()

                # if xx >= 0.9 or xx <= 0.5:
                #     f_acqu_new.append(yy)
                # else:
                #     f_acqu_new.append(0)
            # f_acqu_new = np.array(f_acqu_new)[:, np.newaxis]

        cost_x, _ = self.cost_withGradients(x)
        return -(f_acqu_new*self.space.indicator_constraints(x))/cost_x



        # idxs   = x.flatten().argsort()
        # yplt   = f_acqu.flatten()[idxs]
        # xplt   = x.flatten()[idxs]
        # y_newplt   = f_acqu_new[idxs]
        #
        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(xplt, yplt)
        # plt.plot(xplt, y_newplt)
        #
        # plt.show()
        # exit()
        # return -(f_acqu*self.space.indicator_constraints(x))/cost_x


    def acquisition_function_withGradients(self, x):
        """
        Takes an acquisition and it gradient and weights it so the domain and cost are taken into account.
        """
        f_acqu,df_acqu = self._compute_acq_withGradients(x)
        cost_x, cost_grad_x = self.cost_withGradients(x)
        f_acq_cost = f_acqu/cost_x
        df_acq_cost = (df_acqu*cost_x - f_acqu*cost_grad_x)/(cost_x**2)
        return -f_acq_cost*self.space.indicator_constraints(x), -df_acq_cost*self.space.indicator_constraints(x)

    def optimize(self, duplicate_manager=None):
        """
        Optimizes the acquisition function (uses a flag from the model to use gradients or not).
        """
        if not self.analytical_gradient_acq:
            out = self.optimizer.optimize(f=self.acquisition_function, duplicate_manager=duplicate_manager)
        else:
            out = self.optimizer.optimize(f=self.acquisition_function, f_df=self.acquisition_function_withGradients, duplicate_manager=duplicate_manager)
        return out

    def _compute_acq(self,x):

        raise NotImplementedError('')

    def _compute_acq_withGradients(self, x):

        raise NotImplementedError('')
