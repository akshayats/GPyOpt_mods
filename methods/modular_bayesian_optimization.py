# Copyright (c) 2016, the GPyOpt Authors
# Licensed under the BSD 3-clause license (see LICENSE.txt)

from ..core.bo import BO

class ModularBayesianOptimization(BO):

    """
    Modular Bayesian optimization. This class wraps the optimization loop around the different handlers.

    :param model: GPyOpt model class.
    :param space: GPyOpt space class.
    :param objective: GPyOpt objective class.
    :param acquisition: GPyOpt acquisition class.
    :param evaluator: GPyOpt evaluator class.
    :param X_init: 2d numpy array containing the initial inputs (one per row) of the model.
    :param Y_init: 2d numpy array containing the initial outputs (one per row) of the model.
    :param cost: GPyOpt cost class (defaut, none).
    :param normalize_Y: wheter to normalize the outputs before performing any optimization (default, True).
    :param model_update_interval: interval of collected observations after which the model is updated (default, 1).
    """

    def __init__(self, model, space, objective, acquisition, evaluator, X_init, Y_init=None, cost = None, normalize_Y = True, model_update_interval = 1):

        self.initial_iter = True
        
        # --- Create optimization space
        super(ModularBayesianOptimization ,self).__init__(  model                  = model, 
                                                            space                  = space, 
                                                            objective              = objective, 
                                                            acquisition            = acquisition, 
                                                            evaluator              = evaluator,
                                                            X_init                 = X_init, 
                                                            Y_init                 = Y_init,
                                                            cost                   = cost,
                                                            normalize_Y            = normalize_Y, 
                                                            model_update_interval  = model_update_interval)

        # --- Initilize everyting
        self.run_optimization(0)