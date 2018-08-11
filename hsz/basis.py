# -*- coding: utf-8 -*-
"""
@author: Alex Morgan, UNIVERSITY COLLEGE LONDON.
"""
import attr
import numpy as np
from .State import State

class Basis(object):
    """ Class to represent the a basis of States.
    """
    
    @attr.s()
    class Params(object):
        """ attrs class to represent the basis parameters.
        """
        basis_type = attr.ib(converter=str)
        @basis_type.validator
        def check_basis_type(self, attribute, value):
            if not value in ['ml', 'mj']:
                raise ValueError("basis_type must be either ML or MJ.")
        
        n_min = attr.ib(converter=int)
        @n_min.validator
        def check_n_min(self, attribute, value):
            if not value > 0:
                raise ValueError("n_min must be a positive integer.")
                
        n_max = attr.ib(converter=int)
        @n_max.validator
        def check_n_max(self, attribute, value):
            if not value >= self.n_min:
                raise ValueError("n_max must be larger than or equal to n_min.")
  
        S = attr.ib(default=None)
        @S.validator
        def check_S(self, attribute, value):
            if self.basis_type == 'ml' and not value in [0,1]:
                raise ValueError("S must be either 0, 1, or None (representation both).")
            if self.basis_type == 'mj' and not value in [0,1,None]:
                raise ValueError("S must be either 0, 1 when the basis type is ML.")
                
        L_max = attr.ib(default=None)
        M = attr.ib(default=None)        
        M_max = attr.ib(default=None)
        
        def __str__(self):
            """ To String method
            """
            return 'n={}-{}__L_max={}__S={}__M={}__M_max={}__basis_type={}'.format(
                self.n_min, self.n_max, self.L_max, self.S, self.M, self.M_max, self.basis_type)
        
    def __init__(self, states, basis_type, n_min, n_max, S, L_max=None, M=None, M_max=None):
        self.states = states
        self.params = self.Params(basis_type, n_min, n_max, S, L_max, M, M_max)
        
    def __str__(self):
        """ To String method
        """
        return str(self.params)
    
def basis_states(basis_type, n_min, n_max, **kwargs):
    if basis_type == 'ml':
        return basis_states_ML(basis_type, n_min, n_max, **kwargs)
    elif basis_type == 'mj':
        return basis_states_MJ(basis_type, n_min, n_max, **kwargs)
    else:
        raise Exception("Basis type '{}' not recognised".format(basis_type))

def basis_states_ML(basis_type, n_min, n_max, **kwargs):
    """ Generate the basis set: a list of instances of the attrs class State that 
        satisfy the given ranges of quantum numbers.  By default, all possible 
        states in the range of n_min to n_max are returned.
        
        args:
            n_min             Minimum value of the principal quantum number.

            n_max             Maximum value of the principal quantum number.
        
        kwargs:
            L_max = None      Maximum value of the orbital angular momentum quantum number.
                              If L_max is None 0 < L < n.

            S = None          Value of the total spin quanum number. If S is None S = [0, 1].

            ML = None         Value of the projection of the total angular momentum
                              quantum number. If MJ is None -J <= MJ <= J.

            ML_max = None     Maximum of the absolute value of the projection of the
                              total angular momentum quantum number. If MJ_max and MJ
                              are None -J <= MJ <= J.
    """
    L_max  = kwargs.get('L_max', None)
    S      = kwargs.get('S', None)
    ML     = kwargs.get('ML', None)
    ML_max = kwargs.get('ML_max', None)
    if S == None:
        raise Exception('S must be either 0 or 1 when the basis type is ML.')
    states = []
    n_rng  = np.arange(n_min, n_max + 1, dtype='int')
    # loop over n range
    for n in n_rng:
        # Don't add if n==1 and S==1
        if not(n==1 and S==1):
            if L_max is not None:
                _L_max = min(L_max, n - 1)
            else:
                _L_max = n - 1
            L_rng = np.arange(0, _L_max + 1, dtype='int')
            # loop over L range
            for L in L_rng:
                if L == 0:
                    J = S
                else:
                    J = L
                # loop over ML range
                if ML is None:
                    for _ML in np.arange(-L, L + 1):
                        if (ML_max is None) or (abs(_ML) <= ML_max):
                            states.append(State(n, L, S, J, _ML))
                elif -L <= ML <= L:
                    states.append(State(n, L, S, J, ML))
    return Basis(states, basis_type, n_min, n_max, S, L_max, ML, ML_max)

def basis_states_MJ(basis_type, n_min, n_max, **kwargs):
    """ Generate the basis set: a list of instances of the attrs class State that 
        satisfy the given ranges of quantum numbers.  By default, all possible 
        states in the range of n_min to n_max are returned.
        
        args:
            n_min             Minimum value of the principal quantum number.

            n_max             Maximum value of the principal quantum number.
        
        kwargs:
            L_max = None      Maximum value of the orbital angular momentum quantum number.
                              If L_max is None 0 < L < n.

            S = None          Value of the total spin quanum number. If S is None S = [0, 1].

            MJ = None         Value of the projection of the total angular momentum
                              quantum number. If MJ is None -J <= MJ <= J.

            MJ_max = None     Maximum of the absolute value of the projection of the
                              total angular momentum quantum number. If MJ_max and MJ
                              are None -J <= MJ <= J.
    """
    L_max  = kwargs.get('L_max', None)
    S      = kwargs.get('S', None)
    MJ     = kwargs.get('MJ', None)
    MJ_max = kwargs.get('MJ_max', None)
    states = []
    n_rng  = np.arange(n_min, n_max + 1, dtype='int')
    # loop over n range
    for n in n_rng:
        if L_max is not None:
            _L_max = min(L_max, n - 1)
        else:
            _L_max = n - 1
        L_rng = np.arange(0, _L_max + 1, dtype='int')
        # loop over L range
        for L in L_rng:
            if S is None:
                # singlet and triplet states
                S_vals = [0, 1]
            else:
                S_vals = [S]
            for _S in S_vals:
                # find all J vals and MJ substates
                if L == 0:
                    J = _S
                    if MJ is None:
                        for _MJ in np.arange(-J, J + 1):
                            if MJ_max is None or abs(_MJ) <= MJ_max:
                                states.append(State(n, L, _S, J, _MJ))
                    elif -J <= MJ <= J:
                        states.append(State(n, L, _S, J, MJ))
                elif _S == 0:
                    J = L
                    if MJ is None:
                        for _MJ in np.arange(-J, J + 1):
                            if MJ_max is None or abs(_MJ) <= MJ_max:
                                states.append(State(n, L, _S, J, _MJ))
                    elif -J <= MJ <= J:
                        states.append(State(n, L, _S, J, MJ))
                else:
                    for J in [L + _S, L, L - _S]:
                        if MJ is None:
                            for _MJ in np.arange(-J, J + 1):
                                if MJ_max is None or abs(_MJ) <= MJ_max:
                                    states.append(State(n, L, _S, J, _MJ))
                        elif -J <= MJ <= J:
                            states.append(State(n, L, _S, J, MJ))   
    return Basis(states, basis_type, n_min, n_max, S, L_max, MJ, MJ_max)