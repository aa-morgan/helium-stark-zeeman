# -*- coding: utf-8 -*-
"""
Created on Wed 04 Jul 2017

@author: Alex Morgan, UNIVERSITY COLLEGE LONDON.
"""
from .numerov import radial_overlap
import numpy as np
import os.path
from tqdm import trange
from sympy.physics.wigner import clebsch_gordan, wigner_3j, wigner_6j
from .constants import *
import time

class InteractionMatrix:
    """
    """
    def __init__(self, matrix_type, basis, **kwargs):
        self.type = matrix_type.lower()
        self.basis = basis
        self.num_states = len(self.basis.states)
        self.matrix = None
        self.populate_interaction_matrix(**kwargs)
    
    def populate_interaction_matrix(self, **kwargs):
        """ Populate interaction matrix.
        """
        tqdm_kwargs = dict([(x.replace('tqdm_', ''), kwargs[x]) for x in kwargs.keys() if 'tqdm_' in x])
        cache = kwargs.get('cache_matrices', True)
        if self.matrix is None or cache is False:
            if kwargs.get('load_matrices', False) and \
               self.check_matrix(**kwargs):
                self.matrix = self.load_matrix(**kwargs)['matrix']
            else:
                self.matrix = np.zeros([self.num_states, self.num_states])
                for i in trange(self.num_states, desc='Calculating '+self.type+' terms', **tqdm_kwargs):
                    # off-diagonal elements only
                    for j in range(i, self.num_states):
                        self.matrix[i][j] = self.interaction_term(self.basis.states[i], self.basis.states[j], **kwargs)
                        # assume matrix is symmetric
                        self.matrix[j][i] = self.matrix[i][j]
                if kwargs.get('save_matrices', False):
                    self.save_matrix(**kwargs)  
        else:
            print("Using cached '{}' matrix".format(self.type))
            
    def interaction_term(self, state_1, state_2, **kwargs):
        """ Calculate interaction term
        """
        if self.type == 'stark':
            return stark_interaction(state_1, state_2, self.basis.params.basis_type, **kwargs)
        elif self.type == 'zeeman':
            return zeeman_interaction(state_1, state_2, self.basis.params.basis_type, **kwargs)
        else:
            raise Exception("Interaction term '{}' is not recognised!".format(self.type))  
    
    def save_matrix(self, **kwargs):
        filename = str(self)
        if self.type == 'stark':
            filename += '_angle={}'.format(kwargs.get('field_angle', 0.0))
        save_dir = os.path.join('.', kwargs.get('matrices_dir', ''))
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)
        date = time.strftime("%b %d %Y %H:%M:%S", time.gmtime(time.time()))
        np.savez_compressed(os.path.join(save_dir, filename), 
                            matrix=self.matrix, date=date, params=self.basis.params)
        print("Saved '{}' matrix to, ".format(self.type))
        print('\t', os.path.join(save_dir, '{}.{}'.format(filename, 'npz')))

    def load_matrix(self, **kwargs):
        filename = str(self)
        if self.type == 'stark':
            filename += '_angle={}'.format(kwargs.get('field_angle', 0.0))
        filename += '.npz'
        load_dir = os.path.join('.', kwargs.get('matrices_dir', ''))
        mat = np.load(os.path.join(load_dir, filename))
        print("Loaded '{}' matrix from, ".format(self.type))
        print('\t', os.path.join(load_dir, '{}.{}'.format(filename, 'npz')))
        return mat

    def check_matrix(self, **kwargs):
        filename = str(self)
        if self.type == 'stark':
            filename += '_angle={}'.format(kwargs.get('field_angle', 0.0))
        filename += '.npz'
        load_dir = os.path.join('.', kwargs.get('matrices_dir', ''))
        return os.path.isfile(os.path.join(load_dir, filename))
    
    def __str__(self):
        """ To String method
        """
        return '{}__{}'.format(self.type, self.basis)
    
############################
### Stark effect methods ###
############################
    
def stark_interaction(state_1, state_2, basis_type, **kwargs):
    if basis_type == 'ml':
        return stark_interaction_ML(state_1, state_2, **kwargs)
    elif basis_type == 'mj':
        return stark_interaction_MJ(state_1, state_2, **kwargs)
    else:
        raise Exception("Basis type '{}' not recognised".format(basis_type))

def stark_interaction_ML(state_1, state_2, **kwargs):
    """ Stark interaction between states |n1, l1, m> and |n2, l2, m>.
    """
    dL = state_2.L - state_1.L
    dM = state_2.M - state_1.M
    if (abs(dL) == 1) and (abs(dM) <= 1):
        p = kwargs.get('p', 1.0)
        # Stark interaction
        # TODO: Save the radial_overlap matrix, because this would allow fast recomputation for different angles
        return angular_overlap(state_1.L, state_2.L, state_1.M, state_2.M, **kwargs) * \
               radial_overlap(state_1.n_eff, state_1.L, state_2.n_eff, state_2.L, p=p)
    else:
        return 0.0

def stark_interaction_MJ(state_1, state_2, **kwargs):
    stark_method = kwargs.get('stark_method', '3j')
    if stark_method.lower() == '3j':
        return stark_interaction_Wigner_3j(state_1, state_2, **kwargs)
    elif stark_method.lower() == '6j':
        return stark_interaction_Wigner_6j(state_1, state_2, **kwargs)
    else:
        raise Exception("Stark interaction '{}' method not recognised!".format(stark_method))

def stark_interaction_Wigner_3j(state_1, state_2, **kwargs):
    """ Stark interaction between two states.

        <n' l' S' J' MJ'| H_S |n l S J MJ>.
    """     
    field_angle = kwargs.get('field_angle', 0.0)
    if not np.mod(field_angle, 180.0) == 90.0: # parallel fields
        field_orientation = 'parallel'
    elif not np.mod(field_angle, 180.0) == 0.0: # perpendicular fields
        field_orientation = 'perpendicular'
    else:
        raise Exception('Arbitrary angles not yet supported!')

    delta_L = state_1.L - state_2.L
    delta_S = state_1.S - state_2.S
    delta_MJ = state_1.M - state_2.M
    # Projection of spin, cannot change
    if abs(delta_L) == 1 and delta_S == 0 and \
     ((field_orientation=='parallel'      and     delta_MJ  == 0) or \
      (field_orientation=='perpendicular' and abs(delta_MJ) == 1)):
        # For accumulating each element in the ML sum
        sum_ML = []
        # Loop through all combination of ML for each state
        for MS_1 in np.arange(-state_1.S, state_1.S + 1):
            for MS_2 in np.arange(-state_2.S, state_2.S + 1):
                delta_MS = MS_1 - MS_2
                # Change in projection of spin:  0, +/- 1
                if ((field_orientation=='parallel'      and abs(delta_MS) in [0]) or \
                    (field_orientation=='perpendicular' and abs(delta_MS) in [0,1])):
                    ML_1 = state_1.M - MS_1
                    ML_2 = state_2.M - MS_2
                    if (abs(ML_1) <= state_1.L) and (abs(ML_2) <= state_2.L):
                        _angular_overlap = angular_overlap(state_1.L, state_2.L, ML_1, ML_2, **kwargs)
                        if _angular_overlap != 0.0:
                            sum_ML.append(float(clebsch_gordan(state_1.L, state_1.S, state_1.J,
                                          ML_1, state_1.M - ML_1, state_1.M)) * \
                                          float(clebsch_gordan(state_2.L, state_2.S, state_2.J,
                                          ML_2, state_2.M - ML_2, state_2.M)) * \
                                          _angular_overlap)

        # Stark interaction
        return np.sum(sum_ML) * radial_overlap(state_1.n_eff, state_1.L, state_2.n_eff, state_2.L)
    else:
        return 0.0

def stark_interaction_Wigner_6j(state_1, state_2, **kwargs):
    """ Stark interaction between two states.

        <n' l' S' J' MJ'| H_S |n l S J MJ>.
    """     
    field_angle = kwargs.get('field angle', 0.0)
    if not np.mod(field_angle, 180.0) == 90.0: # parallel fields
        q_arr   = [0]
        tau_arr = [1.]
    elif not np.mod(field_angle, 180.0) == 0.0: # perpendicular fields
        q_arr   = [1,-1]
        tau_arr = [(1./2)**0.5, -(1./2)**0.5]
    else:
        raise Exception('Arbitrary angles not yet supported!')

    delta_L = state_1.L - state_2.L
    delta_S = state_1.S - state_2.S
    delta_MJ = state_1.M - state_2.M  
    if abs(delta_L) == 1 and delta_S == 0:
        S = state_1.S
        sum_q = []
        for q, tau in zip(q_arr, tau_arr):
            sum_q.append( (-1.)**(int(state_1.J - state_1.M)) * \
                        wigner_3j(state_1.J, 1, state_2.J, -state_1.M, -q, state_2.M) * \
                        (-1.)**(int(state_1.L + S + state_2.J + 1.)) * \
                        np.sqrt((2.*state_1.J+1.) * (2.*state_2.J+1.)) * \
                        wigner_6j(state_1.J, 1., state_2.J, state_2.L, S, state_1.L) * \
                        (-1.)**state_1.L * np.sqrt((2.*state_1.L+1.) * (2.*state_2.L+1.)) * \
                        wigner_3j(state_1.L, 1, state_2.L, 0, 0, 0) * tau)

        return np.sum(sum_q) * radial_overlap(state_1.n_eff, state_1.L, state_2.n_eff, state_2.L)
    return 0.0 
                
def angular_overlap(L_1, L_2, M_1, M_2, **kwargs):
    angular_overlap_method = kwargs.get('angular_overlap_method', 'analytical')
    if angular_overlap_method == 'analytical':
        return angular_overlap_analytical(L_1, L_2, M_1, M_2, **kwargs)
    elif angular_overlap_method == 'wigner':
        return angular_overlap_wigner(L_1, L_2, M_1, M_2, **kwargs)
    else:
        raise Exception("Angular overlap method '{}' not recognised!".format(angular_overlap_method))
         
def angular_overlap_wigner(L_1, L_2, M_1, M_2, **kwargs):
    field_angle = kwargs.get('field_angle', 0.0)
    if not np.mod(field_angle, 180.0) == 90.0: # parallel fields
        q_arr   = [0]
        tau_arr = [1.]
    elif not np.mod(field_angle, 180.0) == 0.0: # perpendicular fields
        q_arr   = [1,-1]
        tau_arr = [(1./2)**0.5, (1./2)**0.5]
    else:
        raise Exception('Arbitrary angles not yet supported!')
            
    # For accumulating each element in the angular component, q sum
    sum_q = []
    for q, tau in zip(q_arr, tau_arr):
        sum_q.append(tau * float(wigner_3j(L_2, 1, L_1, -M_2, q, M_1)))
    # Calculate the angular overlap term using Wigner-3J symbols
    _angular_overlap = ((2*L_2+1)*(2*L_1+1))**0.5 * \
                          np.sum(sum_q) * \
                          wigner_3j(L_2, 1, L_1, 0, 0, 0)
    return _angular_overlap
                
def angular_overlap_analytical(L_1, L_2, M_1, M_2, **kwargs):
    """ Angular overlap <l1, m| cos(theta) |l2, m>.
        For Stark interaction
    """
    dL = L_2 - L_1
    dM = M_2 - M_1
    L, M = int(L_1), int(M_1)
    field_angle = kwargs.get('field_angle', 0.0)
    frac_para = np.cos(field_angle*(np.pi/180))**2
    frac_perp = np.sin(field_angle*(np.pi/180))**2
    dM_allow = kwargs.get('dM_allow', [0])
    overlap = 0.0
    if not np.mod(field_angle, 180.0) == 90.0:
        if (dM == 0) and (dM in dM_allow):
            if dL == +1:
                overlap += frac_para * (+(((L+1)**2-M**2)/((2*L+3)*(2*L+1)))**0.5)
            elif dL == -1:
                overlap += frac_para * (+((L**2-M**2)/((2*L+1)*(2*L-1)))**0.5)
        elif (dM == +1) and (dM in dM_allow):
            if dL == +1:
                overlap += frac_para * (-((L+M+2)*(L+M+1)/(2*(2*L+3)*(2*L+1)))**0.5)
            elif dL == -1:
                overlap += frac_para * (+((L-M)*(L-M-1)/(2*(2*L+1)*(2*L-1)))**0.5)
        elif (dM == -1) and (dM in dM_allow):
            if dL == +1:
                overlap += frac_para * (+((L-M+2)*(L-M+1)/(2*(2*L+3)*(2*L+1)))**0.5)
            elif dL == -1:
                overlap += frac_para * (-((L+M)*(L+M-1)/(2*(2*L+1)*(2*L-1)))**0.5)

    if not np.mod(field_angle, 180.0) == 0.0:
        if dM == +1:
            if dL == +1:
                overlap += frac_perp * (+(0.5*(-1)**(M-2*L))  * (((L+M+1)*(L+M+2))/((2*L+1)*(2*L+3)))**0.5)
            elif dL == -1:
                overlap += frac_perp * (-(0.5*(-1)**(-M+2*L)) * (((L-M-1)*(L-M))  /((2*L-1)*(2*L+1)))**0.5)
        elif dM == -1:
            if dL == +1:
                overlap += frac_perp * (+(0.5*(-1)**(M-2*L))  * (((L-M+1)*(L-M+2))/((2*L+1)*(2*L+3)))**0.5)
            elif dL == -1:
                overlap += frac_perp * (-(0.5*(-1)**(-M+2*L)) * (((L+M-1)*(L+M))  /((2*L-1)*(2*L+1)))**0.5)
    return overlap

#############################
### Zeeman effect methods ###
#############################

def zeeman_interaction(state_1, state_2, basis_type, **kwargs):
    if basis_type == 'ml':
        return zeeman_interaction_ML(state_1, state_2, **kwargs)
    elif basis_type == 'mj':
        return zeeman_interaction_MJ(state_1, state_2, **kwargs)
    else:
        raise Exception("Basis type '{}' not recognised".format(basis_type))

def zeeman_interaction_ML(state_1, state_2, **kwargs):
    """ Zeeman interaction between two states.
    """
    if state_1 == state_2:
        return state_1.M
    return 0.0

def zeeman_interaction_MJ(state_1, state_2, **kwargs):
    """ Zeeman interaction between two states.
    """
    delta_S = state_2.S - state_1.S
    delta_L = state_2.L - state_1.L
    delta_J = state_2.J - state_1.J
    delta_MJ = state_2.M - state_1.M
    if delta_MJ == 0 and \
       delta_J in [-1, 0, 1] and \
       delta_S == 0 and \
       delta_L == 0:
        L = state_1.L
        MJ = state_1.M
        S = state_1.S
        g_L2 = g_L * (((2 * L + 1) * L * (L + 1))/6)**0.5
        g_S2 = g_s * (((2 * S + 1) * S * (S + 1))/6)**0.5
        return (-1)**(1 - MJ) * ((2 * state_1.J + 1) * (2 * state_2.J + 1))**0.5 * \
            wigner_3j(state_2.J, 1, state_1.J, -MJ, 0, MJ) * 6**0.5 * ( \
            wigner_6j(L, state_2.J, S, state_1.J, L, 1) * \
            (-1)**(state_1.J + state_2.J + L + S) * g_L2 + \
            wigner_6j(state_1.J, state_2.J, 1, S, S, L) * \
            (-1)**(L + S) * g_S2)
    else:
        return 0.0