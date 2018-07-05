# -*- coding: utf-8 -*-
"""
Created on Wed 04 Jul 2017

@author: Alex Morgan, UNIVERSITY COLLEGE LONDON.
"""
from operator import attrgetter
import numpy as np
from tqdm import trange
from .basis import State, Basis, basis_states
from .interaction_matrix import interaction_matrix
from .constants import *

class Hamiltonian(object):
    """ The total Hamiltonian matrix.  Each element of the basis set is an
        instance of the class 'State', which represents |n L S J MJ>.
    """
    def __init__(self, n_min, n_max, L_max=None, S=None, M=None, M_max=None, basis_type='ml'):
        self.n_min = n_min
        self.n_max = n_max
        self.L_max = L_max
        self.S = S
        self.M = M
        self.M_max = M_max
        self.basis = basis_states(basis_type.lower(), n_min, n_max, L_max=L_max, S=S, M=M, M_max=M_max)
        self.sort_basis('E0', inplace=True)
        self.num_states = len(self.basis.states)
        self._h0_matrix = None
        self._stark_matrix = None
        self._zeeman_matrix = None
      
    def sort_basis(self, attribute, inplace=False):
        """ Sort basis on attribute.
        """
        sorted_basis = sorted(self.basis.states, key=attrgetter(attribute))
        if inplace:
            self.basis.states = sorted_basis
        return sorted_basis

    def attrib(self, attribute):
        """ List of given attribute values from all elements in the basis, e.g., J or E0.
        """
        return [getattr(el, attribute) for el in self.basis.states]

    def where(self, attribute, value):
        """ Indexes of where basis.attribute == value.
        """
        arr = self.attrib(attribute)
        return [i for i, x in enumerate(arr) if x == value]
    
    def filename(self):
        return  'n=' + str(self.n_min) + '-' + str(self.n_max) + '_' + \
                'L_max=' + str(self.L_max) + '_' + \
                'S=' + str(self.S) + '_' + \
                'M=' + str(self.M) + '_' + \
                'M_max=' + str(self.M_max)

    def h0_matrix(self, **kwargs):
        """ Unperturbed Hamiltonian.
        """
        cache = kwargs.get('cache_matrices', True)
        if self._h0_matrix is None or cache is False:
            self._h0_matrix = np.diag(self.attrib('E0'))
        return self._h0_matrix
        
    def stark_map(self, Efield, Bfield=0.0, **kwargs):
        """ The eigenvalues of H_0 + H_S + H_Z, for a range of electric fields.
        
            args:
                Efield           dtype: list      units: V / m      

                Bfield=0.0       dtype: float     units: T
            
            kwargs:
                field_angle=0.0    dtype: [float]

                                 specifies the angle between the electric and magnetic fields.
                                 
                eig_vec=False    dtype: bool

                                 returns the eigenvalues and eigenvectors for 
                                 every field value.

            Nb. A large map with eignvectors can take up a LOT of memory.
        """
        if (not kwargs.get('field_angle', 0.0) == 0.0) and \
            ((not self.M == None) or (not self.M_max == None)):
            print('WARNING: If the fields are not parallel then all'+\
                  ' ML sub-manifolds are required for accurate results!')
        
        tqdm_kwargs = dict([(x.replace('tqdm_', ''), kwargs[x]) for x in kwargs.keys() if 'tqdm_' in x])
        get_eig_vec = kwargs.get('eig_vec', False)
        num_fields = len(Efield)
        # initialise output arrays
        eig_val = np.empty((num_fields, self.num_states), dtype=float)
        if get_eig_vec:
            eig_vec = np.empty((num_fields, self.num_states, self.num_states), dtype=float)
        elif get_eig_vec_elements is not None:
            eig_vec_elements = []
        # optional magnetic field
        if Bfield != 0.0:
            Bz = mu_B * Bfield / En_h
            self._zeeman_matrix = interaction_matrix(matrix_type='zeeman', basis=self.basis, **kwargs)
            H_Z = Bz * self._zeeman_matrix.matrix
        else:
            H_Z = 0.0
        # loop over electric field values
        self._stark_matrix = interaction_matrix(matrix_type='stark', basis=self.basis, **kwargs)
        for i in trange(num_fields, desc="diagonalise Hamiltonian", **tqdm_kwargs):
            Fz = Efield[i] * e * a_0 / En_h
            H_S = Fz * self._stark_matrix.matrix / mu_me
            # diagonalise, assuming matrix is Hermitian.
            if get_eig_vec:
                # eigenvalues and eigenvectors
                eig_val[i], eig_vec[i] = np.linalg.eigh(self.h0_matrix(**kwargs) + H_S + H_Z)           
            else:
                # eigenvalues
                eig_val[i] = np.linalg.eigh(self.h0_matrix(**kwargs) + H_S + H_Z)[0]
        # output
        if get_eig_vec:
            return eig_val * En_h/mu_me, eig_vec
        else:
            return eig_val * En_h/mu_me

    def zeeman_map(self, Bfield, Efield=0.0, **kwargs):
        """ The eigenvalues of H_0 + H_S + H_Z, for a range of magnetic fields.
        
            args:
                Bfield           dtype: list      units: T      

                Efield=0.0       dtype: float     units: V / m
            
            kwargs:
                eig_vec=False    dtype: bool

                                 returns the eigenvalues and eigenvectors for 
                                 every field value.
            
            Nb. A large map with eignvectors can take up a LOT of memory.
        """
        tqdm_kwargs = dict([(x.replace('tqdm_', ''), kwargs[x]) for x in kwargs.keys() if 'tqdm_' in x])
        get_eig_vec = kwargs.get('eig_vec', False)
        get_eig_vec_elements = kwargs.get('eig_vec_elements', None)
        num_fields = len(Bfield)
        # initialise output arrays
        eig_val = np.empty((num_fields, self.num_states), dtype=float)
        if get_eig_vec:
            eig_vec = np.empty((num_fields, self.num_states, self.num_states), dtype=float)
        # optional electric field
        if Efield != 0.0:
            Fz = Efield * e * a_0 / En_h
            self._stark_matrix = interaction_matrix(matrix_type='stark', basis=self.basis, **kwargs)
            H_S = Fz * self._stark_matrix.matrix / mu_me
        else:
            H_S = 0.0
        # loop over magnetic field values
        self._zeeman_matrix = interaction_matrix(matrix_type='zeeman', basis=self.basis, **kwargs)
        for i in trange(num_fields, desc="diagonalise Hamiltonian", **tqdm_kwargs):
            Bz = mu_B * Bfield[i] / En_h
            H_Z =  Bz * self._zeeman_matrix.matrix 
            # diagonalise, assuming matrix is Hermitian.
            if get_eig_vec:
                # eigenvalues and eigenvectors
                eig_val[i], eig_vec[i] = np.linalg.eigh(self.h0_matrix(**kwargs) + H_S + H_Z)          
            else:
                # eigenvalues
                eig_val[i] = np.linalg.eigh(self.h0_matrix(**kwargs) + H_S + H_Z)[0]
        # output
        if get_eig_vec:
            return eig_val * En_h, eig_vec
        else:
            return eig_val * En_h