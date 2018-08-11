# coding: utf-8

from hsfs import HamiltonianMatrix
import numpy as np
import os

#Parameters to set. Make sure there is no space between '=' sign, and an '_' before the variable name
_n_min=20
_n_max=21
_S=1

_Efield=0.0
_field_angle=90.0
_Bfield=0.1

_matrices_dir='saved_matrices'

print(('n_min={}, n_max={}, S={}').format(_n_min, _n_max, _S))
ham0 = HamiltonianMatrix(n_min=_n_min, n_max=_n_max, S=_S)
print('Number of basis states:', '%d'%ham0.num_states)

Efield = np.linspace(_Efield, _Efield, 1) # V /cm
print('E={}, B={}, angle={}'.format(_Efield, _Bfield, _field_angle))
sm0 = ham0.stark_map(Efield*1e2, Bfield=_Bfield, 
                     field_angle=_field_angle, 
                     cache_matrices=False,
                     load_matrices=True,
                     save_matrices=True,
                     matrices_dir=_matrices_dir)

# Save Stark Map
filename_ham = ham0.filename()
filename_fields = 'E={}_B={}_angle={}'.format(_Efield, _Bfield, _field_angle)
filename_full = 'starkMap_{}__{}'.format(filename_ham, filename_fields)
filepath = os.path.join('.', _matrices_dir, filename_full)
if not os.path.isdir(os.path.join('.', _matrices_dir)):
    os.mkdir(os.path.join('.', _matrices_dir))
np.savez_compressed(filepath, matrix=sm0)
print('Filepath={}.npz'.format(filepath))