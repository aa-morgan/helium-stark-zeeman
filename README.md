`Helium Stark Zeeman`
===============

Calculate Stark and Zeeman maps for Rydberg helium, in either the `n, L, S, M_L` or `n, L, S, J, M_J` basis, using the
Numerov method.

M. L. Zimmerman et al., Phys. Rev. A, 20, 2251 (1979)
http://dx.doi.org/10.1103/PhysRevA.20.2251

Install
-------

Install using `setuptools`,
```bash
git clone https://github.com/axm108/helium-stark-zeeman
cd helium-stark-zeeman
python setup.py install
```

Basic usage
-------

Import libraries,
```python
from hsz import HamiltonianMatrix
import numpy as np
from scipy.constants import h
```
Instantiate `HamiltonianMatrix` object,
```python
ham = Hamiltonian(n_min=20, n_max=22, S=1, basis_type='ML')
```

Calculate Stark map,
```python
Efield = np.linspace(0.0, 1000.0, 101) # V /m
sm = ham.stark_map(Efield,
                   Bfield=0.1, # Telsa
                   field_angle=0.0, # Degrees
```

Plot Stark map,
```python
plt.plot(Efield, sm / (h*10**9), '-')
plt.xlabel('electric field (V cm$^{-1}$)')
plt.ylabel('energy / $h$ (GHz)')
```

Parameters
-------

#### Class: `HamiltonianMatrix`
| Parameter  | Description | Data type | Required | Default |
| ------------ | ------------ | ------------ | ------------ | ------------ |
| `n_min` | Minimum value of the principle quantum number `n`, to allow in the basis. | `Int` | Yes | N/A |
| `n_max` | Maximum value of the principle quantum number `n`, to allow in the basis. | `Int` | Yes | N/A |
| `L_max` | Maximum value of the orbital angular momentum quantum number to allow in the basis. `None` means no restriction. | `Int` [or `None`] | No | `None` |
| `S` | Value of the total spin orbital angular momentum. [Singlet: `S=0`, Triplet: `S=1`]. `None` means both.  | `Int` [or `None`] | No | `None` |
| `M` | Single value of the azimuthal quantum number to use in the basis. `None` means no restriction. | `Int` [or `None`] | No | `None` |
| `M_max` | Maximum value of the azimuthal quantum number to allow in the basis. `None` means no restriction. | `Int` [or `None`] | No | `None` |
| `basis_type` | Whether to use the `n, L, S, M_L`, or `n, L, S, J, M_J` basis. Specify using `'ML'`, or `'MJ'`, respectively. | `String` | No | `'ML'` |

#### Method: `stark_map`
| Parameter  | Description | Data type | Required | Default |
| ------------ | ------------ | ------------ | ------------ | ------------ |
| `Efield`  | Electric field values in V/cm | `numpy.array` | Yes  | N/A |
| `Bfield`  | Magnetic field values in T  | `Float`  | No  | `0.0` |
| `field_angle` | Angle between electric and magnetic fields | `Float` | No | `0.0` |
| `cache_matrices` | Cache Stark and Zeeman matrices | `Boolean` | No | `True` |
| `load_matrices` | Load Stark and Zeeman matrices from .npz files | `Boolean` | No | `False` |
| `save_matrices` | Save Stark and Zeeman matrices to .npz files | `Boolean` | No | `False` |
| `matrices_dir` | Directory to save Stark and Zeeman matrices | `String` | No | `'./'` |
| `eig_vec` | Return eigenvectors | `Boolean` | No | `False` |
| `tqdm_disable` | Disable `tqdm` output | `Boolean` | No | `False` |

#### Method: `zeeman_map`
| Parameter  | Description | Data type | Required | Default |
| ------------ | ------------ | ------------ | ------------ | ------------ |
| `Bfield`  | Magnetic field values in T | `numpy.array` | Yes  | N/A |
| `Efield`  | Electric field values in V/cm  | `Float`  | No  | `0.0` |
| `field_angle` | Angle between electric and magnetic fields | `Float` | No | `0.0` |
| `cache_matrices` | Cache Stark and Zeeman matrices | `Boolean` | No | `True` |
| `load_matrices` | Load Stark and Zeeman matrices from .npz files | `Boolean` | No | `False` |
| `save_matrices` | Save Stark and Zeeman matrices to .npz files | `Boolean` | No | `False` |
| `matrices_dir` | Directory to save Stark and Zeeman matrices | `String` | No | `'./'` |
| `eig_vec` | Return eigenvectors | `Boolean` | No | `False` |
| `tqdm_disable` | Disable `tqdm` output | `Boolean` | No | `False` |

Version information
-------------------

| Library  | Version |
| ------------ | ------------ |
| `Python`  | 3.6.1 64bit [GCC 4.2.1 Compatible Apple LLVM 6.0 (clang-600.0.57)] |
| `IPython` | 5.3.0 |
| `OS` | Darwin 17.4.0 x86_64 i386 64bit |
| `attr` | 17.4.0 |
| `matplotlib` | 2.0.2 |
| `numba` | 0.35.0 |
| `numpy` | 1.14.3 |
| `scipy` | 1.00.0 |
| `sympy` | 1.0 |
| `tqdm` | 4.15.0 |
| `version_information` | 1.0.3 |