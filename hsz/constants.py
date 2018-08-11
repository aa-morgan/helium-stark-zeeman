# -*- coding: utf-8 -*-
"""
Created on Wed 04 Jul 2017

@author: Alex Morgan, UNIVERSITY COLLEGE LONDON.
"""
import pandas as pd

#CODATA 2014, DOI: 10.1103/RevModPhys.88.035009
c = 299792458.0 ## speed of light in vacuum
h = 6.626070040e-34
hbar = 1.054571800e-34
Ry = 10973731.568508
e = 1.6021766208e-19
m_e = 9.10938356e-31
alpha = 7.2973525664e-3
m_u = 1.660539040e-27
En_h = alpha**2.0 * m_e * c**2.0
a_0 = hbar/ (m_e * c * alpha)
mu_B = e * hbar / (2.0 * m_e)

## helium
A_r_He = 4.002603254130
ionization_He = 1.9831066637e7
mass_He = A_r_He * m_u
mass_He_core = mass_He - m_e + ionization_He * h /c
## reduced electron mass/ m_e
mu_me = mass_He_core / (mass_He_core + m_e)
## reduced electron mass / core mass,
mu_M = m_e / (mass_He_core + m_e)
## Rydberg constant for helium
Ry_He = Ry * mu_me
## g-factors
g_L = 1 - m_e / mass_He_core
g_s = 2.00231930436182
## Hartree energy for helium
En_h_He = 2*Ry_He*h*c;

def constants_info():
    constant_vals = {
        'speed of light in vacuum, $c$': c,
        'Planks constant, $h$': h,
        'Reduced Planks constant, $\hbar$': hbar,
        'Rydberg constant, $R_{\infty}$': Ry,
        'electron charge, $e$': e,
        'fine structure constant': alpha,
        'atomic mass': m_u,
        'Hatree energy': En_h,
        'Bohr radius, $a_0$': a_0,
        'Bohr magneton, $\mu_B$': mu_B,
        'ionization energy of helium': ionization_He,
        'mass of helium': mass_helium,
        'mass of helium (a.u.)': A_r_He,
        'mass of helium core': mass_He_core,
        'Reduced electron mass / electron mass': mu_me,
        'Reduced electron mass / core mass': mu_M,
        'Rydberg constant for helium': Ry_He
    }
    df = pd.DataFrame(list(constant_vals.items()), columns=['Constant', 'Value'])
    df['Value'] = df['Value'].map('{:.14g}'.format)
    return df