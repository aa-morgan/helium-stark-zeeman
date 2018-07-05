# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:08:34 2017

@author: Adam
"""
from math import ceil, log, exp
import numpy as np
from numba import jit

@jit
def wf(n, l, nmax, step=0.005, rmin=0.65):
    """ Use the Numerov method to find the wavefunction for state n*, l, where
        n* = n - delta.

        nmax ensures that wavefunctions from different values of n can be aligned.
    """
    W1 = -0.5 * n**-2.0
    W2 = (l + 0.5)**2.0
    rmax = 2 * nmax * (nmax + 15)
    r_in = n**2.0 - n * (n**2.0 - l*(l + 1.0))**0.5
    step_sq = step**2.0
    # ensure wf arrays will align using nmax
    if n == nmax:
        i = 0
        r_sub2 = rmax
    else:
        i = int(ceil(log(rmax / (2 * n * (n + 15))) / step))
        r_sub2 = rmax * exp(-i*step)
    i += 1

    # initialise
    r_sub1 = rmax * exp(-i*step)
    rvals = [r_sub2, r_sub1]
    g_sub2 = 2.0 * r_sub2**2.0 * (-1.0 / r_sub2 - W1) + W2
    g_sub1 = 2.0 * r_sub1**2.0 * (-1.0 / r_sub1 - W1) + W2
    y_sub2 = 1e-10
    y_sub1 = y_sub2 * (1.0 + step * g_sub2**0.5)
    yvals = [y_sub2, y_sub1]

    # Numerov method
    i += 1
    r = r_sub1
    while r >= rmin:
        ## next step
        r = rmax * exp(-i*step)
        g = 2.0 * r**2.0 * (-1.0 / r - W1) + W2
        y = (y_sub2 * (g_sub2 - (12.0 / step_sq)) + y_sub1 * \
            (10.0 * g_sub1 + (24.0 / step_sq))) / ((12.0 / step_sq) - g)

        ## check for divergence
        if r < r_in:
            dy = abs((y - y_sub1) / y_sub1)
            dr = (r**(-l-1) - r_sub1**(-l-1)) / r_sub1**(-l-1)
            if dy > dr:
                break

        ## store vals
        rvals.append(r)
        yvals.append(y)

        ## next iteration
        r_sub1 = r
        g_sub2 = g_sub1
        g_sub1 = g
        y_sub2 = y_sub1
        y_sub1 = y
        i += 1

    rvals = np.array(rvals)
    yvals = np.array(yvals)
    # normalisation
    yvals = yvals * (np.sum((yvals**2.0) * (rvals**2.0)))**-0.5
    return rvals, yvals

@jit
def find_first(arr, val):
    """ Index of the first occurence of val in arr.
    """
    i = 0
    while i < len(arr):
        if val == arr[i]:
            return i
        i += 1
    raise Exception('val not found in arr')

@jit
def find_last(arr, val):
    """ Index of the last occurence of val in arr.
    """
    i = len(arr) - 1
    while i > 0:
        if val == arr[i]:
            return i
        i -= 1
    raise Exception('val not found in arr')

@jit
def wf_align(r1, y1, r2, y2):
    """ Align two lists pairs (r, y) on r, assuming r array values overlap
        except at head and tail, and that arrays are reverse sorted.
    """
    if r1[0] != r2[0]:
        # trim front end
        if r1[0] > r2[0]:
            idx = find_first(r1, r2[0])
            r1 = r1[idx:]
            y1 = y1[idx:]
        else:
            idx = find_first(r2, r1[0])
            r2 = r2[idx:]
            y2 = y2[idx:]
    if r1[-1] != r2[-1]:
        # trim back end
        if r1[-1] < r2[-1]:
            idx = find_last(r1, r2[-1])
            r1 = r1[:idx + 1]
            y1 = y1[:idx + 1]
        else:
            idx = find_last(r2, r1[-1])
            r2 = r2[:idx + 1]
            y2 = y2[:idx + 1]
    if r1[0] == r2[0] and r1[-1] == r2[-1] and len(r1) == len(r2):
        return r1, y1, r2, y2
    else:
        raise Exception("Failed to align wavefunctions.")

@jit
def wf_overlap(r1, y1, r2, y2, p=1.0):
    """ Find the overlap between two radial wavefunctions (r, y).
    """
    r1, y1, r2, y2 = wf_align(r1, y1, r2, y2)
    return np.sum(y1 * y2 * r1**(2.0 + p))

@jit(cache=True)
def radial_overlap(n1, l1, n2, l2, p=1.0):
    """ Radial overlap for state n1, l1 and n2 l2.
    """
    nmax = max(n1, n2)
    r1, y1 = wf(n1, l1, nmax)
    r2, y2 = wf(n2, l2, nmax)
    return wf_overlap(r1, y1, r2, y2, p)