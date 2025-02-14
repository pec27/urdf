"""
Copyright (C) 2015-2025 Peter Creasey

Module for loading the BRDF library liburdf.so (in C). 
"""
from __future__ import print_function
from numpy.ctypeslib import ndpointer
import ctypes
from numpy import float64, empty, array, int32, float32, exp, empty_like
from os import path
from time import time
import numpy as np
from numpy.random import standard_normal

_liburdf = None

def initlib():
    """ Init the library (if not already loaded) """
    global _liburdf

    if _liburdf is not None:
        return _liburdf

    name = path.join(path.dirname(path.abspath(__file__)), '../build/liburdf.so')
    if not path.exists(name):
        raise Exception('Library '+str(name)+' does not exist. Maybe you forgot to make it?')

    print('Loading libbrdf - C functions for Unitary BRDF calculations', name)
    _liburdf = ctypes.cdll.LoadLibrary(name)

    # Isotropic BRDF using Jacobi polys
    # void brdf_iso_jacobi(const double *xrxs, const double *nu, const int num_pts, const int kmax, double *out)
    func = _liburdf.brdf_iso_jacobi
    func.restype = None
    func.argtypes = [ndpointer(ctypes.c_double), ndpointer(ctypes.c_double),ctypes.c_int, ctypes.c_int, ndpointer(ctypes.c_double)]

    return _liburdf

def brdf_iso_jacobi(xr,yr,xs,ys,sigma,kmax=None, debug=False):
    """
    xr,yr - Projected components of the reflected direction (xr**2 + yr**2 <= 1). May be array-like.
    xx,ys - Projected components of the specular direction (xs**2 + ys**2 <= 1). May be array-like.
    sigma - Standard deviation of the surface gradients. May be array-like
    [k=30]        - Largest term to include (valid from 0...30)
    [debug=False] - Whether to print timing information
    """
    c_kmax = 30
    if kmax is not None:
        if kmax > c_kmax:
            raise Exception('Cannot use kmax='+str(kmax)+' > 30 due to integer overflow')
        c_kmax = kmax

    num_pts = len(xr+xs+yr+ys+sigma)
    nu = empty(num_pts, dtype=float64)
    nu[:] = exp(-4*sigma*sigma)

    xrxs = empty((num_pts,4), dtype=float64)
    xrxs[:,0] = xr
    xrxs[:,1] = yr
    xrxs[:,2] = xs
    xrxs[:,3] = ys

    out = empty_like(nu)

    lib = initlib()
    t0 = time()
    lib.brdf_iso_jacobi(xrxs,nu, num_pts, c_kmax, out)
    t1 = time()
    if debug:
        print('Time taken', t1-t0, 'i.e {:,} evals/sec'.format(int(num_pts/(t1-t0))))
    return out
