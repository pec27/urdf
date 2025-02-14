"""
Copyright (C) 2015-2025 Peter Creasey

Draw an NxN grid of BRDF values

Requires matplotlib
"""
import numpy as np
import pylab as pl
from urdf.lib import brdf_iso_jacobi

N=300

# NxN grid of projected r
yr, xr = (np.mgrid[:N,:N] + 0.5) * (2 / (N-1)) - 1 
# Projected r only within |r| <= 1
idx = np.flatnonzero((yr*yr+xr*xr).ravel()<1.0)
xr = xr.ravel()[idx]
yr = yr.ravel()[idx]

# Specular direction
theta = 30 * np.pi/180.0 # Angle from normal
xs = np.sin(theta)
ys = 0

# Gradient of surface deviations
sigma = 0.2

# BRDF
f = np.zeros(N*N)
f[idx] = brdf_iso_jacobi(xr,yr,xs,ys, sigma, kmax=10)
f = np.reshape(f, (N,N))

# Plot
pl.imshow(f)
pl.show()


