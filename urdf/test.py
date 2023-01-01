from .lib import brdf_iso_jacobi
from numpy import sin,cos,mgrid,flatnonzero, zeros,pi
import numpy as np
from numpy import *

def test_speed():
    
    npts =300
    dx = 1.0/npts
    yr, xr = mgrid[:npts,:npts] * dx 
    r2 = yr*yr+xr*xr
    idx = flatnonzero(r2.ravel()<1.0)
    xr = xr.ravel()[idx]
    yr = yr.ravel()[idx]

    phi_s = pi/4.0
    theta = 75 * pi/180.0
    print('sin(theta)', sin(theta))
    xs = cos(phi_s)*sin(theta)
    ys = sin(phi_s)*sin(theta)

    sig = 0.1

    for kmax in [3,8,20,3]:
        vals = brdf_iso_jacobi(xr,yr,xs,ys,sig, kmax=kmax, debug=True)
        print('Min val', vals.min())
        print('Total', vals.sum()*dx*dx)

def test_vals():
    """
    Some pre-computed values for the isotropic BRDF. If the 
    library starts to give different values for these, we
    probably want to know about it.
    """
    eps = 1e-8
    xr = array([0.5, 0.1])
    yr = array([0.5, 0.5])
    xi = array([0.5, 0.5])
    yi = array([0.5, 0.5])

    print('Testing k=0')
    x = brdf_iso_jacobi(xr,yr,xi,yr,0.1, kmax=0)
    ans = [ 2.03963153,  0.69381664]
    print('ans', repr(x), 'should be', ans)
    assert(abs(x-ans).max()<eps)

    print('Testing k=1')
    x = brdf_iso_jacobi(xr,yr,xi,yr,0.1, kmax=1)
    ans = [ 5.41271927,  0.74768065]
    print('ans', repr(x), 'should be', ans)
    assert(abs(x-ans).max()<eps)

    # Test the k=20
    print('Testing k=20')
    
    x = brdf_iso_jacobi(xr,yr,xi,yr,0.1, kmax=20)

    ans = [8.23095479, 0.28596027]
    print('ans', x, 'should be', ans)
    assert(abs(x-ans).max()<eps)

    
if __name__=='__main__':

    test_vals() # Assert that the library gives some values we think it should
    test_speed()
