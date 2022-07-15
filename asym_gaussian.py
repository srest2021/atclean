#!/usr/bin/env python
'''
A. Rest
'''

from scipy.interpolate import interp1d
from astropy.modeling.functional_models import Gaussian1D
import numpy as np

def mag2flux(mag):
    flux = 10**((mag-8.9)/-2.5)
    flux = flux * 1e6 # convert to micro Jy
    return flux

def skew_gaussian(peak,sig_minus,sig_plus):
    """
    Creates a skewed gaussian
    """
    x = np.arange(-100,100,.01)
    g1 = Gaussian1D(amplitude = peak,stddev=sig_minus)(x)
    g2 = Gaussian1D(amplitude = peak,stddev=sig_plus)(x)

    ind = np.argmin(abs(x))
    g3 = np.copy(g1)
    g3[ind:] = g2[ind:]
    gauss = np.array([x,g3])
    return gauss

def gauss2lc(time,peak_time,sig_minus,sig_plus,app_mag=None,abs_mag=None,dist=None):
    """
    Matches an asymetric gaussian to the time array presented.

    Inputs
    ------
        time : array
            time of the lightcurve
        peak_time : float
            time for the gaussian peak
        sig_minus : float
            sigma for gaussian rise 
        sig_plus : float
            sigma for gaussian fall
        app_mag : float
            peak apparent magnitude value of gaussian
        abs_mag : float
            peak absolute magnitude value of gaussian
        dist : float
            distance in Mpc to the source

    """
    if (abs_mag is not None) & (dist is not None):
        app_mag = abs_mag + 5*np.log10(dist*1e6/10)
    if app_mag is None:
        raise ValueError("Either app_mag or abs_mag + dist must be specified")
    peak_flux = mag2flux(app_mag)
    g = skew_gaussian(peak_flux,sig_minus,sig_plus)
    g[0,:] += peak_time
    
    gauss = interp1d(g[0],g[1],bounds_error=False,fill_value=0)
    gauss = gauss(time)
    return gauss 
    