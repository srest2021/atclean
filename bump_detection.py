#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys,random,re,os
from copy import deepcopy
from scipy.interpolate import interp1d
from astropy.modeling.functional_models import Gaussian1D

from pdastro import pdastrostatsclass
from atlas_lc import atlas_lc

# suppress deprecation warnings
import warnings
warnings.simplefilter('error', RuntimeWarning)
warnings.filterwarnings("ignore")

"""
SETTINGS
"""

# SN TNS name
tnsname = ''

# SN discovery date
#discovery_date = None

# path to directory that contains SN, control, and other light curves
source_dir = ''

# path to directory where generated tables should be stored
tables_dir = f'{source_dir}/{tnsname}/bump_analysis/tables'

# OPTIONAL: path to text file with light curve of simulated eruption to add
erup_filename = None #f'{source_dir}/{tnsname}/bump_analysis/eruption_m3e-07.dat'

# number of control light curves to load
num_controls = 8

# filter of light curve to analyze
filt = 'o'

# MJD bin size in days of light curve to analyze
mjd_bin_size = 1.0

# search for pre-SN bumps with the following gaussian sigmas
gauss_sigmas = [5, 15, 25, 40, 80, 130, 200, 256, 300]

# select sets of gaussian sigmas to simulate
# where each list corresponds to its matching entry in gauss_sigmas
# if using a simulated eruption, add None entry
sim_sigmas = [[2, 5, 20, 40, 80, 120], # 5
              [2, 5, 20, 40, 80, 120], # 15
              [2, 5, 20, 40, 80, 120], # 15
              [5, 20, 40, 80, 110, 150, 200, 250], # 40
              [5, 20, 40, 80, 110, 150, 200, 250], # 80
              [20, 40, 80, 110, 150, 200, 250, 300], # 130
              [20, 40, 80, 110, 150, 200, 250, 300], # 200
              [20, 40, 80, 110, 150, 200, 250, 300], # 256
              [20, 40, 80, 110, 150, 200, 250, 300]] # 300

# OPTIONAL: select FOM limits to calculate efficiencies 
# where each list corresponds to its matching entry in gauss_sigmas
# if using a simulated eruption, add None entry
fom_limits = [[3.0, 5.0], [10.4, 17.4], [15.1, 25.1], [20.0, 33.4], [25.3, 42.2], [29.0, 48.4], [33.1, 55.2]]

# select range of peak apparent magnitudes to simulate
peak_mag_max = 16 # brightest magnitude
peak_mag_min = 23 # faintest magnitude
n_peaks = 20 # number of magnitudes to generate in log space

# number of iterations of random sigma and peak mjd per peak
iterations = 50000 

# define flag that defines bad measurements (to filter out bad days in lc)
# if more than one, use bitwise OR (this symbol: |) to combine them
flags = 0x800000

# add observation seasons' mjd ranges here
valid_mjd_ranges = [[57365,57622],
                    [57762,57983], 
                    [58120,58383], 
                    [58494,58741-41], 
                    [58822+28,59093], 
                    [59184,59445], 
                    [59566,59835], 
                    [59901,60085]]

# skip any messy control light curves (leave empty list [] if not skipping any)
skip_ctrl = []

"""
UTILITY
"""

def AandB(A,B):
    return np.intersect1d(A,B,assume_unique=False)

def AnotB(A,B):
    return np.setdiff1d(A,B)

def AorB(A,B):
    return np.union1d(A,B)

def not_AandB(A,B):
    return np.setxor1d(A,B)

# Adapted from A. Rest
# get all indices of a dataframe
def get_ix(df):
    return df.index.values

# Adapted from A. Rest
def ix_equals(df, colname, val, indices=None):
    if indices is None:
        ix = get_ix(df)
    else:
        ix = indices

    (keep,) = np.where(df.loc[ix,colname].eq(val))
    ix = ix[keep] 
    
    return ix 

# Adapted from A. Rest
def ix_inrange(df, colname, lowlim=None, uplim=None, indices=None, exclude_lowlim=False, exclude_uplim=False):
    if indices is None:
        ix = get_ix(df)
    else:
        ix = indices

    if not(lowlim is None):
        if exclude_lowlim:
            (keep,) = np.where(df.loc[ix,colname].gt(lowlim))
        else:
            (keep,) = np.where(df.loc[ix,colname].ge(lowlim))
        ix = ix[keep]

    if not(uplim is None):
        if exclude_uplim:
            (keep,) = np.where(df.loc[ix,colname].lt(uplim))
        else:
            (keep,) = np.where(df.loc[ix,colname].le(uplim))
        ix = ix[keep]
        
    return ix 

# convert flux to magnitude 
def flux2mag(flux):
    return -2.5 * np.log10(flux) + 23.9

# convert magnitude to flux
def mag2flux(mag):
    return 10 ** ((mag - 23.9) / -2.5)

"""
ASYM GAUSSIAN
Adapted from A. Rest
"""

class Gaussian:
    def __init__(self, sigma, peak_appmag):
        self.peak_appmag = peak_appmag
        self.sigma = sigma
        self.g = self.new_gaussian(mag2flux(peak_appmag), sigma)

    def new_gaussian(self, peak_flux, sigma):
        x = np.arange(-100,100,.01)
        g1 = Gaussian1D(amplitude=peak_flux, stddev=sigma)(x)
        g2 = Gaussian1D(amplitude=peak_flux, stddev=sigma)(x)

        ind = np.argmin(abs(x))
        g3 = np.copy(g1)
        g3[ind:] = g2[ind:]
        gauss = np.array([x,g3])
        return gauss
    
    # get interpolated function of gaussian at peak MJD (peak_mjd) and match to time array (mjds)
    def gauss2fn(self, mjds, peak_mjd):
        g = deepcopy(self.g)
        g[0,:] += peak_mjd
        
        # interpolate gaussian
        fn = interp1d(g[0],g[1],bounds_error=False,fill_value=0)
        fn = fn(mjds)
        return fn 
    
"""
SIMULATED ERUPTION FROM LIGHT CURVE
"""

class Eruption:
    def __init__(self, filename, sigma=2.8):
        self.peak_appmag = None
        self.sigma = sigma

        self.t = None 
        self.load(filename)

    def load(self, filename):
        print(f'Loading eruption lc at {filename}...')

        try:
            self.t = pd.read_table(filename,delim_whitespace=True,header=None)
            self.t = self.t.rename(columns={0: "MJD", 1: "m"})
        except Exception as e:
            raise RuntimeError(f'ERROR: Could not load eruption at {filename}: {str(e)}')
    
        # app mag -> flux
        self.t['uJy'] = self.t['m'].apply(lambda mag: mag2flux(mag))

    # get interpolated function of eruption light curve
    def erup2fn(self, mjds, peak_mjd, peak_appmag):
        #self.peak_mjd = peak_mjd
        self.peak_appmag = peak_appmag

        peak_idx = self.t['m'].idxmin() # get peak appmag

        self.t['MJD'] -= self.t.loc[peak_idx,'MJD'] # put peak appmag at days=0
        self.t['MJD'] += peak_mjd # put peak appmag at days=peak_mjd

        # scale
        self.t['uJy'] *= mag2flux(peak_appmag)/self.t.loc[peak_idx, 'uJy']
        
        # flux -> app mag
        self.t['m'] = self.t['uJy'].apply(lambda flux: flux2mag(flux)) 
        
        # interpolate lc
        fn = interp1d(self.t['MJD'], self.t['uJy'], bounds_error=False, fill_value=0)
        fn = fn(mjds)
        return fn
    
"""
APPLY ROLLING SUM TO LIGHT CURVE 
"""
    
class LightCurve(atlas_lc):
    def __init__(self, tnsname, mjd_bin_size=1.0, discdate=None):
        atlas_lc.__init__(self, tnsname=tnsname, mjd_bin_size=mjd_bin_size, discdate=discdate, is_averaged=True)
        self.sigma = None

    # load a single averaged light curve without overwriting the 'Mask' column (override)
    def load_lc(self, source_dir, filt, control_index=0):
        if control_index == 0:
            filename = f'{source_dir}/{tnsname}/{tnsname}.{filt}.{self.mjd_bin_size:0.2f}days.lc.txt'
        else:
            filename = f'{source_dir}/{tnsname}/controls/{tnsname}_i{control_index:03d}.{filt}.{self.mjd_bin_size:0.2f}days.lc.txt'

        self.lcs[control_index] = pdastrostatsclass()
        self.lcs[control_index].load_spacesep(filename, delim_whitespace=True)

    # load averaged SN and control light curves (override)
    def load(self, source_dir, filt, num_controls):
        print(f'\nLoading averaged SN light curve and {num_controls} averaged control light curves...')
        
        self.load_lc(source_dir, filt)
        for control_index in range(1, num_controls+1):
            self.load_lc(source_dir, filt, control_index=control_index)

    # clean rolling sum columns
    def clean_rolling_sum(self, control_index):
        dropcols = []
        for col in ['__tmp_SN','SNR','SNRsum','SNRsumnorm']:
            if col in self.lcs[control_index].t.columns:
                dropcols.append(col)
        if len(dropcols) > 0:
            self.lcs[control_index].t.drop(columns=dropcols,inplace=True)

    # clean simulation columns
    def clean_simulations(self, lc):
        dropcols = []
        for col in ['__tmp_SN','uJysim','SNRsim','simLC','SNRsimsum']:
            if col in lc.t.columns:
                dropcols.append(col)
        if len(dropcols) > 0:
            lc.t.drop(columns=dropcols,inplace=True)
        return lc
    
    def apply_rolling_sum(self, control_index, sigma, flag=0x800000, verbose=False):
        self.sigma = sigma

        ix = get_ix(self.lcs[control_index].t)
        if len(ix) < 1:
            raise RuntimeError('ERROR: not enough measurements to apply simulated gaussian')
        good_ix = AandB(ix, self.lcs[control_index].ix_unmasked('Mask', flag)) # all good pre-SN indices
        
        self.clean_rolling_sum(control_index)
        self.lcs[control_index].t.loc[ix, 'SNR'] = 0.0
        self.lcs[control_index].t.loc[good_ix,'SNR'] = self.lcs[control_index].t.loc[good_ix,'uJy']/self.lcs[control_index].t.loc[good_ix,'duJy']

        new_gaussian_sigma = round(self.sigma/self.mjd_bin_size)
        windowsize = int(6 * new_gaussian_sigma)
        halfwindowsize = int(windowsize * 0.5) + 1
        if verbose:
            print(f'# Sigma: {self.sigma:0.2f} days; MJD bin size: {self.mjd_bin_size:0.2f} days; sigma: {new_gaussian_sigma:0.2f} bins; window size: {windowsize} bins')

        # calculate the rolling SNR sum
        l = len(self.lcs[control_index].t.loc[ix])
        dataindices = np.array(range(l) + np.full(l, halfwindowsize))
        temp = pd.Series(np.zeros(l + 2*halfwindowsize), name='SNR', dtype=np.float64)
        temp[dataindices] = self.lcs[control_index].t.loc[ix,'SNR']
        SNRsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
        self.lcs[control_index].t.loc[ix,'SNRsum'] = list(SNRsum[dataindices])
        
        # normalize it
        norm_temp = pd.Series(np.zeros(l + 2*halfwindowsize), name='norm', dtype=np.float64)
        norm_temp[np.array(range(l) + np.full(l, halfwindowsize))] = np.ones(l)
        norm_temp_sum = norm_temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
        self.lcs[control_index].t.loc[ix,'SNRsumnorm'] = list(SNRsum.loc[dataindices] / norm_temp_sum.loc[dataindices] * max(norm_temp_sum.loc[dataindices]))

    def apply_rolling_sums(self, sigma, num_controls):
        print(f'Applying rolling sum to {num_controls} control light curves...')
        for control_index in range(1, num_controls+1):
            self.apply_rolling_sum(control_index, sigma)

    def add_simulation(self, control_index, peak_mjd=None, gaussian=None, eruption=None, peak_appmag=None, flag=0x800000, verbose=False):
        if not(gaussian is None) and not(eruption is None):
            raise RuntimeError('ERROR: both gaussian and eruption passed, but only one allowed')
        elif gaussian is None and eruption is None:
            raise RuntimeError('ERROR: must pass either gaussian or eruption')
        elif not(gaussian is None):
            if peak_mjd is None:
                raise RuntimeError('ERROR: peak_mjd required for adding simulated gaussian')
            if not(peak_appmag is None):
                raise RuntimeError('ERROR: peak_appmag erroneously specified for simulated gaussian')
        elif not(eruption is None):
            if peak_mjd is None or peak_appmag is None:
                raise RuntimeError('ERROR: peak_mjd and peak_appmag required for adding simulated eruption')
        
        lc = deepcopy(self.lcs[control_index])
        lc = self.clean_simulations(lc)
        
        ix = get_ix(lc.t)
        if len(ix) < 1:
            raise RuntimeError('ERROR: not enough measurements to apply simulated gaussian')
        good_ix = AandB(ix, lc.ix_unmasked('Mask', flag)) # all good pre-SN indices

        if verbose:
            if not(gaussian is None):
                print(f'# Adding simulated gaussian: peak MJD = {peak_mjd:0.2f} MJD; peak app mag = {gaussian.peak_appmag:0.2f}; sigma = {gaussian.sigma:0.2f} days')
            else:
                print(f'# Adding simulated eruption: peak MJD = {peak_mjd:0.2f} MJD; peak app mag = {peak_appmag:0.2f}; sigma = {eruption.sigma:0.2f} days')

        lc.t.loc[good_ix,'uJysim'] = lc.t.loc[good_ix,'uJy']
        if not(gaussian is None):
            simflux = gaussian.gauss2fn(lc.t.loc[good_ix,'MJD'], peak_mjd) # get simulated gaussian flux
        else:
            simflux = eruption.erup2fn(lc.t.loc[good_ix,'MJD'], peak_mjd, peak_appmag) # get simulated eruption flux
        lc.t.loc[good_ix,'uJysim'] += simflux # add simulated flux to good indices

        # make sure all bad rows have SNRsim = 0.0 so they have no impact on the rolling SNRsum
        lc.t.loc[ix,'SNRsim'] = 0.0
        # include simflux in the SNR
        lc.t.loc[good_ix,'SNRsim'] = lc.t.loc[good_ix,'uJysim']/lc.t.loc[good_ix,'duJy']

        new_gaussian_sigma = round(self.sigma/self.mjd_bin_size)
        windowsize = int(6 * new_gaussian_sigma)
        halfwindowsize = int(windowsize * 0.5) + 1
        if verbose:
            print(f'# Sigma: {self.sigma:0.2f} days; MJD bin size: {self.mjd_bin_size:0.2f} days; new sigma: {new_gaussian_sigma:0.2f} bins; window size: {windowsize} bins')

        # calculate the rolling SNR sum for SNR with simulated flux
        l = len(self.lcs[control_index].t.loc[ix])
        dataindices = np.array(range(l) + np.full(l, halfwindowsize))
        temp = pd.Series(np.zeros(l + 2*halfwindowsize), name='SNRsim', dtype=np.float64)
        temp[dataindices] = lc.t.loc[ix,'SNRsim']
        SNRsimsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
        lc.t.loc[ix,'SNRsimsum'] = list(SNRsimsum.loc[dataindices])

        return lc
        
"""
SIMULATION DETECTION AND EFFICIENCY TABLES
"""

# get simulation detection dictionary (sd) key for a gauss_sigma peak_appmag pair
def sd_key(gauss_sigma, peak_appmag):
    return f'{gauss_sigma}_{peak_appmag:0.2f}'

# simulation detection table for a gauss_sigma peak_appmag pair
class SimDetecTable:
    def __init__(self, gauss_sigma, iterations=None, peak_appmag=None, peak_flux=None):
        self.gauss_sigma = gauss_sigma
        
        if peak_appmag is None and peak_flux is None:
            raise RuntimeError('ERROR: either peak app mag or flux required to construct simulation detection table')
        if peak_appmag is None:
            self.peak_flux = float(peak_flux)
            self.peak_appmag = flux2mag(self.peak_flux)
        else:
            self.peak_appmag = float(peak_appmag)
            self.peak_flux = mag2flux(self.peak_appmag)

        self.t = None

        if not iterations is None:
            self.setup(iterations)

    def get_filename(self, tables_dir):
        return f'{tables_dir}/simdetec_{sd_key(self.gauss_sigma, self.peak_appmag)}.txt'

    def load(self, tables_dir):
        filename = self.get_filename(tables_dir)
        print(f'Loading simulation detection table simdetec_{sd_key(self.gauss_sigma, self.peak_appmag)}.txt...')

        try:
            self.t = pd.read_table(filename,delim_whitespace=True)
        except Exception as e:
            raise RuntimeError(f'ERROR: Could not load simulation detection table at {filename}: {str(e)}')

    def save(self, tables_dir):
        filename = self.get_filename(tables_dir)
        print(f'Saving simulation detection table simdetec_{sd_key(self.gauss_sigma, self.peak_appmag)}.txt...')
        self.t.to_string(filename, index=False)

    def setup(self, iterations):
        self.t = pd.DataFrame(columns=['gauss_sigma', 'peak_appmag', 'peak_flux', 'peak_mjd', 'sim_gauss_sigma', 'sim_erup_sigma', 'max_fom', 'max_fom_mjd'])
        self.t['gauss_sigma'] = np.full(iterations, self.gauss_sigma)
        self.t['peak_appmag'] = np.full(iterations, self.peak_appmag)
        self.t['peak_flux'] = np.full(iterations, self.peak_flux)
    
    def get_efficiency(self, fom_limit, mjd_ranges=None, sim_sigma=None, erup=False):
        if sim_sigma is None: # get efficiency for all sim_sigmas (entire table)
            sim_sigma_ix = get_ix(self.t)
        elif pd.isnull(sim_sigma):
            return np.nan
        elif erup: # simulated eruption
            sim_sigma_ix = np.where(self.t['sim_erup_sigma'] == sim_sigma)[0]
        else: # simulated gaussian
            sim_sigma_ix = np.where(self.t['sim_gauss_sigma'] == sim_sigma)[0]
        if len(sim_sigma_ix) < 1:
            print(f'WARNING: No sim sigma matches for {sim_sigma}, where erup={erup}; returning NaN...')
            return np.nan
        
        if not mjd_ranges is None:
            # get efficiency for simulations only occurring within mjd_ranges
            mjd_ix = []
            for i in range(len(self.t)):
                if in_valid_season(self.t.loc[i,'peak_mjd'], mjd_ranges):
                    mjd_ix.append(i)
            sim_sigma_ix = AandB(sim_sigma_ix, mjd_ix)
        
        detected_ix = ix_inrange(self.t, 'max_fom', lowlim=fom_limit, indices=sim_sigma_ix)
        efficiency = len(detected_ix)/len(sim_sigma_ix) * 100
        return efficiency
    
def load_sd_dict(gauss_sigmas, peak_appmags, tables_dir):
    sd = {}
    for i in range(len(gauss_sigmas)):
        gauss_sigma = gauss_sigmas[i]
        for peak_index in range(len(peak_appmags)):
            peak_appmag = peak_appmags[peak_index]
            key = sd_key(gauss_sigma, peak_appmag)
            sd[key] = SimDetecTable(gauss_sigma=gauss_sigma, peak_appmag=peak_appmag)
            sd[key].load(tables_dir)
    return sd

class EfficiencyTable:
    def __init__(self, gauss_sigmas, peak_appmags, peak_fluxes, sim_sigmas, sim_erup_sigma=None):
        if not(len(gauss_sigmas) == len(sim_sigmas)):
            raise RuntimeError('ERROR: Each entry in gauss_sigmas must have a matching list in sim_sigmas')
        self.gauss_sigmas = gauss_sigmas
        self.sim_sigmas = sim_sigmas
        self.fom_limits = None

        self.peak_appmags = peak_appmags
        self.peak_fluxes = peak_fluxes

        self.t = None
        #self.sd = None

        self.setup(sim_erup_sigma=sim_erup_sigma)

    def setup(self, sim_erup_sigma=None):
        self.t = pd.DataFrame(columns=['gauss_sigma', 'peak_appmag', 'peak_flux', 'sim_gauss_sigma', 'sim_erup_sigma'])

        for i in range(len(self.gauss_sigmas)):
            gauss_sigma = self.gauss_sigmas[i]
            n = len(self.peak_appmags) * len(self.sim_sigmas[i])
            
            df = pd.DataFrame(columns=['gauss_sigma', 'peak_appmag', 'peak_flux', 'sim_gauss_sigma', 'sim_erup_sigma'])
            df['gauss_sigma'] = np.full(n, gauss_sigma)
            df['peak_appmag'] = np.repeat(self.peak_appmags, len(self.sim_sigmas[i]))
            df['peak_flux'] = np.repeat(self.peak_fluxes, len(self.sim_sigmas[i]))
            
            j = 0
            while(j < n):
                for sim_sigma in self.sim_sigmas[i]:
                    if sim_sigma is None:
                        if sim_erup_sigma is None:
                            raise RuntimeError('ERROR: None element in sim_sigmas array, but no sim_erup_sigma passed')
                        df.loc[j, 'sim_erup_sigma'] = sim_erup_sigma
                    else: 
                        df.loc[j, 'sim_gauss_sigma'] = sim_sigma
                    j += 1
            self.t = pd.concat([self.t, df], ignore_index=True)
    
    def load(self, tables_dir, filename):
        print(f'Loading efficiency table {filename}...')
        filename = f'{tables_dir}/{filename}'
        try:
            self.t = pd.read_table(filename,delim_whitespace=True)
        except Exception as e:
            raise RuntimeError(f'ERROR: Could not load efficiency table at {filename}: {str(e)}')

    def save(self, tables_dir, filename=None):
        if filename is None:
            print(f'Saving efficiency table efficiencies.txt...')
            filename = f'{tables_dir}/efficiencies.txt'
        else:
            print(f'Saving efficiency table {filename}...')
            filename = f'{tables_dir}/{filename}'
        self.t.to_string(filename, index=False)

    def reset(self):
        for col in self.t.columns:
            if re.search('^pct_detec_',col):
                self.t.drop(col, axis=1, inplace=True)

        for i in range(len(self.gauss_sigmas)):
            fom_limits = self.fom_limits[self.gauss_sigmas[i]]
            self.t[f'pct_detec_{fom_limits[0]}'] = np.full(len(self.t), np.nan)
            self.t[f'pct_detec_{fom_limits[1]}'] = np.full(len(self.t), np.nan)

    def set_fom_limits(self, fom_limits):
        if not len(self.gauss_sigmas) == len(self.sim_sigmas) or not len(fom_limits) == len(self.gauss_sigmas):
            raise RuntimeError('ERROR: Each entry in gauss_sigmas must have a matching list in fom_limits')
        
        self.fom_limits = {}
        for i in range(len(self.gauss_sigmas)):
            self.fom_limits[self.gauss_sigmas[i]] = fom_limits[i]

    def get_fom_limit(self, gauss_sigma, sigma):
        if sigma == 3:
          return self.fom_limits[gauss_sigma][0]
        elif sigma == 5:
          return self.fom_limits[gauss_sigma][1]
        else:
          raise RuntimeError('ERROR: when obtaining an FOM limit for given gauss_sigma, sigma must be 3 or 5')

    def get_efficiencies(self, sd, mjd_ranges=None):
        if self.fom_limits is None:
            raise RuntimeError('ERROR: fom_limits is None')
        
        for i in range(len(self.t)):
            gauss_sigma = self.t.loc[i,'gauss_sigma']
            
            if pd.isnull(self.t.loc[i,'sim_gauss_sigma']):
                sim_sigma = self.t.loc[i,'sim_erup_sigma']
                erup = True
            else:
                sim_sigma = self.t.loc[i,'sim_gauss_sigma']
                erup = False

            for fom_limit in self.fom_limits[gauss_sigma]:
                self.t.loc[i,f'pct_detec_{fom_limit}'] = \
                    sd[sd_key(gauss_sigma, self.t.loc[i,'peak_appmag'])].get_efficiency(fom_limit, 
                                                                                        mjd_ranges=mjd_ranges,
                                                                                        sim_sigma=sim_sigma, 
                                                                                        erup=erup)

    def get_subset(self, gauss_sigma=None, fom_limit=None, sim_sigma=None, erup=False):
        colnames = ['gauss_sigma', 'peak_appmag', 'peak_flux']
        
        ix = get_ix(self.t)
        if not(gauss_sigma is None):
            ix = ix_equals(self.t, 'gauss_sigma', gauss_sigma) 

        if sim_sigma is None:
            colnames += ['sim_gauss_sigma', 'sim_erup_sigma']
        else: 
            if erup:
                sim_sigma_colname = 'sim_erup_sigma'
            else:
                sim_sigma_colname = 'sim_gauss_sigma'
            colnames.append(sim_sigma_colname)
            ix = ix_equals(self.t, sim_sigma_colname, sim_sigma, indices=ix)

        if not(fom_limit is None):
            try: 
                colnames.append(f'pct_detec_{fom_limit}')
            except Exception as e:
                raise RuntimeError(f'ERROR: no matching FOM limit {fom_limit}: {str(e)}')
        else:
            for col in self.t.columns:
                if re.search('^pct_detec_',col):
                    colnames.append(col)
        
        return self.t.loc[ix,colnames]
    
    def merge(self, other):
        if not isinstance(other, EfficiencyTable):
            raise RuntimeError(f'ERROR: Cannot merge with object: {type(other)}')
        
        self.gauss_sigmas += other.gauss_sigmas
        self.sim_sigmas += other.sim_sigmas
        if not self.fom_limits is None:
            self.fom_limits.update(other.fom_limits)

        self.t = pd.concat([self.t, other.t], ignore_index=True)


"""
MISCELLANEOUS
"""

# check if MJD is within valid MJD season
def in_valid_season(mjd, valid_seasons):
    in_season = False
    for season in valid_seasons:
        if mjd <= season[1] and mjd >= season[0]:
            in_season = True
    return in_season

# generate list of peak fluxes and app mags
def generate_peaks(peak_mag_min, peak_mag_max, n_peaks):
    peak_mags = list(np.linspace(peak_mag_min, peak_mag_max, num=n_peaks))
    peak_fluxes = list(map(mag2flux, peak_mags))

    peak_mags = [round(item, 2) for item in peak_mags]
    peak_fluxes = [round(item, 2) for item in peak_fluxes]

    return peak_mags, peak_fluxes

"""
GENERATE AND SAVE SIMULATED DETECTION AND EFFICIENCY TABLES
"""

if __name__ == "__main__":
    # load SN and control light currves
    lc = LightCurve(tnsname, mjd_bin_size=mjd_bin_size)
    lc.load(source_dir, filt=filt, num_controls=num_controls)

    # load simulated eruption
    if not(erup_filename is None):
        erup = Eruption(erup_filename)

    # check if tables_dir exists; if not, create new
    if not os.path.exists(tables_dir):
        os.mkdir(tables_dir)

    # generate list of peaks
    peak_appmags, peak_fluxes = generate_peaks(peak_mag_min, peak_mag_max, n_peaks)
    
    # print settings
    print(f'\nRolling sum sigmas in days: ', gauss_sigmas)
    sim_erup_sigma=None
    for i in range(len(gauss_sigmas)):
        gauss_sigma = gauss_sigmas[i]
        print(f'For rolling sum sigma {gauss_sigma}, simulating: ')
        for sim_sigma in sim_sigmas[i]:
            if sim_sigma is None:
                sim_erup_sigma = erup.sigma
                print(f'- eruption from file with sigma {erup.sigma}')
            else:
                print(f'- gaussian with sigma {sim_sigma}')
    print(f'\nPeak magnitudes: ', peak_appmags)
    print(f'Peak fluxes (uJy): ', peak_fluxes)
    print(f'Number of iterations per peak: {iterations}')
    print(f'Flag for bad days: {hex(flags)}')
    print(f'\nSimulating eruptions only within the following MJD seasons: ', valid_mjd_ranges)
    valid_ctrls = [i for i in range(1, num_controls+1) if not i in skip_ctrl]
    print(f'Simulating eruptions only within the following control light curves: ', valid_ctrls)
    
    # construct blank dictionary of simdetec tables
    sd = {}

    # construct blank efficiency table
    e = EfficiencyTable(gauss_sigmas, peak_appmags, peak_fluxes, sim_sigmas, sim_erup_sigma=sim_erup_sigma)

    for i in range(len(gauss_sigmas)):
        gauss_sigma = gauss_sigmas[i]
        print(f'\nUsing rolling sum sigma of {gauss_sigma} days...')

        # apply rolling sum to all control light curves
        lc.apply_rolling_sums(gauss_sigma, num_controls)

        for peak_index in range(len(peak_appmags)):
            peak_appmag = peak_appmags[peak_index]
            peak_flux = peak_fluxes[peak_index]

            # construct new simulation detection table
            key = sd_key(gauss_sigma, peak_appmag)
            sd[key] = SimDetecTable(gauss_sigma=gauss_sigma, iterations=iterations, peak_appmag=peak_appmag)

            # construct list of gaussians for each sim_gauss_sigma
            gaussians = []
            for k in range(len(sim_sigmas[i])):
                if not(sim_sigmas[i][k] is None):
                    gaussians.append(Gaussian(sim_sigmas[i][k], peak_appmag))
            
            print(f'Commencing {iterations} iterations for peak app mag {peak_appmag} (peak flux {peak_flux})...')
            j = 0
            while j < iterations:
                # reset row in case previous iteration failed
                sd[key].t.loc[j, ['peak_mjd', 'sim_gauss_sigma', 'sim_erup_sigma', 'max_fom', 'max_fom_mjd']] = np.full(5, np.nan)

                # pick random control light curve
                rand_control_index = random.choice(valid_ctrls)
                
                # select random peak MJD from start of lc to 50 days before discovery date
                peak_mjd = random.randrange(lc.lcs[rand_control_index].t['MJDbin'].iloc[0]-0.5, int(lc.discdate)-50, 1) + 0.5
                # make sure peak MJD is within an observation season; else redraw
                while not in_valid_season(peak_mjd, valid_mjd_ranges):
                    # redraw random peak mjd
                    peak_mjd = random.randrange(lc.lcs[rand_control_index].t['MJDbin'].iloc[0]-0.5, int(lc.discdate)-50, 1) + 0.5
                sd[key].t.loc[j, 'peak_mjd'] = peak_mjd

                # select random sim sigma
                k = random.randrange(0, len(sim_sigmas[i]), 1)
                if sim_sigmas[i][k] is None:
                    # add simulated eruption
                    sim_sigma = erup.sigma
                    sd[key].t.loc[j, 'sim_erup_sigma'] = sim_sigma
                    sim_lc = lc.add_simulation(rand_control_index, peak_mjd, eruption=erup, peak_appmag=peak_appmag)
                else:
                    # add simulated gaussian
                    sim_sigma = gaussians[k].sigma
                    sd[key].t.loc[j, 'sim_gauss_sigma'] = sim_sigma
                    sim_lc = lc.add_simulation(rand_control_index, peak_mjd, gaussians[k])

                # get max FOM from measurements within 1 sigma of the simulated bump
                sigma_ix = sim_lc.ix_inrange(colnames='MJDbin', lowlim=peak_mjd-sim_sigma, uplim=peak_mjd+sim_sigma)
                sigma_ix = sim_lc.ix_not_null(colnames='MJD', indices=sigma_ix)
                if len(sigma_ix) > 0:
                    max_fom_idx = sim_lc.t.loc[sigma_ix,'SNRsimsum'].idxmax()
                    sd[key].t.loc[j, 'max_fom'] = sim_lc.t.loc[max_fom_idx, 'SNRsimsum']
                    sd[key].t.loc[j, 'max_fom_mjd'] = sim_lc.t.loc[max_fom_idx, 'MJD']
                else:
                    # no valid measurements within this range
                    continue

                j += 1

            sd[key].save(tables_dir)

    print('\nSuccess')

    if not(fom_limits is None):
        print(f'\nUsing FOM limits {fom_limits} to calculate efficiencies...')
        e.set_fom_limits(fom_limits)
        e.get_efficiencies(sd)
        print(e.t.to_string())
    e.save(tables_dir)