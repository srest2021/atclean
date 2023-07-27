#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys,copy,random,math
from scipy.interpolate import interp1d
from pdastro import pdastrostatsclass
from asym_gaussian import gauss2lc
import warnings
warnings.simplefilter('error', RuntimeWarning)
warnings.filterwarnings("ignore")

######### SETTINGS ######### 

# SN TNS name
tnsname = '2023ixf'

# SN discovery date
discovery_date = 60063.727257 - 20.0

# path to directory that contains SN, control, and other light curves
source_dir = '/Users/sofiarest/Desktop/Supernovae/data/misc/2023ixf'

# path to directory where generated tables should be stored
tables_dir = source_dir + '/bump_analysis/tables_test'

# OPTIONAL: path to text file with light curve of simulated eruption to add
erup_filename = source_dir + '/bump_analysis/eruption_m3e-07.dat'

# number of control light curves to load
n_controls = 16

# filter of light curve to analyze
filt = 'o'

# MJD bin size in days of light curve to analyze
mjd_bin_size = 1.0

# search for pre-SN bumps with the following gaussian sigmas
gauss_sigmas = [5, 80, 200]

# select sets of gaussian sigmas to simulate, where each list corresponds to its matching entry in gauss_sigmas
sim_gauss_sigmas = [[2, 5, 13], [30, 80, 120], [150, 200, 250, 300]]

# select range of peak apparent magnitudes to simulate
peak_mag_max = 16 # brightest magnitude
peak_mag_min = 22 # faintest magnitude
n_peaks = 20 # number of magnitudes to generate in log space

# number of iterations of random sigma and peak mjd per peak
iterations = 100 #50000 

# define flag that defines bad measurements (to filter out bad days in lc)
# if more than one, use bitwise OR (this symbol: |) to combine them
flags = 0x800000

# add observation seasons' mjd ranges here
valid_mjd_ranges = [[57463,57622], [57762,57983], [58120,58383], [58494,58741], [58822,59093], [59184,59445], [59566,59835], [59901,60085]]

############################

# for storing lc and control lcs 
global lcs
lcs = {}

# for storing info about the lc
global lc_info
lc_info = {}

# for storing per-peak tables
global detected
detected = {}

# final efficiencies table
efficiencies = pdastrostatsclass(columns=['gauss_sigma','peak_flux', 'peak_appmag'])

def AandB(A,B):
    return np.intersect1d(A,B,assume_unique=False)

def AnotB(A,B):
    return np.setdiff1d(A,B)

def AorB(A,B):
    return np.union1d(A,B)

def not_AandB(A,B):
    return np.setxor1d(A,B)

# convert flux to magnitude 
def flux2mag(flux):
    return -2.5 * np.log10(flux) + 23.9

# convert magnitude to flux
def mag2flux(mag):
    return 10 ** ((mag - 23.9) / -2.5)

# get all indices of a dataframe
def get_ix(df):
    return df.index.values

# generate list of peak fluxes and app mags
def generate_peaks(peak_mag_min, peak_mag_max, n_peaks):
    peak_mags = list(np.linspace(peak_mag_min, peak_mag_max, num=20))
    peak_fluxes = list(map(mag2flux, peak_mags))

    peak_mags = [round(item, 2) for item in peak_mags]
    peak_fluxes = [round(item, 2) for item in peak_fluxes]

    return peak_mags, peak_fluxes

# load a single light curve
def load_lc(source_dir, control_index):
    if control_index == 0:
        filename = f'{source_dir}/{lc_info["tnsname"]}.{lc_info["filt"]}.{lc_info["mjdbinsize"]:0.2f}days.lc.txt'
        print(f'\nLoading light curve for SN {tnsname} at {filename} ...')
    else:
        filename = f'{source_dir}/controls/{lc_info["tnsname"]}_i{control_index:03d}.{lc_info["filt"]}.{lc_info["mjdbinsize"]:0.2f}days.lc.txt'
        print(f'Loading control light curve {control_index:03d} at {filename} ...')

    lcs[control_index] = pdastrostatsclass()
    try:
        lcs[control_index].load_spacesep(filename,delim_whitespace=True)
    except Exception as e:
        raise RuntimeError('ERROR: Could not load light curve at %s: %s' % (filename, str(e)))

# load control light curves
def load_control_lcs(source_dir, n_controls):
    for control_index in range(1,n_controls+1):
        load_lc(source_dir, control_index)

# load SN and control light curves
def load_all_lcs(source_dir, n_controls):
    load_lc(source_dir, 0)
    load_control_lcs(source_dir, n_controls)

def save_detected(gauss_sigma, peak):
    print(f'Saving detected table for gauss_sigma {gauss_sigma} and peak flux {peak}...')
    detected[f'{gauss_sigma}_{peak}'].write(filename=f'{tables_dir}/detected_{gauss_sigma}_{peak:0.2f}.txt')

def save_detected_tables(gauss_sigmas, peaks):
    for gauss_sigma in gauss_sigmas: 
        for peak in peaks:
            save_detected(gauss_sigma, peak)

def save_efficiencies(efficiencies):
    print(f'Saving incomplete efficiencies table...')
    efficiencies.write(filename=f'{tables_dir}/efficiencies.txt')

# check if mjd is within valid mjd season
def in_season(mjd, valid_mjd_ranges):
    for mjd_range in valid_mjd_ranges:
        if mjd >= mjd_range[0] and mjd <= mjd_range[1]:
            return True
    return False

# load simulated eruption from light curve
def load_erup(erup_filename):
    print(f'Loading eruption lc at {erup_filename}...')

    erup = pdastrostatsclass()
    try:
        erup.t = pd.read_table(erup_filename,delim_whitespace=True,header=None)
        erup.t = erup.t.rename(columns={0: "MJD", 1: "m"})
    except Exception as e:
        print(f'ERROR: Could not load eruption at {erup_filename}: {str(e)}')
        sys.exit()
    
    # app mag -> flux
    erup.t['uJy'] = 10 ** ((erup.t['m'] - 23.9) / -2.5)

    return erup

def get_erup_bounds(sim_peakmjd):
    lower_bound = sim_peakmjd - lc_info['erup_peakday']
    upper_bound = sim_peakmjd + (len(erup.t) - lc_info['erup_peakday'] - 1)
    return lower_bound, upper_bound

# get interpolated function of a given simulated eruption light curve
def erup2lc(mjds, erup, sim_peakmjd, sim_appmag):
    min_i = erup.t['m'].idxmin() # get peak appmag

    # store peak appmag's days since start of erup
    lc_info['erup_peakday'] = erup.t.loc[min_i,'MJD']

    erup.t['MJD'] -= erup.t.loc[min_i,'MJD'] # put peak appmag at days=0
    erup.t['MJD'] += sim_peakmjd # put peak appmag at days=peak_mjd

    # scale
    erup.t['uJy'] *= mag2flux(sim_appmag)/erup.t.loc[min_i, 'uJy']
    erup.t['m'] = -2.5 * np.log10(erup.t['uJy']) + 23.9
    
    # interpolate lc
    fn = interp1d(erup.t['MJD'], erup.t['uJy'], bounds_error=False, fill_value=0)
    fn = fn(mjds)
    
    return fn

# optionally add a simulated pre-eruption to the loaded light curve
# then apply a gaussian weighted rolling sum to the SNR
def apply_rolling_sum(lc, gaussian_sigma, sim_gauss=False, erup=None, sim_peakmjd=None, sim_appmag=None, sim_sigma=None, print_=False, flag=0x800000):
    ix = get_ix(lc.t)
    if len(ix) < 1:
        raise RuntimeError('ERROR: not enough measurements to apply simulated gaussian')
    good_ix = AandB(ix, lc.ix_unmasked('Mask', flag)) # all good pre-SN indices

    # make sure there are no lingering simulations
    dropcols=[]
    for col in ['uJysim','SNRsim','simLC','SNRsimsum']:
        if col in lc.t.columns:
            dropcols.append(col)
    if len(dropcols) > 0:
        lc.t.drop(columns=dropcols,inplace=True)
    
    lc.t.loc[ix, 'SNR'] = 0.0
    lc.t.loc[good_ix,'SNR'] = lc.t.loc[good_ix,'uJy']/lc.t.loc[good_ix,'duJy']

    # add simulated gaussian
    sim_erup = not(erup is None)
    if sim_gauss and sim_erup:
        raise RuntimeError('ERROR: Cannot add a simulated gaussian and a simulated eruption at the same time')
    if sim_gauss or sim_erup:
        if sim_gauss and print_:
            print(f'# Adding simulated gaussian: peak at MJD {sim_peakmjd:0.2f}; apparent magnitude {sim_appmag:0.2f}; sigma- and sigma+ of {sim_sigma:0.2f} days')
        if sim_erup and print_:
            print(f'# Adding simulated eruption: peak at MJD {sim_peakmjd:0.2f}, apparent magnitude {sim_appmag:0.2f}')
        
        mjds = lc.t.loc[good_ix,'MJD']
        lc.t.loc[good_ix,'uJysim'] = lc.t.loc[good_ix,'uJy']
        lc.t.loc[ix,'simLC'] = 0.0

        # get simulated gaussian flux and add to light curve flux
        if sim_gauss:
            simflux = gauss2lc(mjds, sim_peakmjd, sim_sigma, sim_sigma, app_mag=sim_appmag)
        if sim_erup:
            simflux = erup2lc(mjds, erup, sim_peakmjd, sim_appmag)
        lc.t.loc[good_ix,'uJysim'] += simflux
        # get the simulated lc for all MJDs
        simflux_all = gauss2lc(lc.t.loc[ix,'MJDbin'], sim_peakmjd, sim_sigma, sim_sigma, app_mag=sim_appmag)
        lc.t.loc[ix,'simLC'] += simflux_all

        # make sure all bad rows have SNRsim = 0.0 so they have no impact on the rolling SNRsum
        lc.t.loc[ix,'SNRsim'] = 0.0
        # include simflux in the SNR
        lc.t.loc[good_ix,'SNRsim'] = lc.t.loc[good_ix,'uJysim']/lc.t.loc[good_ix,'duJy']

    new_gaussian_sigma = round(gaussian_sigma/lc_info['mjdbinsize'])
    windowsize = int(6 * new_gaussian_sigma)
    halfwindowsize = int(windowsize * 0.5) + 1
    if print_:
        print(f'# Sigma: {gaussian_sigma:0.2f} days; MJD bin size: {lc_info["mjdbinsize"]:0.2f} days; sigma: {new_gaussian_sigma:0.2f} bins; window size: {windowsize} bins')

    # calculate the rolling SNR sum
    dataindices = np.array(range(len(lc.t.loc[ix])) + np.full(len(lc.t.loc[ix]), halfwindowsize))
    temp = pd.Series(np.zeros(len(lc.t.loc[ix]) + 2*halfwindowsize), name='SNR', dtype=np.float64)
    temp[dataindices] = lc.t.loc[ix,'SNR']
    SNRsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
    lc.t.loc[ix,'SNRsum'] = list(SNRsum[dataindices])
    # normalize it
    norm_temp = pd.Series(np.zeros(len(lc.t.loc[ix]) + 2*halfwindowsize), name='norm', dtype=np.float64)
    norm_temp[np.array(range(len(lc.t.loc[ix])) + np.full(len(lc.t.loc[ix]), halfwindowsize))] = np.ones(len(lc.t.loc[ix]))
    norm_temp_sum = norm_temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
    lc.t.loc[ix,'SNRsumnorm'] = list(SNRsum.loc[dataindices] / norm_temp_sum.loc[dataindices] * max(norm_temp_sum.loc[dataindices]))

    # calculate the rolling SNR sum for SNR with simflux
    if sim_gauss or sim_erup:
        temp = pd.Series(np.zeros(len(lc.t.loc[ix]) + 2*halfwindowsize), name='SNRsim', dtype=np.float64)
        temp[dataindices] = lc.t.loc[ix,'SNRsim']
        SNRsimsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
        lc.t.loc[ix,'SNRsimsum'] = list(SNRsimsum.loc[dataindices])

    return lc

# optionally load simulated eruption from light curve
sim_erup = False
if not(erup_filename is None):
    print('\nUsing optional simulated eruption from light curve, instead of simulated gaussians')
    sim_erup = True
    erup = load_erup(erup_filename)

# fill lc_info and load SN and control lcs
lc_info['tnsname'] = tnsname
lc_info['discdate'] = discovery_date
lc_info['filt'] = filt
lc_info['mjdbinsize'] = mjd_bin_size
load_all_lcs(source_dir, n_controls)

# generate list of peaks
peak_mags, peak_fluxes = generate_peaks(peak_mag_min, peak_mag_max, n_peaks)

# fill in gauss_sigma column of efficiencies table
gauss_sigma_col = np.full(len(peak_fluxes), gauss_sigmas[0])
if len(gauss_sigmas) > 1:
    for gauss_sigma_index in range(1, len(gauss_sigmas)):
        a = np.full(len(peak_fluxes), gauss_sigmas[gauss_sigma_index])
        gauss_sigma_col = np.concatenate((gauss_sigma_col, a), axis=None)
efficiencies.t['gauss_sigma'] = gauss_sigma_col

print(f'\nGaussian sigmas in days: ', gauss_sigmas)
print(f'Peak magnitudes: ', peak_mags)
print(f'Peak fluxes (uJy): ', peak_fluxes)
if not(sim_erup):
    print(f'Possible simulated gaussian sigmas in days: ', sim_gauss_sigmas)
print(f'Number of iterations per peak: {iterations}')
print(f'Flag for bad days: {hex(flags)}')

gauss_sigma_offset = len(peak_fluxes)
for gauss_sigma_index in range(len(gauss_sigmas)):
    gauss_sigma = gauss_sigmas[gauss_sigma_index]
    print(f'\nUsing gaussian sigma of {gauss_sigma} days...')

    for peak_index in range(len(peak_fluxes)):
        peak_flux = peak_fluxes[peak_index]
        peak_appmag = peak_mags[peak_index]

        peak_index += gauss_sigma_index * gauss_sigma_offset
        efficiencies.t.loc[peak_index, 'peak_flux'] = peak_flux
        efficiencies.t.loc[peak_index, 'peak_appmag'] = peak_appmag
        efficiencies.t.loc[peak_index, 'gauss_sigma'] = gauss_sigma

        if sim_erup:
            print(f'Simulating eruptions from light curve with peak flux {peak_flux:0.2f} (peak appmag {peak_appmag:0.2f})')
        else:
            print(f'Simulating gaussians with peak flux {peak_flux:0.2f} (peak appmag {peak_appmag:0.2f})')

        # initialize detected table for the selected peak flux
        detected[f'{gauss_sigma}_{peak_flux}'] = pdastrostatsclass(columns=['gauss_sigma','peak_flux','peak_appmag','peak_mjd','sim_gauss_sigma','max_SNRsimsum','detection_limit','detected'])
        detected[f'{gauss_sigma}_{peak_flux}'].t['gauss_sigma'] = np.full(iterations, gauss_sigma)
        detected[f'{gauss_sigma}_{peak_flux}'].t['peak_flux'] = np.full(iterations, peak_flux)
        detected[f'{gauss_sigma}_{peak_flux}'].t['peak_appmag'] = np.full(iterations, peak_appmag)
        if sim_erup:
            detected[f'{gauss_sigma}_{peak_flux}'].t['sim_gauss_sigma'] = np.full(iterations, np.nan)
        else:
            detected[f'{gauss_sigma}_{peak_flux}'].t['erup_days'] = np.full(iterations, np.nan)
            detected[f'{gauss_sigma}_{peak_flux}'].t['erup_peakday'] = np.full(iterations, np.nan)

        for i in range(iterations):
            # pick random control light curve
            rand_control_index = random.randrange(1, n_controls+1, 1)
            cur_lc = copy.deepcopy(lcs[rand_control_index])

            # select random width in days
            if not(sim_erup):
                sim_gauss_sigma = random.choice(sim_gauss_sigmas[gauss_sigma_index])
                detected[f'{gauss_sigma}_{peak_flux}'].t.loc[i,'sim_gauss_sigma'] = sim_gauss_sigma
            
            # select random peak MJD from start of lc to 50 days before discovery date
            peak_mjd = random.randrange(cur_lc.t['MJDbin'].iloc[0]-0.5, int(discovery_date)-50, 1) + 0.5
            # make sure peak MJD is within an observation season; else redraw
            while not(in_season(peak_mjd, valid_mjd_ranges)):
                # redraw random peak mjd
                peak_mjd = random.randrange(cur_lc.t['MJDbin'].iloc[0]-0.5, int(discovery_date)-50, 1) + 0.5
            detected[f'{gauss_sigma}_{peak_flux}'].t.loc[i,'peak_mjd'] = peak_mjd

            if sim_erup:
                cur_lc = apply_rolling_sum(cur_lc, 
                                           gauss_sigma, 
                                           erup=copy.deepcopy(erup), 
                                           sim_peakmjd=peak_mjd, 
                                           sim_appmag=peak_appmag, 
                                           flag=flags, print_=False)
                
                # store peak appmag's days since start of erup
                detected[f'{gauss_sigma}_{peak_flux}'].t.loc[i,'erup_days'] = len(erup.t)
                detected[f'{gauss_sigma}_{peak_flux}'].t.loc[i,'erup_peakday'] = lc_info['erup_peakday']
            
            else:
                cur_lc = apply_rolling_sum(cur_lc, 
                                           gauss_sigma, 
                                           sim_gauss=True, 
                                           sim_peakmjd=peak_mjd, 
                                           sim_sigma=sim_gauss_sigma, 
                                           sim_appmag=peak_appmag, 
                                           flag=flags, print_=False)

            # compare max SNRsimsum to detection limit 
            if sim_erup:
                # only calculate max SNRsimsum from measurements where flux > half of peak flux (= fwhm = full width half maximum)
                lower_bound, upper_bound = get_erup_bounds(peak_mjd)
                t_ix = cur_lc.ix_inrange(colnames=['MJD'], lowlim=lower_bound, uplim=upper_bound)
                fwhm = peak_flux/2
                sigma_ix = cur_lc.ix_inrange(colnames=['uJysim'], lowlim=fwhm, indices=t_ix)
            else:
                # only calculate max SNRsimsum from measurements within 1 sigma of the simulated bump
                sigma_ix = cur_lc.ix_inrange(colnames=['MJD'], lowlim=peak_mjd-(sim_gauss_sigma), uplim=peak_mjd+(sim_gauss_sigma))
            if len(sigma_ix > 0):
                max_snrsimsum_ix = cur_lc.t.loc[sigma_ix,'SNRsimsum'].idxmax()
                detected[f'{gauss_sigma}_{peak_flux}'].t.loc[i,'max_SNRsimsum'] = cur_lc.t.loc[max_snrsimsum_ix,'SNRsimsum']
                detected[f'{gauss_sigma}_{peak_flux}'].t.loc[i,'max_SNRsimsum_mjd'] = cur_lc.t.loc[max_snrsimsum_ix,'MJD']
            else:
                # no valid measurements within this mjd range
                detected[f'{gauss_sigma}_{peak_flux}'].t.loc[i,'max_SNRsimsum'] = np.nan
                detected[f'{gauss_sigma}_{peak_flux}'].t.loc[i,'max_SNRsimsum_mjd'] = np.nan
        
        save_detected(gauss_sigma, peak_flux)

print()
efficiencies.write() 
save_efficiencies(efficiencies)