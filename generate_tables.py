#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys,copy,random,math
from pdastro import pdastrostatsclass, AandB, AnotB, AorB
from asym_gaussian import gauss2lc

# suppress deprecation warnings
import warnings
warnings.simplefilter('error', RuntimeWarning)
warnings.filterwarnings("ignore")

# for storing lc and control lcs 
global lcs
lcs = {}

# for storing info about the lc
global lc_info
lc_info = {}

# for storing per-peak tables
global detected
detected = {}

# final table
detected_per_peak = pdastrostatsclass(columns=['gauss_sigma','peak_flux', 'peak_appmag'])

# convert flux to magnitude 
def flux2mag(flux):
    return -2.5 * np.log10(flux) + 23.9

# convert magnitude to flux
def mag2flux(mag):
    return 10 ** ((mag - 23.9) / -2.5)

######### SETTINGS ######### 

# SN TNS name
tnsname = '2017fra'

# SN discovery date
discovery_date = 57939.34 - 20.0

# path to directory that contains SN and control light curves
source_dir = '/Users/sofiarest/Desktop/Supernovae/data/misc/2017fra'

# number of control light curves to load
n_controls = 8

# filter of light curve to analyze
filt = 'o'

# MJD bin size in days of light curve to analyze
mjd_bin_size = 1.0

# search for pre-SN bumps with gaussian sigmas APPROXIMATELY within the following range
gauss_sigma_max = 140
gauss_sigma_min = 2.5
n_gauss_sigmas = 6 # number of gaussian sigmas to generate in log space

# select range of peak apparent magnitudes to simulate
peak_mag_max = 16 #flux2mag(150) # brightest magnitude
peak_mag_min = 22 #flux2mag(2) # faintest magnitude
n_peaks = 20 # number of magnitudes to generate in log space

# select range of gaussian sigmas to simulate APPROXIMATELY within the following range
sim_gauss_sigma_max = 140
sim_gauss_sigma_min = 2.5
n_sim_gauss_sigmas = 6 # number of gaussian sigmas to generate in log space

# number of iterations of random sigma and peak mjd per peak
iterations = 50000 

# define flag that defines bad measurements (to filter out bad days in lc)
# if more than one, use bitwise OR (this symbol: |) to combine them
flags = 0x800000

# add observation seasons' mjd ranges here
valid_mjd_ranges = [[57463,57622], [57762,58047], [58117,58408], [58473,58773], [58846,59131], [59214,59509], [59566,59862], [59958,60085]]

############################

def generate_peaks(peak_mag_min, peak_mag_max, n_peaks):
    peak_mags = list(np.linspace(peak_mag_min, peak_mag_max, num=20))
    peak_fluxes = list(map(mag2flux, peak_mags))

    peak_mags = [round(item, 2) for item in peak_mags]
    peak_fluxes = [round(item, 2) for item in peak_fluxes]

    return peak_mags, peak_fluxes

def generate_gauss_sigmas(gauss_sigma_min, gauss_sigma_max, n_gauss_sigmas):
    log_min = round(math.log2(gauss_sigma_min))
    log_max = round(math.log2(gauss_sigma_max))
    gauss_sigmas = list(np.logspace(log_min, log_max, num=n_gauss_sigmas, base=2))
    gauss_sigmas = [round(item) for item in gauss_sigmas]
    return gauss_sigmas

# optionally add a simulated pre-eruption to the loaded light curve
# then apply a gaussian weighted rolling sum to the SNR
def apply_gaussian(lc, gaussian_sigma, flag=0x800000, sim_gauss=False, sim_peakmjd=None, sim_appmag=None, sim_sigma=None, print_=True):
    if len(lc_info['baseline_ix']) <= 0:
        print('ERROR: Not enough data in baseline flux (before SN start)!')

    ix = lc_info['baseline_ix'] # all pre-SN indices
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
    if sim_gauss:
        if print_:
            print(f'# Adding simulated gaussian: peak at MJD {sim_peakmjd:0.2f}; apparent magnitude {sim_appmag:0.2f}; sigma- and sigma+ of {sim_sigma:0.2f} days')
        mjds = lc.t.loc[good_ix,'MJD']
        lc.t.loc[good_ix,'uJysim'] = lc.t.loc[good_ix,'uJy']
        lc.t.loc[ix,'simLC'] = 0.0

        # get simulated gaussian flux and add to light curve flux
        simflux = gauss2lc(mjds, sim_peakmjd, sim_sigma, sim_sigma, app_mag=sim_appmag)
        lc.t.loc[good_ix,'uJysim'] += simflux
        # get the simulated lc for all MJDs
        simflux_all = gauss2lc(lc.t.loc[ix,'MJDbin'], sim_peakmjd, sim_sigma, sim_sigma, app_mag=sim_appmag)
        lc.t.loc[ix,'simLC'] += simflux_all

        # make sure all bad rows have SNRsim = 0.0 so they have no impact on the rolling SNRsum
        lc.t.loc[ix,'SNRsim'] = 0.0
        # include simflux in the SNR
        lc.t.loc[good_ix,'SNRsim'] = lc.t.loc[good_ix,'uJysim']/lc.t.loc[good_ix,'duJy']

    new_gaussian_sigma = round(gaussian_sigma/lc_info['mjd_bin_size'])
    windowsize = int(6 * new_gaussian_sigma)
    halfwindowsize = int(windowsize * 0.5) + 1
    if print_:
        print(f'# Sigma: {gaussian_sigma:0.2f} days; MJD bin size: {lc_info["mjd_bin_size"]:0.2f} days; sigma: {new_gaussian_sigma:0.2f} bins; window size: {windowsize} bins')

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
    if sim_gauss:
        temp = pd.Series(np.zeros(len(lc.t.loc[ix]) + 2*halfwindowsize), name='SNRsim', dtype=np.float64)
        temp[dataindices] = lc.t.loc[ix,'SNRsim']
        SNRsimsum = temp.rolling(windowsize, center=True, win_type='gaussian').sum(std=new_gaussian_sigma)
        lc.t.loc[ix,'SNRsimsum'] = list(SNRsimsum.loc[dataindices])

    return lc

def load_lc(source_dir, tnsname):
    lcs[0] = pdastrostatsclass()
    filename = f'{source_dir}/{lc_info["tnsname"]}.{lc_info["filt"]}.{lc_info["mjd_bin_size"]:0.2f}days.lc.txt'
    print(f'\nLoading light curve for SN {tnsname} at {filename}...')
    try:
        lcs[0].load_spacesep(filename,delim_whitespace=True)
    except Exception as e:
        print('ERROR: Could not load light curve at %s: %s' % (filename, str(e)))
        sys.exit()
    
    lc_info['baseline_ix'] = lcs[0].ix_inrange(colnames=['MJDbin'], uplim=lc_info['discovery_date']-20, exclude_uplim=True)
    if len(lc_info['baseline_ix'])<=0:
        print('Baseline length is 0! Exiting...')
        sys.exit()
    lc_info['afterdiscdate_ix'] = AnotB(lcs[0].getindices(), lc_info['baseline_ix']) 

def load_control_lcs(source_dir, n_controls): 
    print(f'\nLoading control light curves for SN {tnsname}...')
    for control_index in range(1, n_controls+1):
        lcs[control_index] = pdastrostatsclass()
        filename = f'{source_dir}/controls/{lc_info["tnsname"]}_i{control_index:03d}.{lc_info["filt"]}.{lc_info["mjd_bin_size"]:0.2f}days.lc.txt'
        print(f'# Loading control light curve {control_index:03d} at {filename}...')
        try:
            lcs[control_index].load_spacesep(filename,delim_whitespace=True)
        except Exception as e:
            print(f'ERROR: Could not load contol light curve {control_index:03d} at {filename}: {str(e)}')
            sys.exit()

def load_lcs(source_dir, tnsname, n_controls):
    load_lc(source_dir, tnsname)
    load_control_lcs(source_dir, n_controls)

def save_detected(gauss_sigma, peak):
    print('Saving...')
    detected[f'{gauss_sigma}_{peak}'].write(filename=f'{source_dir}/bump_analysis/tables2/detected_{gauss_sigma}_{peak:0.2f}.txt')

def save_detected_per_peak(detected_per_peak):
    detected_per_peak.write(filename=f'{source_dir}/bump_analysis/tables2/detected_per_peak.txt')

def save_tables(peaks, detected, detected_per_peak):
    for gauss_sigma in gauss_sigmas: 
        for peak in peaks:
            save_detected(gauss_sigma, peak)
    save_detected_per_peak(detected_per_peak)

# check if mjd is within valid mjd season
def in_season(mjd, valid_mjd_ranges):
    for mjd_range in valid_mjd_ranges:
        if mjd >= mjd_range[0] and mjd <= mjd_range[1]:
            return True
    return False

# fill lc_info and load SN and control lcs
lc_info['tnsname'] = tnsname
lc_info['discovery_date'] = discovery_date
lc_info['filt'] = filt
lc_info['mjd_bin_size'] = mjd_bin_size
load_lcs(source_dir, tnsname, n_controls)

# generate list of peaks
peak_mags, peaks = generate_peaks(peak_mag_min, peak_mag_max, n_peaks)
# generate list of gaussian sigmas to look for
gauss_sigmas = generate_gauss_sigmas(gauss_sigma_min, gauss_sigma_max, n_gauss_sigmas)
# generate list of gaussian sigmas to simulate
sim_gauss_sigmas = generate_gauss_sigmas(sim_gauss_sigma_min, sim_gauss_sigma_max, n_sim_gauss_sigmas)

# fill in gauss_sigma column of final table
gauss_sigma_col = np.full(len(peaks), gauss_sigmas[0])
if len(gauss_sigmas) > 1:
    for gauss_sigma_index in range(1, len(gauss_sigmas)):
        a = np.full(len(peaks), gauss_sigmas[gauss_sigma_index])
        gauss_sigma_col = np.concatenate((gauss_sigma_col, a), axis=None)
detected_per_peak.t['gauss_sigma'] = gauss_sigma_col

print(f'\nGaussian sigmas in days: ', gauss_sigmas)
print(f'Peak magnitudes: ', peak_mags)
print(f'Peak fluxes (uJy): ', peaks)
print(f'Possible simulated gaussian sigmas in days: ', sim_gauss_sigmas)
print(f'Number of iterations per peak: {iterations}')
print(f'Flag for bad days: {hex(flags)}')

gauss_sigma_offset = len(peaks)
for gauss_sigma_index in range(len(gauss_sigmas)):
    gauss_sigma = gauss_sigmas[gauss_sigma_index]
    print(f'\nUsing gaussian sigma of {gauss_sigma} days...')

    for peak_index in range(len(peaks)):
        peak = peaks[peak_index]
        peak_appmag = peak_mags[peak_index] #[len(peaks) - 1 - peak_index] #-2.5 * np.log10(peak) + 23.9  # convert peak flux to peak apparent magnitude

        peak_index += gauss_sigma_index * gauss_sigma_offset
        detected_per_peak.t.loc[peak_index, 'peak_flux'] = peak
        detected_per_peak.t.loc[peak_index, 'peak_appmag'] = peak_appmag
        detected_per_peak.t.loc[peak_index, 'gauss_sigma'] = gauss_sigma

        print(f'Simulating gaussians with peak flux {peak:0.2f} (peak appmag {peak_appmag:0.2f})')

        # initialize per-peak table
        detected[f'{gauss_sigma}_{peak}'] = pdastrostatsclass(columns=['gauss_sigma','peak_flux','peak_appmag','peak_mjd','sigma_days','max_SNRsimsum','detection_limit','detected'])
        detected[f'{gauss_sigma}_{peak}'].t['gauss_sigma'] = np.full(iterations, gauss_sigma)
        detected[f'{gauss_sigma}_{peak}'].t['peak_flux'] = np.full(iterations, peak)
        detected[f'{gauss_sigma}_{peak}'].t['peak_appmag'] = np.full(iterations, peak_appmag)
        
        for i in range(iterations):
            # pick random control light curve
            rand_control_index = random.randrange(1, n_controls, 1)
            cur_lc = copy.deepcopy(lcs[rand_control_index])

            # select random width in days
            sim_gauss_sigma = random.choice(sim_gauss_sigmas)
            detected[f'{gauss_sigma}_{peak}'].t.loc[i,'sigma_days'] = sim_gauss_sigma

            # select random peak MJD from start of lc to 50 days before discovery date
            peak_mjd = random.randrange(cur_lc.t['MJDbin'].iloc[0]-0.5, int(discovery_date)-50, 1) + 0.5
            # make sure peak MJD is within an observation season; else redraw
            while not(in_season(peak_mjd, valid_mjd_ranges)):
                # redraw random peak mjd
                peak_mjd = random.randrange(cur_lc.t['MJDbin'].iloc[0]-0.5, int(discovery_date)-50, 1) + 0.5
            detected[f'{gauss_sigma}_{peak}'].t.loc[i,'peak_mjd'] = peak_mjd

            cur_lc = apply_gaussian(cur_lc, gauss_sigma, flag=flags, sim_gauss=True, sim_peakmjd=peak_mjd, sim_sigma=sim_gauss_sigma, sim_appmag=peak_appmag, print_=False)

            # compare max SNRsimsum to detection limit 
            # only calculate max SNRsimsum from measurements within 1 sigma of the simulated bump
            sigma_ix = cur_lc.ix_inrange(colnames=['MJD'], lowlim=peak_mjd-(sim_gauss_sigma), uplim=peak_mjd+(sim_gauss_sigma), indices=lc_info['baseline_ix'])
            if len(sigma_ix > 0):
                max_snrsimsum_ix = cur_lc.t.loc[sigma_ix,'SNRsimsum'].idxmax()
                detected[f'{gauss_sigma}_{peak}'].t.loc[i,'max_SNRsimsum'] = cur_lc.t.loc[max_snrsimsum_ix,'SNRsimsum']
                detected[f'{gauss_sigma}_{peak}'].t.loc[i,'max_SNRsimsum_mjd'] = cur_lc.t.loc[max_snrsimsum_ix,'MJD']
            else:
                # no valid measurements within this mjd range
                detected[f'{gauss_sigma}_{peak}'].t.loc[i,'max_SNRsimsum'] = np.nan
                detected[f'{gauss_sigma}_{peak}'].t.loc[i,'max_SNRsimsum_mjd'] = np.nan

        save_detected(gauss_sigma, peak)

detected_per_peak.write() 
save_detected_per_peak(detected_per_peak)
