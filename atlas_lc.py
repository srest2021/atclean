#!/usr/bin/env python
"""
Author: Sofia Rest
"""

import json, requests, re, time, sys, os, math
from collections import OrderedDict
from astropy.time import Time
from datetime import datetime
import numpy as np
from pdastro import pdastrostatsclass, AandB, AnotB

class atlas_lc:
    def __init__(self, tnsname=None, is_averaged=False, mjd_bin_size=None, discdate=None, ra=None, dec=None):
        self.lcs = {}
        self.num_controls = 0

        self.tnsname = tnsname
        self.is_averaged = is_averaged
        self.mjd_bin_size = mjd_bin_size
        
        self.discdate = discdate
        self.ra = ra
        self.dec = dec

        self.corrected_baseline_ix = None
        self.during_sn_ix = None

        self.dflux_colnames = ['duJy']

    def __str__(self):
        res = f'SN {self.tnsname} light curve'
        if self.is_averaged:
            res += f' (averaged with MJD bin size {self.mjd_bin_size})'
        res += f': RA {self.ra}, Dec {self.dec}, discovery date {self.discdate} MJD'
        return res

    # get RA, Dec, and discovery date information from TNS
    def _get_tns_data(self, api_key, tns_id, bot_name):
        print(f'Obtaining RA, Dec, and/or discovery date from TNS (TNS ID: {tns_id}; TNS bot name: "{bot_name}")...')
        if tns_id == "None" or bot_name == "None":
            raise RuntimeError("# ERROR: Cannot query TNS without TNS ID and bot name! Specify these parameters in params.ini")

        try:
            url = 'https://www.wis-tns.org/api/get/object'
            json_file = OrderedDict([("objname",self.tnsname), ("objid",""), ("photometry","1"), ("spectra","1")])
            data = {'api_key':api_key,'data':json.dumps(json_file)}
            response = requests.post(url, data=data, headers={'User-Agent': 'tns_marker{"tns_id":"%s","type": "bot", "name":"%s"}' % (tns_id, bot_name)})
            json_data = json.loads(response.text,object_pairs_hook=OrderedDict)
            
            if self.ra is None:
                self.ra = json_data['data']['reply']['ra']
            if self.dec is None: 
                self.dec = json_data['data']['reply']['dec']
            discoverydate = json_data['data']['reply']['discoverydate']
        except Exception as e:
            print(json_data['data']['reply'])
            raise RuntimeError('# ERROR in get_tns_data(): '+str(e))

        date = list(discoverydate.partition(' '))[0]
        time = list(discoverydate.partition(' '))[2]
        dateobjects = Time(date+"T"+time, format='isot', scale='utc')
        if self.discdate is None:
            self.discdate = dateobjects.mjd - 20 # make sure no SN flux before discovery date in baseline indices	

    # get baseline indices (any indices before the SN discovery date)
    def get_baseline_ix(self):
        if self.discdate is None:
            raise RuntimeError('ERROR: Cannot get baseline indices because discovery date is None!')
        return self.lcs[0].ix_inrange(colnames=['MJD'],uplim=self.discdate-20,exclude_uplim=True)

    # get a light curve filename for saving/loading
    def get_filename(self, filt, control_index, directory):
        # SN light curve: 				DIRECTORY/2022xxx/2022xxx.o.lc.txt
        # averaged light curve: 		DIRECTORY/2022xxx/2022xxx.o.1.00days.lc.txt
        # control light curve: 			DIRECTORY/2022xxx/controls/2022xxx_i001.o.lc.txt
        # averaged control light curve: DIRECTORY/2022xxx/controls/2022xxx_i001.o.1.00days.lc.txt

        filename = f'{directory}/{self.tnsname}'
        if control_index != 0:
            filename += '/controls'
        filename += f'/{self.tnsname}'
        if control_index != 0:
            filename += f'_i{control_index:03d}'
        filename += f'.{filt}'
        if self.is_averaged:
            filename += f'.{self.mjd_bin_size:0.2f}days'
        filename += '.lc.txt'
        
        #print(f'# Filename: {filename}')
        return filename

    # save a single light curve
    def _save_lc(self, output_dir, control_index=0, filt=None, overwrite=True, keep_empty_bins=True):
        if filt is None:
            # split lc up by filt and save to two separate files
            for filt_ in ['o','c']:
                # if not keeping empty bins in averaged lc, remove all null rows; else keep all
                if self.is_averaged and not keep_empty_bins:
                    ix = AandB(self.lcs[control_index].ix_not_null(colnames=['uJy']), self.lcs[control_index].ix_equal(colnames=['F'],val=filt_))
                else: 
                    ix = self.lcs[control_index].ix_equal(colnames=['F'],val=filt_)
                
                self.lcs[control_index].write(filename=self.get_filename(filt_,control_index,output_dir), indices=ix, overwrite=overwrite)
        else:
            if self.is_averaged and not keep_empty_bins:
                ix = self.lcs[control_index].ix_not_null(colnames=['uJy'])
            else: 
                ix = self.lcs[control_index].getindices()
            self.lcs[control_index].write(filename=self.get_filename(filt,control_index,output_dir), indices=ix, overwrite=overwrite)

    # save SN light curve and, if necessary, control light curves
    def _save(self, output_dir, filt=None, overwrite=True, keep_empty_bins=True):
        if len(self.lcs) < 1:
            print('WARNING: No light curves to save! Skipping...')
        else:
            if self.is_averaged:
                output = f'\nSaving averaged SN light curve and {len(self.lcs)-1} averaged control light curves (keep empty bins: {keep_empty_bins})...'
            else:
                output = f'\nSaving SN light curve and {len(self.lcs)-1} control light curves...'
            print(output)

            for control_index in self.lcs:
                self._save_lc(output_dir, control_index, filt=filt, overwrite=overwrite, keep_empty_bins=keep_empty_bins)

    # load a single light curve
    def load_lc(self, output_dir, filt, is_averaged=False, control_index=0):
        if (len(self.lcs) > 0) and self.is_averaged != is_averaged:
            raise RuntimeError(f'ERROR: cannot load a light curve whose is_averaged status of {is_averaged} does not match previous status of {self.is_averaged}!')
        self.is_averaged = is_averaged

        self.lcs[control_index] = pdastrostatsclass()
        self.lcs[control_index].load_spacesep(self.get_filename(filt, control_index, output_dir), delim_whitespace=True)
        
        # clear previous 'Mask' column
        self.lcs[control_index].t['Mask'] = 0

    # load SN light curve and, if necessary, control light curves for a certain filter
    def _load(self, output_dir, filt, num_controls=0):
        output = f'\nLoading averaged SN light curve and {num_controls} averaged control light curves...' if self.is_averaged else f'\nLoading SN light curve and {num_controls} control light curves...'
        print(output)

        self.num_controls = num_controls

        self.load_lc(output_dir, filt, is_averaged=self.is_averaged)
        for control_index in range(1, num_controls+1):
            self.load_lc(output_dir, filt, is_averaged=self.is_averaged, control_index=control_index)

        self.dflux_colnames = ['duJy'] * (num_controls+1)

    def exists(self, output_dir, filt, control_index=0):
        filename = self.get_filename(filt, control_index, output_dir)
        return os.path.isfile(filename)

    # update given indices of 'Mask' column in the light curve (SN if control index is None) with given flag(s)
    def update_mask_col(self, flag, indices, control_index=0, remove_old=True):
        if remove_old:
            # remove the old mask
            self.lcs[control_index].t['Mask'] = np.bitwise_and(self.lcs[control_index].t['Mask'].astype(int), ~flag)

        if len(indices) > 1:
            flag_arr = np.full(self.lcs[control_index].t.loc[indices,'Mask'].shape, flag)
            self.lcs[control_index].t.loc[indices,'Mask'] = np.bitwise_or(self.lcs[control_index].t.loc[indices,'Mask'].astype(int), flag_arr)
        elif len(indices) == 1:
            self.lcs[control_index].t.loc[indices,'Mask'] = int(self.lcs[control_index].t.loc[indices,'Mask']) | flag

    # get the xth percentile SN flux using given indices
    def get_xth_percentile_flux(self, percentile, indices=None):
        if indices is None or len(indices)==0:
            indices = self.lcs[0].getindices()
        return np.percentile(self.lcs[0].t.loc[indices, 'uJy'], percentile)

    def get_filt_lens(self, control_index=0):
        o_len = len(self.lcs[control_index].ix_equal(colnames=['F'],val='o'))
        c_len = len(self.lcs[control_index].ix_equal(colnames=['F'],val='c'))
        return o_len, c_len
    
    def get_ix(self, control_index=0):
        return self.lcs[control_index].t.index.values

    def get_masked_ix(self, flags, control_index=0):
        flags_ = flags["chisquare"] | flags["uncertainty"] | flags["controls_bad"] | flags["avg_badday"]
        return self.lcs[control_index].ix_masked('Mask',maskval=flags_)

    def get_unmasked_ix(self, flags, control_index=0):
        flags_ = flags["chisquare"] | flags["uncertainty"] | flags["controls_bad"] | flags["avg_badday"]
        return self.lcs[control_index].ix_unmasked('Mask',maskval=flags_)

    def get_pre_SN_ix(self, control_index=0):
        return self.lcs[control_index].ix_inrange('MJD', uplim=self.discdate)
    
    def get_post_SN_ix(self, control_index=0):
        return self.lcs[control_index].ix_inrange('MJD', lowlim=self.discdate)
    
    def prep_for_cleaning(self):
        #print(f'# Dropping extra columns in all light curves...')
        self.drop_extra_columns()
        
        print('# Clearing \'Mask\' column, replacing infs, and calculating flux/dflux in all light curves...')
        for control_index in range(0, self.num_controls+1):
            self.lcs[control_index].t['Mask'] = 0 # clear 'Mask' column
            self.lcs[control_index].t = self.lcs[control_index].t.replace([np.inf, -np.inf], np.nan) # replace infs with NaNs
            self.recalculate_fdf(control_index=control_index)

        self.verify_mjds()

    # drop any added columns from previous iterations
    def drop_extra_columns(self, control_index=0):
        dropcols = []
        for col in ['Noffsetlc', 'uJy/duJy', '__tmp_SN', 'SNR', 'SNRsum', 'SNRsumnorm', 'SNRsim', 'SNRsimsum', 'c2_mean', 'c2_mean_err', 'c2_stdev', 'c2_stdev_err', 'c2_X2norm', 'c2_Ngood', 'c2_Nclip', 'c2_Nmask', 'c2_Nnan', 'c2_abs_stn']:
            if col in self.lcs[control_index].t.columns:
                dropcols.append(col)
        for col in self.lcs[control_index].t.columns:
            if re.search('^c\d_',col): 
                dropcols.append(col)

        # drop any extra columns
        if len(dropcols) > 0: 
            print(f'# Dropping extra columns ({"control light curve "+str(control_index) if control_index > 0 else "SN light curve"}): ',dropcols)
            self.lcs[control_index].t.drop(columns=dropcols,inplace=True)
    
    def recalculate_fdf(self, control_index=0):
        #print(f'# Recalculating flux/dflux column for control light curve {control_index:02d}...')
        self.lcs[control_index].t['uJy/duJy'] = self.lcs[control_index].t['uJy']/self.lcs[control_index].t[self.dflux_colnames[control_index]]
    
    # make sure that for every SN measurement, we have corresponding control light curve measurements at that MJD
    def verify_mjds(self):
        print('# Making sure SN and control light curve MJDs match up exactly...')
        # sort SN lc by MJD
        mjd_sorted_i = self.lcs[0].ix_sort_by_cols('MJD')
        self.lcs[0].t = self.lcs[0].t.loc[mjd_sorted_i]
        sn_sorted = self.lcs[0].t.loc[mjd_sorted_i,'MJD'].to_numpy()

        for control_index in range(1, self.num_controls+1):
            # sort control light curves by MJD
            mjd_sorted_i = self.lcs[control_index].ix_sort_by_cols('MJD')
            control_sorted = self.lcs[control_index].t.loc[mjd_sorted_i,'MJD'].to_numpy()
            
            # compare control light curve to SN light curve and, if out of agreement, fix
            if (len(sn_sorted) != len(control_sorted)) or (np.array_equal(sn_sorted, control_sorted) is False):
                print('## MJDs out of agreement for control light curve %03d, fixing...' % control_index)

                mjds_onlysn = AnotB(sn_sorted, control_sorted)
                mjds_onlycontrol = AnotB(control_sorted, sn_sorted)

                # for the MJDs only in SN, add row with that MJD to control light curve, with all values of other columns NaN
                if len(mjds_onlysn) > 0:
                    #print('### Adding %d NaN rows to control light curve...' % len(mjds_onlysn))
                    for mjd in mjds_onlysn:
                        self.lcs[control_index].newrow({'MJD':mjd,'Mask':0})
                
                # remove indices of rows in control light curve for which there is no MJD in the SN lc
                if len(mjds_onlycontrol) > 0:
                    #print('### Removing %d control light curve row(s) without matching SN row(s)...' % len(mjds_onlycontrol))
                    indices2skip = []
                    for mjd in mjds_onlycontrol:
                        ix = self.lcs[control_index].ix_equal('MJD',mjd)
                        if len(ix)!=1:
                            raise RuntimeError(f'### Couldn\'t find MJD={mjd} in column MJD, but should be there!')
                        indices2skip.extend(ix)
                    indices = AnotB(self.lcs[control_index].getindices(),indices2skip)
                else:
                    indices = self.lcs[control_index].getindices()
                
                ix_sorted = self.lcs[control_index].ix_sort_by_cols('MJD',indices=indices)
                self.lcs[control_index].t = self.lcs[control_index].t.loc[ix_sorted]
                
            self.lcs[control_index].t.reset_index(drop=True, inplace=True)

        print('# Success')

    def _clear_offset(self):
        if 'uJy_offset' in self.lcs[0].t.columns:
            print('# Subtracting previous offset...')
            self.lcs[0].t['uJy'] -= self.lcs[0].t['uJy_offset']
        print('# Setting current offset to 0...')
        self.lcs[0].t['uJy_offset'] = 0

    def _update_offset_col(self, offset, region_ix):
        print(f'## Recording offset {offset:0.2f}...')
        if not 'uJy_offset' in self.lcs[0].t.columns:
            self.lcs[0].t['uJy_offset'] = 0
        self.lcs[0].t.loc[region_ix, 'uJy_offset'] += offset

    def _get_mean(self, region_ix, maskval=None):
        indices = region_ix
        if not maskval is None:
            indices = self.lcs[0].ix_unmasked('Mask', maskval=maskval, indices=region_ix)
        self.lcs[0].calcaverage_sigmacutloop('uJy', noisecol='duJy', Nsigma=3, indices=indices, median_firstiteration=True)
        mean = self.lcs[0].statparams['mean']
        if mean is None:
            print('## WARNING: Could not converge on mean; taking median instead...')
            return np.median(self.lcs[0].t.loc[region_ix, 'uJy'])
        else:
            return round(mean,2)
        
    def _add_offset(self, offset, region_ix, region_name=''):
        self.lcs[0].t.loc[region_ix, 'uJy'] += offset
        s = f'Corrective flux {offset:0.2f} uJy added to {region_name} region'
        print(f'# {s}')
        self._update_offset_col(offset, region_ix)
        return s
    
    def _manual_template_correction(self, region1_offset=None, region2_offset=None, region3_offset=None):
        print('# Proceeding with manual template correction...')
        output = []

        self._clear_offset()

        t1, t2 = 58417, 58882
        region1_ix = self.lcs[0].ix_inrange('MJD', uplim=t1)
        region2_ix = self.lcs[0].ix_inrange('MJD', lowlim=t1, uplim=t2)
        region3_ix = self.lcs[0].ix_inrange('MJD', lowlim=t2)

        if not region1_offset is None:
            output.append(self._add_offset(region1_offset, region1_ix, region_name="first"))
        if not region2_offset is None:
            output.append(self._add_offset(region2_offset, region2_ix, region_name="second"))
        if not region3_offset is None:
            output.append(self._add_offset(region3_offset, region3_ix, region_name="third"))

        return output
    
    def _template_correction(self, maskval=None):
        print('# Proceeding with automatic template correction...')
        output = []

        self._clear_offset()
        
        t1, t2 = 58417, 58882
        region1_ix = self.lcs[0].ix_inrange('MJD', uplim=t1)
        region2_ix = self.lcs[0].ix_inrange('MJD', lowlim=t1, uplim=t2)
        region3_ix = self.lcs[0].ix_inrange('MJD', lowlim=t2)

        print(f'# Calculating offset for Region 1...')
        region1_mean = self._get_mean(region1_ix[-40:], maskval=maskval) # last 40 measurements before t1
        region2a_mean = self._get_mean(region2_ix[:40], maskval=maskval) # first 40 measurements after t1
        region1_offset = region2a_mean - region1_mean
        
        print(f'# Calculating offset for Region 3...')
        region2b_mean = self._get_mean(region2_ix[-40:], maskval=maskval) # last 40 measurements before t2
        region3_mean = self._get_mean(region3_ix[:40], maskval=maskval) # first 40 measurements after t2
        region3_offset = region2b_mean - region3_mean
        
        output.append(self._add_offset(region1_offset, region1_ix, region_name="first"))
        output.append(self._add_offset(region3_offset, region3_ix, region_name="third"))

        ix = self.get_ix()
        if self.discdate > 57600:
            global_ix = ix[:40]
        else:
            global_ix = ix[-40:]
        global_mean = self._get_mean(global_ix, maskval=maskval)
        global_offset = -global_mean
        output.append(self._add_offset(global_offset, ix, region_name="global"))

        return output
    
    def template_correction(self, maskval=None, region1_offset=None, region2_offset=None, region3_offset=None):
        print('\nCorrecting light curve flux due to template changes...')
        if self.discdate is None:
            raise RuntimeError(f'ERROR: discovery date cannot be set to None')

        if not(region1_offset is None and region2_offset is None and region3_offset is None):
            return '\n'.join(self._manual_template_correction(region1_offset=region1_offset, 
                                                              region2_offset=region2_offset, 
                                                              region3_offset=region3_offset))
        else:
            return '\n'.join(self._template_correction(maskval=maskval))