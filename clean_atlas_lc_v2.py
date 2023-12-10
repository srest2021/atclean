#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import sys, argparse, configparser, re, os
from copy import deepcopy
import pandas as pd
import numpy as np

from pdastro import pdastrostatsclass, AandB, AnotB
from atlas_lc import atlas_lc

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

class CleanAtlasLightCurve(atlas_lc):
	def __init__(self, filt, cfg=None, **kwargs):
		atlas_lc.__init__(self, **kwargs)
		self.cfg = cfg
		self.filt = filt

	def get_tns_data(self):
		self._get_tns_data(self.cfg['api_key'], self.cfg['tns_id'], self.cfg['bot_name'])

	def set_tns_data(self, ra, dec, discdate):
		self.ra = ra
		self.dec = dec
		self.discdate = discdate

	def save(self):
		for control_index in range(self.num_controls+1):
			self.drop_extra_columns(control_index=control_index)
		self._save(self.cfg['output_dir'], filt=self.filt)
		print('Success')
	
	def load(self):
		self._load(self.cfg['output_dir'], self.filt, num_controls=self.cfg['num_controls'])
		print('Success')

	def _apply_uncert_cut(self, flag, cut):
		print(f'\nApplying uncertainty cut of {cut:0.2f}...')
		for control_index in range(self.num_controls+1):
			ix = self.get_ix(control_index=control_index)
			kept_ix = self.lcs[control_index].ix_inrange(colnames=['duJy'],uplim=cut)
			cut_ix = AnotB(ix, kept_ix)
			self.update_mask_col(flag, cut_ix, control_index=control_index)
		
			if control_index == 0:
				num_cut = len(cut_ix)
				percent_cut = 100 * num_cut/len(ix)
				output = f'Total percent of data flagged ({hex(flag)}): {percent_cut:0.2f}%'
				print(f'# {output}')
		return output

	def apply_uncert_cut(self, flag):
		return self._apply_uncert_cut(flag, self.cfg['uncert_cut'])
	
	def _get_median_dflux(self, indices=None, control_index=0):
		if indices is None:
			indices = self.get_ix(control_index=control_index)
		return np.median(self.lcs[control_index].t.loc[indices, 'duJy'])

	def _get_flux_stdev(self, indices=None, control_index=0):
		self.lcs[control_index].calcaverage_sigmacutloop('uJy', indices=indices, Nsigma=3.0, median_firstiteration=True)
		return self.lcs[control_index].statparams['stdev']

	def _get_sigma_extra(self, median_dflux, stdev):
		return max(0, np.sqrt(stdev**2 - median_dflux**2))

	def _get_uncert_est_stats(self, uncert_flag, prelim_x2_cut):
		stats = pd.DataFrame(columns=['control_index', 'median_dflux', 'stdev', 'sigma_extra'])
		stats['control_index'] = list(range(1, self.num_controls+1))
		stats.set_index('control_index', inplace=True)

		for control_index in range(1, self.num_controls+1):
			dflux_clean_ix = self.lcs[control_index].ix_unmasked('Mask', maskval=uncert_flag)
			x2_clean_ix = self.lcs[control_index].ix_inrange(colnames=['chi/N'], uplim=prelim_x2_cut, exclude_uplim=True)
			clean_ix = AandB(dflux_clean_ix, x2_clean_ix)

			median_dflux = self._get_median_dflux(indices=clean_ix, control_index=control_index)
			
			stdev = self._get_flux_stdev(indices=clean_ix, control_index=control_index)
			if stdev is None:
				print(f'# WARNING: Could not get flux std dev using clean indices; retrying without preliminary chi-square cut of {prelim_x2_cut}...')
				stdev = self._get_flux_stdev(indices=dflux_clean_ix, control_index=control_index)
				if stdev is None:
					print('# WARNING: Could not get flux std dev using clean indices; retrying with all indices...')
					stdev = self._get_flux_stdev(control_index=control_index)

			sigma_extra = self._get_sigma_extra(median_dflux, stdev)

			stats.loc[control_index, 'median_dflux'] = median_dflux
			stats.loc[control_index, 'stdev'] = stdev
			stats.loc[control_index, 'sigma_extra'] = sigma_extra
		return stats
	
	def _get_final_sigma_extra(self, stats):
		return np.median(stats['sigma_extra'])
	
	def _add_noise(self, sigma_extra):
		self.dflux_colnames = ['duJy_new'] * (self.num_controls+1)
		for control_index in range(self.num_controls+1):
			self.lcs[control_index].t['duJy_new'] = np.sqrt(self.lcs[control_index].t['duJy']*self.lcs[control_index].t['duJy'] + sigma_extra**2)
			self.recalculate_fdf(control_index=control_index)

	def _apply_uncert_est(self, uncert_flag, prelim_x2_cut):
		print('\nApplying true uncertainty estimation...')

		stats = self._get_uncert_est_stats(uncert_flag, prelim_x2_cut)
		final_sigma_extra = self._get_final_sigma_extra(stats)
		#print(f'# Final sigma extra: {final_sigma_extra:0.2f}')
		
		sigma_typical_old = np.median(stats['median_dflux'])
		sigma_typical_new = np.sqrt(final_sigma_extra**2 + sigma_typical_old**2)
		percent_greater = 100 * ((sigma_typical_new - sigma_typical_old)/sigma_typical_old)
		s1 = f'We increase the typical uncertainties from {sigma_typical_old:0.2f} to {sigma_typical_new:0.2f} by adding an additional systematic uncertainty of {final_sigma_extra:0.2f} in quadrature'
		print(f'# {s1}')
		print(f'# New typical uncertainty is {percent_greater:0.2f}% greater than old typical uncertainty')

		if percent_greater >= 10:
			output = f'{s1}.\nThe extra noise was added to the uncertainties of the SN light curve and copied to the "duJy_new" column.'
			print(f'# Proceeding with true uncertainties estimation...')
			print(f'# Calculating new uncertainties in \'duJy_new\' column for each light curve...') 
			self._add_noise(final_sigma_extra)
			print('Success')
		else:
			s2 = 'True uncertainties estimation not needed'
			output = f'{s2}.'
			print(f'# {s2}; skipping procedure')

		return output#, stats
	
	def apply_uncert_est(self, uncert_flag):
		#output, stats = self._apply_uncert_est(uncert_flag, self.cfg['prelim_x2_cut'])
		return self._apply_uncert_est(uncert_flag, self.cfg['prelim_x2_cut']) #output
	
	def _get_all_controls(self):
		controls = [deepcopy(self.lcs[control_index].t) for control_index in self.lcs if control_index > 0]
		all_controls = pdastrostatsclass()
		all_controls.t = pd.concat(controls, ignore_index=True)
		return all_controls
	
	def _get_goodbad_ix(self, lc, ix, stn_cut):
		good_ix = lc.ix_inrange(colnames=['uJy/duJy'], lowlim=-stn_cut, uplim=stn_cut, indices=ix)
		bad_ix = AnotB(ix, good_ix)
		return good_ix, bad_ix
	
	def _get_keptcut_ix(self, lc, ix, cut):
		kept_ix = lc.ix_inrange(colnames=['chi/N'], uplim=cut, indices=ix)
		cut_ix = AnotB(ix, kept_ix)
		return kept_ix, cut_ix
	
	def _get_limcuts_data(self, lc, cut, ix, good_ix, bad_ix, kept_ix=None, cut_ix=None):
		if kept_ix is None or cut_ix is None:
			kept_ix, cut_ix = self._get_keptcut_ix(lc, ix, cut)
		out = {}
		out['PSF Chi-Square Cut'] = cut
		out['N'] = len(ix)
		out['Ngood'] = len(good_ix)
		out['Nbad'] = len(bad_ix)
		out['Nkept'] = len(kept_ix)
		out['Ncut'] = len(cut_ix)
		out['Ngood,kept'] = len(AandB(good_ix,kept_ix))
		out['Ngood,cut'] = len(AandB(good_ix,cut_ix))
		out['Nbad,kept'] = len(AandB(bad_ix,kept_ix))
		out['Nbad,cut'] = len(AandB(bad_ix,cut_ix))
		out['Pgood,kept'] = 100*len(AandB(good_ix,kept_ix))/len(ix)
		out['Pgood,cut'] = 100*len(AandB(good_ix,cut_ix))/len(ix)
		out['Pbad,kept'] = 100*len(AandB(bad_ix,kept_ix))/len(ix)
		out['Pbad,cut'] = 100*len(AandB(bad_ix,cut_ix))/len(ix)
		out['Ngood,kept/Ngood'] = 100*len(AandB(good_ix,kept_ix))/len(good_ix)
		out['Ploss'] = 100*len(AandB(good_ix,cut_ix))/len(good_ix)
		out['Pcontamination'] = 100*len(AandB(bad_ix,kept_ix))/len( kept_ix)
		return out
	
	def _get_limcuts_table(self, lc, ix, good_ix, bad_ix, cut_start, cut_stop, cut_step):
		print(f"# Calculating loss and contamination for chi-square cuts from {cut_start} to {cut_stop}...")

		limcuts = pd.DataFrame(columns=['PSF Chi-Square Cut', 'N', 'Ngood', 'Nbad', 'Nkept', 'Ncut', 'Ngood,kept', 'Ngood,cut', 'Nbad,kept', 'Nbad,cut',
										'Pgood,kept', 'Pgood,cut', 'Pbad,kept', 'Pbad,cut', 'Ngood,kept/Ngood', 'Ploss', 'Pcontamination'])
		
		# for different x2 cuts decreasing from 50
		for cut in range(cut_start, cut_stop+1, cut_step):
			kept_ix, cut_ix = self._get_keptcut_ix(lc, ix, cut)
			percent_kept = 100*(len(kept_ix)/len(ix))
			if percent_kept < 10:
				# less than 10% of measurements kept, so no chi-square cuts beyond this point are valid
				continue
			data = self._get_limcuts_data(lc, cut, ix, good_ix, bad_ix, kept_ix=kept_ix, cut_ix=cut_ix)
			limcuts = pd.concat([limcuts, pd.DataFrame([data])],ignore_index=True)
		return limcuts
	
	def _get_limcuts(self, limcuts, loss_lim, contam_lim):
		contam_lim_cut = None
		loss_lim_cut = None
		contam_case = None
		loss_case = None

		sortby_loss = limcuts.iloc[(limcuts['Ploss']).argsort()].reset_index()
		min_loss = sortby_loss.loc[0,'Ploss']
		max_loss = sortby_loss.loc[len(sortby_loss)-1,'Ploss']
		# if all loss below lim, loss_lim_cut is min cut
		if min_loss < loss_lim and max_loss < loss_lim:
			loss_case = 'below lim'
			loss_lim_cut = limcuts.loc[0,'PSF Chi-Square Cut']
		else:
			# else if all loss above lim, loss_lim_cut is min cut with min% loss
			if min_loss > loss_lim and max_loss > loss_lim:
				loss_case = 'above lim'
				a = np.where(limcuts['Ploss'] == min_loss)[0]
				b = limcuts.iloc[a]
				c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
				loss_lim_cut = c.loc[0,'PSF Chi-Square Cut']
			# else if loss crosses lim at some point, loss_lim_cut is min cut with max% loss <= loss_lim
			else:
				loss_case = 'crosses lim'
				valid_cuts = sortby_loss[sortby_loss['Ploss'] <= loss_lim]
				a = np.where(limcuts['Ploss'] == valid_cuts.loc[len(valid_cuts)-1,'Ploss'])[0]
				# sort by cuts
				b = limcuts.iloc[a]
				c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
				# get midpoint of loss1 and loss2 (two points on either side of lim)
				loss1_i = np.where(limcuts['PSF Chi-Square Cut'] == c.loc[0,'PSF Chi-Square Cut'])[0][0]
				if limcuts.loc[loss1_i,'Ploss'] == loss_lim:
					loss_lim_cut = limcuts.loc[loss1_i,'PSF Chi-Square Cut']
				else:
					loss2_i = loss1_i - 1
					x = np.array([limcuts.loc[loss1_i,'PSF Chi-Square Cut'], limcuts.loc[loss2_i,'PSF Chi-Square Cut']])
					contam_y = np.array([limcuts.loc[loss1_i,'Pcontamination'], limcuts.loc[loss2_i,'Pcontamination']])
					loss_y = np.array([limcuts.loc[loss1_i,'Ploss'], limcuts.loc[loss2_i,'Ploss']])
					contam_line = np.polyfit(x,contam_y,1)
					loss_line = np.polyfit(x,loss_y,1)
					loss_lim_cut = (loss_lim-loss_line[1])/loss_line[0]

		sortby_contam = limcuts.iloc[(limcuts['Pcontamination']).argsort()].reset_index()
		min_contam = sortby_contam.loc[0,'Pcontamination']
		max_contam = sortby_contam.loc[len(sortby_contam)-1,'Pcontamination']
		# if all contam below lim, contam_lim_cut is max cut
		if min_contam < contam_lim and max_contam < contam_lim:
			contam_case = 'below lim'
			contam_lim_cut = limcuts.loc[len(limcuts)-1,'PSF Chi-Square Cut']
		else:
			# else if all contam above lim, contam_lim_cut is max cut with min% contam
			if min_contam > contam_lim and max_contam > contam_lim:
				contam_case = 'above lim'
				a = np.where(limcuts['Pcontamination'] == min_contam)[0]
				b = limcuts.iloc[a]
				c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
				contam_lim_cut = c.loc[len(c)-1,'PSF Chi-Square Cut']
			# else if contam crosses lim at some point, contam_lim_cut is max cut with max% contam <= contam_lim
			else:
				contam_case = 'crosses lim'
				valid_cuts = sortby_contam[sortby_contam['Pcontamination'] <= contam_lim]
				a = np.where(limcuts['Pcontamination'] == valid_cuts.loc[len(valid_cuts)-1,'Pcontamination'])[0]
				# sort by cuts
				b = limcuts.iloc[a]
				c = b.iloc[(b['PSF Chi-Square Cut']).argsort()].reset_index()
				# get midpoint of contam1 and contam2 (two points on either side of lim)
				contam1_i = np.where(limcuts['PSF Chi-Square Cut'] == c.loc[len(c)-1,'PSF Chi-Square Cut'])[0][0]
				if limcuts.loc[contam1_i,'Pcontamination'] == contam_lim:
					contam_lim_cut = limcuts.loc[contam1_i,'PSF Chi-Square Cut']
				else:
					contam2_i = contam1_i + 1
					x = np.array([limcuts.loc[contam1_i,'PSF Chi-Square Cut'], limcuts.loc[contam2_i,'PSF Chi-Square Cut']])
					contam_y = np.array([limcuts.loc[contam1_i,'Pcontamination'], limcuts.loc[contam2_i,'Pcontamination']])
					loss_y = np.array([limcuts.loc[contam1_i,'Ploss'], limcuts.loc[contam2_i,'Ploss']])
					contam_line = np.polyfit(x,contam_y,1)
					loss_line = np.polyfit(x,loss_y,1)
					contam_lim_cut = (contam_lim-contam_line[1])/contam_line[0]

		return contam_lim_cut, loss_lim_cut, contam_case, loss_case
	
	def _choose_btwn_cuts(self, contam_lim_cut, loss_lim_cut, contam_case, loss_case, loss_lim, contam_lim, cut_start, lim_to_prioritize):
		# case 1 and 1: final_cut = 3
		# case 1 and 2: take limit of case 2
		# case 1 and 3: take limit of case 3
		# case 2 and 2: print lims don't work
		# case 2 and 3: choose_btwn_lim_cuts
		# case 3 and 3: choose_btwn_lim_cuts

		case1 = loss_case == 'below lim' or contam_case == 'below lim'
		case2 = loss_case == 'above lim' or contam_case == 'above lim'
		case3 = loss_case == 'crosses lim' or contam_case == 'crosses lim'

		final_cut = None
		if case1 and not case2 and not case3: # 1 and 1
			print('# Valid chi-square cut range from %0.2f to %0.2f! Setting to 3...' % (loss_lim_cut, contam_lim_cut))
			final_cut = cut_start
		elif case1: # 1
			if case2: # and 2
				if loss_case == 'above lim':
					print('# WARNING: contam_lim_cut <= %0.2f falls below limit %0.2f%%, but loss_lim_cut >= %0.2f falls above limit %0.2f%%! Setting to %0.2f...' % (contam_lim_cut, contam_lim, loss_lim_cut, loss_lim, loss_lim_cut))
					final_cut = loss_lim_cut
				else:
					print('# WARNING: loss_lim_cut <= %0.2f falls below limit %0.2f%%, but contam_lim_cut >= %0.2f falls above limit %0.2f%%! Setting to %0.2f...' % (loss_lim_cut, loss_lim, contam_lim_cut, contam_lim, contam_lim_cut))
					final_cut = contam_lim_cut
			else: # and 3
				if loss_case == 'crosses lim':
					print('# Contam_lim_cut <= %0.2f falls below limit %0.2f%% and loss_lim_cut >= %0.2f crosses limit %0.2f%%, setting to %0.2f...' % (contam_lim_cut, contam_lim, loss_lim_cut, loss_lim, loss_lim_cut))
					final_cut = loss_lim_cut
				else:
					print('Loss_lim_cut <= %0.2f falls below limit %0.2f%% and contam_lim_cut >= %0.2f crosses limit %0.2f%%, setting to %0.2f...' % (loss_lim_cut, loss_lim, contam_lim_cut, contam_lim, contam_lim_cut))
					final_cut = contam_lim_cut
		elif case2 and not case3: # 2 and 2
			print('ERROR: chi-square loss_lim_cut >= %0.2f and contam_lim_cut <= %0.2f both fall above limits %0.2f%% and %0.2f%%! Try setting less strict limits. Setting final cut to nan.' % (loss_lim_cut, contam_lim_cut, loss_lim, contam_lim))
			final_cut = np.nan
		else: # 2 and 3 or 3 and 3
			if loss_lim_cut > contam_lim_cut:
				print('# WARNING: chi-square loss_lim_cut >= %0.2f and contam_lim_cut <= %0.2f do not overlap! ' % (loss_lim_cut, contam_lim_cut))
				if lim_to_prioritize == 'contam_lim':
					print('# Prioritizing %s and setting to %0.2f...' % (lim_to_prioritize, contam_lim_cut))
					final_cut = contam_lim_cut
				else:
					print('# Prioritizing %s and setting to %0.2f... ' % (lim_to_prioritize, loss_lim_cut))
					final_cut = loss_lim_cut
			else:
				print('# Valid chi-square cut range from %0.2f to %0.2f! ' % (loss_lim_cut, contam_lim_cut))
				if lim_to_prioritize == 'contam_lim':
					print('# Prioritizing %s and setting to %0.2f... ' % (lim_to_prioritize, loss_lim_cut))
					final_cut = loss_lim_cut
				else:
					print('# Prioritizing %s and setting to %0.2f... ' % (lim_to_prioritize, contam_lim_cut))
					final_cut = contam_lim_cut
		return final_cut
	
	def _apply_x2_cut(self, cut, data, flag):
		print(f'\nApplying chi-square cut of {cut:0.2f}...')
		for control_index in range(self.num_controls+1):
			ix = self.get_ix(control_index=control_index)
			kept_ix = self.lcs[control_index].ix_inrange(colnames=['chi/N'],uplim=cut)
			cut_ix = AnotB(ix, kept_ix)
			self.update_mask_col(flag, cut_ix, control_index=control_index)
		
			if control_index == 0:
				num_cut = len(cut_ix)
				percent_cut = 100 * num_cut/len(ix)
				output1 = f'Chi-square cut {cut:0.2f} selected with {data["Pcontamination"]:0.2f}% contamination and {data["Ploss"]:0.2f}% loss'
				output2 = f'Total percent of data flagged ({hex(flag)}): {percent_cut:0.2f}%'
				print(f'# {output1}\n# {output2}')
				print('Success')
		return f'{output1}\n{output2}'
	
	def _do_x2_cut(self, flag, stn_cut, cut_start, cut_stop, cut_step, loss_lim, contam_lim, lim_to_prioritize, use_preSN_lc=False):
		print('\nCalculating chi-square cut...')

		if lim_to_prioritize != 'loss' and lim_to_prioritize != 'contamination':
			raise RuntimeError("ERROR: limit_to_prioritize must be 'loss' or 'contamination'!")
		
		print(f'# Contamination limit: {contam_lim:0.2f}%; loss limit: {loss_lim:0.2f}%')

		if use_preSN_lc:
			lc_temp = deepcopy(self.lcs[0])
			ix = lc_temp.ix_inrange('MJD', uplim=self.discdate)
		else:
			lc_temp = self._get_all_controls()
			ix = lc_temp.t.index.values
		good_ix, bad_ix = self._get_goodbad_ix(lc_temp, ix, stn_cut)
		
		limcuts = self._get_limcuts_table(lc_temp, ix, good_ix, bad_ix, cut_start, cut_stop, cut_step)
		if limcuts.empty:
			output = f'No cuts kept more than 10% of measurements--chi-square cut not applicable for this SN!'
			print(f'# {output}')
			return output

		contam_lim_cut, loss_lim_cut, contam_case, loss_case = self._get_limcuts(limcuts, loss_lim, contam_lim)
		contam_data = self._get_limcuts_data(lc_temp, contam_lim_cut, ix, good_ix, bad_ix)
		loss_data = self._get_limcuts_data(lc_temp, loss_lim_cut, ix, good_ix, bad_ix)

		print('# Contamination cut according to given contamination limit, with %0.2f%% contamination and %0.2f%% loss: %0.2f' % (contam_data['Pcontamination'], contam_data['Ploss'], contam_lim_cut))
		if contam_case == 'above lim':
			print('# WARNING: Contamination cut not possible with contamination <= contam_lim %0.1f!' % contam_lim)
		print('# Loss cut according to given loss limit, with %0.2f%% contamination and %0.2f%% loss: %0.2f' % (loss_data['Pcontamination'], loss_data['Ploss'], loss_lim_cut))
		if loss_case == 'above lim':
			print('# WARNING: Loss cut not possible with loss <= loss_lim %0.2f!' % loss_lim)

		final_cut = self._choose_btwn_cuts(contam_lim_cut, loss_lim_cut, contam_case, loss_case, loss_lim, contam_lim, cut_start, lim_to_prioritize)
		data = self._get_limcuts_data(lc_temp, final_cut, ix, good_ix, bad_ix)

		if np.isnan(final_cut):
			output = 'Final suggested chi-square cut could not be determined'
			print(f'# ERROR: {output}. We suggest rethinking your contamination and loss limits.')
			return output
		
		return self._apply_x2_cut(final_cut, data, flag)
	
	def apply_x2_cut(self, flag):
		return self._do_x2_cut(flag, 
							   self.cfg['x2_cut_params']['stn_cut'], 
							   self.cfg['x2_cut_params']['cut_start'], 
							   self.cfg['x2_cut_params']['cut_stop'], 
							   self.cfg['x2_cut_params']['cut_step'],
							   self.cfg['x2_cut_params']['loss_lim'],
							   self.cfg['x2_cut_params']['contam_lim'],
							   self.cfg['x2_cut_params']['lim_to_prioritize'],
							   self.cfg['x2_cut_params']['use_preSN_lc'])
	
	def _get_control_stats(self, flags):
		print('# Calculating control light curve statistics...')

		len_mjd = len(self.lcs[0].t['MJD'])

		# construct arrays for control lc data
		uJy = np.full((self.num_controls, len_mjd), np.nan)
		duJy = np.full((self.num_controls, len_mjd), np.nan)
		Mask = np.full((self.num_controls, len_mjd), 0, dtype=np.int32)
		
		for control_index in range(1, self.num_controls+1):
			if (len(self.lcs[control_index].t) != len_mjd) or (np.array_equal(self.lcs[0].t['MJD'], self.lcs[control_index].t['MJD']) is False):
				raise RuntimeError(f'ERROR: SN lc not equal to control lc for control_index {control_index}! Rerun or debug verify_mjds().')
			else:
				uJy[control_index-1,:] = self.lcs[control_index].t['uJy']
				duJy[control_index-1,:] = self.lcs[control_index].t[self.dflux_colnames[control_index]]
				Mask[control_index-1,:] = self.lcs[control_index].t['Mask']

		c2_param2columnmapping = self.lcs[0].intializecols4statparams(prefix='c2_',format4outvals='{:.2f}',skipparams=['converged','i'])

		for index in range(uJy.shape[-1]):
			pda4MJD = pdastrostatsclass()
			pda4MJD.t['uJy'] = uJy[0:,index]
			pda4MJD.t[self.dflux_colnames[0]] = duJy[0:,index]
			pda4MJD.t['Mask'] = np.bitwise_and(Mask[0:,index], flags['chisquare']|flags['uncertainty'])
			
			pda4MJD.calcaverage_sigmacutloop('uJy',
											noisecol=self.dflux_colnames[0],
											maskcol='Mask',
											maskval=(flags['chisquare']|flags['uncertainty']),
											verbose=1, Nsigma=3.0, median_firstiteration=True)
			self.lcs[0].statresults2table(pda4MJD.statparams, c2_param2columnmapping, destindex=index)
	
	def _apply_controls_cut(self, flags, x2_max, stn_max, Nclip_max, Ngood_min):
		print(f'\nApplying control light curve cut...')

		self._get_control_stats(flags)

		self.lcs[0].t['c2_abs_stn'] = self.lcs[0].t['c2_mean'] / self.lcs[0].t['c2_mean_err']

		# flag measurements according to given bounds
		flag_x2_ix = self.lcs[0].ix_inrange(colnames=['c2_X2norm'], lowlim=x2_max, exclude_lowlim=True)
		flag_stn_ix = self.lcs[0].ix_inrange(colnames=['c2_abs_stn'], lowlim=stn_max, exclude_lowlim=True)
		flag_nclip_ix = self.lcs[0].ix_inrange(colnames=['c2_Nclip'], lowlim=Nclip_max, exclude_lowlim=True)
		flag_ngood_ix = self.lcs[0].ix_inrange(colnames=['c2_Ngood'], uplim=Ngood_min, exclude_uplim=True)
		self.update_mask_col(flags['controls_x2'], flag_x2_ix)
		self.update_mask_col(flags['controls_stn'], flag_stn_ix)
		self.update_mask_col(flags['controls_Nclip'], flag_nclip_ix)
		self.update_mask_col(flags['controls_Ngood'], flag_ngood_ix)

		# update mask column with control light curve cut on any measurements flagged according to given bounds
		zero_Nclip_ix = self.lcs[0].ix_equal('c2_Nclip', 0)
		unmasked_ix = self.lcs[0].ix_unmasked('Mask', maskval=flags['controls_x2']|flags['controls_stn']|flags['controls_Nclip']|flags['controls_Ngood'])
		self.update_mask_col(flags['controls_questionable'], AnotB(unmasked_ix, zero_Nclip_ix))
		self.update_mask_col(flags['controls_bad'], AnotB(self.get_ix(),unmasked_ix))

		# copy over SN's control cut flags to control light curve 'Mask' column
		flags_arr = np.full(self.lcs[0].t['Mask'].shape, (flags['controls_bad']|
														  flags['controls_questionable']|
														  flags['controls_x2']|
														  flags['controls_stn']|
														  flags['controls_Nclip']|
														  flags['controls_Ngood']))
		flags_to_copy = np.bitwise_and(self.lcs[0].t['Mask'], flags_arr)
		for control_index in range(1,self.num_controls+1):
			self.lcs[control_index].t['Mask'] = self.lcs[control_index].t['Mask'].astype(np.int32)
			if len(self.lcs[control_index].t) < 1:
				continue
			elif len(self.lcs[control_index].t) == 1:
				self.lcs[control_index].t.loc[0,'Mask']= int(self.lcs[control_index].t.loc[0,'Mask']) | flags_to_copy
			else:
				self.lcs[control_index].t['Mask'] = np.bitwise_or(self.lcs[control_index].t['Mask'], flags_to_copy)

		len_ix = len(self.get_ix())
		s = []
		s.append('Percent of data above x2_max bound (%s): %0.2f%%' 
		   % (hex(flags['controls_x2']), 100 * len(self.lcs[0].ix_masked('Mask',maskval=flags['controls_x2'])) / len_ix))
		s.append('Percent of data above stn_max bound (%s): %0.2f%%' 
		   % (hex(flags['controls_stn']), 100 * len(self.lcs[0].ix_masked('Mask',maskval=flags['controls_stn'])) / len_ix))
		s.append('Percent of data above Nclip_max bound (%s): %0.2f%%' 
		   % (hex(flags['controls_Nclip']), 100 * len(self.lcs[0].ix_masked('Mask',maskval=flags['controls_Nclip'])) / len_ix))
		s.append('Percent of data below Ngood_min bound (%s): %0.2f%%' 
		   % (hex(flags['controls_Ngood']), 100 * len(self.lcs[0].ix_masked('Mask',maskval=flags['controls_Ngood'])) / len_ix))
		s.append('Total percent of data flagged as questionable (not masked with control light curve flags but Nclip > 0) (%s): %0.2f%%' 
		   % (hex(flags['controls_questionable']), 100 * len(self.lcs[0].ix_masked('Mask',maskval=flags['controls_questionable'])) / len_ix))
		s.append('Total percent of data flagged as bad (%s): %0.2f%%' 
		   % (hex(flags['controls_bad']), 100 * len(self.lcs[0].ix_masked('Mask',maskval=flags['controls_bad'])) / len_ix))
		print('# Control light curve cut results:')
		for out in s:
			print(f'## {out}')
		output = '\n'.join(s)
		
		print('Success')
		return output
	
	def apply_controls_cut(self, flags):
		return self._apply_controls_cut(flags, 
								  		self.cfg['controls_cut_params']['x2_max'],
										self.cfg['controls_cut_params']['stn_max'],
										self.cfg['controls_cut_params']['Nclip_max'],
										self.cfg['controls_cut_params']['Ngood_min'])
	
	def _average(self, avglc, flags, Nclip_max, Ngood_min, x2_max, control_index=0, mjd_bin_size=1.0, flux2mag_sigmalimit=3.0):
		if control_index == 0:
			print(f'Now averaging SN light curve...')
		else:
			print(f'Now averaging control light curve {control_index:03d}...')

		avglc.lcs[control_index] = pdastrostatsclass(columns=['MJD','MJDbin','uJy','duJy','stdev','x2','Nclip','Ngood','Nexcluded','Mask'],hexcols=['Mask'])

		mjd = int(np.amin(self.lcs[control_index].t['MJD']))
		mjd_max = int(np.amax(self.lcs[control_index].t['MJD']))+1

		#good_ix = self.lcs[control_index].ix_unmasked('Mask', maskval=flags['chisquare']|flags['uncertainty']|flags['controls_bad'])

		while mjd <= mjd_max:
			range_ix = self.lcs[control_index].ix_inrange(colnames=['MJD'], lowlim=mjd, uplim=mjd+mjd_bin_size, exclude_uplim=True)
			range_good_ix = self.lcs[control_index].ix_unmasked('Mask', 
													   maskval=flags['chisquare']|flags['uncertainty']|flags['controls_bad'], 
													   indices=range_ix) #AandB(good_ix,range_ix)

			# add new row to averaged light curve
			new_row = {'MJDbin':mjd+0.5*mjd_bin_size, 'Nclip':0, 'Ngood':0, 'Nexcluded':len(range_ix)-len(range_good_ix), 'Mask':0}
			avglc_index = avglc.lcs[control_index].newrow(new_row)
			
			# if no measurements present, flag or skip over day
			if len(range_ix) < 1:
				avglc.update_mask_col(flags['avg_badday'], [avglc_index], control_index=control_index, remove_old=False)
				mjd += mjd_bin_size
				continue
			
			# if no good measurements, average values anyway and flag
			if len(range_good_ix) < 1:
				# average flux
				self.lcs[control_index].calcaverage_sigmacutloop('uJy', noisecol=self.dflux_colnames[control_index], indices=range_ix, Nsigma=3.0, median_firstiteration=True)
				fluxstatparams = deepcopy(self.lcs[control_index].statparams)
				
				# get average mjd
				#indices=fluxstatparams['ix_good']
				self.lcs[control_index].calcaverage_sigmacutloop('MJD', indices=range_ix, Nsigma=0, median_firstiteration=False)
				avg_mjd = self.lcs[control_index].statparams['mean']

				# add row and flag
				avglc.lcs[control_index].add2row(avglc_index, {'MJD':avg_mjd, 
															   'uJy':fluxstatparams['mean'] if not fluxstatparams['mean'] is None else np.nan, 
															   'duJy':fluxstatparams['mean_err'] if not fluxstatparams['mean_err'] is None else np.nan, 
															   'stdev':fluxstatparams['stdev'] if not fluxstatparams['stdev'] is None else np.nan,
															   'x2':fluxstatparams['X2norm'] if not fluxstatparams['X2norm'] is None else np.nan,
															   'Nclip':fluxstatparams['Nclip'] if not fluxstatparams['Nclip'] is None else np.nan,
															   'Ngood':fluxstatparams['Ngood'] if not fluxstatparams['Ngood'] is None else np.nan,
															   'Mask':0})
				self.update_mask_col(flags['avg_badday'], range_ix, control_index=control_index, remove_old=False)
				avglc.update_mask_col(flags['avg_badday'], [avglc_index], control_index=control_index, remove_old=False)

				mjd += mjd_bin_size
				continue
			
			# average good measurements
			self.lcs[control_index].calcaverage_sigmacutloop('uJy', noisecol=self.dflux_colnames[control_index], indices=range_good_ix, Nsigma=3.0, median_firstiteration=True)
			fluxstatparams = deepcopy(self.lcs[control_index].statparams)

			if fluxstatparams['mean'] is None or len(fluxstatparams['ix_good']) < 1:
				self.update_mask_col(flags['avg_badday'], range_ix, control_index=control_index, remove_old=False)
				avglc.update_mask_col(flags['avg_badday'], [avglc_index], control_index=control_index, remove_old=False)
				mjd += mjd_bin_size
				continue

			# get average mjd
			# TO DO: SHOULD NOISECOL HERE BE DUJY OR NONE??
			self.lcs[control_index].calcaverage_sigmacutloop('MJD', noisecol=self.dflux_colnames[control_index], indices=fluxstatparams['ix_good'], Nsigma=0, median_firstiteration=False)
			avg_mjd = self.lcs[control_index].statparams['mean']

			# add row to averaged light curve
			avglc.lcs[control_index].add2row(avglc_index, {'MJD':avg_mjd, 
														   'uJy':fluxstatparams['mean'], 
														   'duJy':fluxstatparams['mean_err'], 
														   'stdev':fluxstatparams['stdev'],
														   'x2':fluxstatparams['X2norm'],
														   'Nclip':fluxstatparams['Nclip'],
														   'Ngood':fluxstatparams['Ngood'],
														   'Mask':0})
			
			# flag clipped measurements in lc
			if len(fluxstatparams['ix_clip']) > 0:
				self.update_mask_col(flags['avg_ixclip'], fluxstatparams['ix_clip'], control_index=control_index, remove_old=False)
			
			# if small number within this bin, flag measurements
			if len(range_good_ix) < 3:
				self.update_mask_col(flags['avg_smallnum'], range_ix, control_index=control_index, remove_old=False) # TO DO: CHANGE TO RANGE_IX??
				avglc.update_mask_col(flags['avg_smallnum'], [avglc_index], control_index=control_index, remove_old=False)
			# else check sigmacut bounds and flag
			else:
				is_bad = False
				if fluxstatparams['Ngood'] < Ngood_min:
					is_bad = True
				if fluxstatparams['Nclip'] > Nclip_max:
					is_bad = True
				if not(fluxstatparams['X2norm'] is None) and fluxstatparams['X2norm'] > x2_max:
					is_bad = True
				if is_bad:
					self.update_mask_col(flags['avg_badday'], range_ix, control_index=control_index, remove_old=False)
					avglc.update_mask_col(flags['avg_badday'], [avglc_index], control_index=control_index, remove_old=False)

			mjd += mjd_bin_size
		
		avglc.lcs[control_index].flux2mag('uJy','duJy','m','dm', zpt=23.9, upperlim_Nsigma=flux2mag_sigmalimit)
		avglc.drop_extra_columns(control_index=control_index)
		for col in ['Nclip','Ngood','Nexcluded','Mask']: 
			avglc.lcs[control_index].t[col] = avglc.lcs[control_index].t[col].astype(np.int32)
		return avglc

	def _apply_averaging(self, flags, Nclip_max, Ngood_min, x2_max, mjd_bin_size=1.0, flux2mag_sigmalimit=3.0):
		print('\nAveraging light curve(s)...')
		avglc = CleanAtlasLightCurve(self.filt, 
							   		 tnsname=self.tnsname, 
									 is_averaged=True, 
									 mjd_bin_size=mjd_bin_size, 
									 discdate=self.discdate)
		print('# Parameters: MJD bin size = %0.1f day(s), Nclip_max = %d, Ngood_min = %d, x2_max = %0.2f... ' % (mjd_bin_size, Nclip_max, Ngood_min, x2_max))

		for control_index in range(self.num_controls+1):
			avglc = self._average(avglc, flags, Nclip_max, Ngood_min, x2_max, control_index=control_index, mjd_bin_size=mjd_bin_size, flux2mag_sigmalimit=flux2mag_sigmalimit)

		s = 'Total percent of binned data flagged (%s): %0.2f%%' % (flags['avg_badday'], 100 * len(avglc.lcs[0].ix_masked('Mask',maskval=flags['avg_badday'])) / len(avglc.lcs[0].t))
		print(f'# {s}')
		output = f'{s}.'
		print('Success')
		return avglc, output

	def apply_averaging(self, flags):
		return self._apply_averaging(flags,
							   		 self.cfg['averaging_params']['Nclip_max'],
									 self.cfg['averaging_params']['Ngood_min'],
									 self.cfg['averaging_params']['x2_max'],
									 self.cfg['averaging_params']['mjd_bin_size'],
									 self.cfg['flux2mag_sigmalimit'])

class CleaningLoop():
	def __init__(self, args, settings):
		self.tnsnames = args.tnsnames
		self.settings = settings
		self.snlist = None

		self.lc_objs = {}

		# flags for each cut
		self.flags = {'chisquare':0x1,

					  'uncertainty':0x2,

					  'controls_bad':0x400000,
					  'controls_questionable':0x80000,
					  'controls_x2':0x100,
					  'controls_stn':0x200,
					  'controls_Nclip':0x400,
					  'controls_Ngood':0x800,

					  'avg_badday':0x800000,
					  'avg_ixclip':0x1000,
					  'avg_smallnum':0x2000}

		#self.apply_template_correction = args.template_correction
		self.apply_uncert_est = args.uncert_est
		self.apply_uncert_cut = args.uncert_cut
		self.apply_x2_cut = args.x2_cut
		self.apply_controls_cut = args.controls_cut
		self.apply_averaging = args.averaging

		self.plot = args.plot

	def load_snlist(self):
		if self.settings['snlist_filename'] != 'None':
			filename = f'{self.settings["output_dir"]}/{self.settings["snlist_filename"]}'
			print(f'\nLoading SN list at {filename}...')
			try:
				self.snlist = pd.read_table(filename, delim_whitespace=True)
				print('Success')
			except Exception as e:
				raise RuntimeError(f'# Could not load SN list at {filename}: {str(e)}')
		else:
			print(f'\nSkipping SN list at {self.settings["snlist_filename"]}...')

	def get_lc_data(self, tnsname, filt):
		k = f'{tnsname}_{filt}'
		# check if data exists in snlist
		if not self.snlist is None:
			row = self.snlist[self.snlist['tnsname'] == tnsname]
			if len(row) < 1:
				print()
				self.lc_objs[k].get_tns_data()
			else:
				print(f'\nSetting RA, Dec, and discovery date using {self.settings["snlist_filename"]}...')
				if len(row) > 1:
					# use first row
					row = row[0]
				index = row.index[0]
				self.lc_objs[k].set_tns_data(self.snlist.loc[index,'ra'], self.snlist.loc[index,'dec'], self.snlist.loc[index,'discovery_date'])
		else: # else get TNS data
			print()
			self.lc_objs[k].get_tns_data()
		print(self.lc_objs[k])

	def begin_readme(self, tnsname):
		f = open(f'{self.settings["output_dir"]}/{tnsname}/README.md','w+')
		f.write(f"# SN {tnsname} Light Curve Cleaning and Averaging")
		f.write(f'\n\nThe ATLAS SN light curves are separated by filter (orange and cyan) and labelled as such in the file name. Averaged light curves contain an additional number in the file name that represents the MJD bin size used. Control light curves are located in the "controls" subdirectory and follow the same naming scheme, only with their control index added after the SN name.')
		
		f.write(f'\n\nThe following details the file names for each of the light curve versions:')
		f.write(f'\n\t- SN light curves: {tnsname}.o.lc.txt and {tnsname}.c.lc.txt')
		if self.apply_averaging:
			f.write(f'\n\t- Averaged light curves: {tnsname}.o.{self.settings["averaging_params"]["mjd_bin_size"]:0.2f}days.lc.txt and {tnsname}.c.{self.settings["averaging_params"]["mjd_bin_size"]:0.2f}days.lc.txt')
		if self.apply_controls_cut:
			f.write(f'\n\t- Control light curves, where X=001,...,{self.settings["num_controls"]:03d}: {tnsname}_iX.o.lc.txt and {tnsname}_iX.c.lc.txt')

		f.write(f'\n\nThe following summarizes the hex values in the "Mask" column of each light curve for each cut applied (see below sections for more information on each cut): ')
		if self.apply_uncert_cut:
			f.write(f'\n\t- Uncertainty cut: {hex(self.flags["uncert"])}')
		if self.apply_x2_cut:
			f.write(f'\n\t- Chi-square cut: {hex(self.flags["x2"])}')
		if self.apply_controls_cut:
			f.write(f'\n\t- Control light curve cut: {hex(self.flags["controls_bad"])}')
		if self.apply_averaging:
			f.write(f'\n\t- Bad day (for averaged light curves): {hex(self.flags["avg_badday"])}')

		return f
	
	def add_to_readme(self, f, k, uncert_cut_output=None, uncert_est_output=None, x2_cut_output=None, controls_cut_output=None, averaging_output=None):
		f.write(f'\n\n## FILTER: {self.lc_objs[k].filt}')

		if self.apply_uncert_cut:
			f.write(f'\n\n### Uncertainty cut\n')
			f.write(uncert_cut_output)
		
		if self.apply_uncert_est:
			f.write(f'\n\n### True uncertainties estimation\n')
			f.write(uncert_est_output)
			
		if self.apply_x2_cut:
			f.write(f'\n\n### Chi-square cut\n')
			f.write(x2_cut_output)

		if self.apply_controls_cut:
			f.write(f'\n\n### Control light curve cut\n')
			f.write(controls_cut_output)

		f.write(f'\n\nAfter the cuts are applied, the light curves are resaved with the new "Mask" column.')
		percent_flagged = len(self.lc_objs[k].get_masked_ix(self.flags)) * 100 / len(self.lc_objs[k].lcs[0].t)
		f.write(f'\nTotal percent of data flagged as bad ({hex(self.flags["uncert"]|self.flags["x2"]|self.flags["controls_bad"]|self.flags["avg_badday"])}): {percent_flagged:0.2f}')
		
		if self.apply_averaging:
			f.write(f'\n\n### Averaging cleaned light curves\n')
			f.write(averaging_output)
			f.write(f'\nThe averaged light curves are then saved in a new file with the MJD bin size added to the filename.')

		return f

	def loop(self):
		self.load_snlist()

		for obj_index in range(len(self.tnsnames)):
			tnsname = self.tnsnames[obj_index]
			print(f'\nCOMMENCING LOOP FOR SN {tnsname}')

			f = self.begin_readme(tnsname)

			for filt in ['o', 'c']:
				k = f'{tnsname}_{filt}' # key for lc_objs dict using TNS name and current filter

				print(f'\nSETTING FILTER: {filt}')
				self.lc_objs[k] = CleanAtlasLightCurve(filt, cfg=self.settings, tnsname=tnsname)
				self.get_lc_data(tnsname, filt)

				self.lc_objs[k].load()

				print('\nPreparing for cleaning...')
				self.lc_objs[k].prep_for_cleaning()

				#template_correction_output = None
				#if self.apply_template_correction:
					#template_correction_output = self.lc_objs[k].apply_template_correction()

				uncert_cut_output = None
				if self.apply_uncert_cut:
					uncert_cut_output = self.lc_objs[k].apply_uncert_cut(self.flags['uncertainty'])
				
				uncert_est_output = None
				if self.apply_uncert_est:
					uncert_est_output = self.lc_objs[k].apply_uncert_est(self.flags['uncertainty'])

				x2_cut_output = None
				if self.apply_x2_cut:
					x2_cut_output = self.lc_objs[k].apply_x2_cut(self.flags['chisquare'])

				controls_cut_output = None
				if self.apply_controls_cut:
					controls_cut_output = self.lc_objs[k].apply_controls_cut(self.flags)

				averaging_output = None
				if self.apply_averaging:
					avglc, averaging_output = self.lc_objs[k].apply_averaging(self.flags)
					
					if self.settings['overwrite']:
						avglc._save(self.settings['output_dir'], filt=filt)

				if self.settings['overwrite']:
					self.lc_objs[k].save()

				f = self.add_to_readme(f, k,
						   			   uncert_cut_output=uncert_cut_output,
						   			   uncert_est_output=uncert_est_output,
									   x2_cut_output=x2_cut_output,
									   controls_cut_output=controls_cut_output,
									   averaging_output=averaging_output)

				if self.plot:
					print('WARNING: Plotting not implemented yet! Skipping...')

			f.close()

# define command line arguments
def define_args(parser=None, usage=None, conflict_handler='resolve'):
	if parser is None:
			parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
		
	parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
	parser.add_argument('--num_controls', type=int, default=None, help='Number of control light curves to load and clean')
	
	#parser.add_argument('-t', '--template_correction', default=False, action='store_true', help='apply automatic ATLAS template change correction')
	parser.add_argument('-e', '--uncert_est', default=False, action='store_true', help='apply true uncertainty estimation')
	parser.add_argument('-u', '--uncert_cut', default=False, action='store_true', help='apply uncertainty cut')
	parser.add_argument('-x', '--x2_cut', default=False, action='store_true', help='apply chi-square cut')
	parser.add_argument('-c', '--controls_cut', default=False, action='store_true', help='apply control light curve cut')
	parser.add_argument('-g', '--averaging', default=False, action='store_true', help='average light curves and cut bad days')
	parser.add_argument('-m', '--mjd_bin_size', type=float, default=None, help='MJD bin size in days for averaging')

	parser.add_argument('-p', '--plot', default=False, action='store_true', help='plot each cut and save into PDF file')
	parser.add_argument('--xlim_lower', type=float, default=None, help='if plotting, manually set lower x axis limit to a certain MJD')
	parser.add_argument('--xlim_upper', type=float, default=None, help='if plotting, manually set upper x axis limit to a certain MJD')
	parser.add_argument('--ylim_lower', type=float, default=None, help='if plotting, manually set lower y axis limit to a certain uJy')
	parser.add_argument('--ylim_upper', type=float, default=None, help='if plotting, manually set upper y axis limit to a certain uJy')

	parser.add_argument('-f','--cfg_filename', default='settings.ini', type=str, help='file name of ini file with settings for this class')
	parser.add_argument('-o','--overwrite', default=False, action='store_true', help='don\'t overwrite existing file with same file name')
		
	return parser

# load config file
def load_cfg(cfg_filename):
	print('LOADING SETTINGS FROM CONFIG FILE AND CMD ARGUMENTS...')
	cfg = configparser.ConfigParser()
	try:
		print(f'Loading config file at {cfg_filename}')
		cfg.read(cfg_filename)
	except Exception as e:
		raise RuntimeError(f'ERROR: Could not load config file at {cfg_filename}: {e}')
	return cfg

def load_settings():
	args = define_args().parse_args()
	cfg = load_cfg(args.cfg_filename)

	settings = {
		'output_dir': cfg['general']['output_dir'],
		'api_key': cfg['tns_cred']['api_key'],
		'tns_id': cfg['tns_cred']['tns_id'],
		'bot_name': cfg['tns_cred']['bot_name'],
		'snlist_filename': cfg['general']['snlist_filename'],
		'overwrite': args.overwrite,
		'num_controls': int(cfg['general']['num_controls']) if args.num_controls is None else args.num_controls,
		'flux2mag_sigmalimit': float(cfg['general']['flux2mag_sigmalimit'])
	}

	cleaning = CleaningLoop(args, settings)
	
	if cleaning.apply_uncert_est:
		cleaning.settings['prelim_x2_cut'] = bool(cfg['uncert_est']['prelim_x2_cut'])

	if cleaning.apply_uncert_cut:
		cleaning.settings['uncert_cut'] = float(cfg['uncert_cut']['cut'])

	if cleaning.apply_x2_cut:
		cleaning.settings['x2_cut'] = None
		if cfg['x2_cut']['override_cut'].isdigit():
			override_cut = float(cfg['x2_cut']['override_cut'])
			cleaning.settings['x2_cut'] = override_cut
		else:
			cleaning.settings['x2_cut_params'] = {
				'stn_cut': float(cfg['x2_cut']['stn_bound']),
				'cut_start': int(cfg['x2_cut']['min_cut']),
				'cut_stop': int(cfg['x2_cut']['max_cut']),
				'cut_step': int(cfg['x2_cut']['cut_step']),
				'contam_lim': float(cfg['x2_cut']['contamination_limit']),
				'loss_lim': float(cfg['x2_cut']['loss_limit']),
				'lim_to_prioritize': cfg['x2_cut']['limit_to_prioritize'],
				'use_preSN_lc': bool(cfg['x2_cut']['use_preSN_lc'])
			}
	
	if cleaning.apply_controls_cut:
		cleaning.settings['controls_cut_params'] = {
			'x2_max': float(cfg['controls_cut']['x2_max']),
			'stn_max': float(cfg['controls_cut']['stn_max']),
			'Nclip_max': int(cfg['controls_cut']['Nclip_max']),
			'Ngood_min': int(cfg['controls_cut']['Ngood_min'])
		}
	
	if cleaning.apply_averaging:
		cleaning.settings['averaging_params'] = {
			'mjd_bin_size': float(cfg['averaging']['mjd_bin_size']),
			'x2_max': float(cfg['averaging']['x2_max']),
			'Nclip_max': int(cfg['averaging']['Nclip_max']),
			'Ngood_min': int(cfg['averaging']['Ngood_min'])
		}
	
	if cleaning.plot:
		cleaning.settings['plot_params'] = {
			'xlim_lower': args.xlim_lower,
			'xlim_upper': args.xlim_upper,
			'ylim_lower': args.ylim_lower,
			'ylim_upper': args.ylim_upper
		}

	print(f'Success: {cleaning.settings}')
	return cleaning

if __name__ == "__main__":
	cleaning = load_settings()
	cleaning.loop()