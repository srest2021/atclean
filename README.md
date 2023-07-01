# Interactive Jupyter Notebooks and Python Scripts for Cleaning ATLAS Light Curves

## Table of Contents
- [Jupyter Notebooks](#jupyter-notebooks)
	- [`clean_atlas_lc.v2.ipynb`](#clean_atlas_lcv2ipynb) (for one SN, apply all cuts - chi-squares, uncertainties, control light curves - and average light curve)
- [Python Scripts](#python-scripts)
    - [Quick setup in `parans.ini`](#quick-setup-in-paramsini)
    - [`download_atlas_lc.py`](#download_atlas_lcpy) (for one or more SNe, download light curves and, optionally, control light curves)
    - [`clean_atlas_lc.py`](#clean_atlas_lcpy) (for one or more SNe, apply all cuts - chi-squares, uncertainties, control light curves - and average light curves)

## Jupyter Notebooks

### `clean_atlas_lc.v2.ipynb`
#### (applies all cuts--chi-squares, uncertainties, and control light curves--and applies averaging)
Using an already downloaded light curve, determine the best chi-square cut, apply the chi-square cut and an uncertainty cut, and average the light curve with bad day flagging. Then, save the light curve with the flags.

Control light curve cut currently in progress of being implemented and added to this notebook (this cut also requires pre-downloaded control light curves). To easily download control light curves in order to load them into this notebook, see the **`download_atlas_lc.py`** section to run this script.

## Python Scripts

### Quick setup in `params.ini`
Open the config file `params.ini` and replace the following fields with your information.
1. Replace `[ATLAS credentials]` `username` with your ATLAS username. You will be prompted for your ATLAS password when/if you run `download_atlas_lc.py`.
2. If given a TNS API key, the script will automatically fetch an object's RA, Dec, and discovery date from TNS. 
	- If you have a key, set `[TNS credentials]` `api_key` to your TNS API key. Then, set `[TNS credentials]` `tns_id` to the TNS ID of your bot and `[TNS credentials]` `bot_name` to the name of your bot.
	- If you don't have a key, you have two options:
		1. Obtain your own key from TNS. A key is obtained automatically once you create a bot, which you can do [here](https://www.wis-tns.org/bots) (log in, then click '+Add bot'). 
		2. Manually add this information to a table in a text file titled `snlist.txt`. (You can change this file's name in `[Input/output settings]` `snlist_filename`.) This text file is automatically generated inside the output directory after a single run of the script and stores infomation about SNe from previous iterations of the code; however, you can also edit/add in your own SN TNS names, coordinates, etc. It should include six columns (`tnsname`, `ra`, `dec`, `discovery_date`, `closebright_ra`, and `closebright_dec`), and empty cells should be marked as `NaN`. 

			Example `snlist.txt` file:
			```
			tnsname           ra          dec  discovery_date  closebright_ra  closebright_dec
			2023bee 8:56:11.6303 -3:19:32.095    59956.749940             NaN              NaN
			2017fra 15:31:51.530 +37:24:44.71    57939.340000             NaN              NaN
			2022oqm 15:09:08.211 +52:32:05.14    59751.190000             NaN              NaN
			2021foa 13:17:12.354 -17:15:25.77    59268.450000             NaN              NaN
			```
3. Replace `[Input/output settings]` `output_dir` with the directory address in which the light curve files and the `snlist.txt` file will be stored.
4. You can also change the sigma limit when converting flux to magnitude (magnitudes are limits when dmagnitudes are NaN). If you intend to download control light curves, you can change the radius of the circle pattern of locations around the SN and the total number of control light curves.

### `download_atlas_lc.py` 
#### (downloads SN light curve and, optionally, control light curves)
This script allows you to download ATLAS light curve(s) using [ATLAS's REST API](https://fallingstar-data.com/forcedphot/apiguide/) and [TNS's API](https://www.wis-tns.org/content/tns-getting-started) (to optionally fetch RA, Dec, and discovery date information for the SN). 

**Arguments** (will override default config file settings if specified):
- First provide TNS name(s) of object(s) to download
- `--coords`: manually input comma-separated RA and Dec coordinates of a single SN to override querying TNS and SN list text file
- `--discdate`: manually input discovery date in MJD of a single SN to override querying TNS and SN list text file
- `-c` or `--controls`: download control light curves in addition to SN light curve
- `--closebright`: manually input comma-separated RA and Dec coordinates of a close bright object interfering with the target SN's light curve and use as the center of the circle of control light curves
	- If you are downloading multiple SNe at once and have multiple sets of coordinates, they must be manually input to the SN list text file instead of this argument, as they cannot be automatically fetched from TNS.
- `--ctrl_coords`: manually input file name of txt file located within `output_dir` that contains control light curve coordinates; this txt file need only contain two space-separated columns titled "ra" and "dec" with the RA and Dec coordinates of the controls
- `-u` or `--username`: override default username given in `params.ini` config file
- `-a` or `--tns_api_key`: override default TNS API key given in `params.ini` config file
- `-f ` or `--cfg_filename`: provide a different config file filename (default is `params.ini`)
- `-l` or `--lookbacktime_days`: specify a lookback time in days (if not specified, script will download full light curve)
- `--dont_overwrite`: don't overwrite existing light curves with the same filename

**Example commands**:
- `download_atlas_lc.py 2019vxm` - downloads full SN 2019vxm light curve using ATLAS password 'XXX'
- `download_atlas_lc.py 2019vxm -l 100` - downloads SN 2019vxm light curve with a lookback time of 100 days
- `download_atlas_lc.py 2019vxm -c` downloads full SN 2019vxm and control light curves
- `download_atlas_lc.py 2019vxm 2020lse -c` downloads full SN and control light curves for SN 2019vxm AND SN 2020lse

### `clean_atlas_lc.py`
#### (applies all cuts - chi-squares, uncertainties, control light curves - and averages light curves)
Using the default settings in `params.ini`, load previously downloaded light curves and apply any of the chi-square, uncertainty, and control light curve cuts, average the light curves and flag bad days in both original and averaged light curves, then save both original and averaged light curves with the updated 'Mask' columns.

The **uncertainty cut** is a static procedure currently set at a constant value of 160. To change, set the `cut` field in `params.ini` to a different value.

We also attempt to **account for an extra noise source** in the data by estimating the true typical uncertainty, deriving the additional systematic uncertainty, and lastly **applying this extra noise to a new uncertainty column**. This new uncertainty column will be used in the cuts following this section. This procedure can be turned off or back on in `params.ini` through the `estimate_true_uncertainties` field. Here is the exact procedure we use:
1. Keep the uncertainty cut at 160 and apply a preliminary chi-square cut at 20. Filter out any measurements flagged by these two cuts.
2. Calculate the true typical uncertainty `sigma_true_typical` by taking a 3σ cut of the unflagged baseline flux and getting the standard deviation.
3. If `sigma_true_typical` is greater than the median uncertainty of the unflagged baseline flux, proceed with estimating the extra noise to add. Otherwise, keep the old chi-square cut and skip this procedure. 
4. Calculate the extra noise source `sigma_extra` using the following formula, where the median uncertainty `median_duJy` is taken from the unflagged baseline flux:
    - `sigma_extra^2 = sigma_true_typical^2 - median_duJy^2`
5. Apply `sigma_extra` to the existing uncertainty `duJy` and save into a new uncertainty column `duJy_new` using the following formula:
    - `duJy_new = sqrt(duJy^2 + sigma_extra^2)`
6. Repeat for each control light curve. For cuts following this procedure, use the new uncertainty column with the extra noise added instead of the old uncertainty column.

The **chi-square cut** procedure may be dynamic (default) or static. In order to apply a static cut at a constant value, set the `override_cut` parameter in the `Chi-square cut settings` section to that value; otherwise, leave set at `None` to apply the dynamic cut. More in-depth explanation of each parameter, its meaning, and overall procedures is located in **`clean_atlas_lc.v2.ipynb`**.

The **control light curve cut** uses a set of quality control light curves to determine the reliability of each SN measurement. Since we know that control light curve flux must be consistent with 0, any lack of consistency may indicate something wrong with the SN measurement at this epoch. We thus examine each SN epoch and its corresponding control light curve measurements at that epoch, apply a 3-sigma-clipped average, calculate statistics, and then cut bad epochs based on those returned statistics. We cut any measurements in the SN light curve for the given epoch for which statistics fulfill any of the following criteria (can be changed in `params.ini` under the correct section):
- A returned chi-square > 2.5 (to change, set field `x2_max`)
- A returned abs(flux/dflux) > 3.0 (to change, set field `stn_max`)
- Number of measurements averaged < 2 (to change, set field `Nclip_max`)
- Number of measurements clipped > 4 (to change, set field `Ngood_min`)

Note that this cut may not greatly affect certain SNe depending on the quality of the light curve. Its main purpose is to account for inconsistent flux in the case of systematic interference from bright objects, etc. that also affect the area around the SN. Therefore, normal SN light curves will usually see <1%-2% of data flagged as bad in this cut.

Our goal with the **averaging** procedure is to identify and cut out bad days by taking a 3σ-clipped average of each day. For each day, we calculate the 3σ-clipped average of any SN measurements falling within that day and use that average as our flux for that day. Because the ATLAS survey takes about 4 exposures every 2 days, we usually average together approximately 4 measurements per epoch (can be changed in `params.ini` by setting field `mjd_bin_size` to desired number of days). However, out of these 4 exposures, only measurements not cut in the previous methods are averaged in the 3σ-clipped average cut. (The exception to this statement would be the case that all 4 measurements are cut in previous methods; in this case, they are averaged anyway and flagged as a bad day.) Then we cut any measurements in the SN light curve for the given epoch for which statistics fulfill any of the following criteria (can be changed in `params.ini` under the correct section): 
- A returned chi-square > 4.0 (to change, set field `x2_max`)
- Number of measurements averaged < 2 (to change, set field `Nclip_max`)
- Number of measurements clipped > 1 (to change, set field `Ngood_min`)

For this part of the cleaning, we still need to improve the cutting at the peak of the SN (important epochs are sometimes cut, maybe due to fast rise, etc.).

The last optional step of this procedure is to **check for any pre-SN eruptions** (with specified size of `gaussian_sigma` days in `params.ini`) in the SN light curve. We apply a rolling gaussian weighted sum to the SN's flux/dflux ratio in order to amplify these possible precursor bumps, then plot. We can also apply this rolling sum to the control light curves (controlled by field `apply_to_controls` in `params.ini`) in order to establish the detection limit for this SN. Optionally, we can insert simulated precursor bump(s) into the SN light curve at specified MJDs, apparent magnitudes, and sigmas in order to test their amplification by the rolling sum. 

**Arguments** (will override default config file settings if specified):
- First provide TNS name(s) of object(s) to clean
- `-x` or `--chisquares`: apply chi-square cut
- `-u` or `--uncertainties`: apply uncertainty cut
- `-c` or `--controls`: apply control light curve cut
- `-g` or `--average`: average the light curve, cut bad days, and save as new file
	- `-m` or `--mjd_bin_size`: set MJD bin size in days
- `-b` or `--detect_bumps`: apply a rolling gaussian weighted sum to the SN's flux/dflux in order to amplify possible precursor bumps
	- If you add this argument, the script will additionally apply the previously specified cuts and averaging to the control light curves, as these are needed to identify the detection limit for that particular SN. So **if you have not downloaded control light curves for this SN**, you must change the `apply_to_controls` field in `params.ini` to `False`. Additionally, **if you have not added the `-g` averaging argument to the command, this procedure will not execute**. 
	- `--sim_gaussian`: comma-separated peakMJD list, peak_appmag, gaussian_sigma
		- This argument will allow you to add a gaussian at peakMJD with a peak apparent magnitude of peak_appmag and a sigma of gaussian_sigma in days.
		- This argument will allow you to simulate a pre-SN eruption within your light curve and analyze whether or not the rolling gaussian weighted sum successfully amplifies it. 
- `-p` or `--plot`: saves a PDF file of plots depicting the SN light curve, control light curves if necessary, and which measurements are flagged in each cut
	- You can use the arguments `--xlim_lower`, `--xlim_upper`, `--ylim_lower`, and `--ylim_upper` to set the x and y axis limits of the plots manually.
- `-f ` or `--cfg_filename`: provide a different config file filename (default is `params.ini`)
- `--dont_overwrite`: don't overwrite existing light curves with the same filename
- `-a` or `--tns_api_key`: override default TNS API key given in `params.ini` config file

**Example commands**:
- `clean_atlas_lc.py 2019vxm -x -u -c -g -p` - applies chi-square, uncertainty, and control light curve cuts to SN 2019vxm, averages the SN light curves and saves the averaged light curves, then saves plots of these cuts into PDF
- `clean_atlas_lc.py 2019vxm -x` - applies ONLY chi-square cut to SN 2019vxm