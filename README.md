# Interactive Jupyter Notebooks and Python Scripts for Cleaning ATLAS Light Curves

## Jupyter Notebooks

### `clean_atlas_lc.ipynb`
#### (applies all cuts--chi-squares, uncertainties, and control light curves--and applies averaging)
Using an already downloaded light curve, determine the best chi-square cut, apply the chi-square cut and an uncertainty cut, and average the light curve with bad day flagging. Then, save the light curve with the flags.

Control light curve cut currently in progress of being implemented and added to this notebook (this cut also requires pre-downloaded control light curves). To easily download control light curves in order to load them into this notebook, see the **`download_atlas_lc.py`** section to run this script.

## Python Scripts

### Quick setup in `atlas_lc_settings.ini`
Open the config file `atlas_lc_settings.ini` and replace the following fields with your information.
1. Replace `[ATLAS credentials]` `username` with your ATLAS username. You will be prompted for your ATLAS password when/if you run `download_atlas_lc.py`.
2. Replace `[TNS credentials]` `api_key` with your TNS API key (or ask Sofia to send you hers).
3. Replace `[Input/output settings]` `output_dir` with the directory address in which the light curve files will be stored.

### `download_atlas_lc.py` 
#### (downloads SN light curve and, optionally, control light curves)
Download an ATLAS light curve using ATLAS's REST API (more information located here: https://fallingstar-data.com/forcedphot/apiguide/) and TNS's API. 

Configure default settings for downloading and saving light curves in **`atlas_lc_settings.ini`**. You must add your ATLAS credentials and TNS API key to this file for the script to work properly. Then, set the proper output directory for the files. You can also change the sigma limit when converting flux to magnitude (magnitudes are limits when dmagnitudes are NaN). If you intend to download control light curves, you can change the radius of the circle pattern of locations around the SN and the total number of control light curves.

Arguments (will override default config file settings if specified):
- First provide TNS name(s) of object(s) to download
- `-c` or `--controls`: download control light curves in addition to SN light curve
- `-b` or `--closebright`: RA and Dec coordinates separated by comma of a close bright object interfering with the target object's light curve
	- These coordinates will become the center of the circle of control light curves.
	- Note that you may only specify these for 1 object, so it is recommended to only download 1 object's group of light curves when using this command.
	- If the script crashes in the middle of downloading control light curves, you can use the `--start_from` argument to specify a certain control index to start downloading from. For example, if you've already downloaded 3 control light curves and you get an error downloading control light curve 004, you can simply rerun the script with the argument `--start_from 4`.
- `-u` or `--username`: override default username given in `atlas_lc_settings.ini` config file
- `-a` or `--tns_api_key`: override default TNS API key given in `atlas_lc_settings.ini` config file
- `-f ` or `--cfg_filename`: provide a different config file filename (default is `atlas_lc_settings.ini`)
- `-l` or `--lookbacktime_days`: specify a lookback time in days (if not specified, script will download full light curve)
- `--dont_overwrite`: don't overwrite existing light curves with the same filename

Example commands:
- `download_atlas_lc.py 2019vxm` - downloads full SN 2019vxm light curve using ATLAS password 'XXX'
- `download_atlas_lc.py 2019vxm -l 100` - downloads SN 2019vxm light curve with a lookback time of 100 days
- `download_atlas_lc.py 2019vxm -c` downloads full SN 2019vxm and control light curves
- `download_atlas_lc.py 2019vxm 2020lse -c` downloads full SN and control light curves for SN 2019vxm AND SN 2020lse

### `clean_atlas_lc.py`
#### (applies all cuts - chi-squares, uncertainties, control light curves - and averages light curves)
Using the default settings in `atlas_lc_settings.ini`, load previously downloaded light curves and apply any of the chi-square, uncertainty, and control light curve cuts, average the light curves and flag bad days in both original and averaged light curves, then save both original and averaged light curves with the updated 'Mask' columns.

The chi-square cut procedure may be dynamic (default) or static. In order to apply a static cut at a constant value, set the `override_cut` parameter in the `Chi-square cut settings` section to that value; otherwise, leave set at `None` to apply the dynamic cut. More in-depth explanation of each parameter, its meaning, and overall procedures is located in **`clean_atlas_lc.ipynb`**.

The uncertainty cut is a static procedure currently set at a constant value of 160. To change, set the `cut` parameter in the `Uncertainty cut settings` section to a different value.

The control light curve cut examines each SN epoch and its corresponding control light curve measurements at that epoch, applies a 3-sigma-clipped average, calculates statistics, and then cuts bad epochs based on those returned statistics. The four parameters in these sections define the bounds on these returned statistics if a SN measurement is to be kept.

Our goal with the averaging procedure is to identify and cut out bad days by taking a 3σ-clipped average of each day. For each day, we calculate the 3σ-clipped average of any SN measurements falling within that day and use that average as our flux for that day. Because the ATLAS survey takes about 4 exposures every 2 days, we usually average together approximately 4 measurements per epoch (can be changed in `atlas_lc_settings.ini` by setting variable `mjd_bin_size` to desired number of days). However, out of these 4 exposures, only measurements not cut in the previous methods are averaged in the 3σ-clipped average cut. (The exception to this statement would be the case that all 4 measurements are cut in previous methods; in this case, they are averaged anyway and flagged as a bad day.) 

Then we cut any measurements in the SN light curve for the given epoch for which statistics fulfill any of the following criteria (can be changed in `atlas_lc_settings.ini`: 
- A returned chi-square > 4.0 (change variable `x2_max`)
- Number of measurements averaged < 2 (change variable `Nclip_max`)
- Number of measurements clipped > 1 (change variable `Ngood_min`)

For this part of the cleaning, we still need to improve the cutting at the peak of the SN (important epochs are sometimes cut, maybe due to fast rise, etc.).

Arguments (will override default config file settings if specified):
- First provide TNS name(s) of object(s) to clean
- `-x` or `--chisquares`: apply chi-square cut
- `-u` or `--uncertainties`: apply uncertainty cut
- `-c` or `--controls`: apply control light curve cut
- `-g` or `--average`: average the light curve, cut bad days, and save as new file
	- `-m` or `--mjd_bin_size`: set MJD bin size in days
- `-p` or `--plot`: saves a PDF file of plots depicting the SN light curve, control light curves if necessary, and which measurements are flagged in each cut
	- You can use the arguments `--xlim_lower`, `--xlim_upper`, `-ylim_lower`, and `--ylim_upper` to set the x and y axis limits of the plots manually.
- `-f ` or `--cfg_filename`: provide a different config file filename (default is `atlas_lc_settings.ini`)
- `--dont_overwrite`: don't overwrite existing light curves with the same filename
- `-a` or `--tns_api_key`: override default TNS API key given in `atlas_lc_settings.ini` config file

Example commands:
- `clean_atlas_lc.py 2019vxm -x -u -c -a -p` - applies chi-square, uncertainty, and control light curve cuts to SN 2019vxm, averages the SN light curves and saves the averaged light curves, then saves plots of these cuts into PDF
- `clean_atlas_lc.py 2019vxm -x` - applies ONLY chi-square cut to SN 2019vxm