# ATClean
### Interactive Jupyter Notebooks and Python Scripts for Cleaning ATLAS Light Curves

## Table of Contents
- [Jupyter Notebooks](#jupyter-notebooks)
	- [`clean_atlas_lc.v4.ipynb`](#clean_atlas_lcv4ipynb) (for one SN, apply all cuts - chi-squares, uncertainties, control light curves - and average light curve)
- [Python Scripts](#python-scripts)
    - [Quick setup in `settings.ini`](#quick-setup-in-settingsini)
    - [`download_atlas_lc.py`](#download_atlas_lcpy) (for one or more SNe, download light curves and, optionally, control light curves)
    - [`clean_atlas_lc_v2.py`](#clean_atlas_lcv2py) (for one or more SNe, apply all cuts - chi-squares, uncertainties, control light curves - and average light curves)

## Jupyter Notebooks

### `clean_atlas_lc.v4.ipynb`
#### Estimates true uncertainties, applies all cuts (chi-squares, uncertainties, control light curves), averages light curves, and corrects for ATLAS template changes.
Using previously downloaded SN and control light curves:
- Apply uncertainty cut
- Estimate true uncertainties
- Apply chi-square cut 
- Apply control light curve cut
- Average cleaned light curves
- Optionally correct for ATLAS template changes
- Save both original and averaged light curves with the updated 'Mask' columns

Example notebooks for example SNe are located in the `/extern` folder of this repository.

To easily download control light curves in order to load them into this notebook, see the **`download_atlas_lc.py`** section to run this script.

## Python Scripts

### Quick setup in `settings.ini`
Open the config file `settings.ini` and replace the following fields with your information.
1. Replace `[atlas_cred]` `username` with your ATLAS username. You will be prompted for your ATLAS password when/if you run `download_atlas_lc.py`.
2. If given a TNS API key, the script will automatically fetch an object's RA, Dec, and discovery date from TNS. 
	- If you have a key, set `[tns_cred]` `api_key` to your TNS API key. Then, set `[tns_cred]` `tns_id` to the TNS ID of your bot and `[tns_cred]` `bot_name` to the name of your bot.
	- If you don't have a key, you have two options:
		1. Obtain your own key from TNS. A key is obtained automatically once you create a bot, which you can do [here](https://www.wis-tns.org/bots) (log in, then click '+Add bot'). 
		2. Manually add this information to a table in a text file titled `sninfo.txt`. (You can change this file's name in `[general]` `sninfo_filename`.) This text file is automatically generated inside the output directory after a single run of the script and stores infomation about SNe from previous iterations of the code; however, you can also edit/add in your own SN TNS names, coordinates, etc. It should include six columns (`tnsname`, `ra`, `dec`, `discovery_date`, `closebright_ra`, and `closebright_dec`), and empty cells should be marked as `NaN`. 

			Example `sninfo.txt` file (also located in `extern/sninfo.txt`:
			```
			tnsname           ra          dec  discovery_date  closebright_ra  closebright_dec
			2023bee 8:56:11.6303 -3:19:32.095    59956.749940             NaN              NaN
			2017fra 15:31:51.530 +37:24:44.71    57939.340000             NaN              NaN
			2022oqm 15:09:08.211 +52:32:05.14    59751.190000             NaN              NaN
			2021foa 13:17:12.354 -17:15:25.77    59268.450000             NaN              NaN
			GRB220514B  246.5583144  61.03944515    60100.000000             NaN              NaN
			GRB200111A  99.29999197  37.07916637    60100.000000             NaN              NaN
			```
3. Replace `[general]` `output_dir` with the directory address in which the light curve files and the `sninfo.txt` file will be stored.
4. You can also change the sigma limit when converting flux to magnitude (magnitudes are limits when dmagnitudes are NaN). If you intend to download control light curves, you can change the radius of the circle pattern of locations around the SN and the total number of control light curves.

### `download_atlas_lc.py` 
#### Downloads SN light curve and, optionally, control light curves.
This script allows you to download ATLAS light curve(s) using [ATLAS's REST API](https://fallingstar-data.com/forcedphot/apiguide/) and [TNS's API](https://www.wis-tns.org/content/tns-getting-started) (to optionally fetch RA, Dec, and discovery date information for the SN). 

In order to change the number of control light curves downloaded, replace `[general]` `num_controls` with the number of desired control light curves.

**Arguments** (will override default config file settings if specified):
- First provide TNS name(s) of object(s) to download
- `--coords`: manually input comma-separated RA and Dec coordinates of a single SN to override querying TNS and SN list text file
- `--discdate`: manually input discovery date in MJD of a single SN to override querying TNS and SN list text file
- `-c` or `--controls`: download control light curves in addition to SN light curve
- `--closebright`: manually input comma-separated RA and Dec coordinates of a close bright object interfering with the target SN's light curve and use as the center of the circle of control light curves
	- If you are downloading multiple SNe at once and have multiple sets of coordinates, they must be manually input to the SN list text file instead of this argument, as they cannot be automatically fetched from TNS.
- `--ctrl_coords`: manually input file name of txt file located within `output_dir` that contains control light curve coordinates; this txt file need only contain two space-separated columns titled "ra" and "dec" with the RA and Dec coordinates of the controls (example file located in `extern/ctrl_coords.txt`)
- `-u` or `--username`: override default username given in `settings.ini` config file
- `-a` or `--tns_api_key`: override default TNS API key given in `settings.ini` config file
- `-f ` or `--cfg_filename`: provide a different config file filename (default is `settings.ini`)
- `-l` or `--lookbacktime_days`: specify a lookback time in days (if not specified, script will download full light curve)
- `--dont_overwrite`: don't overwrite existing light curves with the same filename

You can easily download light curves for a single SN without TNS credentials or `sninfo.txt` straight from the command line:

`./download_atlas_lc.py 2019vxm --coords 10:41:02.190,-27:05:00.42 --discdate 58985.264`

As well as its control light curves by adding -c:

`./download_atlas_lc.py 2019vxm -c --coords 10:41:02.190,-27:05:00.42 --discdate 58985.264`

**Other example commands**:
- `download_atlas_lc.py 2019vxm` - downloads full SN 2019vxm light curve using ATLAS password 'XXX'
- `download_atlas_lc.py 2019vxm -l 100` - downloads SN 2019vxm light curve with a lookback time of 100 days
- `download_atlas_lc.py 2019vxm -c` downloads full SN 2019vxm and control light curves
- `download_atlas_lc.py 2019vxm 2020lse -c` downloads full SN and control light curves for SN 2019vxm AND SN 2020lse

### `clean_atlas_lc_v2.py`
#### Estimates true uncertainties, applies all cuts (chi-squares, uncertainties, control light curves), and averages light curves.
Using the default settings in `settings.ini` and previously downloaded SN and control light curves:
- Apply uncertainty cut (`-u`)
- Estimate true uncertainties (`-e`)
- Apply chi-square cut (`-x`)
- Apply control light curve cut (`-c`)
- Average cleaned light curves (`-g`)
- Optionally correct for ATLAS template changes (**NOT CURRENTLY FUNCTIONAL**)
- Save both original and averaged light curves with the updated 'Mask' columns

**Arguments** (will override default config file settings if specified):
- First provide TNS name(s) of object(s) to clean
- `--num_controls`: set number of control light curves to load and clean
- `-u` or `--uncert_cut`: apply uncertainty cut
- `-e` or `--uncert_est`: apply uncertainties estimation
- `-x` or `--x2_cut`: apply chi-square cut
- `-c` or `--controls_cut`: apply control light curve cut
- `-g` or `--averaging`: average the light curve, cut bad days, and save as new file
	- `-m` or `--mjd_bin_size`: set MJD bin size in days
- `-p` or `--plot`: saves a PDF file of plots depicting the SN light curve, control light curves if necessary, and which measurements are flagged in each cut
	- You can use the arguments `--xlim_lower`, `--xlim_upper`, `--ylim_lower`, and `--ylim_upper` to set the x and y axis limits of the plots manually.
- `-f ` or `--cfg_filename`: provide a different config file filename (default is `settings.ini`)
- `-o` or `--overwrite`: overwrite existing light curves with the same filename

**Example commands**:
- `_clean_atlas_lc_v2.py 2019vxm -x -u -c -g -p -o` - applies chi-square, uncertainty, and control light curve cuts to SN 2019vxm and saves the light curves, averages the SN light curves and saves the averaged light curves, then saves plots of these cuts into PDF
- `clean_atlas_lc_v2.py 2019vxm -x -o` - applies ONLY chi-square cut to SN 2019vxm, then saves the light curves

#### Uncertainty cut (`-u`)
The **uncertainty cut** is a static procedure currently set at a default value of 160. To change, set the `[uncert_cut]` `cut` field in `settings.ini`.

#### True uncertainties estimation (`-e`)
We also attempt to **account for an extra noise source** in the data by estimating the true typical uncertainty, deriving the additional systematic uncertainty, and lastly **applying this extra noise to a new uncertainty column**. This new uncertainty column will be used in the cuts following this section. Here is the exact procedure we use:
1. Keep the previously applied uncertainty cut and apply a preliminary chi-square cut at 20 (default value; to change, set the `uncert_est` `prelim_x2_cut` field in `settings.ini`). Filter out any measurements flagged by these two cuts.
2.  Calculate the extra noise source for each control light curve using the following formula. The median uncertainty, `median_∂µJy`, is taken from the unflagged baseline flux. `σ_true_typical` is calculated by applying a 3-σ cut of the measurements cleaned in step 1, then getting the standard deviation.
    - `σ_extra^2 = σ_true_typical^2 - median_∂µJy^2`
3. Calculate the final extra noise source by taking the median of all `σ_extra`.
4. Decide whether or not to recommend addition of the extra noise source. First, get `σ_typical_old` by taking the median of the control light curves' `median_∂µJy`. Next, get `σ_typical_new` using the following formula:
    - `σ_typical_new^2 = σ_extra^2 + σ_typical_old^2`

    If `σ_typical_new` is 10% greater than `σ_typical_old`, recommend addition of the extra noise.
5. Apply the extra noise source to the existing uncertainty using the following formula:
    - `∂µJy_new^2 = ∂µJy_old^2 + σ_extra^2`
6. For cuts following this procedure, use the new uncertainty column with the extra noise added instead of the old uncertainty column.

#### Chi-square cut (`-x`)
The **chi-square cut** is a static procedure currently set at a default value of 5. To change, set the `[x2_cut]` `cut` field in `settings.ini`.

We use two factors, <strong>contamination</strong> and <strong>loss</strong>, to analyze the effectiveness of a PSF chi-square cut for the target SN, with flux/dflux as the deciding factor of what constitutes a good measurement vs. a bad measurement. 
- When calculating contamination and loss, we decide what will determine a good measurement vs. a bad measurement using a factor outside of the chi-square values. Our chosen factor is the absolute value of flux (µJy) divided by dflux (dµJy). The recommended boundary is a value of 3, such that any measurements with |µJy/dµJy| <= 3 are regarded as "good" measurements, and any measurements with |µJy/dµJy| > 3 are regarded as "bad" measurements. You can set this boundary to a different number by setting the `[x2_cut]` `stn_bound` field.
- We aim to separate good measurements from bad using the calculated chi-square cut by minimizing as much loss *and* contamination as possible. 
- We define contamination $C$ for a certain chi-square cut to be the number of bad kept measurements over the total number of kept measurements.
    - $C = N_{bad, kept}/N_{kept}$
- We define loss $L$ for a certain chi-square cut to be the number of good cut measurements over the total number of good measurements.
    - $L = N_{good,cut}/N_{good}$
- In order to provide informative output on contamination and loss, we calculate these measures for a range of possible chi-square cuts and optionally plot them (use `-p` argument to do so). 
    - We set the upper and lower bounds of a range of possible cuts for which to calculate $C$ and $L$. We start at a low value of 3 (to change, set field `[x2_cut]` `cut_start`) and end at 50 (to change, set field `[x2_cut]` `cut_stop`) with a step size of 1 (to change, set field `[x2_cut]` `cut_step`). For chi-square cuts falling on or between `cut_start` and `cut_stop` in increments of `cut_step`, we can begin to calculate contamination and loss percentages.
    - Since we can assume that the expected value of the control light curve flux is 0, we use these measurements by default to calculate and plot contamination and loss for the range of possible cuts. However, you can use the pre-SN flux instead by setting the `[x2_cut]` `use_preSN_lc` field.
- We output the calculated contamination and loss for the applied chi-square cut in a table with each SN's TNS name and filter in the output directory.

We set our default chi-square cut to 5, and defer overriding of that cut for a particular SN to the user, given the optional informative plot on alternative cuts with respect to contamination and loss. Again, the user can override this cut by changing the `[x2_cut]` `cut` field and rerunning the script.

<div class="alert alert-block alert-warning">
<b>Warning:</b> For very bright SNe, the chi-square values may increase during the SN even for good measurements due to imperfection in PSF fitting. Therefore, we recommend that the user double-check the chi-square values or the output plots to verify that the cut is working as intended, and override the cut with a custom value if needed.
</div>

#### Control light curve cut (`-c`)
The **control light curve cut** uses a set of quality control light curves to determine the reliability of each SN measurement. Since we know that control light curve flux must be consistent with 0, any lack of consistency may indicate something wrong with the SN measurement at this epoch. We thus examine each SN epoch and its corresponding control light curve measurements at that epoch, apply a 3-sigma-clipped average, calculate statistics, and then cut bad epochs based on those returned statistics. We cut any measurements in the SN light curve for the given epoch for which statistics fulfill any of the following criteria (fields can be changed in `settings.ini`):
- A returned chi-square > 2.5 (to change, set field `[controls_cut]` `x2_max`)
- A returned abs(flux/dflux) > 3.0 (to change, set field `[controls_cut]` `stn_max`)
- Number of measurements averaged < 2 (to change, set field `[controls_cut]` `Nclip_max`)
- Number of measurements clipped > 4 (to change, set field `[controls_cut]` `Ngood_min`)

Note that this cut may not greatly affect certain SNe depending on the quality of the light curve. Its main purpose is to account for inconsistent flux in the case of systematic interference from bright objects, etc. that also affect the area around the SN. Therefore, normal SN light curves will usually see <1%-2% of data flagged as bad in this cut.

#### Averaging and cutting bad days (`-g`)
Our goal with the **averaging** procedure is to identify and cut out bad days by taking a 3σ-clipped average of each day. For each day, we calculate the 3σ-clipped average of any SN measurements falling within that day and use that average as our flux for that day. Because the ATLAS survey takes about 4 exposures every 2 days, we usually average together approximately 4 measurements per epoch (can be changed in `settings.ini` by setting field `[averaging]` `mjd_bin_size` to desired number of days). However, out of these 4 exposures, only measurements not cut in the previous methods are averaged in the 3σ-clipped average cut. (The exception to this statement would be the case that all 4 measurements are cut in previous methods; in this case, they are averaged anyway and flagged as a bad day.) Then we cut any measurements in the SN light curve for the given epoch for which statistics fulfill any of the following criteria (can be changed in `settings.ini` under `[averaging]`): 
- A returned chi-square > 4.0 (to change, set field `[averaging]` `x2_max`)
- Number of measurements averaged < 2 (to change, set field `[averaging]` `Nclip_max`)
- Number of measurements clipped > 1 (to change, set field `[averaging]` `Ngood_min`)

For this part of the cleaning, we still need to improve the cutting at the peak of the SN (important epochs are sometimes cut, maybe due to fast rise, etc.).