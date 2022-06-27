# Interactive Jupyter Notebooks and ATLAS Downloading Script for Cleaning ATLAS Light Curves

## Jupyter Notebooks

### `chisquarecut_master.ipynb` 
#### (applies chi-square-cut)
Using an already downloaded light curve, determine the best chi-square cut for that particular object. Then, save the light curve with the appropriate chi-square cut flag.

### `clean_atlas_lc.ipynb` 
#### (applies all cuts - chi-squares, uncertainties, control light curves, and averaging)
Using an already downloaded light curve, determine the best chi-square cut, apply the chi-square cut and an uncertainty cut, and average the light curve. Then, save the light curve with the flags.

Control light curve cut currently in progress of being implemented and added to this notebook (this cut also requires pre-downloaded control light curves). To easily download control light curves in order to load them into this notebook, see the **`download_atlas_lc.py`** section to run this script.

## Python Scripts

### `download_atlas_lc.py` 
#### (downloads SN light curve and, optionally, control light curves)
Download an ATLAS light curve using ATLAS's REST API (more information located here: https://fallingstar-data.com/forcedphot/apiguide/) and TNS's API. 

Configure default settings for downloading and saving light curves in **`atlaslc.ini`**. You must add your ATLAS credentials and TNS API key to this file for the script to work properly. Then, set the proper output directory for the files. You can also change the default flux and dflux column names as well as the sigma limit when converting flux to magnitude (magnitudes are limits when dmagnitudes are NaN). If you intend to download control light curves, you can change the radius of the circle pattern of locations around the SN and the total number of control light curves.

Arguments (will override default config file settings if specified):
- First provide TNS name(s) of object(s) to download
- `-c` or `--controls`: download control light curves in addition to SN light curve
- `-b` or `--closebright`: RA and Dec coordinates separated by comma of a close bright object interfering with the target object's light curve
	- These coordinates will become the center of the circle of control light curves
	- Note that you may only specify these for 1 object, so it is recommended to only download 1 object's group of light curves when using this command
- `-u` or `--username`: override default username given in `atlaslc.ini` config file
- `-p` or `--password`: REQUIRED! provide ATLAS password
- `-a` or `--tns_api_key`: override default TNS API key given in `atlaslc.ini` config file
- `-f ` or `--cfg_filename`: provide a different config file filename (default is `atlaslc.ini`)
- `-l` or `--lookbacktime_days`: specify a lookback time in days (if not specified, script will download full light curve)
- `-o` or `--dont_overwrite`: don't overwrite existing light curves with the same filename

Example commands:
- `download_atlas_lc.py 2019vxm -p 'XXX'` - downloads full SN 2019vxm light curve using ATLAS password 'XXX'
- `download_atlas_lc.py 2019vxm -l 100 -p 'XXX'` - downloads SN 2019vxm light curve with a lookback time of 100 days
- `download_atlas_lc.py 2019vxm -c -p 'XXX'` downloads full SN 2019vxm and control light curves
- `download_atlas_lc.py 2019vxm 2020lse -c -p 'XXX'` downloads full SN and control light curves for SN 2019vxm AND SN 2020lse