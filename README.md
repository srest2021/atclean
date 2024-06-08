# ATClean

### Interactive Jupyter Notebooks and Python Scripts for Cleaning and Analyzing ATLAS Light Curves

View our paper [here](https://arxiv.org/abs/2405.03747) for details on our cleaning and binning method, as well as our pre-SN outburst detection process!

## Table of Contents

<!-- - [Python Scripts](#python-scripts)
    - [Setup in `config.ini`](#setup-in-configini): Set the default configuration for the `convert.py`, `download.py`, and `clean.py` scripts.
    - [`download.py`](#downloadpy): Download one or more SNe and their control light curves from the ATLAS forced photometry server.
    - [`convert.py`](#convertpy) (**WIP**): Convert a non-ATLAS light curve into an ATClean-readable format, so that it may be run through any of the following scripts.
    - [`clean.py`](#cleanpy): Apply one or more default and/or custom cuts and binning to one or more SNe and their control light curves. 
    - [`generate_sim_tables.py`](#generate_sim_tablespy): Part of our pre-SN outburst detection analysis. Generate tables of simulations (SimTables) by specifying the type of model and possible parameter values.
    - [`generate_detec_tables.py`](#generate_detec_tablespy): Part of our pre-SN outburst detection analysis. For each row in each SimTable, add the simulation to a random control light curve and record its max FOM and MJD, then update the rows and save as SimDetecTables. Optionally calculate efficiencies using specified FOM detection limits.

- [Jupyter Notebooks](#jupyter-notebooks)
    - [`clean.ipynb`](#cleanipynb) (**WIP**): An in-depth walkthrough of our cleaning and binning process for a single SN and its control light curves. 
    - [`atlas_lc_template_correction.ipynb`](#atlas_lc_template_correctionipynb) (**WIP**): A standalone walkthrough of our ATLAS template change correction.
    - [`simdetec_analysis.ipynb`](#simdetec_analysisipynb) (**WIP**): Part of our pre-SN outburst detection analysis. An in-depth walkthrough analysis of the generated SimDetecTables and efficiencies for a given SN and its control light curves.

- [Dependencies](#dependencies) -->

## Python Scripts

### Setup in `config.ini`

Open the `config.ini` file, which contains configuration for input/output directory paths, `convert.py`, `download.py`, and `clean.py`. The following toggles describe each field section by section. Bolded fields denote fields that we urge the user to change before attempting to use the scripts.

Note that these configurations may also be overridden by command line arguments.

#### Input/output directory paths: `dir` config section
- `raw_input`: This parameter sets the path to the directory containing raw TESS or other non-ATLAS light curves that are not in ATClean-readable format. These files are typically in their original downloaded format and need to be processed or converted by `convert.py` before they can be used by the ATClean software.

- **`atclean_input`**: This parameter specifies the path to the directory where light curves that are in ATClean-readable format are stored. These files have either been directly downloaded from the ATLAS server by `download.py` or converted from the raw format by `convert.py`.

- **`output`**: This parameter designates the path to the directory where all output data, including cleaned and binned light curves, plots, efficiency tables, and other results, will be saved. This directory serves as the main repository for the results produced by the ATClean pipeline.

- `sninfo_filename`: This parameter provides the name of the SN info file located inside the output directory. This space-separated `.txt` file contains essential information about SNe, including the TNS name and optionally the RA, Dec, and MJD0 (MJD at which the transient begins). This file may be provided manually with the correct column names (`tnsname`, `ra`, `dec`, and `mjd0`, with blank fields denoted by `NaN`) or generated and updated automatically if TNS credentials are provided.

#### Credentials for TNS API and ATLAS server: `credentials` config section

ATLAS credentials:

- **`atlas_username`**: The username required to download the ATLAS data from the server.

- **`atlas_password`**: The password required to download the ATLAS data from the server.

TNS credentials (these are optional, but if not provided, the user must provide a SN info file and set the `sninfo_filename` field):

- `tns_api_key`: The API key for accessing the TNS API.

- `tns_id`: The user ID for the TNS account.

- `tns_bot_name`: The name of the bot used for retrieving data from the TNS.

To set up a TNS bot and get an API key, login to TNS, then navigate [here](https://www.wis-tns.org/bots) and click on "Add bot".

### `download.py`

This script allows you to download ATLAS light curve(s) using [ATLAS's REST API](https://fallingstar-data.com/forcedphot/apiguide/) and [TNS's API](https://www.wis-tns.org/content/tns-getting-started) (to optionally fetch RA, Dec, and discovery date information for the SN). To download light curves from the ATLAS server, make sure to add your ATLAS username and password to the `credentials` section in `config.ini`. 

This script allows the user to download a single SN or a batch of SNe, as well as their optional control light curves. 
- To download control light curves, use the `-c` argument.
- To specify the number of control light curves to download, verify that the `num_controls` field in `config.ini` is set to the correct number, or use the `--num_controls` argument. 
- To specify the radius in arcseconds of the circle pattern of control light curves from the center location, verify that the `radius` field in `config.ini` is set to the correct number, or use the `--radius` argument. 

We additionally allow the user to change the center location of the control light curve circle pattern. By default, we use the SN location. However, this location can be changed via the `--closebright` argument. 

We allow the user to either specify certain RA and Dec coordinates, or to query the TNS API and automatically retrieve and save the coordinates. 
- To manually specify RA, Dec, and MJD0 for a single SN, use the `--coords` and `--mjd0` arguments. 
- To manually specify RA, Dec, and MJD0 for one or more SNe, you must provide a SN info file in the output directory. 
    - The file must be space-separated and include at least the following columns: `tnsname`, `ra`, `dec`, and `mjd0`. Any blank or unknown fields should be denoted by `NaN`. 
    - If TNS credentials are provided, the script will query TNS for any blank or unknown fields and update the SN info file with the missing information.
- To automatically retrieve RA, Dec, and MJD0 from TNS, simply provide your credentials in `config.ini`. The SN info file will be automatically generated and maintained inside the output directory. 

Lastly, we give users the chance to either specify the control light curve coordinates through a control coordinates table, or have them calculated automatically as a circle pattern.
- To specify all control light curve coordinates, you must provide a control coordinates file name using the `--ctrl_coords` argument.  
    - The file must be space-separated and include at least the following columns: `ra` and `dec`.
- To automatically calculate the control light curve coordinates, simply do not use the `--ctrl_coords` argument.

#### `download` config section in `config.ini`

- `flux2mag_sigmalimit`: The sigma limit used when converting flux to magnitude. Magnitudes are set as limits when their uncertainties are `NaN`.

- `num_controls`: The number of control light curves to be downloaded and used for analysis.

- `radius`: The radius in arcseconds for the circle pattern of control light curves around a center location (default is SN location). 

- `closebright_min_dist`: The minimum distance in arcseconds from the SN location to a control light curve location. This distance is used when the center of the circle pattern is set to a nearby bright object, and helps avoid any control locations landing on top of or too close to the SN.

#### Arguments
Arguments will override default config file settings if specified.
- First provide TNS name(s) of the object(s) to download. 
- `--sninfo_file`: Specifies the SN info file name. 
    - Type: str
    - Default: `None` (i.e., the `sninfo_filename` field in `config.ini`)
    - Usage: `--sninfo_file sninfo.txt`
- `--config_file`: Specifies the file name of the .ini file with settings for this script.
    - Type: str
    - Default: `config.ini`
    - Usage: `--config_file config.ini`
- `-l`, `--lookbacktime`: Sets the lookback time in days. 
    - Type: int
    - Default: `None` (i.e., look back as far as possible)
    - Usage: `-l 100` or `--lookbacktime 100`
- `--max_mjd`: Specifies the maximum MJD to download data up to.
    - Type: float
    - Default: `None`
    - Usage: `--max_mjd 59500.0`
- `-o`, `--overwrite`: If specified, existing files with the same name will be overwritten.
    - Type: bool
    - Default: `False`
    - Usage: `-o` or `--overwrite`
- `--coords`: comma-separated RA and Dec of the SN light curve to download.
    - Type: str
    - Default: `None` (i.e., reference the SN info file or query TNS)
    - Usage: `--coords 10.684,41.269`
- `--mjd0`: The start date of the SN in MJD.
    - Type: float
    - Default: `None` (i.e., reference the SN info file or query TNS)
    - Usage: `--mjd0 58800.0`
- `-c`, `--controls`: If specified, control light curves will be downloaded in addition to the SN light curve.
    - Type: bool
    - Default: `False`
    - Usage: `-c` or `--controls`
- `-n`, `--num_controls`: The number of control light curves to download per SN.
    - Type: int
    - Default: `None` (i.e., the `num_controls` field in `config.ini`)
    - Usage: `-n 5` or `--num_controls 5`
- `-r`, `--radius`: The radius of the control light curve circle pattern around the center location (by default the SN location), in arcseconds.
    - Type: float
    - Default: `None`
    - Usage: `-r 20.0` or `--radius 20.0`
- `--ctrl_coords`: Specifies the file name of the control coordinates table.
    - Type: str
    - Default: `None` (i.e., calculate the coordinates of each control light curve)
    - Usage: `--ctrl_coords ctrl_coords.txt`
- `--closebright`: Comma-separated RA and Dec of a nearby bright object interfering with the SN light curve. This object becomes the center of the control light curve circle pattern.
    - Type: str
    - Default: `None` (i.e., use the SN location as the center of the circle pattern)
    - Usage: `--closebright 10.684,41.269`

#### Example commands
- Download a single SN: `./download.py 2020lse -o`
- Download a batch of SNe: `./download.py 2020lse 2019vxm 2023ixf -o`
- Specify coordinates and MJD0: `./download.py 2020lse --coords 10:41:02.190,-27:05:00.42 --mjd0 58985.264 -o`
- Specify SN info table: `./download.py 2020lse 2019vxm 2023ixf --sninfo_file my_custom_SN_info.txt -o`
- Download only the last 10 days of data: `./download.py 2020lse -l 10 -o`
- Download control light curves: `./download.py 2020lse -c -o`
- Specify specify radius and number of control light curves: `./download.py 2020lse -c --radius 34 --num_controls 16 -o`
- Specify control light curve coordinates `./download.py 2020lse -c --ctrl_coords /path/to/control_coordinates.txt -o`
- Change the center location of the control light curve circle pattern: `./download.py 2020lse -c --closebright 10:41:02.290,-27:05:00.52 -o`

### `convert.py`

**WIP**

#### `convert` config section in `config.ini`

- `mjd_column_name`: The name of the column in the raw input data that contains MJD values.

- `flux_column_name`: The name of the column in the raw input data that contains flux values.

- `uncertainty_column_name`: The name of the column in the raw input data that contains flux error or uncertainty values.

- `chisquare_column_name` (optional, can be set to `None`): The name of the column in the raw input data that contains chi-square values.

- `filter_column_name` (optional, can be set to `None` for one filter): The name of the column in the raw input data that contains filter values.

- `filters`: A comma-separated list of filters as they appear in the filter column. If only one filter is used, provide a short identifier for that filter to be used in the filenames (for example, `tess` for TESS light curves). This parameter helps in distinguishing data from different filters or surveys.

#### Arguments
**WIP**

#### Example commands
**WIP**

### `clean.py`

**WIP**
<details>
<summary>uncert_est section (settings for true uncertainties estimation)</summary>

#### `uncert_est` section (settings for true uncertainties estimation in `clean.py`)

- `temp_x2_max_value`: A temporary, very high chi-square cut value used to eliminate the most egregious outliers from the data. This is an initial step in uncertainty estimation to ensure grossly incorrect data points are removed.

</details>

<details>
<summary>uncert_cut section (settings for uncertainty cut)</summary>

#### `uncert_cut` section (settings for the uncertainty cut in `clean.py`)

- `max_value`: The maximum allowable value for the uncertainties (`duJy` column). Measurements with uncertainties above this threshold will be flagged.

- `flag`: The flag value *in hex* assigned to measurements that exceed the maximum allowable uncertainty.

</details>

<details>
<summary>x2_cut section (settings for chi-square cut)</summary>

#### `x2_cut` section (settings for the chi-square cut in `clean.py`)

- `max_value`: The maximum allowable value for the PSF chi-squares (`chi/N` column). Measurements with a chi-square above this threshold will be flagged.

- `flag`: The flag value *in hex* assigned to measurements that exceed the maximum allowable chi-square.

We additionally calculate contamination and loss for a range of possible chi-square cuts in order to examine the effectiveness of the chosen cut:

- `stn_bound`: $|\frac{\text{flux}}{\text{dflux}}|$ bound that determines a "good" measurement vs. "bad" measurement.

- `min_cut`: The minimum chi-square value, inclusive, for which to calculate contamination and loss. 

- `min_cut`: The maximum chi-square value, inclusive, for which to calculate contamination and loss. 

- `cut_step`: The step size for the `[min_cut, max_cut]` range.

- `use_pre_mjd0_lc`: If set to `True`, we use the pre-MJD0 light curve to calculate contamination and loss. If set to `False` *(recommended)*, we instead use the control light curves.

</details>

<details>
<summary>controls_cut section (settings for control light curve cut)</summary>

#### `controls_cut` section (settings for the control light curve cut in `clean.py`)

- `bad_flag`: The flag value *in hex* assigned to measurements identified as bad. 

- `questionable_flag`: The flag value *in hex* assigned to measurements identified as questionable.

We perform a $3\sigma$-clipped average on each epoch across control light curves. The following criteria are used on the calculated statistics of each epoch, *not the statistics of the individual measurements,* to identify SN measurements in that epoch as "bad" or "questionable". Measurements violating *all* of the following thresholds will be flagged as "bad".

- `x2_max`: The chi-square threshold for an epoch.

    - `x2_flag`: The flag value *in hex* assigned to epochs with chi-square values exceeding `x2_max`.

- `stn_max`: The $|\frac{\text{flux}}{\text{dflux}}|$ threshold for an epoch.

    - `stn_flag`: The flag value *in hex* assigned to epochs with $|\frac{\text{flux}}{\text{dflux}}|$ exceeding `stn_max`.

- `Nclip_max`: The threshold of clipped control measurements for an epoch. 

    - `Nclip_flag`: The flag value *in hex* assigned to epochs with the number of clipped control measurements exceeding `Nclip_max`.

- `Ngood_min`: The threshold of good control measurements for an epoch. 

    - `Ngood_flag`: The flag value *in hex* assigned to epochs with the number of good control measurements falling below `Ngood_min`.

</details>

<details>
<summary>Custom cuts</summary>

#### Custom cuts section (additional static cuts that can be applied to any column of the light curve in `clean.py`)

Custom cuts run during cleaning allow you to define additional filtering criteria based on specific column values. We provide an example template for specifying custom cuts below. 

```
[example_cut]
    column: duJy
    max_value: 160
    min_value: None
    flag: 0x1000000
```

Note that the section title of the cut (in the above example, `[example_cut]`) must end in `_cut` in order to be properly parsed by `clean.py`.

- `column`: The name of the column in the data that the custom cut will be applied to.

- `max_value`: The maximum allowable value for the specified column. Measurements with values exceeding this threshold will be flagged. This parameter can be set to `None` if no upper limit is required.

- `min_value`: The minimum allowable value for the specified column. Measurements with values falling below this threshold will be flagged. This parameter can be set to `None` if no lower limit is required.

- `flag`: The flag value *in hex* assigned to measurements that do not meet the defined criteria. Flag values of `0x1000000` and above are available for use in custom cuts.

</details>

<details>
<summary>averaging section (settings for averaging and the bad day cut)</summary>

#### `averaging` section (settings for averaging light curves and applying the bad day cut)

For both the SN and control locations, we bin the light curve and perform a $3\sigma$-clipped average on the unflagged measurements in each bin. We use the calculated average and its error as flux and uncertainty values in the averaged light curves.

- `bad_flag`: The flag value *in hex* assigned to measurements identified as "bad".

- `ixclip_flag`: The flag value *in hex* assigned to bins where one or more measurements have been clipped during the $3\sigma$-clipped average.

- `smallnum_flag`: The flag value *in hex* assigned to bins with 2 or less unflagged measurements in a bin. 

- `mjd_bin_size`: The bin size in days.

The following criteria are used on the calculated statistics of each bin, *not the statistics of the individual measurements,* to identify measurements in that bin as "bad". Measurements violating *all* of the following thresholds will be flagged as "bad".

- `x2_max`: The chi-square threshold for a bin. 

- `Nclip_max`: The threshold of clipped measurements for a bin.

- `Ngood_min`: The threshold of good measurements for a bin. 

</details>

#### Arguments
**WIP**

#### Example commands
**WIP**

### `generate_sim_tables.py`

**WIP**

#### Configuration file: `simulation_settings.json`
**WIP**

#### Arguments
**WIP**

#### Example commands
**WIP**

### `generate_detec_tables.py`

**WIP**

#### Configuration file: `detection_settings.json`
**WIP**

#### Arguments
**WIP**

#### Example commands
**WIP**

## Jupyter Notebooks

### `clean.ipynb`
**WIP**

### `atlas_lc_template_correction.ipynb`
**WIP**

### `simdetec_analysis.ipynb`
**WIP**

## Dependencies
- Python 3.11.6 or higher
- typing
- os
- sys
- requests
- argparse
- configparser
- math
- pandas
- numpy
- getpass
- astropy
- abc
- re
- time
- json
- io
- copy
- itertools
- scipy