# ATClean

### Interactive Jupyter Notebooks and Python Scripts for Cleaning and Analyzing ATLAS Light Curves

View our paper [here](https://arxiv.org/abs/2405.03747) for details on our cleaning and binning method, as well as our pre-SN outburst detection process!

## Table of Contents

- [Python Scripts](#python-scripts)
    - [Setup in `config.ini`](#setup-in-configini)
    - [`download.py`](#downloadpy): Download one or more SNe and their control light curves from the ATLAS forced photometry server.
    - [`convert.py`](#convertpy) (**WIP**): Convert a non-ATLAS light curve into an ATClean-readable format, so that it may be run through any of the following scripts.
    - [`clean.py`](#cleanpy): Apply one or more default and/or custom cuts and binning to one or more SNe and their control light curves. 
    - [`generate_sim_tables.py`](#generate_sim_tablespy): Part of our pre-SN outburst detection analysis. Generate tables of simulations (SimTables) by specifying the type of model and possible parameter values.
    - [`generate_detec_tables.py`](#generate_detec_tablespy): Part of our pre-SN outburst detection analysis. For each row in each SimTable, add the simulation to a random control light curve and record its max FOM and MJD, then update the rows and save as SimDetecTables. Optionally calculate efficiencies using specified FOM detection limits.

- [Jupyter Notebooks](#jupyter-notebooks)
    - [`clean.ipynb`](#cleanipynb) (**WIP**): An in-depth walkthrough of our cleaning and binning process for a single SN and its control light curves. 
    - [`atlas_lc_template_correction.ipynb`](#atlas_lc_template_correctionipynb): A standalone walkthrough of our ATLAS template change correction.
    - [`simdetec_analysis.ipynb`](#simdetec_analysisipynb) (**WIP**): Part of our pre-SN outburst detection analysis. An in-depth walkthrough analysis of the generated SimDetecTables and efficiencies for a given SN and its control light curves.

- [Dependencies](#dependencies)

## Python Scripts

### Setup in `config.ini`

Open the `config.ini` file, which contains configuration for input/output directory paths, `convert.py`, `download.py`, and `clean.py`. The following toggles describe each field section by section. Bolded fields denote fields that we urge the user to change before attempting to use the scripts.

<details>
<summary>dir section</summary>

#### `dir` section (input/output directory paths)
- `raw_input`: This parameter sets the path to the directory containing raw TESS or other non-ATLAS light curves that are not in ATClean-readable format. These files are typically in their original downloaded format and need to be processed or converted by `convert.py` before they can be used by the ATClean software.

- **`atclean_input`**: This parameter specifies the path to the directory where light curves that are in ATClean-readable format are stored. These files have either been directly downloaded from the ATLAS server by `download.py` or converted from the raw format by `convert.py`.

- **`output`**: This parameter designates the path to the directory where all output data, including cleaned and binned light curves, plots, efficiency tables, and other results, will be saved. This directory serves as the main repository for the results produced by the ATClean pipeline.

- `sninfo_filename`: This parameter provides the name of the SN info file located inside the output directory. This space-separated `.txt` file contains essential information about SNe, including the TNS name and optionally the RA, Dec, and MJD0 (MJD at which the transient begins). This file may be provided manually with the correct column names (`tnsname`, `ra`, `dec`, and `mjd0`, with blank fields denoted by `NaN`) or generated and updated automatically if TNS credentials are provided.

</details>

<details>
<summary>convert section</summary>

#### `convert` section (settings for `convert.py`)

- `mjd_column_name`: The name of the column in the raw input data that contains MJD values.

- `flux_column_name`: The name of the column in the raw input data that contains flux values.

- `uncertainty_column_name`: The name of the column in the raw input data that contains flux error or uncertainty values.

- `chisquare_column_name` (optional, can be set to `None`): The name of the column in the raw input data that contains chi-square values.

- `filter_column_name` (optional, can be set to `None` for one filter): The name of the column in the raw input data that contains filter values.

- `filters`: A comma-separated list of filters as they appear in the filter column. If only one filter is used, provide a short identifier for that filter to be used in the filenames (for example, `tess` for TESS light curves). This parameter helps in distinguishing data from different filters or surveys.

</details>

<details>
<summary>credentials section (for TNS API and ATLAS server)</summary>

#### `credentials` section (for TNS API and ATLAS server)

ATLAS credentials:

- **`atlas_username`**: The username required to download the ATLAS data from the server.

- **`atlas_password`**: The password required to download the ATLAS data from the server.

TNS credentials (these are optional, but if not provided, the user must provide a SN info file and set the `sninfo_filename` field):

- `tns_api_key`: The API key for accessing the TNS API.

- `tns_id`: The user ID for the TNS account.

- `tns_bot_name`: The name of the bot used for retrieving data from the TNS.

To set up a TNS bot and get an API key, login to TNS, then navigate [here](https://www.wis-tns.org/bots) and click on "Add bot".

</details>

<details>
<summary>download section</summary>

#### `download` section (settings for `download.py`)

- `flux2mag_sigmalimit`: The sigma limit used when converting flux to magnitude. Magnitudes are set as limits when their uncertainties are `NaN`.

- `num_controls`: The number of control light curves to be downloaded and used for analysis.

- `radius`: The radius in arcseconds for the circle pattern of control light curves around a center location (default is SN location). 

- `closebright_min_dist`: The minimum distance in arcseconds from the SN location to a control light curve location. This distance is used when the center of the circle pattern is set to a nearby bright object, and helps avoid any control locations landing on top of or too close to the SN.

</details>

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

Note that the section title of the cut (in the above example, `[example_cut]`) must end in `_cut` in order to be properly parsed by `clean.py`.)

- `column`: The name of the column in the data that the custom cut will be applied to.

- `max_value`: The maximum allowable value for the specified column. Measurements with values exceeding this threshold will be flagged. This parameter can be set to `None` if no upper limit is required.

- `min_value`: The minimum allowable value for the specified column. Measurements with values falling below this threshold will be flagged. This parameter can be set to `None` if no lower limit is required.

- `flag`: The flag value assigned to measurements that do not meet the defined criteria. Flag values of `0x1000000` and above are available for use in custom cuts.

</details>

<details>
<summary>averaging section (settings for averaging and the bad day cut)</summary>

#### `averaging` section (settings for averaging light curves and applying the bad day cut)



</details>