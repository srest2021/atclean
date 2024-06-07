# ATClean
### Interactive Jupyter Notebooks and Python Scripts for Cleaning and Analyzing ATLAS Light Curves

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

Open the `config.ini` file, which contains configuration for input/output directory paths, `convert.py`, `download.py`, and `clean.py`. The following describes each field section by section. Bolded fields denote fields that we highly recommend the user change before attempting to use the scripts.

#### `dir` section
- `raw_input`: This parameter sets the path to the directory containing raw TESS or other non-ATLAS light curves that are not in ATClean-readable format. These files are typically in their original downloaded format and need to be processed or converted by `convert.py` before they can be used by the ATClean software.
- **`atclean_input`**: This parameter specifies the path to the directory where light curves that are in ATClean-readable format are stored. These files have either been directly downloaded from the ATLAS server by `download.py` or converted from the raw format by `convert.py`.
- **`output`**: This parameter designates the path to the directory where all output data, including cleaned and binned light curves, plots, efficiency tables, and other results, will be saved. This directory serves as the main repository for the results produced by the ATClean pipeline.
- `sninfo_filename`: This parameter provides the name of the SN info file located inside the output directory. This space-separated `.txt` file contains essential information about SNe, including the TNS name and optionally the RA, Dec, and MJD0 (MJD at which the transient begins). This file may be provided manually with the correct column names (`tnsname`, `ra`, `dec`, and `mjd0`, with blank fields denoted by `NaN`) or updated automatically if TNS credentials are provided.