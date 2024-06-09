#!/usr/bin/env python

from abc import ABC, abstractmethod
import itertools
import math
import os
import random
import argparse, re
from copy import deepcopy
import sys
from typing import Dict, List, Self, Tuple
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import root
from astropy.modeling.functional_models import Gaussian1D

from download import make_dir_if_not_exists
from pdastro import pdastrostatsclass
from generate_sim_table import (
    SimTable,
    SimTables,
    load_json_config,
    parse_params,
    GAUSSIAN_MODEL_NAME,
    ASYMMETRIC_GAUSSIAN_MODEL_NAME,
    print_progress_bar,
)
from lightcurve import SimDetecLightCurve, SimDetecSupernova, Simulation


NON_PARAM_COLNAMES = [
    "sigma_kern",
    "peak_appmag",
    "peak_flux",
    "filename",
    "model_name",
    "mjd_colname",
    "mag_colname",
    "flux_colname",
    "control_index",
    "max_fom",
    "max_fom_mjd",
]


# convert flux to magnitude
def flux2mag(flux: float):
    return -2.5 * np.log10(flux) + 23.9


# convert magnitude to flux
def mag2flux(mag: float):
    return 10 ** ((mag - 23.9) / -2.5)


class AsymmetricGaussian(Simulation):
    def __init__(self, model_name: str = ASYMMETRIC_GAUSSIAN_MODEL_NAME, **kwargs):
        """
        Initialize an asymmetric Gaussian simulation object.

        :param sigma_plus: Sigma or kernel size of one half of the Gaussian.
        :param sigma_minus: Sigma or kernel size of the other half of the Gaussian.
        :param peak_appmag: Peak apparent magnitude of the Gaussian.
        :param model_name: Name of the Gaussian model in the config file.
        """
        Simulation.__init__(self, model_name=model_name, **kwargs)
        self.g = None
        self.sigma_plus: float = None
        self.sigma_minus: float = None

    def new(self, sigma_plus: float, sigma_minus: float, peak_appmag: float):
        self.sigma_plus = sigma_plus
        self.sigma_minus = sigma_minus
        self.peak_appmag = peak_appmag

        peak_flux = mag2flux(peak_appmag)
        x = np.arange(-100, 100, 0.01)
        g1 = Gaussian1D(amplitude=peak_flux, stddev=sigma_minus)(x)
        g2 = Gaussian1D(amplitude=peak_flux, stddev=sigma_plus)(x)

        ind = np.argmin(abs(x))
        g3 = np.copy(g1)
        g3[ind:] = g2[ind:]

        self.g = np.array([x, g3])

    def get_sim_flux(
        self,
        mjds,
        peak_appmag: float,
        sigma_sim_plus: float = None,
        sigma_sim_minus: float = None,
        peak_mjd: float = None,
        **kwargs,
    ):
        """
        Get the interpolated function of the Gaussian at a given peak MJD and match it to the given time array.

        :param mjds: Time array of MJDs.
        :param peak_appmag: Desired peak apparent magnitude of the Gaussian.
        :param sigma: Desired sigma or kernel size of the Gaussian.
        :param peak_mjd: MJD at which the Gaussian should reach its peak apparent magnitude.

        :return: The simulated flux corresponding to the given time array.
        """
        if sigma_sim_plus is None or sigma_sim_minus is None:
            raise RuntimeError(
                "ERROR: sim_sigma_plus and sim_sigma_minus required to get flux of simulated asymmetric Gaussian."
            )
        if peak_mjd is None:
            raise RuntimeError(
                "ERROR: Peak MJD required to get flux of simulated asymmetric Gaussian."
            )

        self.new(sigma_sim_plus, sigma_sim_minus, peak_appmag)

        g = deepcopy(self.g)
        g[0, :] += peak_mjd

        fn = interp1d(g[0], g[1], bounds_error=False, fill_value=0)
        sim_flux = fn(mjds)
        return sim_flux

    def __str__(self):
        return (
            super().__str__()
            + f", sigma plus = {self.sigma_plus}, sigma_minus = {self.sigma_minus}"
        )


class Gaussian(AsymmetricGaussian):
    def __init__(self, model_name: str = GAUSSIAN_MODEL_NAME, **kwargs):
        """
        Initialize a Gaussian simulation object.

        :param sigma: Sigma or kernel size of the Gaussian.
        :param peak_appmag: Peak apparent magnitude of the Gaussian.
        :param model_name: Name of the Gaussian model in the config file.
        """
        AsymmetricGaussian.__init__(self, model_name=model_name, **kwargs)

    def get_sim_flux(
        self,
        mjds,
        peak_appmag: float,
        sigma_sim: float = None,
        peak_mjd: float = None,
        **kwargs,
    ):
        return super().get_sim_flux(
            mjds,
            peak_appmag,
            sigma_sim_plus=sigma_sim,
            sigma_sim_minus=sigma_sim,
            peak_mjd=peak_mjd,
            **kwargs,
        )

    def __str__(self):
        return Simulation().__str__() + f", sigma = {self.sigma_plus}"


class Model(Simulation):
    def __init__(
        self,
        filename: str,
        mjd_colname: str = False,
        mag_colname: str = False,
        flux_colname: str = False,
        model_name: str = "pre_SN_outburst",
        **kwargs,
    ):
        """
        Initialize a model simulation object.

        :param filename: File name of the model to load.
        :param mjd_colname: MJD column name in the model file (None if present but no column name; False if not present).
        :param mag_colname: Magnitude column name in the model file (None if present but no column name; False if not present).
        :param flux_colname: Flux column name in the model file (None if present but no column name; False if not present).
        :param model_name: Name of the model assigned in the config file.
        """
        Simulation.__init__(self, model_name=model_name, **kwargs)
        self.t = None

        self.load(
            filename,
            mjd_colname=mjd_colname,
            mag_colname=mag_colname,
            flux_colname=flux_colname,
        )

    def load(
        self,
        filename: str,
        mjd_colname: str = False,
        mag_colname: str = False,
        flux_colname: str = False,
        verbose: bool = False,
    ):
        """
        Load the model from a file into a DataFrame.
        Discern which column is which using the given column names.
        If necessary, create any missing MJD, magnitude, or flux columns.

        :param filename: File name of the model to load.
        :param mjd_colname: MJD column name in the model file (None if present but no column name; False if not present).
        :param mag_colname: Magnitude column name in the model file (None if present but no column name; False if not present).
        :param flux_colname: Flux column name in the model file (None if present but no column name; False if not present).
        """
        if verbose:
            print(f"\nLoading model at {filename}...")

        if mag_colname is False and flux_colname is False:
            raise RuntimeError(
                f"ERROR: Model must have either mag or flux column. Please set one or both fields to null or the correct column name."
            )

        try:
            header = "infer"
            if not mjd_colname and not mag_colname and not flux_colname:
                # all three column names are null or false
                header = None
            self.t = pd.read_table(filename, delim_whitespace=True, header=header)
        except Exception as e:
            raise RuntimeError(f"ERROR: Could not load model at {filename}: {str(e)}")

        if mjd_colname is False:
            # create MJD column and make it the first column
            columns = ["MJD"] + self.t.columns
            self.t["MJD"] = range(len(self.t))
            self.t = self.t[columns]
        else:
            # rename column to "MJD"
            self.t.rename(
                columns={0 if mjd_colname is None else mjd_colname: "MJD"}, inplace=True
            )

        if mag_colname is False:
            # flux column must be present
            # rename flux column to "uJy"
            self.t.rename(
                columns={1 if flux_colname is None else flux_colname: "uJy"},
                inplace=True,
            )
            # create mag column
            self.t["m"] = self.t["uJy"].apply(lambda flux: flux2mag(flux))
        else:
            # rename mag column to "m"
            self.t.rename(
                columns={1 if mag_colname is None else mag_colname: "m"}, inplace=True
            )
            if flux_colname is False:
                # create flux column
                self.t["uJy"] = self.t["m"].apply(lambda mag: mag2flux(mag))
            else:
                # rename flux column to "uJy"
                self.t.rename(
                    columns={2 if flux_colname is None else flux_colname: "uJy"},
                    inplace=True,
                )

        if verbose:
            print(self.t[["MJD", "m", "uJy"]].head().to_string())
            print("Success")

    def get_sim_flux(self, mjds, peak_appmag: float, peak_mjd: float = None, **kwargs):
        """
        Get the interpolated function of the model at a given peak MJD and peak apparent magnitude and match it to the given time array.

        :param mjds: Time array of MJDs.
        :param peak_appmag: Desired peak apparent magnitude of the model.
        :param peak_mjd: MJD at which the model should reach its peak apparent magnitude.

        :return: The simulated flux array corresponding to the given time array.
        """
        self.peak_appmag = peak_appmag
        if peak_mjd is None:
            raise RuntimeError("ERROR: Peak MJD required to construct simulated model.")

        # get original peak appmag index
        peak_idx = self.t["m"].idxmin()

        # scale flux to the desired peak appmag
        self.t["uJy"] *= mag2flux(peak_appmag) / self.t.loc[peak_idx, "uJy"]

        # recalulate appmag column
        self.t["m"] = self.t["uJy"].apply(lambda flux: flux2mag(flux))

        # put peak appmag at peak_mjd
        self.t["MJD"] -= self.t.loc[peak_idx, "MJD"]
        self.t["MJD"] += peak_mjd

        # interpolate lc and match to time array
        fn = interp1d(self.t["MJD"], self.t["uJy"], bounds_error=False, fill_value=0)
        sim_flux = fn(mjds)
        return sim_flux

    def __str__(self):
        return super().__str__()


class SimDetecTable(SimTable):
    def __init__(self, sigma_kern: float, peak_appmag: float, **kwargs):
        SimTable.__init__(self, peak_appmag, **kwargs)
        self.sigma_kern: float = sigma_kern

    def get_params_at_index(self, index: int, time_colname: str = None) -> Dict:
        """
        Get a dictionary of the parameter column-value pairs of the Simulation object at a certain row.
        Any known non-parameter column names will be skipped.

        :param index: Index of the table from which to get the parameter column-value pairs.
        :param time_colname: Any peak MJD, MJD0, or time-related parameter name that denotes where to inject the Simulation. If provided, will be additionally skipped when getting the column-value pairs.
        """
        if time_colname is None:
            skip_params = NON_PARAM_COLNAMES
        else:
            skip_params = NON_PARAM_COLNAMES.append(time_colname)

        colnames = []
        for colname in self.t.columns:
            if not colname in skip_params:
                colnames.append(colname)
        return dict(self.t.loc[index, colnames])

    def update_row_at_index(self, index: int, data: Dict):
        """
        Update a certain row of the table.

        :param index: Index of the table to update.
        :param data: Dictionary of column-value pairs.
        """
        # self.t.loc[index, data.keys()] = np.array(list(data.values()))
        for key, value in data.items():
            self.t.at[index, key] = value

    def get_detec_filename(self, model_name: str, detec_tables_dir: str) -> str:
        """
        Get the filename of the SimDetecTable.

        :param model_name: Name of the model of which the SimDetecTable contains simulations.
        :param detec_tables_dir: Directory where the SimDetecTable is located.
        """
        return f"{detec_tables_dir}/simdetec_{model_name}_{self.sigma_kern}_{self.peak_appmag:0.2f}.txt"

    def load_detec_table(self, model_name: str, detec_tables_dir: str):
        """
        Load an existing SimDetecTable.

        :param model_name: Name of the model of which the SimDetecTable contains simulations.
        :param detec_tables_dir: Directory where the SimDetecTable is located.
        """
        filename = self.get_detec_filename(model_name, detec_tables_dir)
        try:
            self.load_spacesep(filename, delim_whitespace=True)
        except Exception as e:
            raise RuntimeError(
                f"ERROR: Could not load SimDetecTable at {filename}: {str(e)}"
            )

    def load_from_sim_table(self, model_name: str, sim_tables_dir: str):
        """
        Load an existing SimTable and turn it into a SimDetecTable.

        :param model_name: Name of the model of which the SimTable contains simulations.
        :param sim_tables_dir: Directory where the SimTable is located.
        """
        super().load_sim_table(model_name, sim_tables_dir)
        self.t["sigma_kern"] = self.sigma_kern

    def save_detec_table(self, model_name: str, detec_tables_dir: str):
        """
        Save the current SimDetecTable.

        :param model_name: Name of the model of which the SimDetecTable contains simulations.
        :param detec_tables_dir: Directory where the SimDetecTable should be saved.
        """
        filename = self.get_detec_filename(model_name, detec_tables_dir)
        self.write(filename=filename, overwrite=True, index=False)

    def get_efficiency(self, fom_limit: float, **params):
        """
        Get the efficiency where columns match all the given values and are within all the given ranges.

        :param params: Arbitrary number of pairs of column = value, column = range, or column = list of ranges.
        Example usage for columns A, B, C: self.get_efficiency(10.0, A=2, B=[5, 6], C=[[1, 2], [3, 4]])

        :return: Efficiency of the rows that match the criteria.
        """
        col_ix = self.getindices()
        for column, value in params.items():
            if isinstance(value, list):
                if all(isinstance(v, list) for v in value):
                    mask = self.t.loc[col_ix, column].apply(
                        lambda x: any(a <= x <= b for a, b in value)
                    )
                    col_ix = self.t.index[mask].tolist()
                else:
                    col_ix = self.ix_inrange(
                        colnames=column, lowlim=value[0], uplim=value[1], indices=col_ix
                    )
            else:
                col_ix = self.ix_equal(column, value, indices=col_ix)

        detected_ix = self.ix_inrange("max_fom", lowlim=fom_limit, indices=col_ix)
        efficiency = 100 * len(detected_ix) / len(col_ix)
        return efficiency


class SimDetecTables(SimTables):
    def __init__(self, peak_appmags: List, model_name: str, sigma_kerns: List):
        SimTables.__init__(self, peak_appmags, model_name)
        self.sigma_kerns: List = sigma_kerns
        self.d: Dict[float, Dict[float, SimDetecTable]] = {}

    def get_table(self, sigma_kern: float, peak_appmag: float):
        return self.d[sigma_kern][peak_appmag]

    def update_row_at_index(
        self, sigma_kern: float, peak_appmag: float, index: int, data: Dict
    ):
        self.d[sigma_kern][peak_appmag].update_row_at_index(index, data)

    def get_efficiency(
        self, sigma_kern: float, peak_appmag: float, fom_limit: float, **params
    ):
        return self.d[sigma_kern][peak_appmag].get_efficiency(fom_limit, **params)

    def save_detec_table(
        self, sigma_kern: float, peak_appmag: float, detec_tables_dir: str
    ):
        self.d[sigma_kern][peak_appmag].save_detec_table(
            self.model_name, detec_tables_dir
        )

    def save_all(self, detec_tables_dir: str):
        print(f"\nSaving SimDetecTables in directory: {detec_tables_dir}")
        make_dir_if_not_exists(detec_tables_dir)
        for sigma_kern in self.d.keys():
            for table in self.d[sigma_kern].values():
                table.save_detec_table(self.model_name, detec_tables_dir)
        print("Success")

    def load_all_from_sim_tables(self, sim_tables_dir: str):
        """
        Load existing SimTables and turn them into SimDetecTables.

        :param sim_tables_dir: Directory where the SimTables are located.
        """
        print(
            f"\nConstructing SimDetecTables from existing SimTables in directory: {sim_tables_dir}"
        )
        for sigma_kern in self.sigma_kerns:
            self.d[sigma_kern] = {}
            for peak_appmag in self.peak_appmags:
                self.d[sigma_kern][peak_appmag] = SimDetecTable(sigma_kern, peak_appmag)
                self.d[sigma_kern][peak_appmag].load_from_sim_table(
                    self.model_name, sim_tables_dir
                )
        print("Success")

    def load_all(self, detec_tables_dir: str):
        """
        Load existing SimDetecTables.

        :param detec_tables_dir: Directory where the SimDetecTables are located.
        """
        print(f"\nLoading SimDetecTables from directory: {detec_tables_dir}")
        for sigma_kern in self.sigma_kerns:
            self.d[sigma_kern] = {}
            for peak_appmag in self.peak_appmags:
                self.d[sigma_kern][peak_appmag] = SimDetecTable(sigma_kern, peak_appmag)
                self.d[sigma_kern][peak_appmag].load_detec_table(
                    self.model_name, detec_tables_dir
                )
        print("Success")


class EfficiencyTable(pdastrostatsclass):
    def __init__(
        self,
        sigma_kerns: List,
        peak_appmags: List,
        params: Dict[str, List],
        **kwargs,
    ):
        """
        Initialize an EfficiencyTable.

        :param sigma_kerns: List of detection algorithm kernel sizes.
        :param peak_appmags: List of possible simulation peak apparent magnitudes.
        :param params: Dictionary of parameter names and lists of possible values.
        """
        pdastrostatsclass.__init__(self, **kwargs)

        self.sigma_kerns: List = sigma_kerns
        self.peak_appmags: List = peak_appmags
        self.peak_fluxes: List = list(map(mag2flux, peak_appmags))

        self.params: Dict[str, List] = params
        if "peak_appmag" in self.params.keys():
            del self.params["peak_appmag"]

    def setup(self, time_colname: str):
        """
        Set up the table columns for sigma_kerns, peak_appmags, and other known parameter values.
        The time column name, if known, will be skipped when constructing columns.

        :param time_colname: Any peak MJD, MJD0, or time-related parameter name that denotes where to inject the Simulation.
        """
        all_params = dict(
            {
                "sigma_kern": self.sigma_kerns,
                "peak_appmag": self.peak_appmags,
            },
            **self.params,
        )

        if time_colname in all_params.keys():
            del all_params[time_colname]

        print(f"\nSetting up efficiency table with columns: {all_params.keys()}")

        keys, values = zip(*(all_params).items())
        combinations = list(itertools.product(*values))
        self.t = pd.DataFrame(combinations, columns=keys)

        self.t["peak_flux"] = self.t["peak_appmag"].apply(lambda mag: mag2flux(mag))

        col_order = ["sigma_kern", "peak_appmag", "peak_flux"]
        other_cols = [col for col in self.t.columns if col not in col_order]
        col_order += other_cols
        self.t = self.t[col_order]

    # create dictionary of FOM limits, with sigma_kerns as the keys
    def set_fom_limits(self, fom_limits: List | Dict):
        """
        Create dictionary of FOM limits, with sigma_kerns as the keys.
        """
        if fom_limits is None:
            return None

        if isinstance(fom_limits, list):
            res = {}
            for i in range(len(self.sigma_kerns)):
                res[self.sigma_kerns[i]] = fom_limits[i]
        elif len(fom_limits) != len(self.sigma_kerns):
            raise RuntimeError(
                "ERROR: Each entry in sigma_kerns must have a matching list in fom_limits"
            )
        else:
            res = fom_limits
        return res

    def get_efficiencies(
        self,
        sd: SimDetecTables,
        fom_limits: List | Dict[float, List[float]],
        time_colname: str,
        **kwargs,
    ):
        """
        For each row in the efficiency table, compute efficiencies for the FOM limits corresponding to the given sigma_kern.

        :param sd: SimDetecTables object that contains simulation information, max FOM, and other data needed to run the detection algorithm. Each SimDetecTable corresponds to one sigma_kern x peak_appmag combination.
        :param fom_limits: Dictionary with sigma_kerns as keys and lists of FOM limits as values.
        :param time_colname: Any peak MJD, MJD0, or time-related parameter name that denotes where to inject the Simulation.
        """

        fom_limits = self.set_fom_limits(fom_limits)

        l = len(self.t)
        print("Calculating efficiencies...")
        print_progress_bar(0, l, prefix="Progress:", suffix="Complete", length=50)
        for i in range(l):
            sigma_kern = self.t.loc[i, "sigma_kern"]
            peak_appmag = self.t.loc[i, "peak_appmag"]
            sim_detec_table = sd.get_table(sigma_kern, peak_appmag)

            if kwargs:
                params = kwargs
            else:
                params = sim_detec_table.get_params_at_index(
                    i, time_colname=time_colname
                )

            for fom_limit in fom_limits[sigma_kern]:
                try:
                    efficiency = sd.get_efficiency(
                        sigma_kern, peak_appmag, fom_limit, **params
                    )
                except Exception as e:
                    raise RuntimeError(
                        f"ERROR: Could not calculate efficiency for sigma_kern={sigma_kern}, peak_appmag={peak_appmag:0.2f}, fom_limit={fom_limit:0.2f}: {str(e)}"
                    )
                self.t.loc[i, f"pct_detec_{fom_limit:0.2f}"] = efficiency

            print_progress_bar(
                i + 1, l, prefix="Progress:", suffix="Complete", length=50
            )

    def get_subset(self, fom_limits: List = None, **kwargs):
        """
        Get a subset of the table where the columns match the given values.

        :param fom_limits: List of FOM limits to get columns for. Set to None for all FOM limit columns.
        :param kwargs: Arbitrary number of pairs of column = value, column = range, or column = list of ranges.
        Example usage: self.get_efficiency(10.0, sigma_kern=2, sigma_sim=[5, 6], peak_mjd=[[58000, 58100], [58200, 58300]])

        :return: Efficiency of the rows that match the criteria.
        """
        ix = self.getindices()
        colnames: List = ["sigma_kern", "peak_appmag"] + self.params.keys()

        if not fom_limits is None:
            for col in self.t.columns:
                for fom_limit in fom_limits:
                    if re.search(f"^pct_detec_{fom_limit:0.2f}", col):
                        colnames.append(col)
        else:
            for col in self.t.columns:
                if re.search("^pct_detec_", col):
                    colnames.append(col)

        for column, value in kwargs.items():
            if isinstance(value, list):
                if all(isinstance(v, list) for v in value):
                    mask = self.t.loc[ix, column].apply(
                        lambda x: any(a <= x <= b for a, b in value)
                    )
                    ix = self.t.index[mask].tolist()
                else:
                    ix = self.ix_inrange(
                        colnames=column, lowlim=value[0], uplim=value[1], indices=ix
                    )
            else:
                ix = self.ix_equal(column, value, indices=ix)

        return self.t.loc[ix, colnames]

    def reset_table(self):
        """
        Remove any previously calculated efficiency columns.
        """
        for col in self.t.columns:
            if re.search("^pct_detec_", col):
                self.t.drop(col, axis=1, inplace=True)

        for i in range(len(self.sigma_kerns)):
            fom_limits = self.fom_limits[self.sigma_kerns[i]]
            for fom_limit in fom_limits:
                self.t[f"pct_detec_{fom_limit:0.2f}"] = np.full(len(self.t), np.nan)

    def merge_tables(self, other: Self):
        """
        Add table content, sigma_kerns, and fom_limits from another EfficiencyTable.
        """
        if not isinstance(other, EfficiencyTable):
            raise RuntimeError(
                f"ERROR: Cannot merge EfficiencyTable with object type: {type(other)}"
            )

        self.sigma_kerns += other.sigma_kerns
        if not self.fom_limits is None:
            self.fom_limits.update(other.fom_limits)

        self.t = pd.concat([self.t, other.t], ignore_index=True)

    def load(self, detec_tables_dir: str, model_name: str):
        filename = f"{detec_tables_dir}/efficiencies_{model_name}.txt"
        print(f"Loading efficiency table at {filename}...")
        try:
            self.load_spacesep(filename, delim_whitespace=True)
        except Exception as e:
            raise RuntimeError(
                f"ERROR: Could not load efficiency table at {filename}: {str(e)}"
            )

    def save(self, detec_tables_dir: str, model_name: str):
        filename = f"{detec_tables_dir}/efficiencies_{model_name}.txt"
        print(f"Saving efficiency table as {filename}...")
        self.write(filename=filename, overwrite=True, index=False)

    def __str__(self):
        return self.t.to_string()


# TODO: documentation
class SimDetecLoop(ABC):
    def __init__(self, sigma_kerns: List, **kwargs):
        self.sigma_kerns: List = sigma_kerns
        self.peak_appmags: List = None
        self.peak_fluxes: List = None

        self.sn: SimDetecSupernova = None
        self.e: EfficiencyTable = None
        self.sd: SimDetecTables = None

    @abstractmethod
    def set_peak_mags_and_fluxes(
        self,
        model_name: str = None,
        sim_tables_dir: str | None = None,
        detec_tables_dir: str | None = None,
        **kwargs,
    ):
        """
        Set the lists of peak apparent magnitudes and fluxes.
        Default behavior searches the given directory for SimTable or SimDetecTable filenames with matching model names.

        :param model_name: Name of the model for which to search for possible peak apparent magnitudes.
        :param sim_tables_dir: Directory where the SimTables are located. If given, search for SimTables.
        :param detec_tables_dir: Directory where the SimDetecTables are located. If given, search for SimDetecTables.
        """
        if model_name is None or (sim_tables_dir is None and detec_tables_dir is None):
            raise RuntimeError(
                "ERROR: Please either provide a model name and SimTables or SimDetecTables directory, or overwrite this function with your own."
            )
        if not sim_tables_dir is None and not detec_tables_dir is None:
            raise RuntimeError(
                "ERROR: Please provide either a SimTables directory or a SimDetecTables directory, not both."
            )

        self.peak_appmags = set()
        if detec_tables_dir is None:
            pattern = re.compile(
                f"sim_{re.escape(model_name)}_([0-9]*\.[0-9]{{2}})\.txt"
            )
        else:
            pattern = re.compile(
                f"simdetec_{re.escape(model_name)}_" + r"\d+_(\d+\.\d{2})\.txt$"
            )

        filenames = (
            os.listdir(sim_tables_dir)
            if detec_tables_dir is None
            else os.listdir(detec_tables_dir)
        )
        for filename in filenames:
            match = pattern.match(filename)
            if match:
                peak_appmag = float(match.group(1))
                self.peak_appmags.add(peak_appmag)

        self.peak_appmags = list(self.peak_appmags)
        self.peak_appmags.sort()

        self.peak_fluxes = list(map(mag2flux, self.peak_appmags))

    @abstractmethod
    def load_sn(
        self,
        data_dir: str,
        tnsname: str,
        num_controls: int,
        mjdbinsize: float = 1.0,
        filt: str = "o",
    ):
        """
        Load the averaged SN and its control light curves.

        :param data_dir: Directory where the SN folder is located.
        :param tnsname: TNS name of the SN to load.
        :param num_controls: Number of averaged control light curves to load.
        :param mjdbinsize: MJD bin size of the averaged light curves to load.
        :param filt: Filter of the averaged light curves to load.
        """
        self.sn = SimDetecSupernova(tnsname, mjdbinsize=mjdbinsize, filt=filt)
        self.sn.load_all(data_dir, num_controls=num_controls)
        self.sn.remove_rolling_sums()
        self.sn.remove_simulations()

    @abstractmethod
    def load_sd(
        self,
        model_name: str,
        sim_tables_dir: str | None = None,
        detec_tables_dir: str | None = None,
    ):
        """
        Either load existing SimTables and construct SimDetecTables out of them, or load existing SimDetecTables.

        :param model_name: Name of the model whose tables will be loaded.
        :param sim_tables_dir: Directory where the SimTables are located. If given, load SimTables.
        :param detec_tables_dir: Directory where the SimDetecTables are located. If given, load SimDetecTables.
        """
        self.sd = SimDetecTables(self.peak_appmags, model_name, self.sigma_kerns)
        if sim_tables_dir is None and detec_tables_dir is None:
            raise RuntimeError(
                "ERROR: Please either provide a SimTables or SimDetecTables directory."
            )
        if not sim_tables_dir is None and not detec_tables_dir is None:
            raise RuntimeError(
                "ERROR: Please provide either a SimTables directory or a SimDetecTables directory, not both."
            )

        if detec_tables_dir is None:
            self.sd.load_all_from_sim_tables(sim_tables_dir)
        else:
            self.sd.load_all(detec_tables_dir)

    @abstractmethod
    def load_sim(
        self,
        table_row: Dict,
    ) -> Simulation:
        """
        Construct and return a Simulation object given a row from a SimTable or SimDetecTable.
        """
        model_name = table_row["model_name"]
        filename = (
            None
            if not isinstance(table_row["filename"], str)
            else table_row["filename"]
        )

        def get_col_val(colname, table_row):
            if colname in table_row:
                return None if np.isnan(table_row[colname]) else table_row[colname]
            else:
                return False

        mjd_colname = get_col_val("mjd_colname", table_row)
        mag_colname = get_col_val("mag_colname", table_row)
        flux_colname = get_col_val("flux_colname", table_row)

        if args.model_name == GAUSSIAN_MODEL_NAME:
            print("Using Gaussian simulations")
            sim = Gaussian()
        elif args.model_name == ASYMMETRIC_GAUSSIAN_MODEL_NAME:
            print("Using asymmetric Gaussian simulations")
            sim = AsymmetricGaussian()
        else:
            print(
                f'Using "{model_name}" simulations with MJD column {mjd_colname}, mag column {mag_colname}, and flux column {flux_colname} at filename: {filename}'
            )
            sim = Model(
                filename=filename,
                mjd_colname=mjd_colname,
                mag_colname=mag_colname,
                flux_colname=flux_colname,
                model_name=model_name,
            )
        return sim

    @abstractmethod
    def get_max_fom_indices(
        self,
        sim_lc: SimDetecLightCurve,
        **kwargs,
    ):
        """
        From a light curve with a Simulation injected, return the indices within which to search for the max FOM.

        :param sim_lc: SimDetecLightCurve with a Simulation injected.
        """
        pass

    @abstractmethod
    def update_sd_row(
        self,
        sigma_kern: float,
        peak_appmag: float,
        index: int,
        control_index: int,
        max_fom: float,
        max_fom_mjd: float,
    ):
        """
        Update a certain row of a SimDetecTable with info about an injected Simulation, i.e., the index of the control light curve it was added to, the max FOM, and the max FOM MJD.

        :param sigma_kern: Sigma of the desired SimDetecTable.
        :param peak_appmag: Peak apparent magnitude of the desired SimDetecTable.
        :param index: Index of the row to update.
        :param control_index: Control light curve index to which the Simulation was added.
        :param max_fom: Max FOM of the simulated flux.
        :param max_fom_mjd: MJD of the max FOM of the simulated flux.
        """
        data = {
            "control_index": control_index,
            "max_fom": max_fom,
            "max_fom_mjd": max_fom_mjd,
        }
        self.sd.update_row_at_index(sigma_kern, peak_appmag, index, data)

    @abstractmethod
    def calculate_efficiencies(
        self,
        fom_limits: List | Dict[float, List[float]],
        params: Dict[str, List],
        detec_tables_dir: str,
        model_name: str,
        time_param_name: str,
        **kwargs,
    ):
        """
        Construct and save an EfficiencyTable that contains efficiencies for every combination of a Simulation's sigma_kern, peak_appmag, and other parameters EXCEPT the time parameter.

        :param fom_limits: Dictionary of sigma_kerns and their corresponding FOM limits.
        :param params: Dictionary of parameter names and possible values.
        :param detec_tables_dir: Directory where the EfficiencyTable should be saved.
        :param model_name: Name of the model for which to calculate efficiencies.
        :param time_param_name: Any peak MJD, MJD0, or time-related parameter name that denotes where to inject the Simulation.
        """
        self.e = EfficiencyTable(self.sigma_kerns, self.peak_appmags, params)
        self.e.setup(time_param_name)
        self.e.get_efficiencies(self.sd, fom_limits, time_param_name)
        self.e.save(detec_tables_dir, model_name)

    @abstractmethod
    def loop(
        self,
        valid_control_ix: List,
        detec_tables_dir: str,
        flag=0x800000,
        **kwargs,
    ):
        """
        Loop over each possible sigma_kern, then each possible peak_appmag, then each row in that corresponding SimDetecTable.
        For each row, inject a Simulation with the specified parameters into a random control light curve.
        Update the row with information about where it was injected, what/where its max FOM is, etc.

        :param valid_control_ix: List of indices of control light curves which may be randomly selected to have a Simulation injected.
        :param detec_tables_dir: Directory where the SimDetecTables should be saved.
        :param flag: Flag that denotes bad days in the averaged light curves.
        """
        pass


class AtlasSimDetecLoop(SimDetecLoop):
    def __init__(self, sigma_kerns: List, **kwargs):
        super().__init__(sigma_kerns, **kwargs)

    def set_peak_mags_and_fluxes(self, model_name=None, sim_tables_dir=None, **kwargs):
        return super().set_peak_mags_and_fluxes(
            model_name=model_name, sim_tables_dir=sim_tables_dir, **kwargs
        )

    def load_sn(self, data_dir, tnsname, num_controls, mjdbinsize=1, filt="o"):
        return super().load_sn(data_dir, tnsname, num_controls, mjdbinsize, filt)

    def load_sd(self, model_name, sim_tables_dir=None, detec_tables_dir=None):
        return super().load_sd(
            model_name, sim_tables_dir=sim_tables_dir, detec_tables_dir=detec_tables_dir
        )

    def load_sim(self, data: Dict) -> Simulation:
        return super().load_sim(data)

    def get_max_fom_indices(
        self, sim_lc: SimDetecLightCurve, peak_mjd=None, sigma_sim=None, **kwargs
    ):
        """
        Get indices of MJD within 1 sigma of the peak MJD.
        For Gaussians, use the sigma provided.
        For Charlie's model, use the manually calculated value of 2.8.
        """
        if peak_mjd is None:
            raise RuntimeError("ERROR: A peak MJD is required to find the max FOM.")
        if sigma_sim is None:
            # replace with default sigma sim for Charlie's model
            sigma_sim = 2.8

        # measurements within 1 sigma of the peak MJD
        indices = sim_lc.ix_inrange(
            colnames="MJDbin", lowlim=peak_mjd - sigma_sim, uplim=peak_mjd + sigma_sim
        )
        return indices

    def update_sd_row(
        self, sigma_kern, peak_appmag, index, control_index, max_fom, max_fom_mjd
    ):
        return super().update_sd_row(
            sigma_kern, peak_appmag, index, control_index, max_fom, max_fom_mjd
        )

    def calculate_efficiencies(
        self,
        fom_limits: List | Dict[float, List[float]],
        params: Dict[str, List],
        detec_tables_dir: str,
        model_name: str,
        time_param_name: str,
        **kwargs,
    ):
        return super().calculate_efficiencies(
            fom_limits, params, detec_tables_dir, model_name, time_param_name, **kwargs
        )

    def loop(
        self,
        valid_control_ix: List,
        detec_tables_dir: str,
        flag=0x800000,
        **kwargs,
    ):
        # loop through each rolling sum kernel size
        for sigma_kern in self.sigma_kerns:
            print(f"\n\tUsing rolling sum kernel size sigma_kern={sigma_kern} days...")
            self.sn.apply_rolling_sums(sigma_kern, flag=flag)

            # loop through each possible peak apparent magnitude
            for peak_appmag in self.peak_appmags:
                sim_detec_table = self.sd.get_table(sigma_kern, peak_appmag)
                print(
                    f"\nCommencing {len(sim_detec_table.t)} simulations for peak app mag {peak_appmag} (peak flux {mag2flux(peak_appmag):0.2f} uJy)..."
                )

                # load the Simulation object based on the data in the first row
                # we assume here that every row adds the same type of model
                sim = self.load_sim(dict(sim_detec_table.t.loc[0, :]))

                for i in range(len(sim_detec_table.t)):
                    # pick random control light curve
                    rand_control_index = random.choice(valid_control_ix)

                    # add the simulated flux to the chosen control light curve
                    params = sim_detec_table.get_params_at_index(i)
                    sim_lc = self.sn.avg_lcs[rand_control_index].add_simulation(
                        sim, peak_appmag, flag=flag, remove_old=True, **params
                    )

                    # get the max simulated FOM within certain indices of the light curve
                    indices = self.get_max_fom_indices(sim_lc, **params)
                    max_fom_mjd, max_fom = sim_lc.get_max_fom(indices=indices)

                    # update the corresponding row in the SimDetecTable
                    self.update_sd_row(
                        sigma_kern,
                        peak_appmag,
                        i,
                        rand_control_index,
                        max_fom,
                        max_fom_mjd,
                    )

                self.sd.save_detec_table(sigma_kern, peak_appmag, detec_tables_dir)
                print("Success")

        print("\nFinished generating all SimDetecTables")


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)
    parser.add_argument("model_name", type=str, help="name of model to use")
    parser.add_argument(
        "-d",
        "--detec_config_file",
        default="detection_settings.json",
        type=str,
        help="file name of JSON file with SimDetecTable generation and efficiency calculation settings",
    )
    parser.add_argument(
        "-f",
        "--sim_config_file",
        default="simulation_settings.json",
        type=str,
        help="file name of JSON file with model information and SimTable generation settings",
    )
    parser.add_argument(
        "--skip_generate",
        default=False,
        action="store_true",
        help="skip SimDetecTable generation and directly load them instead; use -e argument to calculate efficiencies on loaded SimDetecTables",
    )
    parser.add_argument(
        "-e",
        "--efficiencies",
        default=False,
        action="store_true",
        help="calculate efficiencies using FOM limits",
    )

    return parser


if __name__ == "__main__":
    args = define_args().parse_args()
    detec_config = load_json_config(args.detec_config_file)
    sim_config = load_json_config(args.sim_config_file)

    if " " in args.model_name:
        raise RuntimeError("ERROR: Model name cannot have spaces.")
    try:
        model_settings = sim_config[args.model_name]
    except Exception as e:
        raise RuntimeError(
            f"ERROR: Could not find model {args.model_name} in model config file: {str(e)}"
        )

    sn_info = detec_config["sn_info"]
    data_dir = detec_config["data_dir"]
    sim_tables_dir = detec_config["sim_tables_dir"]
    detec_tables_dir = detec_config["detec_tables_dir"]
    sigma_kerns = [obj["sigma_kern"] for obj in detec_config["sigma_kerns"]]
    valid_control_ix = [
        i
        for i in range(1, sn_info["num_controls"] + 1)
        if not i in detec_config["skip_control_ix"]
    ]

    simdetec = AtlasSimDetecLoop(sigma_kerns)
    simdetec.load_sn(
        data_dir,
        sn_info["tnsname"],
        sn_info["num_controls"],
        sn_info["mjd_bin_size"],
        sn_info["filt"],
    )

    if args.skip_generate:
        simdetec.set_peak_mags_and_fluxes(
            model_name=args.model_name, detec_tables_dir=detec_tables_dir
        )
        simdetec.load_sd(args.model_name, detec_tables_dir=detec_tables_dir)
    else:
        simdetec.set_peak_mags_and_fluxes(
            model_name=args.model_name, sim_tables_dir=sim_tables_dir
        )
        simdetec.load_sd(args.model_name, sim_tables_dir=sim_tables_dir)

        simdetec.loop(valid_control_ix, detec_tables_dir, flag=sn_info["badday_flag"])

    if args.efficiencies:
        parsed_params = parse_params(model_settings)
        fom_limits = [obj["fom_limits"] for obj in detec_config["sigma_kerns"]]
        simdetec.calculate_efficiencies(
            fom_limits,
            parsed_params,
            detec_tables_dir,
            args.model_name,
            model_settings["time_parameter_name"],
        )
