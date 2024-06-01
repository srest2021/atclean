#!/usr/bin/env python

from abc import ABC, abstractmethod
import argparse, re
from copy import deepcopy
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from astropy.modeling.functional_models import Gaussian1D

from download import make_dir_if_not_exists
from pdastro import pdastrostatsclass
from generate_sim_table import SimTable, SimTables, load_json_config
from lightcurve import SimDetecLightCurve, SimDetecSupernova, Simulation


# convert flux to magnitude
def flux2mag(flux: float):
    return -2.5 * np.log10(flux) + 23.9


# convert magnitude to flux
def mag2flux(mag: float):
    return 10 ** ((mag - 23.9) / -2.5)


class Gaussian(Simulation):
    def __init__(self, sigma, peak_appmag=None, model_name="gaussian", **kwargs):
        """
        Initialize a Gaussian simulation object.

        :param sigma: Sigma or kernel size of the Gaussian.
        :param peak_appmag: Peak apparent magnitude of the Gaussian.
        :param model_name: Name of the Gaussian model in the config file.
        """
        Simulation.__init__(self, model_name=model_name, **kwargs)
        self.sigma = sigma
        if peak_appmag is None:
            raise RuntimeError("ERROR: Gaussian must have a peak_appmag.")
        self.peak_appmag = peak_appmag
        self.g = self.new_gaussian(mag2flux(peak_appmag), sigma)

    def new_gaussian(self, peak_flux, sigma):
        x = np.arange(-100, 100, 0.01)
        g1 = Gaussian1D(amplitude=peak_flux, stddev=sigma)(x)
        g2 = Gaussian1D(amplitude=peak_flux, stddev=sigma)(x)

        ind = np.argmin(abs(x))
        g3 = np.copy(g1)
        g3[ind:] = g2[ind:]
        gauss = np.array([x, g3])
        return gauss

    # get interpolated function of gaussian at peak MJD (peak_mjd)
    # and match to time array (mjds)
    def get_sim_flux(self, mjds, peak_mjd=None, **kwargs):
        """
        Get the interpolated function of the Gaussian at a given peak MJD and match it to the given time array.

        :param mjds: Time array of MJDs.
        :param peak_mjd: MJD at which the Gaussian should reach its peak apparent magnitude.

        :return: The simulated flux corresponding to the given time array.
        """
        if peak_mjd is None:
            raise RuntimeError(
                "ERROR: Peak MJD required to get flux of simulated Gaussian."
            )

        g = deepcopy(self.g)
        g[0, :] += peak_mjd

        fn = interp1d(g[0], g[1], bounds_error=False, fill_value=0)
        sim_flux = fn(mjds)
        return sim_flux

    def __str__(self):
        return (
            super().__str__()
            + f": peak appmag = {self.peak_appmag:0.2f}, sigma = {self.sigma}"
        )


class Model(Simulation):
    def __init__(
        self,
        filename,
        mjd_colname=False,
        mag_colname=False,
        flux_colname=False,
        model_name="pre_SN_outburst",
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
        self.peak_appmag = None
        self.t = None

        self.load(
            filename,
            mjd_colname=mjd_colname,
            mag_colname=mag_colname,
            flux_colname=flux_colname,
        )

    def load(
        self,
        filename,
        mjd_colname=False,
        mag_colname=False,
        flux_colname=False,
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

        print(self.t[["MJD", "m", "uJy"]].head().to_string())
        print("Success")

    # get interpolated function of model at peak MJD (peak_mjd)
    # and match to time array (mjds)
    def get_sim_flux(self, mjds, peak_mjd=None, peak_appmag=None, **kwargs):
        """
        Get the interpolated function of the model at a given peak MJD and peak apparent magnitude and match it to the given time array.

        :param mjds: Time array of MJDs.
        :param peak_mjd: MJD at which the model should reach its peak apparent magnitude.
        :param peak_appmag: Desired peak apparent magnitude of the model.

        :return: The simulated flux corresponding to the given time array.
        """
        if peak_mjd is None:
            raise RuntimeError("ERROR: Peak MJD required to construct simulated model.")

        # get original peak appmag index
        peak_idx = self.t["m"].idxmin()

        if not peak_appmag is None:
            self.peak_appmag = peak_appmag

            # scale flux to the desired peak appmag
            self.t["uJy"] *= mag2flux(peak_appmag) / self.t.loc[peak_idx, "uJy"]

            # recalulate appmag column
            self.t["m"] = self.t["uJy"].apply(lambda flux: flux2mag(flux))
        else:
            self.peak_appmag = self.t.loc[peak_idx, "m"]

        # put peak appmag at peak_mjd
        self.t["MJD"] -= self.t.loc[peak_idx, "MJD"]
        self.t["MJD"] += peak_mjd

        # interpolate lc and match to time array
        fn = interp1d(self.t["MJD"], self.t["uJy"], bounds_error=False, fill_value=0)
        sim_flux = fn(mjds)
        return sim_flux

    def __str__(self):
        return super().__str__() + f": peak appmag = {self.peak_appmag:0.2f}"


class SimDetecTable(SimTable):
    def __init__(self, sigma_kern, peak_appmag, **kwargs):
        SimTable.__init__(peak_appmag, **kwargs)
        self.sigma_kern = sigma_kern

    def update_row_at_index(self, index, data: Dict):
        """
        Update a certain row of the table.

        :param index: Index of the table to update.
        :param data: Dictionary of column-value pairs.
        """
        # self.t.loc[index, data.keys()] = np.array(list(data.values()))
        for key, value in data.items():
            self.t.at[index, key] = value

    def get_filename(self, model_name, tables_dir):
        return f"{tables_dir}/simdetec_{model_name}_{self.sigma_kern}_{self.peak_appmag:0.2f}.txt"

    def load_from_sim_table(self, model_name, sim_tables_dir):
        """
        Load an existing SimTable and turn it into a SimDetecTable.

        :param model_name: Name of the model assigned in the config file.
        :param sim_tables_dir: Directory where the SimTable is located.
        """
        super().load(model_name, sim_tables_dir)
        self.t["sigma_kern"] = self.sigma_kern

    def get_efficiency(self, fom_limit, **kwargs):
        """
        Get the efficiency where columns match all the given values and are within all the given ranges.

        :param kwargs: Arbitrary number of pairs of column = value, column = range, or column = list of ranges.
        Example usage for columns A, B, C: self.get_efficiency(10.0, A=2, B=[5, 6], C=[[1, 2], [3, 4]])

        :return: Efficiency of the rows that match the criteria.
        """
        col_ix = self.getindices()
        for column, value in kwargs.items():
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
        SimTables.__init__(peak_appmags, model_name)
        self.sigma_kerns = sigma_kerns
        self.d: Dict[float, Dict[float, SimDetecTable]] = {}

    def get_efficiency(self, sigma_kern, peak_appmag, fom_limit, **kwargs):
        return self.d[sigma_kern][peak_appmag].get_efficiency(fom_limit, **kwargs)

    def save_all(self, tables_dir):
        print(f"\nSaving SimDetecTables in directory: {tables_dir}")
        make_dir_if_not_exists(tables_dir)
        for sigma_kern in self.d.keys():
            for table in self.d[sigma_kern].values():
                table.save(self.model_name, tables_dir)
        print("Success")

    def load_all_from_sim_tables(self, sim_tables_dir):
        """
        Load existing SimTables and turn them into SimDetecTables.

        :param sim_tables_dir: Directory where the SimTable is located.
        """
        print(
            f"\nConstructing SimDetecTables from existing SimTables in directory: {sim_tables_dir}"
        )
        for sigma_kern in self.sigma_kerns:
            for peak_appmag in self.peak_appmags:
                self.d[sigma_kern][peak_appmag] = SimDetecTable(sigma_kern, peak_appmag)
                self.d[sigma_kern][peak_appmag].load_from_sim_table(
                    self.model_name, sim_tables_dir
                )

    def load_all(self, tables_dir):
        """
        Load existing SimDetecTables.

        :param sim_tables_dir: Directory where the SimDetecTables are located.
        """
        print(f"\nLoading SimDetecTables from directory: {tables_dir}")
        for sigma_kern in self.sigma_kerns:
            for peak_appmag in self.peak_appmags:
                self.d[sigma_kern][peak_appmag] = SimDetecTable(sigma_kern, peak_appmag)
                self.d[sigma_kern][peak_appmag].load(self.model_name, tables_dir)


class EfficiencyTable(pdastrostatsclass):
    def __init__(
        self,
        sigma_kerns,
        peak_appmags,
        params,
        fom_limits=None,
        **kwargs,
    ):
        pdastrostatsclass.__init__(self, **kwargs)

        self.sigma_kerns = sigma_kerns
        self.fom_limits = self.get_fom_limits(fom_limits)

        self.peak_appmags = peak_appmags
        self.peak_fluxes = list(map(mag2flux, peak_appmags))

        self.setup(params)

    def setup(self, params):
        # TODO
        self.params: Dict[str, List] = params

    # create dictionary of FOM limits, with sigma_kerns as the keys
    def get_fom_limits(self, fom_limits):
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

    def get_efficiencies(self):
        # TODO
        pass

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

    # remove previously calculated efficiency columns
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

    def merge_tables(self, other):
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

    def load(self, tables_dir, filename="efficiencies.txt"):
        filename = f"{tables_dir}/{filename}"
        print(f"Loading efficiency table at {filename}...")
        try:
            self.load_spacesep(filename, delim_whitespace=True)
        except Exception as e:
            raise RuntimeError(
                f"ERROR: Could not load efficiency table at {filename}: {str(e)}"
            )

    def save(self, tables_dir, filename="efficiencies.txt"):
        filename = f"{tables_dir}/{filename}"
        print(f"Saving efficiency table as {filename}...")
        self.write(filename=filename, overwrite=True, index=False)

    def __str__(self):
        return self.t.to_string()


# TODO: documentation
class SimDetecLoop(ABC):
    def __init__(self, output_dir: str, tables_dir: str, sigma_kerns: List, **kwargs):
        self.output_dir = output_dir
        self.tables_dir = tables_dir
        self.sigma_kerns = sigma_kerns

        self.sn: SimDetecSupernova = None
        # self.e: EfficiencyTable = None
        # self.sd: SimDetecTables = None

    @abstractmethod
    def set_peak_mags_and_fluxes(
        self,
        peak_mag_min: float = 23.0,
        peak_mag_max: float = 16.0,
        n_peaks: int = 20,
        **kwargs,
    ):
        self.peak_appmags = list(np.linspace(peak_mag_min, peak_mag_max, num=n_peaks))
        self.peak_fluxes = list(map(mag2flux, self.peak_appmags))

    @abstractmethod
    def load_sn(self, data_dir, tnsname, num_controls, mjdbinsize=1.0, filt="o"):
        self.sn.load_all(data_dir, num_controls=num_controls)
        self.sn.remove_rolling_sums()
        self.sn.remove_simulations()

    @abstractmethod
    def modify_light_curve(self, control_index=0):
        pass

    @abstractmethod
    def get_max_fom(self, sim_lc: SimDetecLightCurve, **kwargs) -> Tuple[float, float]:
        pass

    @abstractmethod
    def add_row_to_sd(
        self,
        sigma_kern,
        peak_appmag,
        index,
        sim: Simulation,
        control_index,
        max_fom,
        max_fom_mjd,
    ):
        pass
        # sim_params = sim.get_row_params()
        # other_params = {
        # 	'control_index': control_index,
        # 	'max_fom': max_fom,
        # 	'max_fom_mjd': max_fom_mjd
        # }
        # row = other_params.update(sim_params)
        # self.sd.update_row_at_index(sigma_kern, peak_appmag, index, row)

    @abstractmethod
    def calculate_efficiencies(self):
        pass

    @abstractmethod
    def loop(self):
        pass


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

    parser.add_argument(
        "-f",
        "--config_file",
        default="detection_settings.json",
        type=str,
        help="file name of JSON file with general settings",
    )
    parser.add_argument(
        "-f",
        "--sim_config_file",
        default="simulation_settings.json",
        type=str,
        help="file name of JSON file with model settings",
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
    config = load_json_config(args.config_file)
    sim_config = load_json_config(args.sim_config_file)

    sn_info = config["sn_info"]

    data_dir = config["data_dir"]
    sim_tables_dir = config["sim_tables_dir"]
    detec_tables_dir = config["detec_tables_dir"]

    sigma_kerns = [obj["sigma_kern"] for obj in config["sigma_kerns"]]
    skip_control_ix = config["skip_control_ix"]
    fom_limits = None
    if args.efficiencies:
        fom_limits = [obj["fom_limits"] for obj in config["sigma_kerns"]]

    # model = Model(
    #     filename=sim_config["charlie_model"]["filename"],
    #     mjd_colname=sim_config["charlie_model"]["mjd_column_name"],
    #     mag_colname=sim_config["charlie_model"]["mag_column_name"],
    #     flux_colname=sim_config["charlie_model"]["flux_column_name"],
    # )
