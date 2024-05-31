#!/usr/bin/env python

from abc import ABC, abstractmethod
import argparse
from copy import deepcopy
from typing import List, Tuple
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from astropy.modeling.functional_models import Gaussian1D

from generate_sim_table import load_json_config
from lightcurve import SimDetecLightCurve, SimDetecSupernova


# convert flux to magnitude
def flux2mag(flux: float):
    return -2.5 * np.log10(flux) + 23.9


# convert magnitude to flux
def mag2flux(mag: float):
    return 10 ** ((mag - 23.9) / -2.5)


class Simulation(ABC):
    def __init__(self, model_name=None, **kwargs):
        """
        Initialize the Simulation object.
        """
        self.model_name = model_name

    @abstractmethod
    def get_sim_flux(self, mjds, **kwargs):
        """
        Compute the simulated flux for given MJDs.

        :param mjds: List or array of MJDs.

        :return: An array of flux values corresponding to the input MJDs.
        """
        pass

    def __str__(self):
        return f'Simulation with model name "{self.model_name}"'


class Gaussian(Simulation):
    def __init__(self, sigma, peak_appmag=None, model_name="gaussian", **kwargs):
        Simulation.__init__(self, model_name=model_name, **kwargs)
        self.sigma = sigma
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
        print(f"Loading model light curve at {filename}...")

        try:
            self.t = pd.read_table(filename, delim_whitespace=True, header=None)
        except Exception as e:
            raise RuntimeError(f"ERROR: Could not load model at {filename}: {str(e)}")

        if mjd_colname is False:
            # create MJD column
            columns = ["MJD"] + self.t.columns
            self.t["MJD"] = range(len(self.t))
            self.t = self.t[columns]
        else:
            # rename column to "MJD"
            self.t.rename(
                columns={0 if mjd_colname is None else mjd_colname: "MJD"}, inplace=True
            )

        if mag_colname is False:
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

    # get interpolated function of model at peak MJD (peak_mjd)
    # and match to time array (mjds)
    def get_sim_flux(self, mjds, peak_mjd=None, peak_appmag=None, **kwargs):
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


# TODO: documentation
class SimDetecLoop(ABC):
    def __init__(self, output_dir: str, tables_dir: str, sigma_kerns: List, **kwargs):
        self.output_dir = output_dir
        self.tables_dir = tables_dir
        self.sigma_kerns = sigma_kerns

        self.sn: SimDetecSupernova = None
        self.e: EfficiencyTable = None
        self.sd: SimDetecTables = None

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
