#!/usr/bin/env python

import argparse
import os
import sys
from typing import Dict, List, Optional

from clean import parse_config_cuts
from download import (
    Credentials,
)
from lightcurve import (
    AveragedSupernova,
    Supernova,
)
from utils import (
    BadDayCut,
    ChiSquareCut,
    ControlLightCurveCut,
    CustomLogger,
    CutList,
    PresetColumnNames,
    SnInfoTable,
    UncertaintyCut,
    find_all_filts,
    get_allowed_presets,
    get_mjd0_from_tns,
    load_config,
    load_preset_column_names_from_config,
    make_dir_if_not_exists,
)
from plot import PlotLimits, PlotPdf


class PlotLoop:
    def __init__(
        self,
        colnames: PresetColumnNames,
        input_dir: str,
        output_dir: str,
        credentials: Credentials,
        sninfo_filename: Optional[str] = None,
        overwrite: bool = False,
    ):
        self.logger = CustomLogger()

        self.colnames = colnames
        self.sn: Optional[Supernova] = None
        self.avg_sn: Optional[AveragedSupernova] = None
        self.cut_list: Optional[CutList] = None
        self.p: Optional[PlotPdf] = None

        self.credentials: Credentials = credentials
        self.input_dir: str = input_dir
        self.output_dir: str = output_dir
        self.overwrite: bool = overwrite

        self.sninfo: SnInfoTable = SnInfoTable(
            self.output_dir, filename=sninfo_filename
        )

    def load_sn(
        self,
        tnsname: str,
        mjd0: float,
        filt: str,
        num_controls: int = 0,
        cleaned: bool = False,
    ):
        self.sn = Supernova(self.colnames, tnsname=tnsname, mjd0=mjd0, filt=filt)
        try:
            self.sn.load_all(
                self.output_dir, num_controls=num_controls, cleaned=cleaned
            )
        except Exception as e:
            raise RuntimeError(f"Could not load light curves: {str(e)}")

    def load_avg_sn(self, tnsname: str, mjd0: float, filt: str, mjdbinsize: float):
        self.avg_sn = AveragedSupernova(
            self.colnames,
            tnsname=tnsname,
            mjd0=mjd0,
            filt=filt,
            mjdbinsize=mjdbinsize,
        )
        try:
            self.avg_sn.load_all(output_dir, num_controls=num_controls)
        except Exception as e:
            raise RuntimeError(f"Could not load light curves: {str(e)}")

    def plot_lcs(
        self,
        tnsname: str,
        mjd0: float,
        filt: str,
        num_controls: int = 0,
        plot_uncert_est: bool = False,
        custom_lims: PlotLimits | None = None,
    ):
        self.logger.subheader(f"Plotting filter: {filt}")

        # load the cleaned SN and control light curves
        self.load_sn(tnsname, mjd0, filt, num_controls=num_controls, cleaned=True)

        if self.cut_list.has(BadDayCut.name()):
            # load averaged light curves
            self.load_avg_sn(
                tnsname,
                mjd0,
                filt,
                self.cut_list.get(BadDayCut.name()).mjd_bin_size,
            )

        # initialize PDF of diagnostic plots
        self.p = PlotPdf(f"{self.output_dir}/{tnsname}", tnsname, filt=filt)

        # plot original SN light curve and control light curves
        self.p.plot_SN(
            self.sn,
            custom_lims=custom_lims,
            plot_controls=True,
            plot_template_changes=True,
        )

        self.p.plot_all_controls(self.sn, custom_lims=custom_lims, include_sn=True)

        uncert_cut = self.cut_list.get(UncertaintyCut.name())
        if not uncert_cut is None:
            # plot uncertainty cut
            self.p.plot_cut(
                self.sn,
                uncert_cut.flag,
                custom_lims=custom_lims,
                title=UncertaintyCut.name(),
            )

        if plot_uncert_est:
            # plot true uncertainties estimation
            self.p.plot_uncert_est(self.sn)

        x2_cut = self.cut_list.get(ChiSquareCut.name())
        if not x2_cut is None:
            # plot chi-square cut
            self.p.plot_cut(
                self.sn, x2_cut.flag, custom_lims=custom_lims, title=ChiSquareCut.name()
            )

        controls_cut = self.cut_list.get(ControlLightCurveCut.name())
        if not controls_cut is None:
            # plot control light curve cut
            self.p.plot_cut(
                self.sn,
                controls_cut.flag,
                custom_lims=custom_lims,
                title=ControlLightCurveCut.name(),
            )

        custom_cuts = self.cut_list.get_custom_cuts()
        for cut in custom_cuts.values():
            # plot custom cut
            self.p.plot_cut(
                self.sn, cut.flag, custom_lims=custom_lims, title=cut.name()
            )

        # plot cleaned light curve using all previous cuts
        previous_flags = self.cut_list.get_previous_flags(BadDayCut.name())
        self.p.plot_cut(
            self.sn, previous_flags, custom_lims=custom_lims, title="All previous cuts"
        )
        self.p.plot_cleaned_SN(
            self.sn,
            previous_flags,
            custom_lims=custom_lims,
            plot_controls=True,
            plot_flagged=False,
        )

        badday_cut = self.cut_list.get(BadDayCut.name())
        if not badday_cut is None:
            # plot bad day cut
            self.p.plot_cut(
                self.avg_sn,
                badday_cut.flag,
                custom_lims=custom_lims,
                title=BadDayCut.name(),
            )

            # plot averaged light curves
            self.p.plot_averaged_SN(
                self.avg_sn,
                badday_cut.flag,
                custom_lims=custom_lims,
                plot_controls=True,
                plot_flagged=False,
            )

        # save the plots
        if not self.overwrite and os.path.exists(self.p.filename):
            self.logger.warning(
                f"Overwrite set to {self.overwrite} and file already exists at {self.p.filename}; skipping saving",
                dots=True,
            )
        else:
            self.p.save_pdf()

    def loop(
        self,
        tnsnames: List[str],
        cut_list: CutList,
        num_controls: int = 0,
        mjd0=None,
        filters: Optional[List[str]] = None,
        plot_uncert_est: bool = False,
        lims: Optional[PlotLimits] = None,
    ):
        self.cut_list = cut_list

        for obj_index in range(len(tnsnames)):
            tnsname = tnsnames[obj_index]
            self.logger.header(f"Plotting light curves for {tnsname}")

            make_dir_if_not_exists(f"{output_dir}/{tnsname}")

            if filters is None:
                filters = find_all_filts(self.output_dir, tnsname)

            if mjd0 is None:
                mjd0, coords = get_mjd0_from_tns(tnsname, self.sninfo, self.credentials)
                if not coords is None:
                    self.logger.body(f"Setting MJD0 to {mjd0}", newline=True)
                    self.sninfo.update_row(tnsname, coords=coords, mjd0=mjd0)
            else:
                self.logger.body(f"Setting MJD0 to {mjd0}", newline=True)

            for filt in filters:
                self.plot_lcs(
                    tnsname,
                    mjd0,
                    filt,
                    num_controls=num_controls,
                    plot_uncert_est=plot_uncert_est,
                    custom_lims=lims,
                )


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

    parser.add_argument(
        "tnsnames", nargs="+", help="TNS names of the transients to clean"
    )
    parser.add_argument(
        "-p",
        "--preset",
        type=str,
        default="atlas",
        help="preset name from config file (ex. atlas, rubin, tess)",
    )
    parser.add_argument(
        "--config_file",
        default="config.ini",
        type=str,
        help="file name of .ini file with settings for this class",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        default=False,
        action="store_true",
        help="overwrite existing file with same file name",
    )
    parser.add_argument(
        "--filters",
        type=str,
        default=None,
        help="comma-separated list of filters to plot",
    )
    parser.add_argument(
        "--num_controls",
        type=int,
        default=None,
        help="number of control light curves to load and plot",
    )
    parser.add_argument(
        "--mjd0", type=float, default=None, help="transient start date in MJD"
    )

    # possible cuts
    parser.add_argument(
        "-e",
        "--uncert_est",
        default=False,
        action="store_true",
        help="plot true uncertainty estimation",
    )
    parser.add_argument(
        "-u",
        "--uncert_cut",
        default=False,
        action="store_true",
        help="plot uncertainty cut",
    )
    parser.add_argument(
        "-x",
        "--x2_cut",
        default=False,
        action="store_true",
        help="plot chi-square cut",
    )
    parser.add_argument(
        "-c",
        "--controls_cut",
        default=False,
        action="store_true",
        help="plot control light curve cut",
    )
    parser.add_argument(
        "-g",
        "--averaging",
        default=False,
        action="store_true",
        help="plot averaged light curves and bad day cut",
    )
    parser.add_argument(
        "-m",
        "--mjd_bin_size",
        type=float,
        default=None,
        help="MJD bin size in days for averaging",
    )
    parser.add_argument(
        "--custom_cuts",
        default=False,
        action="store_true",
        help="scan config file for custom cuts and plot them",
    )

    # x and y limits
    parser.add_argument(
        "--xlim_lower",
        default=None,
        type=float,
        help="Lower x limit",
    )
    parser.add_argument(
        "--xlim_upper",
        default=None,
        type=float,
        help="Upper x limit",
    )
    parser.add_argument(
        "--ylim_lower",
        default=None,
        type=float,
        help="Lower y limit",
    )
    parser.add_argument(
        "--ylim_upper",
        default=None,
        type=float,
        help="Upper y limit",
    )

    return parser


if __name__ == "__main__":
    logger = CustomLogger()

    args = define_args().parse_args()
    config = load_config(args.config_file)

    if len(args.tnsnames) < 1:
        raise RuntimeError("Please specify at least one TNS name to plot.")
    if len(args.tnsnames) > 1 and not args.mjd0 is None:
        raise RuntimeError(f"Cannot specify one MJD0 {args.mjd0} for a batch of SNe.")
    logger.info(f"List of transients to plot: {args.tnsnames}", newline=True)

    colnames = load_preset_column_names_from_config(args.preset, config)

    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]
    sninfo_filename = config["dir"]["sninfo_filename"]
    make_dir_if_not_exists(input_dir)
    make_dir_if_not_exists(output_dir)
    print()
    logger.info(f"ATClean input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Overwrite existing files: {args.overwrite}")
    logger.info(f"Filters: {args.filters}")
    if args.mjd0:
        logger.info(f"MJD0: {args.mjd0}")
    num_controls = (
        args.num_controls
        if not args.num_controls is None
        else int(config["download"]["num_controls"])
    )
    logger.info(f"Number of control light curves to load and plot: {num_controls}")

    creds = Credentials(
        config["credentials"]["atlas_username"],
        config["credentials"]["atlas_password"],
        config["credentials"]["tns_api_key"],
        config["credentials"]["tns_id"],
        config["credentials"]["tns_bot_name"],
    )
    print()
    logger.secret(f"TNS ID: {creds.tns_id}")
    logger.secret(f"TNS bot name: {creds.tns_bot_name}")

    lims = PlotLimits(
        xlower=args.xlim_lower,
        xupper=args.xlim_upper,
        ylower=args.ylim_lower,
        yupper=args.ylim_upper,
    )
    if not lims.is_empty():
        logger.info(f"{lims}")

    cut_list = parse_config_cuts(args, config, colnames)

    print()
    plotloop = PlotLoop(
        colnames,
        input_dir,
        output_dir,
        creds,
        sninfo_filename=sninfo_filename,
        overwrite=args.overwrite,
    )

    plotloop.loop(
        args.tnsnames,
        cut_list,
        num_controls=num_controls,
        mjd0=args.mjd0,
        filters=args.filters,
        plot_uncert_est=args.uncert_est,
        lims=lims if not lims.is_empty() else None,
    )
