import argparse
from typing import Dict, List

from clean import parse_config_cuts, parse_config_filters
from download import Credentials, load_config, make_dir_if_not_exists
from lightcurve import (
    AveragedSupernova,
    CutList,
    SnInfoTable,
    Supernova,
    get_mjd0,
    get_tns_coords_from_json,
    get_tns_mjd0_from_json,
    query_tns,
)
from plot import PlotPdf


class PlotLoop:
    def __init__(
        self,
        input_dir: str,
        output_dir: str,
        credentials: Credentials,
        sninfo_filename: str = None,
        overwrite: bool = False,
    ):
        self.sn: Supernova = None
        self.avg_sn: AveragedSupernova = None
        self.cut_list: CutList = None
        self.p: PlotPdf = None

        self.credentials: Credentials = credentials
        self.input_dir: str = input_dir
        self.output_dir: str = output_dir
        self.overwrite: bool = overwrite

        self.sninfo: SnInfoTable = SnInfoTable(
            self.output_dir, filename=sninfo_filename
        )

    def plot_lcs(self, tnsname: str, mjd0: float, filt: str, num_controls: int = 0):
        # TODO
        pass

    def loop(
        self,
        tnsnames: List[str],
        num_controls: int = 0,
        mjd0=None,
        filters: List[str] = ["o", "c"],
        cut_list: CutList = None,
    ):
        self.cut_list = cut_list

        for obj_index in range(len(tnsnames)):
            tnsname = tnsnames[obj_index]
            print(f"\n\tPLOTTING LIGHT CURVES FOR: SN {tnsname}")

            make_dir_if_not_exists(f"{output_dir}/{tnsname}")

            if mjd0 is None:
                mjd0, coords = get_mjd0(tnsname, self.sninfo, self.credentials)
                if not coords is None:
                    print(f"Setting MJD0 to {mjd0}")
                    self.sninfo.update_row(tnsname, coords=coords, mjd0=mjd0)
                    print("Success")
            else:
                print(f"\nSetting MJD0 to {mjd0}")

            for filt in filters:
                self.plot_lcs(
                    tnsname,
                    mjd0,
                    filt,
                    num_controls=num_controls,
                )


# define command line arguments
def define_args(parser=None, usage=None, conflict_handler="resolve"):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage, conflict_handler=conflict_handler)

    parser.add_argument(
        "tnsnames", nargs="+", help="TNS names of the transients to clean"
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
        "-t",
        "--template_correction",
        default=False,
        action="store_true",
        help="plot ATLAS template change correction",
    )
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
        "--custom_cuts",
        default=False,
        action="store_true",
        help="scan config file for custom cuts and plot them",
    )

    return parser


if __name__ == "__main__":
    args = define_args().parse_args()
    config = load_config(args.config_file)

    # TODO

    if len(args.tnsnames) < 1:
        raise RuntimeError("ERROR: Please specify at least one TNS name to plot.")
    if len(args.tnsnames) > 1 and not args.mjd0 is None:
        raise RuntimeError(
            f"ERROR: Cannot specify one MJD0 {args.mjd0} for a batch of SNe."
        )
    print(f"\nList of transients to plot: {args.tnsnames}")

    input_dir = config["dir"]["atclean_input"]
    output_dir = config["dir"]["output"]
    sninfo_filename = config["dir"]["sninfo_filename"]
    make_dir_if_not_exists(input_dir)
    make_dir_if_not_exists(output_dir)
    print(f"\nATClean input directory: {input_dir}")
    print(f"Output directory: {output_dir}")

    print(f'TNS ID: {config["credentials"]["tns_id"]}')
    print(f'TNS bot name: {config["credentials"]["tns_bot_name"]}')

    print(f"Overwrite existing files: {args.overwrite}")
    filters = parse_config_filters(args, config)
    print(f"Filters: {filters}")
    if args.mjd0:
        print(f"MJD0: {args.mjd0}")
    num_controls = (
        args.num_controls
        if args.num_controls
        else int(config["download"]["num_controls"])
    )
    print(f"Number of control light curves to load and plot: {num_controls}")

    # TODO: fix
    cut_list = parse_config_cuts(args, config)

    print()
    credentials = Credentials(
        config["credentials"]["atlas_username"],
        config["credentials"]["atlas_password"],
        config["credentials"]["tns_api_key"],
        config["credentials"]["tns_id"],
        config["credentials"]["tns_bot_name"],
    )
    plotloop = PlotLoop(
        input_dir,
        output_dir,
        credentials,
        sninfo_filename=sninfo_filename,
        overwrite=args.overwrite,
    )

    plotloop.loop(
        args.tnsnames,
        cut_list=cut_list,
        num_controls=num_controls,
        mjd0=args.mjd0,
        filters=filters,
    )
