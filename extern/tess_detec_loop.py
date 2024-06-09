import random
from typing import Dict, List
from generate_detec_table import SimDetecLoop, define_args
from generate_sim_table import load_json_config, parse_params
from lightcurve import SimDetecSupernova, SimDetecLightCurve, Simulation


class TessSimDetecLoop(SimDetecLoop):
    def __init__(self, sigma_kerns: List, **kwargs):
        super().__init__(sigma_kerns, **kwargs)

    def set_peak_mags_and_fluxes(
        self,
        model_name: str = None,
        sim_tables_dir: str | None = None,
        detec_tables_dir: str | None = None,
        **kwargs,
    ):
        return super().set_peak_mags_and_fluxes(
            model_name, sim_tables_dir, detec_tables_dir, **kwargs
        )

    def load_sn(
        self,
        data_dir: str,
        tnsname: str,
        num_controls: int,
        mjdbinsize: float = 1,
        filt: str = "o",
    ):
        return super().load_sn(data_dir, tnsname, num_controls, mjdbinsize, filt)

    def load_sd(
        self,
        model_name: str,
        sim_tables_dir: str | None = None,
        detec_tables_dir: str | None = None,
    ):
        return super().load_sd(model_name, sim_tables_dir, detec_tables_dir)

    def load_sim(self, table_row: Dict) -> Simulation:
        return super().load_sim(table_row)

    def get_max_fom_indices(self, sim_lc: SimDetecLightCurve, **kwargs):
        return super().get_max_fom_indices(sim_lc, **kwargs)

    def update_sd_row(
        self,
        sigma_kern: float,
        peak_appmag: float,
        index: int,
        control_index: int,
        max_fom: float,
        max_fom_mjd: float,
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
        self, valid_control_ix: List, detec_tables_dir: str, flag=0x800000, **kwargs
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

    simdetec = TessSimDetecLoop(sigma_kerns)
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
