##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
import os
import shutil
import glob
import json
import numpy
from matplotlib import pyplot
from .config import Config
from .messaging import Messaging


##########################################################
class Plotting:
    def __init__(self, config: Config):
        self.config = config
        self.loaded_files = []

    def get_file_path(self, case_name: str, variant_name: str) -> str:
        # name
        name = f"{case_name}-{variant_name}"
        results_pattern = self.config.results_pattern

        # calc
        fname = f"result-{name}.json"
        avail = glob.glob(f"{results_pattern}/*/*/{fname}")

        # has none
        if avail is None or len(avail) == 0:
            return None

        # we can sort by date and take the last one
        avail.sort()
        fpath = avail[-1]

        # ok
        return fpath

    def load_case_variant_data(self, case_name: str, variant_name: str) -> dict:
        # get path
        fpath = self.get_file_path(case_name, variant_name)

        # if not found skip
        if fpath is None:
            return None

        # progress
        Messaging.step(f"Loading {fpath}")

        # keep track
        self.loaded_files.append(fpath)

        # load it
        with open(fpath, "r") as fp:
            # load
            return json.load(fp)

    def load_variants(self, data: dict, case_name: str) -> bool:
        # extract vars
        variant_names = self.config.config["variants"].keys()

        # loop on all variants
        for variant_name in variant_names:
            # load json
            case_variant_data = self.load_case_variant_data(case_name, variant_name)

            # merge
            if case_variant_data is not None:
                data["variants"].append(variant_name)
                data["means"].append(case_variant_data["results"][0]["mean"])
                data["median"].append(case_variant_data["results"][0]["median"])
                data["stddev"].append(case_variant_data["results"][0]["stddev"])
                data["min"].append(case_variant_data["results"][0]["min"])
                data["max"].append(case_variant_data["results"][0]["max"])
                data["q1"].append(
                    numpy.quantile(case_variant_data["results"][0]["times"], 0.25)
                )
                data["q3"].append(
                    numpy.quantile(case_variant_data["results"][0]["times"], 0.75)
                )

        # ok
        return True

    def build_data(self):
        # extract vars
        # case_names = self.config.case_names
        case_names = self.config.config["cases"].keys()

        # to fill
        data = {}

        # loop
        for name in case_names:
            # build
            case_data = {
                "variants": [],
                "means": [],
                "median": [],
                "min": [],
                "max": [],
                "q1": [],
                "q3": [],
                "stddev": [],
            }

            # load
            self.load_variants(case_data, name)

            # attach
            if len(case_data["variants"]) != 0:
                data[name] = case_data

        # return
        return data

    def keep_info_of_previous_runs(self):
        # vars
        results = self.config.results

        # dump list of loaded paths
        with open(f"{results}/loaded-result-files.json", "w+") as fp:
            json.dump(self.loaded_files, fp=fp, indent="\t")

        # copy the loaded files localy to get a full view
        target = f"{results}/imported_previous_runs"
        os.makedirs(target, exist_ok=True)
        for file in self.loaded_files:
            # only the non local ones
            if results not in file:
                target_path = os.path.join(target, os.path.basename(file))
                shutil.copyfile(file, target_path)

    def plot(self):
        # info
        Messaging.section("Plotting performances")

        # extract
        cases = self.config.config["cases"]

        # get data
        data = self.build_data()

        # plot
        for case_name, case_data in data.items():
            # progress
            Messaging.step(f"Plotting {case_name}")

            # extract title
            case_config = cases[case_name]
            results = self.config.results
            title = case_config.get("title", case_name)

            # start plot
            fig, ax = pyplot.subplots()


            # build graph
            x_pos = numpy.arange(len(case_data["variants"]))
            # ax.bar(x_pos, entry['median'], yerr=[entry['min'],entry['max']], align='center', alpha=0.5, ecolor='black', capsize=10)
            ax.bar(
                x_pos,
                case_data["median"],
                yerr=case_data["stddev"],
                align="center",
                alpha=0.5,
                ecolor="black",
                capsize=10,
            )
            ax.set_ylabel("Average runtime (seconds)")
            ax.set_xticks(x_pos)
            ax.set_xticklabels(case_data["variants"], rotation=45, ha="right")
            ax.set_title(f"CROCO - {title}")
            ax.yaxis.grid(True)

            # Save the figure and show
            pyplot.tight_layout()
            os.makedirs(f"{results}/{case_name}", exist_ok=True)
            pyplot.savefig(f"{results}/{case_name}/plot-{case_name}.png")
            pyplot.savefig(f"{results}/{case_name}/plot-{case_name}.svg")
            pyplot.close()

            # keep track of the previous run used results
            if not self.config.no_previous:
                self.keep_info_of_previous_runs()
