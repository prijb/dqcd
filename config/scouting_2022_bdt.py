from analysis_tools import ObjectCollection, Category, Process, Dataset, Feature, Systematic
from analysis_tools.utils import DotDict
from analysis_tools.utils import join_root_selection as jrs
from plotting_tools import Label
from collections import OrderedDict

from config.ul_2018 import Config as ultralegacy_config

class Config(ultralegacy_config):

    def add_regions(self, **kwargs):
        regions = [
            Category("loose_bdt", "Loose bdt region", selection="{{bdt}} > 0.45"),
            Category("tight_bdt", "Tight bdt region", selection="{{bdt}} > 0.99"),

            Category("loose_bdt_scenarioA", "Loose bdt (A) region", selection="{{bdt_scenarioA}} > 0.65"),
            Category("tight_bdt_scenarioA", "Tight bdt (A) region", selection="{{bdt_scenarioA}} > 0.98"),

            #Category("loose_bdt_scenarioB1", "Loose bdt (B1) region", selection="{{bdt_scenarioB1}} > 0.6"),
            #Category("tight_bdt_scenarioB1", "Tight bdt (B1) region", selection="{{bdt_scenarioB1}} > 0.98"),

            #Category("loose_bdt_scenarioB2", "Loose bdt (B2) region", selection="{{bdt_scenarioB2}} > 0.55"),
            #Category("tight_bdt_scenarioB2", "Tight bdt (B2) region", selection="{{bdt_scenarioB2}} > 0.92"),

            #Category("loose_bdt_scenarioC", "Loose bdt (C) region", selection="{{bdt_scenarioC}} > 0.55"),
            #Category("tight_bdt_scenarioC", "Tight bdt (C) region", selection="{{bdt_scenarioC}} > 0.8"),
        ]
        return ObjectCollection(regions)

    def add_datasets(self):
        sample_path = "/vols/cms/pb4918/qcd_scouting_new/"

        samples = {
            "qcd_1000toInf": ("qcd_1000toInf"),
            "qcd_120to170": ("qcd_120to170"),
            "qcd_15to20": ("qcd_15to20"),
            "qcd_170to300": ("qcd_170to300"),
            "qcd_20to30": ("qcd_20to30"),
            "qcd_300to470": ("qcd_300to470"),
            "qcd_30to50": ("qcd_30to50"),
            "qcd_470to600": ("qcd_470to600"),
            "qcd_50to80": ("qcd_50to80"),
            "qcd_600to800": ("qcd_600to800"),
            "qcd_800to1000": ("qcd_800to1000"),
            "qcd_80to120": ("qcd_80to120"),
        }

        #170/300, 470/600, and 1000/Inf still using 2018 xs
        xs = {
            "qcd_15to20": 2982000,
            "qcd_20to30": 2679000,
            "qcd_30to50": 1497000,
            "qcd_50to80": 402900,
            "qcd_80to120": 96200,
            "qcd_120to170": 22980,
            "qcd_170to300": 7055,
            "qcd_300to470": 699.1,
            "qcd_470to600": 59.24,
            "qcd_600to800": 21.37,
            "qcd_800to1000": 3.913,
            "qcd_1000toInf": 1.078,
        }

        #Adding datasets (skip till 101 because those are used for training BDT)
        datasets = [
            Dataset("qcd_1000toInf",
                folder=sample_path + samples["qcd_1000toInf"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_1000toInf"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_1000toInf"),
                xs=xs["qcd_1000toInf"],
                merging={
                    "base": 10,
                },
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_120to170",
                folder=sample_path + samples["qcd_120to170"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_120to170"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_120to170"),
                xs=xs["qcd_120to170"],
                merging={
                    "base": 10,
                },
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_15to20",
                folder=sample_path + samples["qcd_15to20"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_15to20"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_15to20"),
                xs=xs["qcd_15to20"],
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_170to300",
                folder=sample_path + samples["qcd_170to300"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_170to300"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_170to300"),
                xs=xs["qcd_170to300"],
                merging={
                    "base": 10,
                },
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_20to30",
                folder=sample_path + samples["qcd_20to30"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_20to30"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_20to30"),
                xs=xs["qcd_20to30"],
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_300to470",
                folder=sample_path + samples["qcd_300to470"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_300to470"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_300to470"),
                xs=xs["qcd_300to470"],
                merging={
                    "base": 10,
                },
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_30to50",
                folder=sample_path + samples["qcd_30to50"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_30to50"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_30to50"),
                xs=xs["qcd_30to50"],
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_470to600",
                folder=sample_path + samples["qcd_470to600"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_470to600"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_470to600"),
                xs=xs["qcd_470to600"],
                merging={
                    "base": 10,
                },
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_50to80",
                folder=sample_path + samples["qcd_50to80"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_50to80"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_50to80"),
                xs=xs["qcd_50to80"],
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_600to800",
                folder=sample_path + samples["qcd_600to800"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_600to800"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_600to800"),
                xs=xs["qcd_600to800"],
                merging={
                    "base": 10,
                },
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_80to120",
                folder=sample_path + samples["qcd_80to120"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_80to120"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_80to120"),
                xs=xs["qcd_80to120"],
                merging={
                    "base": 10,
                },
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("qcd_800to1000",
                folder=sample_path + samples["qcd_800to1000"],
                skipFiles=["{}/output_{}.root".format(
                    sample_path + samples["qcd_800to1000"], i)
                    for i in range(1, 101)],
                process=self.processes.get("qcd_800to1000"),
                xs=xs["qcd_800to1000"],
                merging={
                    "base": 10,
                    "singlev_cat3": 10
                },
                skipped_files_must_be_in_dataset=False,
            ),

            Dataset("Scouting22F",
                folder=[
                    "/vols/cms/pb4918/" + "StoreNTuple/Scouting22F/Nanotronv2/",
                ],
                process=self.processes.get("data"),
                merging={
                    "base": 20,
                },
                tags=["ul"],
                runPeriod="F",
            ),

            Dataset("scenarioA_mpi_4_mA_1p33_ctau_10",
                folder = "/vols/cms/pb4918/" + "StoreNTuple/DQCD_v2/Run3/Scouting/scenarioA_mpi_4_mA_1p33_ctau_10/",
                #dataset = "/scenarioA_mpi_4_mA_1p33_ctau_10/jleonhol-nanotron-3b50327cf5b3a9483d26e0670720126c/USER",
                process=self.processes.get("scenarioA_mpi_4_mA_1p33_ctau_10"),
                check_empty=False,
                skipFiles=["{}/output_{}.root".format(
                    "/vols/cms/pb4918/" + "StoreNTuple/DQCD_v2/Run3/Scouting/scenarioA_mpi_4_mA_1p33_ctau_10", i)
                    for i in range(1, 61)],
                #skipFiles=[f"/store/user/jleonhol/samples/nanotron/scenarioA_mpi_4_mA_1p33_ctau_10/nanotron/231124_165003/0000/nano_{i}.root"
                #    for i in range(1, 21)],
                #prefix="gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms",
                xs=0.522,
            ),
        ]
        return ObjectCollection(datasets)

    def add_weights(self):
        weights = DotDict()
        weights.default = "1"

        #weights.total_events_weights = ["puWeight"]
        # weights.total_events_weights = ["genWeight"]
        weights.total_events_weights = ["1"]

        #weights.base = ["puWeight", "PUjetID_SF", "idWeight", "trigSF"]  # others needed
        weights.base = ["1"]  # others needed

        for category in self.categories:
            weights[category.name] = weights.base

        return weights

    # other methods

# config = Config("base", year=2018, ecm=13, lumi_pb=59741)
config = Config("base", year=2022, ecm=13.6, lumi_pb=35200, isUL=True)

