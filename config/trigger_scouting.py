from analysis_tools import ObjectCollection, Category, Process, Dataset, Feature, Systematic
from analysis_tools.utils import DotDict
from analysis_tools.utils import join_root_selection as jrs
from plotting_tools import Label
from collections import OrderedDict

from config.legacy_2018 import Config as legacy_config

class Config(legacy_config):

    def add_datasets(self):
        datasets = [
            Dataset("Scouting2022F",
                dataset = "/ScoutingPFRun3/ppradeep-2022F-3588c69a11675da582eb4e73ebc9784b/USER",
                process = self.processes.get("Scouting2022F"),
                check_empty=False,
            ),
            Dataset("DileptonMinBias2022",
                dataset = "/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/ppradeep-2022-8c2f5715331d2ea6d01c115b141b34b0/USER",
                process = self.processes.get("DileptonMinBias2022"),
                check_empty=False,
            ),
        ]
        return ObjectCollection(datasets)

config = Config("base", year=2022, ecm=13.6, lumi_pb=35200, isUL=True)