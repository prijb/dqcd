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
                dataset = "/ScoutingPFRun3/ppradeep-2022F-f0e3b442cad690715385ee64027967ff/USER",
                process = self.processes.get("Scouting2022F"),
                check_empty=False,
            ),
        ]
        return ObjectCollection(datasets)

config = Config("base", year=2022, ecm=13.6, lumi_pb=35200, isUL=True)