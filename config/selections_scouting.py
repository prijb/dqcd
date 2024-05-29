from analysis_tools import ObjectCollection, Category, Process, Dataset, Feature, Systematic
from analysis_tools.utils import DotDict
from analysis_tools.utils import join_root_selection as jrs
from plotting_tools import Label
from collections import OrderedDict

from config.legacy_2018 import Config as legacy_config

class Config(legacy_config):

    def add_datasets(self):
        eos_path = "/eos/user/j/jleonhol/scouting_2022/"
        sample_path = "/vols/cms/pb4918/StoreNTuple/SnTScouting/Data/"
        datasets = [
            Dataset("Scouting",
                folder = sample_path,
                process = self.processes.get("Scouting2022F"),
                file_pattern = "output_Data(.*).root",
                check_empty=False,
            ),
            #Dataset("Scouting",
            #    prefix = "eosuser.cern.ch",
            #    folder = eos_path,
            #    process = self.processes.get("Scouting2022F"),
            #    file_pattern = "output_Data(.*).root",
            #    check_empty=False,
            #),
        ]
        return ObjectCollection(datasets)

config = Config("base", year=2022, ecm=13.6, lumi_pb=35200, isUL=True)