from analysis_tools import ObjectCollection, Category, Process, Dataset, Feature, Systematic
from analysis_tools.utils import DotDict
from analysis_tools.utils import join_root_selection as jrs
from plotting_tools import Label
from collections import OrderedDict

from config.legacy_2018 import Config as legacy_config

class Config(legacy_config):

    def add_datasets(self):
        self.tree_name = "tout"
        eos_path = "/eos/user/j/jleonhol/scouting_2022/"
        sample_path = "/vols/cms/pb4918/StoreNTuple/SnTScouting/Data/"
        datasets = [
            Dataset("Scouting2022BLocal",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("Scouting2022B"),
                file_pattern = "output_DataB(.*).root",
                check_empty=False,
            ),
            Dataset("Scouting2022CLocal",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("Scouting2022C"),
                file_pattern = "output_DataC(.*).root",
                check_empty=False,
            ),
            Dataset("Scouting2022DLocal",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("Scouting2022D"),
                file_pattern = "output_DataD(.*).root",
                check_empty=False,
            ),
            Dataset("Scouting2022ELocal",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("Scouting2022E"),
                file_pattern = "output_DataE(.*).root",
                check_empty=False,
            ),
            Dataset("Scouting2022G",
                folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput/",
                process = self.processes.get("Scouting2022G"),
                file_pattern = "output_DataG(.*).root",
                prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty=False,
            ),
            Dataset("Scouting2022GLocal",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("Scouting2022G"),
                file_pattern = "output_DataG(.*).root",
                check_empty=False,
            ),
            Dataset("Scouting2022F",
                folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput/",
                process = self.processes.get("Scouting2022F"),
                file_pattern = "output_DataF(.*).root",
                prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty=False,
            ),
            Dataset("Scouting2022FLocal",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("Scouting2022F"),
                file_pattern = "output_DataF(.*).root",
                check_empty=False,
            ),
            Dataset("DileptonMinBias",
                folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput/",
                process = self.processes.get("DileptonMinBias"),
                file_pattern = "output_DileptonMinBias(.*).root",
                prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty=False,
            ),
            Dataset("DileptonMinBiasLocal",
                #folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutput/",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputAllTrig/",
                process = self.processes.get("DileptonMinBias"),
                file_pattern = "output_DileptonMinBias(.*).root",
                check_empty=False,
            )
        ]
        return ObjectCollection(datasets)
    
    def add_categories(self):
        categories = [
            Category("base", label="base", selection="evtn >= 0"),
            Category("cat1", "Selection, cat. 1", selection="""(DimuonBestLxy > 0.0) && (DimuonBestLxy < 0.2)"""),
            Category("cat2", "Selection, cat. 2", selection="""(DimuonBestLxy > 0.2) && (DimuonBestLxy < 1.0)"""),
            Category("cat3", "Selection, cat. 3", selection="""(DimuonBestLxy > 1.0) && (DimuonBestLxy < 2.4)"""),
            Category("cat4", "Selection, cat. 4", selection="""(DimuonBestLxy > 2.4) && (DimuonBestLxy < 3.1)"""),
            Category("cat5", "Selection, cat. 5", selection="""(DimuonBestLxy > 3.1) && (DimuonBestLxy < 7.0)"""),
            Category("cat6", "Selection, cat. 6", selection="""(DimuonBestLxy > 7.0) && (DimuonBestLxy < 11.0)"""),
            Category("cat7", "Selection, cat. 7", selection="""(DimuonBestLxy > 11.0) && (DimuonBestLxy < 16.0)"""),
            Category("cat8", "Selection, cat. 8", selection="""(DimuonBestLxy > 16.0) && (DimuonBestLxy < 70.0)""")
        ]
        return ObjectCollection(categories)

    def add_features(self):
        features = [
            Feature("JPsi_mass_denom", "DimuonBestMass",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"),
            Feature("JPsi_mass_num", "DimuonBestMass",
                selection="DimuonBestPass==1",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV")
        ]
        return ObjectCollection(features)

    

config = Config("base", year=2022, ecm=13.6, lumi_pb=35200, isUL=True)