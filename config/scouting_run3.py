from analysis_tools import ObjectCollection, Category, Process, Dataset, Feature, Systematic
from analysis_tools.utils import DotDict
from analysis_tools.utils import join_root_selection as jrs
from plotting_tools import Label
from collections import OrderedDict

from config.legacy_2018 import Config as legacy_config

class Config(legacy_config):
    def add_datasets(self):
        self.tree_name = "tout"
        datasets = [
            #2022
            Dataset("Scouting2022B",
                #folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput/",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("data"),
                file_pattern = "output_DataB(.*).root",
                #prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty = False,
                runPeriod = "2022"
            ),
            Dataset("Scouting2022C",
                #folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput/",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("data"),
                file_pattern = "output_DataC(.*).root",
                #prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty = False,
                runPeriod = "2022"
            ),
            Dataset("Scouting2022D",
                #folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput/",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("data"),
                file_pattern = "output_DataD(.*).root",
                #prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty = False,
                runPeriod = "2022"
            ),
            Dataset("Scouting2022E",
                #folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput/",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("data"),
                file_pattern = "output_DataE(.*).root",
                #prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty = False,
                runPeriod = "2022"
            ),
            Dataset("Scouting2022F",
                #folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput/",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("data"),
                file_pattern = "output_DataF(.*).root",
                #prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty = False,
                runPeriod = "2022"
            ),
            Dataset("Scouting2022G",
                #folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput/",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutputv2/",
                process = self.processes.get("data"),
                file_pattern = "output_DataG(.*).root",
                #prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty = False,
                runPeriod = "2022"
            ),
            #2023
            Dataset("Scouting2023C",
                #folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput2023/",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutput2023/",
                process = self.processes.get("data"),
                file_pattern = "output_DataC(.*).root",
                #prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty = False,
                runPeriod = "2023"
            ),
            Dataset("Scouting2023D",
                #folder = "/store/user/ppradeep/Run3ScoutingOutput/LooperOutput2023/",
                folder = "/vols/cms/pb4918/StoreNTuple/SnTScouting/LooperOutput2023/",
                process = self.processes.get("data"),
                file_pattern = "output_DataD(.*).root",
                #prefix = "redirector.t2.ucsd.edu:1095/",
                check_empty = False,
                runPeriod = "2023"
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
            Feature("DimuonMass", "DimuonBestMass",
                binning=(15000, 0.0, 150.0),
                x_title=Label("Dimuon mass"),
                units="GeV"),
            Feature("DimuonMassPass", "DimuonBestMass",
                selection="DimuonBestPass==1",
                binning=(15000, 0.0, 150.0),
                x_title=Label("Dimuon mass"),
                units="GeV"),
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