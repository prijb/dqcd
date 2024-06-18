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
            Category("presel_cat1", "Preselection, cat. 1", selection="""DimuonLxyDenom[(DimuonLxyDenom > 0.0) && (DimuonLxyDenom < 0.2)].size() > 0"""),
            Category("presel_cat4", "Preselection, cat. 4", selection="""DimuonLxyDenom[(DimuonLxyDenom > 2.4) && (DimuonLxyDenom < 3.1)].size() > 0"""),
            Category("sel_cat1", "Selection, cat. 1", selection="""(DimuonLxyPassBest > 0.0) && (DimuonLxyPassBest < 0.2)"""),
            Category("sel_cat4", "Selection, cat. 4", selection="""(DimuonLxyPassBest > 2.4) && (DimuonLxyPassBest < 3.1)""")
        ]
        return ObjectCollection(categories)

    def add_features(self):
        features = [
            Feature("JPsi_mass_denom", "DimuonMassDenom",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"),
            Feature("JPsi_mass_num", "DimuonMassPassBest",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"),
            # Add lxy binned JPsi mass
            Feature("JPsi_mass_denom_lxy1", "DimuonMassDenom[DimuonLxyDenom > 0.0 && DimuonLxyDenom < 0.2]",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"),
            Feature("JPsi_mass_num_lxy1", "DimuonMassPassBest",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                selection="DimuonLxyPassBest > 0.0 && DimuonLxyPassBest < 0.2",
                units="GeV"),
            Feature("JPsi_mass_denom_lxy4", "DimuonMassDenom[DimuonLxyDenom > 2.4 && DimuonLxyDenom < 3.1]",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"),
            Feature("JPsi_mass_num_lxy4", "DimuonMassPassBest",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                selection="DimuonLxyPassBest > 2.4 && DimuonLxyPassBest < 3.1",
                units="GeV"),
        ]
        return ObjectCollection(features)

    

config = Config("base", year=2022, ecm=13.6, lumi_pb=35200, isUL=True)