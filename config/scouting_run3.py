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
            Category("cat1", "Selection, cat. 1", selection="""(EventDimuonsBestLxy > 0.0) && (EventDimuonsBestLxy < 0.2)"""),
            Category("cat2", "Selection, cat. 2", selection="""(EventDimuonsBestLxy > 0.2) && (EventDimuonsBestLxy < 1.0)"""),
            Category("cat3", "Selection, cat. 3", selection="""(EventDimuonsBestLxy > 1.0) && (EventDimuonsBestLxy < 2.4)"""),
            Category("cat4", "Selection, cat. 4", selection="""(EventDimuonsBestLxy > 2.4) && (EventDimuonsBestLxy < 3.1)"""),
            Category("cat5", "Selection, cat. 5", selection="""(EventDimuonsBestLxy > 3.1) && (EventDimuonsBestLxy < 7.0)"""),
            Category("cat6", "Selection, cat. 6", selection="""(EventDimuonsBestLxy > 7.0) && (EventDimuonsBestLxy < 11.0)"""),
            Category("cat7", "Selection, cat. 7", selection="""(EventDimuonsBestLxy > 11.0) && (EventDimuonsBestLxy < 16.0)"""),
            Category("cat8", "Selection, cat. 8", selection="""(EventDimuonsBestLxy > 16.0) && (EventDimuonsBestLxy < 70.0)""")
        ]
        return ObjectCollection(categories)

    def add_features(self):
        features = [
            #Dimuon displacement
            Feature("DimuonLxy", "EventDimuonsBestLxy",
                binning=(1200, 0.0, 120.0),
                x_title=Label("Dimuon Lxy"),
                units="cm"
            ),
            Feature("DimuonLxyPass", "EventDimuonsBestLxy",
                selection="EventDimuonsBestPass==true",
                binning=(1200, 0.0, 120.0),
                x_title=Label("Dimuon Lxy"),
                units="cm"
            ),
            Feature("DimuonLxyPassMaterialVeto", "EventDimuonsBestLxy",
                selection="EventDimuonsBestPassMaterialVeto==true",
                binning=(1200, 0.0, 120.0),
                x_title=Label("Dimuon Lxy"),
                units="cm"
            ),
            #Dimuon mass
            Feature("DimuonMass", "EventDimuonsBestMass",
                binning=(15000, 0.0, 150.0),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("DimuonMassPass", "EventDimuonsBestMass",
                selection="EventDimuonsBestPass==true",
                binning=(15000, 0.0, 150.0),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("DimuonMassPassMinusDisplacement", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusDisplacement==true",
                binning=(15000, 0.0, 150.0),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("DimuonMassPassMinusMaterialVeto", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusMaterialVeto==true",
                binning=(15000, 0.0, 150.0),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("DimuonMassPassMinusExcessHits", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusExcessHits==true",
                binning=(15000, 0.0, 150.0),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("DimuonMassPassMinusPileupVeto", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusPileupVeto==true",
                binning=(15000, 0.0, 150.0),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("DimuonMassPassMinusKinematics", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusKinematics==true",
                binning=(15000, 0.0, 150.0),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            #JPsi
            Feature("JpsiMass", "EventDimuonsBestMass",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("JpsiMassPass", "EventDimuonsBestMass",
                selection="EventDimuonsBestPass==true",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("JpsiMassPassMinusDisplacement", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusDisplacement==true",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("JpsiMassPassMinusMaterialVeto", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusMaterialVeto==true",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("JpsiMassPassMinusExcessHits", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusExcessHits==true",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("JpsiMassPassMinusPileupVeto", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusPileupVeto==true",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"
            ),
            Feature("JpsiMassPassMinusKinematics", "EventDimuonsBestMass",
                selection="EventDimuonsBestPassMinusKinematics==true",
                binning=(30, 2.8, 3.4),
                x_title=Label("Dimuon mass"),
                units="GeV"
            )
        ]
        return ObjectCollection(features)

    

config = Config("base", year=2022, ecm=13.6, lumi_pb=35200, isUL=True)