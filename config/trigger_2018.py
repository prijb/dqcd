from analysis_tools import ObjectCollection, Category, Process, Dataset, Feature, Systematic
from analysis_tools.utils import DotDict
from analysis_tools.utils import join_root_selection as jrs
from plotting_tools import Label
from collections import OrderedDict

from config.legacy_2018 import Config as legacy_config


class Config(legacy_config):

    def add_datasets(self):

        datasets = [
            Dataset("LambdaBToJpsiLambda",
                dataset = "/LambdaBToJpsiLambda_JpsiToMuMu_TuneCP5_13TeV-pythia8-evtgen/jleonhol-JPsiMC-f7d89a2f103706ffd8ab9007e324774d/USER", 
                process=self.processes.get("LambdaBToJpsiLambda"),
                check_empty=False,
            ),
            Dataset("egamma",
                dataset = "/EGamma/jleonhol-egamma-6838dccfc09128c984377350a24eba4a/USER",
                process=self.processes.get("egamma"),
                check_empty=False,
            ),
            Dataset("egamma_eraA",
                dataset = "/EGamma/jleonhol-egammaA-6838dccfc09128c984377350a24eba4a/USER",
                process=self.processes.get("egamma"),
                check_empty=False,
            ),
            Dataset("egamma_eraB",
                dataset = "/EGamma/jleonhol-egammaB-6838dccfc09128c984377350a24eba4a/USER",
                process=self.processes.get("egamma"),
                check_empty=False,
            ),
            Dataset("egamma_eraC",
                dataset = "/EGamma/jleonhol-egammaC-6838dccfc09128c984377350a24eba4a/USER" ,
                process=self.processes.get("egamma"),
                check_empty=False,
            ),
            Dataset("SingleMuon",
                dataset = "/SingleMuon/jleonhol-SingleMuon-9e237a9f56e9d058139fa37ffdb183f0/USER" ,
                process=self.processes.get("SingleMuon"),
                check_empty=False,
            ),
            Dataset("BToJpsiJPsiToMuMu",
                dataset = "/BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen/jleonhol-BJPsiMC-f7d89a2f103706ffd8ab9007e324774d/USER",
                process=self.processes.get("BToJpsiJPsiToMuMu"),
                check_empty=False,
            ),
            Dataset("BuToJpsiK",
                dataset = "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/jleonhol-nanotronv2-571c6e4dc467acb2f3a7892cb8ebd34e/USER",
                process=self.processes.get("BuToJpsiK"),
                check_empty=False,
            ),
            Dataset("BuToJpsiK_inclusive_muonSV",
                dataset = "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/jleonhol-nanotronv2_v3-209802c58913bedf9ef81b726e08129a/USER",
                #dataset = "/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/jleonhol-nanotronv2_slimmedmuons-c3e54d33997269b36a54d7557d26fd32/USER",
                process=self.processes.get("BuToJpsiK"),
                prefix="gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms",
                check_empty=False,
            ),
            Dataset("data_2018d_bph1",
                dataset="/ParkingBPH1/jleonhol-nanotronv2-205145b8a3c6bd3ea858a0dbe549c313/USER",
                process=self.processes.get("data"),
                merging={
                    "base": 20,
                },
                tags=["ul"],
                runPeriod="D",
                prefix="gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms",
            ),
        ]
        return ObjectCollection(datasets)

    # other methods

# config = Config("base", year=2018, ecm=13, lumi_pb=59741)
config = Config("base", year=2018, ecm=13, lumi_pb=33600, isUL=True)
