import os

from Corrections.JME.PUjetID_SF import PUjetID_SFRDFProducer
from analysis_tools.utils import import_root
import correctionlib

ROOT = import_root()
correctionlib.register_pyroot_binding()

class SnTDimuonsReweightProducer():
    def __init__(self, *args, **kwargs):
        
        filename = "/vols/cms/pb4918/dqcd_analysis/Oct24/dqcd/data/NUM_Data_DEN_MC_subleading_pt_dimuon_eta_2022_syst.json"

        self.runPeriod = kwargs.pop("runPeriod")
        self.isMC = kwargs.pop("isMC")

        if self.isMC:
            if "/libCorretionsWrapper.so" not in ROOT.gSystem.GetLibraries():
                ROOT.gInterpreter.Load("libCorrectionsWrapper.so")

            ROOT.gInterpreter.Declare(os.path.expandvars(
                '#include "$CMSSW_BASE/src/Base/Modules/interface/correctionWrapper.h"'))
            ROOT.gInterpreter.ProcessLine(
                f'auto corr_kinematic = MyCorrections("{os.path.expandvars(filename)}", '
                    '"NUM_Data_DEN_MC_sublead_pt_dimuon_eta_2022_syst");'
            )

            if not os.getenv("_SnTKinematicReweight"):
                os.environ["_SnTKinematicReweight"] = "SnTKinematicReweight"
                ROOT.gInterpreter.Declare("""
                    using Vfloat = ROOT::RVec<float>;
                    using Vint = ROOT::RVec<int>;
                    using Vbool = ROOT::RVec<bool>;
                    //Use variables from the dimuon selection module
                    float get_snt_kinematic_reweight(float EventDimuonsBestSubleadingPt, float EventDimuonsBestEta, std::string syst){
                        return corr_kinematic.eval({EventDimuonsBestSubleadingPt, EventDimuonsBestEta, syst});
                    }
                """)
    
    def run(self, df):
        if self.isMC:
            branches = ['kinematic', 'kinematic_up', 'kinematic_down']
            for branch_name, syst in zip(branches, ["sf", "systup", "systdown"]):
                df = df.Define(branch_name, """get_snt_kinematic_reweight(EventDimuonsBestSubleadingPt, EventDimuonsBestEta, "%s")""" % syst)
        else:
            branches = []
        return df, branches

def SnTDimuonsReweight(**kwargs):
    return lambda: SnTDimuonsReweightProducer(**kwargs)