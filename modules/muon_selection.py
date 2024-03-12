from analysis_tools.utils import import_root
ROOT = import_root()

class DQCDMuonSelectionRDFProducer():
    def __init__(self, *args, **kwargs):
        self.year = kwargs.pop("year", 2018)

        ROOT.gInterpreter.Declare("""
            #include "DataFormats/Math/interface/deltaR.h"
            using Vint = const ROOT::RVec<int>&;
            using Vfloat = const ROOT::RVec<float>&;
            ROOT::RVec<int> match_col1_col2(Vfloat pt1, Vfloat eta1, Vfloat phi1,
                    Vfloat pt2, Vfloat eta2, Vfloat phi2, float max_dpt, float max_dr) {
                ROOT::RVec<int> matching(eta1.size(), -1);
                for (auto i = 0; i < eta1.size(); i++) {
                    float min_dR = 999;
                    int min_dR_index = -1;
                    for (auto j = 0; j < eta2.size(); j++) {
                        auto dR = reco::deltaR(eta1[i], phi1[i], eta2[j], phi2[j]);
                        if (dR < max_dr && dR < min_dR && fabs((pt2[j] / pt1[i]) - 1) < max_dpt) {
                            min_dR = dR;
                            min_dR_index = j;
                        }
                    }
                    matching[i] = min_dR_index;
                }
                return matching;
            }
            std::vector<ROOT::RVec<int>> get_leading_elems(Vfloat vec) {
                ROOT::RVec<int> leading(vec.size(), 0);
                ROOT::RVec<int> subleading(vec.size(), 0);
                int lead_index = -1;
                int sublead_index = -1;
                float lead_value = -999.;
                float sublead_value = -999.;
                for (size_t i = 0; i < vec.size(); i++) {
                    if (vec[i] > lead_value) {
                        sublead_value = lead_value;
                        lead_value = vec[i];
                        sublead_index = lead_index;
                        sublead_index = i;
                    } else if (vec[i] > sublead_value) {
                        sublead_value = vec[i];
                        sublead_index = i;
                    }
                }
                if (lead_index != -1) {
                    leading[lead_index] = 1;
                }
                if (sublead_index != -1) {
                    subleading[sublead_index] = 1;
                }
                return {leading, subleading};
            }
        """)

    def run(self, df):
        df = df.Define("MuonBPark_isLooseMuon", """(MuonBPark_looseId == 1) &&
            (MuonBPark_pt > 3.) && (abs(MuonBPark_eta) < 2.5)""")
        df = df.Define("MuonBPark_isMuonWithEtaAndPtReq", """(MuonBPark_isLooseMuon == 1) &&
            (MuonBPark_pt > 5.) && (abs(MuonBPark_eta) < 2.4)""")
        df = df.Define("MuonBPark_isTriggeringMuon", """(MuonBPark_isLooseMuon == 1) &&
            (MuonBPark_pt > 9.) && (abs(MuonBPark_eta) < 1.5 && abs(MuonBPark_sip3d) > 6.)""")

        # filtering
        df = df.Filter("MuonBPark_pt[MuonBPark_isLooseMuon == 1].size() > 0", ">= 1 loose muon")
        df = df.Filter("MuonBPark_pt[MuonBPark_isMuonWithEtaAndPtReq == 1].size() > 0", ">= 1 muon with pt and eta req")
        df = df.Filter("MuonBPark_pt[MuonBPark_isTriggeringMuon == 1].size() > 0", ">= 1 triggering muon")

        # trigger flag
        if self.year == 2018:
            df = df.Define("DisplacedMuonTrigger_flag", " || ".join([
                "HLT_Mu9_IP6_part0",
                "HLT_Mu9_IP6_part1",
                "HLT_Mu9_IP6_part2",
                "HLT_Mu9_IP6_part3",
                "HLT_Mu9_IP6_part4"
            ]))
        df = df.Filter("DisplacedMuonTrigger_flag > 0", "Pass trigger")

        # cpf candidates
        #df = df.Define("cpf_pt", "sqrt(cpf_px * cpf_px + cpf_py * cpf_py)")
        df = df.Define("cpf_p", "sqrt(cpf_px * cpf_px + cpf_py * cpf_py + cpf_pz * cpf_pz)")
        df = df.Define("cpf_eta", "atanh(cpf_pz/cpf_p)")
        df = df.Define("cpf_phi", "atan2(cpf_py, cpf_px)")
        df = df.Define("cpf_mu_dR", "atan2(cpf_py, cpf_px)")
        # df = df.Define("MuonBPark_cpf_match", """match_col1_col2(
        #     MuonBPark_pt, MuonBPark_eta, MuonBPark_phi,
        #     cpf_pt, cpf_eta, cpf_phi,
        #     0.1, 0.02)""")

        df = df.Define("leading", "get_leading_elems(MuonBPark_pt)").Define(
            "MuonBPark_isLeading", "leading[0]").Define("MuonBPark_isSubleading", "leading[1]")

        # match to trigger muons
        # trigger matched
        df = df.Define("MuonBPark_trigger_matched", """(MuonBPark_isTriggeringMuon > 0) &&
            (MuonBPark_isTriggering > 0) && (MuonBPark_fired_HLT_Mu9_IP6 > 0)""")
        df = df.Filter("MuonBPark_pt[MuonBPark_trigger_matched > 0].size() > 0", ">= 1 trigger-matched muon")

        # tighter eta and pt reqs
        df = df.Define("MuonBPark_isMuonWithTighterEtaAndPtReq", """MuonBPark_isLooseMuon == 1 &&
            MuonBPark_pt > 10. && abs(MuonBPark_eta) < 1.5""")

        return df, ["MuonBPark_isLooseMuon", "MuonBPark_isTriggeringMuon",
            "MuonBPark_isMuonWithEtaAndPtReq",
            #"MuonBPark_isMuonWithEtaAndPtReq", "MuonBPark_cpf_match",
            "MuonBPark_isLeading", "MuonBPark_isSubleading",
            "MuonBPark_trigger_matched", "MuonBPark_isMuonWithTighterEtaAndPtReq"]

#Scouting version
class DQCDMuonSelectionScoutingRDFProducer():
    def __init__(self, *args, **kwargs):
        self.year = kwargs.pop("year", 2022)

        ROOT.gInterpreter.Declare("""
            #include "DataFormats/Math/interface/deltaR.h"
            using Vint = const ROOT::RVec<int>&;
            using Vfloat = const ROOT::RVec<float>&;
            ROOT::RVec<int> match_col1_col2(Vfloat pt1, Vfloat eta1, Vfloat phi1,
                    Vfloat pt2, Vfloat eta2, Vfloat phi2, float max_dpt, float max_dr) {
                ROOT::RVec<int> matching(eta1.size(), -1);
                for (auto i = 0; i < eta1.size(); i++) {
                    float min_dR = 999;
                    int min_dR_index = -1;
                    for (auto j = 0; j < eta2.size(); j++) {
                        auto dR = reco::deltaR(eta1[i], phi1[i], eta2[j], phi2[j]);
                        if (dR < max_dr && dR < min_dR && fabs((pt2[j] / pt1[i]) - 1) < max_dpt) {
                            min_dR = dR;
                            min_dR_index = j;
                        }
                    }
                    matching[i] = min_dR_index;
                }
                return matching;
            }
            std::vector<ROOT::RVec<int>> get_leading_elems(Vfloat vec) {
                ROOT::RVec<int> leading(vec.size(), 0);
                ROOT::RVec<int> subleading(vec.size(), 0);
                int lead_index = -1;
                int sublead_index = -1;
                float lead_value = -999.;
                float sublead_value = -999.;
                for (size_t i = 0; i < vec.size(); i++) {
                    if (vec[i] > lead_value) {
                        sublead_value = lead_value;
                        lead_value = vec[i];
                        sublead_index = lead_index;
                        sublead_index = i;
                    } else if (vec[i] > sublead_value) {
                        sublead_value = vec[i];
                        sublead_index = i;
                    }
                }
                if (lead_index != -1) {
                    leading[lead_index] = 1;
                }
                if (sublead_index != -1) {
                    subleading[sublead_index] = 1;
                }
                return {leading, subleading};
            }
        """)

    def run(self, df):
        df = df.Define("Muon_isMuonWithEtaAndPtReq", """(Muon_pt > 3.) && (abs(Muon_eta) < 2.4)""")

        # filtering to have two muons passing this cut
        df = df.Filter("Muon_pt[Muon_isMuonWithEtaAndPtReq == 1].size() > 1", ">= 2 pT>3 and |eta|<2.4 muons")

        # trigger flag
        if (self.year == 2022):
            df = df.Define("ScoutingMuonTrigger_flag", " || ".join([
                "L1_DoubleMu_15_7",
                "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7",
                "L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18",
                "L1_DoubleMu4_SQ_OS_dR_Max1p2",
                "L1_DoubleMu4p5_SQ_OS_dR_Max1p2"
            ]))
        
        elif (self.year == 2023):
            df = df.Define("ScoutingMuonTrigger_flag", "DST_Run3_DoubleMu3_PFScoutingPixelTracking")

        df = df.Filter("ScoutingMuonTrigger_flag > 0", "Pass Scouting triggers")


        df = df.Define("leading", "get_leading_elems(Muon_pt)").Define(
            "Muon_isLeading", "leading[0]").Define("Muon_isSubleading", "leading[1]")

        #No trigger matching involved since they are trigger muons already

        # tighter eta and pt reqs
        df = df.Define("Muon_isMuonWithTighterEtaAndPtReq", """Muon_pt > 8. && abs(Muon_eta) < 2.0""")

        return df, ["Muon_isMuonWithEtaAndPtReq",
            "Muon_isLeading", "Muon_isSubleading",
            "Muon_isMuonWithTighterEtaAndPtReq"]

def DQCDMuonSelectionRDF(*args, **kwargs):
    return lambda: DQCDMuonSelectionRDFProducer(*args, **kwargs)

def DQCDMuonSelectionScoutingRDF(*args, **kwargs):
    return lambda: DQCDMuonSelectionScoutingRDFProducer(*args, **kwargs)

