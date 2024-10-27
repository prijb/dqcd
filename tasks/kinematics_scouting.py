import law
import os
from array import array

from cmt.base_tasks.base import DatasetTaskWithCategory, HTCondorWorkflow, SGEWorkflow, InputData

from analysis_tools.utils import create_file_dir, import_root

class KinematicsScoutingSnT(DatasetTaskWithCategory, law.LocalWorkflow, HTCondorWorkflow, SGEWorkflow):
    def create_branch_map(self):
        return len(self.dataset.get_files(
            os.path.expandvars("$CMT_TMP_DIR/%s/" % self.config_name), add_prefix=False,
            check_empty=True))

    def workflow_requires(self):
        return {"data": InputData.req(self)}

    def requires(self):
        return {"data": InputData.req(self, file_index=self.branch)}

    #def output(self):
    #    return self.local_target(f"data_{self.addendum}{self.branch}.root")

    def output(self):
        return self.local_target("histos_%s.root" % self.branch)

    def add_to_root(self, root):
        root.gInterpreter.Declare("""
        #include "DataFormats/Math/interface/deltaR.h"
        using Vbool = ROOT::RVec<bool>;
        using Vint = ROOT::RVec<int>;
        using Vfloat = ROOT::RVec<float>;
                            
        //Vertex struct
        struct event_vertices{
            //PV information (copies of the same number)
            Vfloat PV_x;
            Vfloat PV_y;
            Vfloat PV_z;
            //Muon information
            Vint mu1idx;
            Vint mu2idx;
            Vfloat mu1pt;
            Vfloat mu2pt;
            Vfloat mu1eta;
            Vfloat mu2eta;
            Vfloat mu1phi;
            Vfloat mu2phi;
            Vfloat mu1phiCorr;
            Vfloat mu2phiCorr;
            Vfloat mu1dxyCorr;
            Vfloat mu2dxyCorr;
            Vfloat mu1dxye;
            Vfloat mu2dxye;
            //SV information
            Vbool passMaterialVeto;
            Vint excesshits;
            //Vint assocSVOverlap;
            Vfloat px; 
            Vfloat py;
            Vfloat pz; 
            Vfloat pt;
            Vfloat x;
            Vfloat y;
            Vfloat z;
            Vfloat mass; 
            Vfloat lxy;
            Vfloat l3d;
            Vfloat prob;
            Vfloat subleading_pt;
            Vfloat dimuon_eta;
            //Gets the highest chi2 vertex index
            int bestvtxidx;
        };

        //Denominator definition (SVs where both the SV and the dimuon pass the preselection and nothing else)
        auto getDenomDimuons(Vint Muon_ch, Vfloat Muon_pt, Vfloat Muon_eta, Vfloat Muon_phi, Vfloat Muon_phiCorr,
        Vfloat Muon_dxyCorr, Vfloat Muon_dxye, Vbool Muon_selected,
        Vint Muon_bestAssocSVOverlapIdx, Vint Muon_bestAssocSVIdx, Vint Muon_nhitsbeforesv,
        float PV_x, float PV_y, float PV_z,
        Vfloat SV_x, Vfloat SV_y, Vfloat SV_z,
        Vfloat SV_xe, Vfloat SV_ye, Vfloat SV_ze,
        Vfloat SV_minDistanceFromDet_x, Vfloat SV_minDistanceFromDet_y, Vfloat SV_minDistanceFromDet_z, Vbool SV_onModuleWithinUnc,
        Vfloat SV_lxy, Vfloat SV_l3d, Vfloat SV_prob, Vbool SV_selected, ROOT::RVec<uint> SV_index){
        
            event_vertices dimuons;
            //Make sure that the same muon is not used for two different vertices
            std::vector<int> muon_sel_idx{}; 

            //Loop over the muons
            for(unsigned i=0; i<Muon_pt.size(); i++){
                //Basic cuts
                if(Muon_bestAssocSVOverlapIdx.at(i) > -1) continue;
                if(std::find(muon_sel_idx.begin(), muon_sel_idx.end(), i) != muon_sel_idx.end()) continue;
                if(Muon_bestAssocSVIdx.at(i) < 0) continue;
                if(!Muon_selected.at(i)) continue;

                for(unsigned j=i+1; j<Muon_pt.size(); j++){
                    if(Muon_bestAssocSVOverlapIdx.at(j) > -1) continue; 
                    if(std::find(muon_sel_idx.begin(), muon_sel_idx.end(), j) != muon_sel_idx.end()) continue;
                    if(Muon_bestAssocSVIdx.at(j) < 0) continue;
                    if(!Muon_selected.at(j)) continue;

                    bool same_vertex = (Muon_bestAssocSVIdx.at(i) == Muon_bestAssocSVIdx.at(j));
                    bool is_neutral = ((Muon_ch.at(i) + Muon_ch.at(j)) == 0);

                    //Skip if they're not neutral or not from the same vertex
                    if(!is_neutral || !same_vertex) continue;

                    //Find the position of the common vertex
                    int common_vertex;
                    for(int vpos=0; vpos<SV_index.size(); vpos++){
                        if(SV_index.at(vpos) == Muon_bestAssocSVIdx.at(i)){
                            common_vertex = vpos;
                            break;
                        }
                    }

                    //Skip the vertex if it fails preselection
                    if(!SV_selected.at(common_vertex)) continue;

                    //Get the material veto info
                    bool pass_material_veto = true;
                    if (SV_onModuleWithinUnc.at(common_vertex)) pass_material_veto = false;
                    if((TMath::Abs(SV_minDistanceFromDet_x.at(common_vertex)) < 0.81) && (TMath::Abs(SV_minDistanceFromDet_y.at(common_vertex)) < 3.24) && (TMath::Abs(SV_minDistanceFromDet_z.at(common_vertex)) < 0.0145)) pass_material_veto = false;

                    //Count the excess hits from the two muons
                    int excess_hits_dimuon = Muon_nhitsbeforesv.at(i) + Muon_nhitsbeforesv.at(j);

                    //Get the momentum of the system
                    ROOT::Math::PtEtaPhiMVector Muon_i_corrected(Muon_pt.at(i), Muon_eta.at(i), Muon_phiCorr.at(i), 0.10566);
                    ROOT::Math::PtEtaPhiMVector Muon_j_corrected(Muon_pt.at(j), Muon_eta.at(j), Muon_phiCorr.at(j), 0.10566);
                    ROOT::Math::PtEtaPhiMVector dimuon_corrected = Muon_i_corrected + Muon_j_corrected;


                    //This vertex will now be stored as the denominator
                    dimuons.PV_x.push_back(PV_x);
                    dimuons.PV_y.push_back(PV_y);
                    dimuons.PV_z.push_back(PV_z);
                    dimuons.mu1idx.push_back(i);
                    dimuons.mu2idx.push_back(j);
                    dimuons.mu1pt.push_back(Muon_pt.at(i));
                    dimuons.mu2pt.push_back(Muon_pt.at(j));
                    dimuons.mu1eta.push_back(Muon_eta.at(i));
                    dimuons.mu2eta.push_back(Muon_eta.at(j));
                    dimuons.mu1phi.push_back(Muon_phi.at(i));
                    dimuons.mu2phi.push_back(Muon_phi.at(j));
                    dimuons.mu1phiCorr.push_back(Muon_phiCorr.at(i));
                    dimuons.mu2phiCorr.push_back(Muon_phiCorr.at(j));
                    dimuons.mu1dxyCorr.push_back(Muon_dxyCorr.at(i));
                    dimuons.mu2dxyCorr.push_back(Muon_dxyCorr.at(j));
                    dimuons.mu1dxye.push_back(Muon_dxye.at(i));
                    dimuons.mu2dxye.push_back(Muon_dxye.at(j));
                    dimuons.passMaterialVeto.push_back(pass_material_veto);
                    dimuons.excesshits.push_back(excess_hits_dimuon);
                    dimuons.px.push_back(dimuon_corrected.px());
                    dimuons.py.push_back(dimuon_corrected.py());
                    dimuons.pz.push_back(dimuon_corrected.pz());
                    dimuons.pt.push_back(dimuon_corrected.pt());
                    dimuons.x.push_back(SV_x.at(common_vertex));
                    dimuons.y.push_back(SV_y.at(common_vertex));
                    dimuons.z.push_back(SV_z.at(common_vertex));
                    dimuons.mass.push_back(dimuon_corrected.M());
                    dimuons.lxy.push_back(SV_lxy.at(common_vertex));
                    dimuons.l3d.push_back(SV_l3d.at(common_vertex));
                    dimuons.prob.push_back(SV_prob.at(common_vertex));
                                  
                    //Get the subleading pt and dimuon eta
                    float subleading_pt = (Muon_pt.at(i) > Muon_pt.at(j)) ? Muon_pt.at(j) : Muon_pt.at(i);
                    TVector3 sv_position(SV_x.at(common_vertex), SV_y.at(common_vertex), SV_z.at(common_vertex));
                    float sv_eta = sv_position.Eta();
                    dimuons.subleading_pt.push_back(subleading_pt);
                    dimuons.dimuon_eta.push_back(sv_eta);
                }
            }

            //Iterate to get best vertex
            float best_prob = -1.;
            for(unsigned i=0; i<dimuons.prob.size(); i++){
                if(dimuons.prob.at(i) > best_prob){
                    best_prob = dimuons.prob.at(i);
                    dimuons.bestvtxidx = i;
                }
            }
            return dimuons;
        }            
        """)

    def run(self):
        ROOT = import_root()
        self.add_to_root(ROOT)

        df = ROOT.RDataFrame("tout", self.input()["data"][0].path)

        mu_trigs = None

        if self.dataset.runPeriod == "2022":
            mu_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"

        elif self.dataset.runPeriod == "2023":
            mu_trigs = "(Run3_DoubleMu3_PFScouting==true)"

        else:
            print("Run period not provided, assuming 2022")
            mu_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"

        df = df.Define("MuTrigger", mu_trigs)
        df = df.Filter("MuTrigger==true")

        #Make the secondary vertices
        df = df.Define("EventDimuons", "getDenomDimuons(Muon_ch, Muon_pt, Muon_eta, Muon_phi, Muon_phiCorr, Muon_dxyCorr, Muon_dxye, Muon_selected, Muon_bestAssocSVOverlapIdx, Muon_bestAssocSVIdx, Muon_nhitsbeforesv, PV_x, PV_y, PV_z, SV_x, SV_y, SV_z, SV_xe, SV_ye, SV_ze, SV_minDistanceFromDet_x, SV_minDistanceFromDet_y, SV_minDistanceFromDet_z, SV_onModuleWithinUnc, SV_lxy, SV_l3d, SV_prob, SV_selected, SV_index)")
        df = df.Define("nDimuons", "EventDimuons.mass.size()")
        #Get the best vertex and corresponding quantities
        df = df.Filter("nDimuons>0").Define("EventDimuonsBestIdx", "EventDimuons.bestvtxidx")
        df = df.Define("EventDimuonsBestMass", "EventDimuons.mass.at(EventDimuonsBestIdx)")
        df = df.Define("EventDimuonsBestSubleadingPt", "EventDimuons.subleading_pt.at(EventDimuonsBestIdx)").Define("EventDimuonsBestEta", "EventDimuons.dimuon_eta.at(EventDimuonsBestIdx)")

        #Define the bins
        pt_bins_2d = array('d', [0, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40, 45, 50])
        eta_bins_2d = array('d', [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

        histos = {}
        histos["h_dimuon_subleading_pt"] = df.Histo1D(("h_dimuon_subleading_pt", "; Subleading muon pT (GeV); Events/2 GeV", 50, 0, 50), "EventDimuonsBestSubleadingPt")
        histos["h_dimuon_eta"] = df.Histo1D(("h_dimuon_eta", "; Dimuon eta; Events/0.1", 50, -2.5, 2.5), "EventDimuonsBestEta")
        histos["h_sublead_pT_dimuon_eta"] = df.Histo2D(("h_sublead_pT_dimuon_eta", "; Subleading muon pT (GeV); Dimuon eta; Events", len(pt_bins_2d)-1, pt_bins_2d, len(eta_bins_2d)-1, eta_bins_2d), "EventDimuonsBestSubleadingPt", "EventDimuonsBestEta")

        histo_file = ROOT.TFile.Open(create_file_dir(self.output().path), "RECREATE")
        for histo in histos.values():
            histo.Write()
        histo_file.Close()
