import law
import os
from array import array

from cmt.base_tasks.base import DatasetTaskWithCategory, HTCondorWorkflow, SGEWorkflow, InputData

from analysis_tools.utils import create_file_dir, import_root

class DimuonKinematicsSnT(DatasetTaskWithCategory, law.LocalWorkflow, HTCondorWorkflow, SGEWorkflow):
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
        #include "TLorentzVector.h"
        #include "TVector3.h"
        #include "DataFormats/Math/interface/deltaR.h"
        using Vbool = ROOT::RVec<bool>;
        using Vint = ROOT::RVec<int>;
        using Vfloat = ROOT::RVec<float>;
                            
        //Vertex struct
        struct event_vertices{
            Vint mu1idx;
            Vint mu2idx;
            Vint mu3idx;
            Vint mu4idx;
            Vint passveto;
            Vint excesshits;
            Vint category;
            Vint assocSVOverlap;
            Vfloat px; 
            Vfloat py;
            Vfloat pz; 
            Vfloat x;
            Vfloat y;
            Vfloat z;
            Vfloat mass; 
            Vfloat lxy;
            Vfloat l3d;
            Vfloat prob;
            Vfloat subleading_pt;
            Vfloat dimuon_eta;
            int bestvtxidx;
        };

        //Denominator definition (all OS dimuons with a common vertex)
        auto getSimpleDimuons(Vint Muon_ch, Vfloat Muon_pt, Vfloat Muon_eta, Vfloat Muon_phi, Vfloat Muon_phiCorr,
        Vfloat Muon_dxyCorr, Vfloat Muon_dxye, Vint Muon_bestAssocSVIdx,
        float PV_x, float PV_y, float PV_z,
        Vfloat SV_x, Vfloat SV_y, Vfloat SV_z,
        Vfloat SV_xe, Vfloat SV_ye, Vfloat SV_ze,
        Vfloat SV_lxy, Vfloat SV_l3d, Vfloat SV_prob, ROOT::RVec<uint> SV_index){

            event_vertices dimuons;
            std::vector<int> muon_sel_idx{};

            //Iterate over muons
            for(unsigned i=0; i<Muon_pt.size(); i++){
                for(unsigned j=i+1; j<Muon_pt.size(); j++){
                    bool same_vertex = (Muon_bestAssocSVIdx.at(i) == Muon_bestAssocSVIdx.at(j));
                    bool is_neutral = ((Muon_ch.at(i) + Muon_ch.at(j)) == 0); 

                    if(same_vertex && is_neutral){
                        int pass_veto = 1;
                        int common_vertex;
                        for(int vpos=0; vpos<SV_index.size(); vpos++){
                            if(SV_index.at(vpos) == Muon_bestAssocSVIdx.at(i)){
                                common_vertex = vpos;
                                break;
                            }
                        }
                        //Muon lorentz vectors
                        ROOT::Math::PtEtaPhiMVector Muon_i_corr(Muon_pt.at(i), Muon_eta.at(i), Muon_phiCorr.at(i), 0.10566);
                        ROOT::Math::PtEtaPhiMVector Muon_i_uncorr(Muon_pt.at(i), Muon_eta.at(i), Muon_phi.at(i), 0.10566);
                        ROOT::Math::PtEtaPhiMVector Muon_j_corr(Muon_pt.at(j), Muon_eta.at(j), Muon_phiCorr.at(j), 0.10566);
                        ROOT::Math::PtEtaPhiMVector Muon_j_uncorr(Muon_pt.at(j), Muon_eta.at(j), Muon_phi.at(j), 0.10566);

                        //SV vector
                        ROOT::Math::XYZVector SV_vec(SV_x.at(common_vertex) - PV_x, SV_y.at(common_vertex) - PV_y, SV_z.at(common_vertex) - PV_z);
                        //Position vector
                        TVector3 SV_pos(SV_x.at(common_vertex), SV_y.at(common_vertex), SV_z.at(common_vertex));

                        //Dimuon vector
                        ROOT::Math::PtEtaPhiMVector dimuon_corr = Muon_i_corr + Muon_j_corr;
                        ROOT::Math::PtEtaPhiMVector dimuon_uncorr = Muon_i_uncorr + Muon_j_uncorr;

                        dimuons.mu1idx.push_back(i);
                        dimuons.mu2idx.push_back(j);
                        dimuons.mu3idx.push_back(-1);
                        dimuons.mu4idx.push_back(-1);
                        dimuons.passveto.push_back(0);
                        dimuons.excesshits.push_back(0);
                        dimuons.category.push_back(0);
                        dimuons.assocSVOverlap.push_back(0);
                        dimuons.px.push_back(dimuon_corr.px());
                        dimuons.py.push_back(dimuon_corr.py());
                        dimuons.pz.push_back(dimuon_corr.pz());
                        dimuons.x.push_back(SV_x.at(common_vertex));
                        dimuons.y.push_back(SV_y.at(common_vertex));
                        dimuons.z.push_back(SV_z.at(common_vertex));
                        dimuons.mass.push_back(dimuon_corr.mass());
                        dimuons.lxy.push_back(SV_lxy.at(common_vertex));
                        dimuons.l3d.push_back(SV_l3d.at(common_vertex));
                        dimuons.prob.push_back(SV_prob.at(common_vertex));
                                  
                        float subleading_pt = (Muon_pt.at(i) > Muon_pt.at(j)) ? Muon_pt.at(j) : Muon_pt.at(i);
                        float dimuon_eta = SV_pos.Eta();
                                  
                        dimuons.subleading_pt.push_back(subleading_pt);
                        dimuons.dimuon_eta.push_back(dimuon_eta);
                        
                        

                    }  
                }        
            }
            //Iterate over dimuons to get best probability in system
            float best_prob = -1.;
            for(int i=0; i<dimuons.prob.size(); i++){
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

        muscouting_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"

        df = df.Define("MuTrigger", muscouting_trigs)
        df = df.Filter("MuTrigger == true")

        df = df.Define("nSV", "SV_x.size()").Define("nSVOverlap", "SVOverlap_x.size()")
        df = df.Define("EventDimuons", "getSimpleDimuons(Muon_ch, Muon_pt, Muon_eta, Muon_phi, Muon_phiCorr, Muon_dxyCorr, Muon_dxye, Muon_bestAssocSVIdx, PV_x, PV_y, PV_z, SV_x, SV_y, SV_z, SV_xe, SV_ye, SV_ze, SV_lxy, SV_l3d, SV_prob, SV_index)")
        df = df.Define("DimuonLxy", "EventDimuons.lxy").Define("DimuonMass", "EventDimuons.mass").Define("DimuonSubleadingPt", "EventDimuons.subleading_pt").Define("DimuonEta", "EventDimuons.dimuon_eta").Define("DimuonBestVtxIdx", "EventDimuons.bestvtxidx")


        #Define dimuon variables
        df = df.Filter("DimuonMass.size() > 0")
        df = df.Define("DimuonMassBest", "DimuonMass.at(EventDimuons.bestvtxidx)")
        df = df.Define("DimuonSubleadingPtBest", "DimuonSubleadingPt.at(EventDimuons.bestvtxidx)")
        df = df.Define("DimuonEtaBest", "DimuonEta.at(EventDimuons.bestvtxidx)")
        df = df.Define("DimuonLxyBest", "DimuonLxy.at(EventDimuons.bestvtxidx)")

        #Overall kinematics
        pt_bins_2d_kine = array('d',[0.0, 3.0, 5.0, 10.0, 15.0, 30.0, 50.0])
        eta_bins_2d_kine = array('d',[-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
        lxy_bins_2d = array('d',[0.0, 0.2, 1.0, 2.4, 3.1, 7.0, 11.0, 16.0, 70.0])

        histos = {}
        
        #Kinematics for best vertex
        histos["h_dimuon_mass"] = df.Histo1D(("h_dimuon_mass", "; Dimuon mass (GeV); Events/0.01 GeV", 15000, 0, 150), "DimuonMassBest")
        histos["h_dimuon_subleading_pt"] = df.Histo1D(("h_dimuon_subleading_pt", "; Subleading muon pT (GeV); Events/2 GeV", 100, 0, 50), "DimuonSubleadingPtBest")
        histos["h_dimuon_eta"] = df.Histo1D(("h_dimuon_eta", "; Dimuon eta; Events/0.1", 50, -2.5, 2.5), "DimuonEtaBest")
        histos["h_sublead_pT_dimuon_eta"] = df.Histo2D(("h_sublead_pT_dimuon_eta", "; Subleading muon pT (GeV); Dimuon eta; Events", len(pt_bins_2d_kine)-1, pt_bins_2d_kine, len(eta_bins_2d_kine)-1, eta_bins_2d_kine), "DimuonSubleadingPtBest", "DimuonEtaBest")

        histo_file = ROOT.TFile.Open(create_file_dir(self.output().path), "RECREATE")
        for histo in histos.values():
            histo.Write()
        histo_file.Close()







    
