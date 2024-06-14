import os
from analysis_tools.utils import import_root
ROOT = import_root()

class SnTMakeDimuonsRDFProducer():
    def __init__(self, *args, **kwargs):
        
        ROOT.gInterpreter.Declare("""
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

                    }  
                }        
            }
            return dimuons;
        }

        //Gets the dimuons for the event
        auto getDimuons(Vint Muon_ch, Vfloat Muon_pt, Vfloat Muon_eta, Vfloat Muon_phi, Vfloat Muon_phiCorr, 
        Vfloat Muon_dxyCorr, Vfloat Muon_dxye, Vbool Muon_selected,
        Vint Muon_bestAssocSVOverlapIdx, Vint Muon_bestAssocSVIdx, Vint Muon_nhitsbeforesv,
        float PV_x, float PV_y, float PV_z,
        Vfloat SV_x, Vfloat SV_y, Vfloat SV_z,
        Vfloat SV_xe, Vfloat SV_ye, Vfloat SV_ze,
        Vfloat SV_minDistanceFromDet_x, Vfloat SV_minDistanceFromDet_y, Vfloat SV_minDistanceFromDet_z, Vbool SV_onModuleWithinUnc,
        Vfloat SV_lxy, Vfloat SV_l3d, Vfloat SV_prob, Vbool SV_selected, ROOT::RVec<uint> SV_index){
            
            event_vertices dimuons;
            std::vector<int> muon_sel_idx{};
            
            //Iterate over muons
            for(unsigned i=0; i<Muon_pt.size(); i++){
                //if(Muon_bestAssocSVOverlapIdx.at(i) > -1) continue; //Using SnT arbitration
                if(std::find(muon_sel_idx.begin(), muon_sel_idx.end(), i) != muon_sel_idx.end()) continue; //Skip muons that have already been used for vertexing
                if(Muon_bestAssocSVIdx.at(i) < 0) continue;
                if(!Muon_selected.at(i)) continue;
                
                for(unsigned j=i+1; j<Muon_pt.size(); j++){
                    //if(Muon_bestAssocSVOverlapIdx.at(j) > -1) continue; //Using SnT arbitration
                    if(std::find(muon_sel_idx.begin(), muon_sel_idx.end(), j) != muon_sel_idx.end()) continue; //Skip muons that have already been used for vertexing
                    if(Muon_bestAssocSVIdx.at(j) < 0) continue;
                    if(!Muon_selected.at(j)) continue;
                    
                    bool same_vertex = (Muon_bestAssocSVIdx.at(i) == Muon_bestAssocSVIdx.at(j));
                    bool is_neutral = ((Muon_ch.at(i) + Muon_ch.at(j)) == 0);
                    
                    //If same vertex and neutral, evaluate vertices
                    if(same_vertex && is_neutral){
                        int pass_veto = 1;
                        
                        int common_vertex;
                        for(int vpos=0; vpos<SV_index.size(); vpos++){
                            if(SV_index.at(vpos) == Muon_bestAssocSVIdx.at(i)){
                                common_vertex = vpos;
                                break;
                            }
                        } 
                        
                        //Basic selections
                        bool material_veto = false;
                        if (SV_onModuleWithinUnc.at(common_vertex)) material_veto = true;
                        if((TMath::Abs(SV_minDistanceFromDet_x.at(common_vertex)) < 0.81) && (TMath::Abs(SV_minDistanceFromDet_y.at(common_vertex)) < 3.24) && (TMath::Abs(SV_minDistanceFromDet_z.at(common_vertex)) < 0.0145)) material_veto = true;
                            
                        if (!SV_selected.at(common_vertex)) continue;
                        if(material_veto) continue;
                            
                        //Excess hits cut
                        bool pass_excess_hits = true;
                        int excess_hits_dimuon = Muon_nhitsbeforesv.at(i) + Muon_nhitsbeforesv.at(j);
                            
                        //New version of the excess hits cut
                        if((SV_lxy.at(common_vertex) < 11.)){
                            if (excess_hits_dimuon > 0) pass_excess_hits = false;
                        }   
                        if((SV_lxy.at(common_vertex) > 11.)&&(SV_lxy.at(common_vertex) < 16.)){
                            if (excess_hits_dimuon > 1) pass_excess_hits = false;
                        }   
                        if((SV_lxy.at(common_vertex) > 16.)){
                            if (excess_hits_dimuon > 2) pass_excess_hits = false;
                        } 
                        
                        //Muon lorentz vectors
                        ROOT::Math::PtEtaPhiMVector Muon_i_corr(Muon_pt.at(i), Muon_eta.at(i), Muon_phiCorr.at(i), 0.10566);
                        ROOT::Math::PtEtaPhiMVector Muon_i_uncorr(Muon_pt.at(i), Muon_eta.at(i), Muon_phi.at(i), 0.10566);
                        ROOT::Math::PtEtaPhiMVector Muon_j_corr(Muon_pt.at(j), Muon_eta.at(j), Muon_phiCorr.at(j), 0.10566);
                        ROOT::Math::PtEtaPhiMVector Muon_j_uncorr(Muon_pt.at(j), Muon_eta.at(j), Muon_phi.at(j), 0.10566);
                        
                        //SV vector
                        ROOT::Math::XYZVector SV_vec(SV_x.at(common_vertex) - PV_x, SV_y.at(common_vertex) - PV_y, SV_z.at(common_vertex) - PV_z);
                            
                        //Dimuon vector
                        ROOT::Math::PtEtaPhiMVector dimuon_corr = Muon_i_corr + Muon_j_corr;
                        ROOT::Math::PtEtaPhiMVector dimuon_uncorr = Muon_i_uncorr + Muon_j_uncorr;
                        
                        //Derived quantities from vectors
                        double Muon_i_uncorr_mag = std::sqrt(std::pow(Muon_i_uncorr.Px(), 2) + std::pow(Muon_i_uncorr.Py(), 2) + std::pow(Muon_i_uncorr.Pz(), 2));
                        double Muon_j_uncorr_mag = std::sqrt(std::pow(Muon_j_uncorr.Px(), 2) + std::pow(Muon_j_uncorr.Py(), 2) + std::pow(Muon_j_uncorr.Pz(), 2));
                        double dimuon_corr_mag = std::sqrt(std::pow(dimuon_corr.Px(), 2) + std::pow(dimuon_corr.Py(), 2) + std::pow(dimuon_corr.Pz(), 2));
                        double dimuon_uncorr_mag = std::sqrt(std::pow(dimuon_uncorr.Px(), 2) + std::pow(dimuon_uncorr.Py(), 2) + std::pow(dimuon_uncorr.Pz(), 2));
                        double SV_vec_mag = std::sqrt(std::pow(SV_vec.X(), 2) + std::pow(SV_vec.Y(), 2) + std::pow(SV_vec.Z(), 2));
                        
                        //Now get angular variables from the derived quantities
                        TVector3 Muon_i_uncorr_3d(Muon_i_uncorr.px(), Muon_i_uncorr.py(), Muon_i_uncorr.pz());
                        TVector3 Muon_j_uncorr_3d(Muon_j_uncorr.px(), Muon_j_uncorr.py(), Muon_j_uncorr.pz());
                        TVector3 dimuon_uncorr_3d(dimuon_uncorr.px(), dimuon_uncorr.py(), dimuon_uncorr.pz());
                        TVector3 SV_vec_3d(SV_vec.x(), SV_vec.y(), SV_vec.z());
                        
                        double dphi = TMath::Abs(dimuon_uncorr_3d.DeltaPhi(SV_vec_3d));
                        double dphi_dimuon = TMath::Abs(Muon_i_uncorr_3d.DeltaPhi(Muon_j_uncorr_3d));
                        double deta_dimuon = TMath::Abs(Muon_i_uncorr.Eta() - Muon_j_uncorr.Eta());
                        double alpha3d = TMath::Abs(dimuon_uncorr_3d.Angle(SV_vec_3d));
                        double alpha3d_dimuon = TMath::Abs(Muon_i_uncorr_3d.Angle(Muon_j_uncorr_3d));

                        //Corrected muon dxy and dxysig 
                        double rel_dxy_i = TMath::Abs((Muon_dxyCorr.at(i)*dimuon_corr.pt())/(SV_lxy.at(common_vertex)*dimuon_corr.mass()));
                        double rel_dxy_j = TMath::Abs((Muon_dxyCorr.at(j)*dimuon_corr.pt())/(SV_lxy.at(common_vertex)*dimuon_corr.mass()));
                        float Muon_dxysig_i = (TMath::Abs(Muon_dxyCorr.at(i)/Muon_dxye.at(i)));
                        float Muon_dxysig_j = (TMath::Abs(Muon_dxyCorr.at(j)/Muon_dxye.at(j)));
                            
                        //All quantities which decide whether the SV passes the veto or not
                        if(Muon_dxysig_i < 2.) pass_veto = 0;
                        if(Muon_dxysig_j < 2.) pass_veto = 0;
                        if(rel_dxy_i < 0.1) pass_veto = 0;
                        if(rel_dxy_j < 0.1) pass_veto = 0;
                        if(!pass_excess_hits) pass_veto = 0;
                        if(TMath::Log10(deta_dimuon/dphi_dimuon) > 1.25) pass_veto = 0;
                        
                        //Angular cuts on SV
                        if(alpha3d_dimuon > (0.9*TMath::Pi())) continue;
                        if(dphi_dimuon > (0.9*TMath::Pi())) continue;
                        if(dphi > (TMath::Pi()/2.)) continue;
                        if(alpha3d > (TMath::Pi()/2.)) continue;

                        //Diagnostic tool to catch dimuons that can also be associated to an overlapping SV
                        int isAssocSVOverlap = 0;
                        if((Muon_bestAssocSVOverlapIdx.at(i) > -1) && (Muon_bestAssocSVOverlapIdx.at(j) > -1)) isAssocSVOverlap = 1;
                        
                        //Fill the event dimuon struct
                        muon_sel_idx.push_back(i);
                        muon_sel_idx.push_back(j);
                        
                        dimuons.mu1idx.push_back(i);
                        dimuons.mu2idx.push_back(j);
                        dimuons.mu3idx.push_back(-1);
                        dimuons.mu4idx.push_back(-1);
                        dimuons.passveto.push_back(pass_veto);
                        dimuons.excesshits.push_back(excess_hits_dimuon);
                        dimuons.category.push_back(0);
                        dimuons.assocSVOverlap.push_back(isAssocSVOverlap);
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
                    }
                }
            }
            return dimuons;                     
        }
        //Function to get best dimuon value
        auto getBestDimuon(Vfloat dimuon_val, Vfloat dimuon_prob){
            float best_val = -1;
            if(dimuon_prob.size() > 0){
                best_val = dimuon_val[ArgMax(dimuon_prob)];
            }
            return best_val;
        }                  
        """)

    def run(self, df):

        muscouting_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"

        df = df.Define("MuTrigger", muscouting_trigs)

        df = df.Filter("passHLT==1")
        df = df.Filter("MuTrigger==true")
        df = df.Define("nSV", "SV_selected.size()").Define("nSVOverlap", "SVOverlap_x.size()")
        df = df.Define("EventDimuonsDenom", "getSimpleDimuons(Muon_ch, Muon_pt, Muon_eta, Muon_phi, Muon_phiCorr, Muon_dxyCorr, Muon_dxye, Muon_bestAssocSVIdx, PV_x, PV_y, PV_z, SV_x, SV_y, SV_z, SV_xe, SV_ye, SV_ze, SV_lxy, SV_l3d, SV_prob, SV_index)")
        df = df.Define("EventDimuons", "getDimuons(Muon_ch, Muon_pt, Muon_eta, Muon_phi, Muon_phiCorr, Muon_dxyCorr, Muon_dxye, Muon_selected, Muon_bestAssocSVOverlapIdx, Muon_bestAssocSVIdx, Muon_nhitsbeforesv, PV_x, PV_y, PV_z, SV_x, SV_y, SV_z, SV_xe, SV_ye, SV_ze, SV_minDistanceFromDet_x, SV_minDistanceFromDet_y, SV_minDistanceFromDet_z, SV_onModuleWithinUnc, SV_lxy, SV_l3d, SV_prob, SV_selected, SV_index)")

        #Define plot variables
        #Denominator
        df = df.Define("DimuonMassDenom", "EventDimuonsDenom.mass").Define("DimuonProbDenom", "EventDimuonsDenom.prob").Define("DimuonLxyDenom", "EventDimuonsDenom.lxy")
        #Numerator
        df = df.Define("DimuonMass", "EventDimuons.mass").Define("DimuonPass", "EventDimuons.passveto").Define("DimuonProb", "EventDimuons.prob").Define("DimuonAssocSVOverlap", "EventDimuons.assocSVOverlap").Define("DimuonLxy", "EventDimuons.lxy")
        df = df.Define("DimuonMassPass", "DimuonMass[DimuonPass==1]").Define("DimuonProbPass", "DimuonProb[DimuonPass==1]").Define("DimuonLxyPass", "DimuonLxy[DimuonPass==1]")
        df = df.Define("DimuonMassPassBest", "getBestDimuon(DimuonMassPass, DimuonProbPass)").Define("DimuonProbPassBest", "getBestDimuon(DimuonProbPass, DimuonProbPass)").Define("DimuonLxyPassBest", "getBestDimuon(DimuonLxyPass, DimuonProbPass)")
        #df = df.Filter("DimuonMassPass.size() > 0")
        #df = df.Define("DimuonMassPassAssocSVOverlap", "DimuonMass[DimuonPass==1 && DimuonAssocSVOverlap==1]").Define("DimuonProbPassAssocSVOverlap", "DimuonProb[DimuonPass==1 && DimuonAssocSVOverlap==1]")
        #df = df.Define("DimuonMaxProbIdx", "ArgMax(DimuonProbPass)")
        #df = df.Define("DimuonMassPassBest", "DimuonMassPass[ArgMax(DimuonProbPass)]").Define("DimuonProbPassBest", "DimuonProbPass[ArgMax(DimuonProbPass)]")

        return df, ["DimuonMassDenom", "DimuonProbDenom", "DimuonLxyDenom", "DimuonMassPass", "DimuonProbPass", "DimuonLxyPass", "DimuonMassPassBest", "DimuonProbPassBest", "DimuonLxyPassBest"]
        #return df, ["DimuonMassDenom", "DimuonProbDenom", "DimuonMassPass", "DimuonProbPass", "DimuonMaxProbIdx", "DimuonMassPassBest", "DimuonProbPassBest", "DimuonMassPassAssocSVOverlap", "DimuonProbPassAssocSVOverlap"]

def SnTMakeDimuonsRDF(*args, **kwargs):
    return lambda: SnTMakeDimuonsRDFProducer(*args, **kwargs)