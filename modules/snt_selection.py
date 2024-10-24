#This code applies SnT dimuon selections on event dimuon objects
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

        //Apply a bunch of selections
                                  
        //Selection for the muon displacement                          
        Vbool muDisplacement(event_vertices dimuons){
            Vbool pass_displacement_cut;
            unsigned n_dimuons = dimuons.mass.size();

            //Loop over dimuons and apply the displacement selection on the muons
            for(unsigned i=0; i<n_dimuons; i++){
                bool pass_cut = true;

                //Get the quantities for the calculation
                float dimuon_pt = dimuons.pt.at(i);
                float dimuon_mass = dimuons.mass.at(i);
                float dimuon_lxy = dimuons.lxy.at(i);
                float muon1_dxy = dimuons.mu1dxyCorr.at(i);
                float muon2_dxy = dimuons.mu2dxyCorr.at(i);
                float muon1_dxye = dimuons.mu1dxye.at(i);
                float muon2_dxye = dimuons.mu2dxye.at(i);

                //dxysig and relative dxysig cut
                double muon1_dxysig = TMath::Abs(muon1_dxy)/muon1_dxye;
                double muon2_dxysig = TMath::Abs(muon2_dxy)/muon2_dxye;
                double muon1_reldxysig = TMath::Abs((muon1_dxy*dimuon_pt)/(dimuon_lxy*dimuon_mass));
                double muon2_reldxysig = TMath::Abs((muon2_dxy*dimuon_pt)/(dimuon_lxy*dimuon_mass));
                                  
                if(muon1_dxysig < 2.0) pass_cut = false;
                if(muon2_dxysig < 2.0) pass_cut = false
                if(muon1_reldxysig < 0.1) pass_cut = false;
                if(muon2_reldxysig < 0.1) pass_cut = false;
                                  
                pass_displacement_cut.push_back(pass_cut);
            }
            return pass_displacement_cut;
        }
                                  
        //Material veto selection
        Vbool materialVeto(event_vertices dimuons){
            Vbool pass_material_veto;
            unsigned n_dimuons = dimuons.mass.size();
                                  
            //Loop over dimuons and apply the material veto selection
            for(unsigned i=0; i<n_dimuons; i++){
                bool pass_cut = dimuons.passMaterialVeto.at(i);
                pass_material_veto.push_back(pass_cut);
            }
            return pass_material_veto;
        }
                                  
        //Excess hits selection
        Vbool excessHits(event_vertices dimuons){
            Vbool pass_excess_hits;
            unsigned n_dimuons = dimuons.mass.size();
                                  
            //Loop over dimuons and apply the excess hits selection
            for(unsigned i=0; i<n_dimuons; i++){
                bool pass_cut = true;
                                  
                //Dependent on displacement
                dimuon_lxy = dimuons.lxy.at(i);
                if(dimuon_lxy < 11.){
                    if(dimuons.excesshits.at(i) > 0) pass_cut = false;                  
                }
                if((dimuon_lxy > 11.) && (dimuon_lxy < 16.)){
                    if(dimuons.excesshits.at(i) > 1) pass_cut = false;
                }
                if(dimuon_lxy > 16.){
                    if(dimuons.excesshits.at(i) > 2) pass_cut = false;
                }
                
                pass_excess_hits.push_back(pass_cut);
        }                       
            return pass_excess_hits;
        }
                                  
        //Pileup veto selection
        Vbool pileupVeto(event_vertices dimuons){
            Vbool pass_pileup_veto;
            unsigned n_dimuons = dimuons.mass.size();
                                  
            //Loop over dimuons and apply the pileup veto selection
            for(unsigned i=0; i<n_dimuons; i++){
                bool pass_cut = true;
                
                //Load required quantities (using uncorrected phi)
                ROOT::Math::PtEtaPhiMVector muon1(dimuons.mu1pt.at(i), dimuons.mu1eta.at(i), dimuons.mu1phi.at(i), 0.10566);
                ROOT::Math::PtEtaPhiMVector muon2(dimuons.mu2pt.at(i), dimuons.mu2eta.at(i), dimuons.mu2phi.at(i), 0.10566);
                //3D vector
                TVector3 muon1_vector(muon1.px(), muon1.py(), muon1.pz());
                TVector3 muon2_vector(muon2.px(), muon2.py(), muon2.pz());
                                  
                //Calculate the angle between the two muons
                double dphi = TMath::Abs(muon1_vector.DeltaPhi(muon2_vector));
                double deta = TMath::Abs(muon1.Eta() - muon2.Eta());
                
                //Pileup veto
                if(TMath::Log10(deta/dphi) > 1.25) pass_cut = false;
                
                pass_pileup_veto.push_back(pass_cut);
            }
            return pass_pileup_veto;
        }


        //Dimuon kinematics selection
        Vbool dimuonKinematics(event_vertices dimuons){
            Vbool pass_dimuon_kinematics;
            unsigned n_dimuons = dimuons.mass.size();
                                  
            //Loop over dimuons and apply the dimuon kinematics selection
            for(unsigned i=0; i<n_dimuons; i++){
                bool pass_cut = true;
                
                //Make vectors (using uncorrected phi)
                ROOT::Math::PtEtaPhiMVector muon1(dimuons.mu1pt.at(i), dimuons.mu1eta.at(i), dimuons.mu1phi.at(i), 0.10566);
                ROOT::Math::PtEtaPhiMVector muon2(dimuons.mu2pt.at(i), dimuons.mu2eta.at(i), dimuons.mu2phi.at(i), 0.10566);
                ROOT::Math::PtEtaPhiMVector dimuon = muon1 + muon2;
                ROOT::Math::XYZVector sv(dimuons.x.at(i) - dimuons.PV_x.at(i), dimuons.y.at(i) - dimuons.PV_y.at(i), dimuons.z.at(i) - dimuons.PV_z.at(i));
                
                //Get necessary 3 vectors
                TVector3 muon1_vector(muon1.px(), muon1.py(), muon1.pz());
                TVector3 muon2_vector(muon2.px(), muon2.py(), muon2.pz());
                TVector3 dimuon_vector(dimuon.px(), dimuon.py(), dimuon.pz());
                TVector3 sv_vector(sv.x(), sv.y(), sv.z());
                                  
                //Now calculate the angular quantities 
                double dphi_dimuon_sv = TMath::Abs(dimuon_vector.DeltaPhi(sv_vector));
                double dphi_dimuon = TMath::Abs(muon1_vector.DeltaPhi(muon2_vector));
                double alpha3d_dimuon_sv = TMath::Abs(dimuon_vector.Angle(sv_vector));
                double alpha3d_dimuon = TMath::Abs(muon1_vector.Angle(muon2_vector));
                                  
                //Apply the cuts
                if(dphi_dimuon_sv > (TMath::Pi()/2.)) pass_cut = false;
                if(dphi_dimuon > (0.9*TMath::Pi())) pass_cut = false;
                if(alpha3d_dimuon_sv > (TMath::Pi()/2.)) pass_cut = false;
                if(alpha3d_dimuon > (0.9*TMath::Pi())) pass_cut = false;

                pass_dimuon_kinematics.push_back(pass_cut);  
            }
            return pass_dimuon_kinematics;
        }
        """)
    
    def run(self, df):

        #Print dataset details (debugging)
        print("Dataset is data:", self.dataset.process.isData)
        print("Dataset run period:", self.dataset.runPeriod)

        #Define the triggers (year dependent)
        eg_trigs = None
        jet_trigs = None
        ortho_trigs = None
        mu_trigs = None

        if self.dataset.runPeriod == "2022":
            eg_trigs = "(L1_DoubleEG_LooseIso16_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso18_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso20_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso22_LooseIso12_er1p5==true)||(L1_SingleLooseIsoEG28er2p1==true)||(L1_SingleLooseIsoEG28er1p5==true)||(L1_SingleLooseIsoEG30er1p5==true)||(L1_SingleIsoEG28er2p1==true)||(L1_SingleIsoEG30er2p1==true)||(L1_SingleIsoEG32er2p1==true)"
            jet_trigs = "(L1_HTT200er==true)||(L1_HTT255er==true)||(L1_HTT280er==true)||(L1_HTT320er==true)||(L1_HTT360er==true)||(L1_HTT400er==true)||(L1_HTT450er==true)||(L1_ETT2000==true)||(L1_SingleJet180==true)||(L1_SingleJet200==true)||(L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5==true)"
            ortho_trigs = eg_trigs + "||" + jet_trigs
            mu_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"
        
        elif self.dataset.runPeriod == "2023":
            eg_trigs = "(DST_Run3_EG30_PFScoutingPixelTracking==true)||(DST_Run3_EG16_EG12_PFScoutingPixelTracking==true)"
            jet_trigs = "(DST_Run3_JetHT_PFScoutingPixelTracking==true)"
            ortho_trigs = eg_trigs + "||" + jet_trigs
            mu_trigs = "(Run3_DoubleMu3_PFScouting==true)"

        else:
            print("Run period not provided, assuming 2022")
            eg_trigs = "(L1_DoubleEG_LooseIso16_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso18_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso20_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso22_LooseIso12_er1p5==true)||(L1_SingleLooseIsoEG28er2p1==true)||(L1_SingleLooseIsoEG28er1p5==true)||(L1_SingleLooseIsoEG30er1p5==true)||(L1_SingleIsoEG28er2p1==true)||(L1_SingleIsoEG30er2p1==true)||(L1_SingleIsoEG32er2p1==true)"
            jet_trigs = "(L1_HTT200er==true)||(L1_HTT255er==true)||(L1_HTT280er==true)||(L1_HTT320er==true)||(L1_HTT360er==true)||(L1_HTT400er==true)||(L1_HTT450er==true)||(L1_ETT2000==true)||(L1_SingleJet180==true)||(L1_SingleJet200==true)||(L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5==true)"
            ortho_trigs = eg_trigs + "||" + jet_trigs
            mu_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"

        df = df.Define("MuTrigger", mu_trigs)
        df = df.Filter("PassHLT==1").Filter("MuTrigger==true")

        #Make the secondary vertices (only preselection applied)
        df = df.Define("EventDimuons", "getDenomDimuons(Muon_charge, Muon_pt, Muon_eta, Muon_phi, Muon_phiCorr, Muon_dxyCorr, Muon_dxye, Muon_selected, Muon_bestAssocSVOverlapIdx, Muon_bestAssocSVIdx, Muon_nhitsbeforesv, PV_x, PV_y, PV_z, SV_x, SV_y, SV_z, SV_xe, SV_ye, SV_ze, SV_minDistanceFromDet_x, SV_minDistanceFromDet_y, SV_minDistanceFromDet_z, SV_onModuleWithinUnc, SV_lxy, SV_l3d, SV_prob, SV_selected, SV_index)")
        df = df.Define("nDimuons", "EventDimuons.mass.size()")
        #Define the cuts
        df = df.Define("EventDimuonsPassDisplacement", "muDisplacement(EventDimuons)")
        df = df.Define("EventDimuonsPassMaterialVeto", "materialVeto(EventDimuons)")
        df = df.Define("EventDimuonsPassExcessHits", "excessHits(EventDimuons)")
        df = df.Define("EventDimuonsPassPileupVeto", "pileupVeto(EventDimuons)")
        df = df.Define("EventDimuonsPassKinematics", "dimuonKinematics(EventDimuons)")
        #Get the best vertex and corresponding quantities
        df = df.Filter("nDimuons>0").Define("EventDimuonsBestIdx", "EventDimuons.bestvtxidx")
        df = df.Define("EventDimuonsBestMass", "EventDimuons.mass.at(EventDimuonsBestIdx)").Define("EventDimuonsBestLxy", "EventDimuons.lxy.at(EventDimuonsBestIdx)")
        #Get flags for whether the best vertex passes each of the cuts
        df = df.Define("EventDimuonsBestPassDisplacement", "EventDimuonsPassDisplacement.at(EventDimuonsBestIdx)").Define("EventDimuonsBestPassMaterialVeto", "EventDimuonsPassMaterialVeto.at(EventDimuonsBestIdx)").Define("EventDimuonsBestPassExcessHits", "EventDimuonsPassExcessHits.at(EventDimuonsBestIdx)").Define("EventDimuonsBestPassPileupVeto", "EventDimuonsPassPileupVeto.at(EventDimuonsBestIdx)").Define("EventDimuonsBestPassKinematics", "EventDimuonsPassKinematics.at(EventDimuonsBestIdx)")
        #Get the reweighting variables
        df = df.Define("EventDimuonsBestSubleadingPt", "EventDimuons.subleading_pt.at(EventDimuonsBestIdx)").Define("EventDimuonsBestEta", "EventDimuons.dimuon_eta.at(EventDimuonsBestIdx)")
        #Define N selections
        df = df.Define("EventDimuonsBestPass", "(EventDimuonsBestPassDisplacement)&&(EventDimuonsBestPassMaterialVeto)&&(EventDimuonsBestPassExcessHits)&&(EventDimuonsBestPassPileupVeto)&&(EventDimuonsBestPassKinematics)")
        #Define N-1 selections
        df = df.Define("EventDimuonsBestPassMinusDisplacement", "(EventDimuonsBestPassMaterialVeto)&&(EventDimuonsBestPassExcessHits)&&(EventDimuonsBestPassPileupVeto)&&(EventDimuonsBestPassKinematics)")
        df = df.Define("EventDimuonsBestPassMinusMaterialVeto", "(EventDimuonsBestPassDisplacement)&&(EventDimuonsBestPassExcessHits)&&(EventDimuonsBestPassPileupVeto)&&(EventDimuonsBestPassKinematics)")
        df = df.Define("EventDimuonsBestPassMinusExcessHits", "(EventDimuonsBestPassDisplacement)&&(EventDimuonsBestPassMaterialVeto)&&(EventDimuonsBestPassPileupVeto)&&(EventDimuonsBestPassKinematics)")
        df = df.Define("EventDimuonsBestPassMinusPileupVeto", "(EventDimuonsBestPassDisplacement)&&(EventDimuonsBestPassMaterialVeto)&&(EventDimuonsBestPassExcessHits)&&(EventDimuonsBestPassKinematics)")
        df = df.Define("EventDimuonsBestPassMinusKinematics", "(EventDimuonsBestPassDisplacement)&&(EventDimuonsBestPassMaterialVeto)&&(EventDimuonsBestPassExcessHits)&&(EventDimuonsBestPassPileupVeto)")
        
        return df, ["nDimuons", "EventDimuonsBestMass", "EventDimuonsBestLxy", "EventDimuonsBestSubleadingPt", "EventDimuonsBestEta", "EventDimuonsBestPassDisplacement", "EventDimuonsBestPassMaterialVeto", "EventDimuonsBestPassExcessHits", "EventDimuonsBestPassPileupVeto", "EventDimuonsBestPassKinematics", "EventDimuonsBestPass", "EventDimuonsBestPassMinusDisplacement", "EventDimuonsBestPassMinusMaterialVeto", "EventDimuonsBestPassMinusExcessHits", "EventDimuonsBestPassMinusPileupVeto", "EventDimuonsBestPassMinusKinematics"]


def SnTMakeDimuonsRDF(*args, **kwargs):
    return lambda: SnTMakeDimuonsRDFProducer(*args, **kwargs)

