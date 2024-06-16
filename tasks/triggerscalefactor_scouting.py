import law
import os
from array import array

from cmt.base_tasks.base import DatasetTaskWithCategory, HTCondorWorkflow, SGEWorkflow, InputData

from analysis_tools.utils import create_file_dir, import_root

class TriggerSFScouting(DatasetTaskWithCategory, law.LocalWorkflow, HTCondorWorkflow, SGEWorkflow):
    def create_branch_map(self):
        return len(self.dataset.get_files(
            os.path.expandvars("$CMT_TMP_DIR/%s/" % self.config_name), add_prefix=False,
            check_empty=True))

    def workflow_requires(self):
        return {"data": InputData.req(self)}

    def requires(self):
        return {"data": InputData.req(self, file_index=self.branch)}

    def output(self):
        return self.local_target(f"data_{self.addendum}{self.branch}.root")

    def output(self):
        return self.local_target("histos_%s.root" % self.branch)

    def add_to_root(self, root):
        root.gInterpreter.Declare("""
            #include "TLorentzVector.h"
            #include "TVector3.h"
            float getdR(float eta1, float phi1, float eta2, float phi2) {
                float deta = eta1 - eta2;
                float dphi = phi1 - phi2;
                float pi = float(TMath::Pi());
                if(dphi > pi) dphi -= 2*pi;
                return std::sqrt(deta*deta + dphi*dphi);
            }
                                
            struct Dimuon_system {
                float subleading_pt;
                float min_dxy;
                float dimuon_lxy;
                float dimuon_dr;
                float dimuon_mass;
                float dimuon_pt;
                float dimuon_eta;
            };

            //Event, by construction, has only one of these SVs, so return the first one        
            int get_sv_idx(ROOT::RVec<int> SV_nMuonMatched){
                for (int i=0; i<SV_nMuonMatched.size(); i++) {
                    if (SV_nMuonMatched[i]==2) return i;
                }
                return -1;
            }

            //Gets the dimuon system for the event          
            Dimuon_system make_dimuon_system(int SV_chosenidx, ROOT::RVec<float> SV_mass, ROOT::RVec<float> SV_lxy, ROOT::RVec<int> SV_mu1idx, ROOT::RVec<int> SV_mu2idx, ROOT::RVec<float> SV_x, ROOT::RVec<float> SV_y, ROOT::RVec<float> SV_z, 
            ROOT::RVec<float> Muon_pt, ROOT::RVec<float> Muon_eta, ROOT::RVec<float> Muon_phi, ROOT::RVec<float> Muon_dxy) {
                Dimuon_system dimuon;
                TLorentzVector mu1, mu2;
                TVector3 SV_vec(SV_x[SV_chosenidx], SV_y[SV_chosenidx], SV_z[SV_chosenidx]);
                mu1.SetPtEtaPhiM(Muon_pt[SV_mu1idx[SV_chosenidx]], Muon_eta[SV_mu1idx[SV_chosenidx]], Muon_phi[SV_mu1idx[SV_chosenidx]], 0.1055);
                mu2.SetPtEtaPhiM(Muon_pt[SV_mu2idx[SV_chosenidx]], Muon_eta[SV_mu2idx[SV_chosenidx]], Muon_phi[SV_mu2idx[SV_chosenidx]], 0.1055);
                TLorentzVector dimuon_vector = mu1 + mu2;
                float dR = getdR(Muon_eta[SV_mu1idx[SV_chosenidx]], Muon_phi[SV_mu1idx[SV_chosenidx]], Muon_eta[SV_mu2idx[SV_chosenidx]], Muon_phi[SV_mu2idx[SV_chosenidx]]);
                dimuon.subleading_pt = std::min(Muon_pt[SV_mu1idx[SV_chosenidx]], Muon_pt[SV_mu2idx[SV_chosenidx]]);
                dimuon.min_dxy = std::min(std::abs(Muon_dxy[SV_mu1idx[SV_chosenidx]]), std::abs(Muon_dxy[SV_mu2idx[SV_chosenidx]]));
                dimuon.dimuon_lxy = TMath::Abs(SV_lxy[SV_chosenidx]);
                dimuon.dimuon_dr = dR;
                dimuon.dimuon_mass = SV_mass[SV_chosenidx];
                dimuon.dimuon_pt = dimuon_vector.Pt();
                //Note: Eta here is position and not momentum eta which is extracted from SV_x, SV_y and SV_z
                dimuon.dimuon_eta = SV_vec.Eta();
                return dimuon;
            }
        """)

    def run(self):
        ROOT = import_root()
        self.add_to_root(ROOT)

        df = ROOT.RDataFrame("Events", self.input()["data"][0].path)

        eg_trigs = "(L1_DoubleEG_LooseIso16_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso18_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso20_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso22_LooseIso12_er1p5==true)||(L1_SingleLooseIsoEG28er2p1==true)||(L1_SingleLooseIsoEG28er1p5==true)||(L1_SingleLooseIsoEG30er1p5==true)||(L1_SingleIsoEG28er2p1==true)||(L1_SingleIsoEG30er2p1==true)||(L1_SingleIsoEG32er2p1==true)"
        jet_trigs = "(L1_HTT200er==true)||(L1_HTT255er==true)||(L1_HTT280er==true)||(L1_HTT320er==true)||(L1_HTT360er==true)||(L1_HTT400er==true)||(L1_HTT450er==true)||(L1_ETT2000==true)||(L1_SingleJet180==true)||(L1_SingleJet200==true)||(L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5==true)"
        ortho_trigs = eg_trigs + "||" + jet_trigs
        muscouting_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"

        #Define orthogonal dataframe
        df_orthogonal = df.Define("OrthoTrigger", ortho_trigs).Define("MuTrigger", muscouting_trigs)
        df_orthogonal = df_orthogonal.Filter("OrthoTrigger==true")

        #Now get the denominator
        #Muon and SV selection
        df_denominator = df_orthogonal.Filter("nMuon==2").Filter("(Muon_pt[0]>3.0)&&(Muon_pt[1]>3.0)").Filter("(std::abs(Muon_eta[0])<2.4)&&(std::abs(Muon_eta[1])<2.4)").Filter("Muon_charge[0]!=Muon_charge[1]")
        df_denominator = df_denominator.Filter("ROOT::VecOps::Sum(SV_nMuonMatched==2)==1")
        #Get the index of the SV that has nMuonMatched==2
        df_denominator = df_denominator.Define("SV_chosenidx", "get_sv_idx(SV_nMuonMatched)")
        df_denominator = df_denominator.Filter("std::abs(SV_xError[SV_chosenidx])<0.05").Filter("std::abs(SV_yError[SV_chosenidx])<0.05").Filter("std::abs(SV_zError[SV_chosenidx])<0.1").Filter("(SV_chi2[SV_chosenidx]/SV_ndof[SV_chosenidx])<5")
        #Set up dimuon quantities
        df_denominator = df_denominator.Filter("(SV_mass[SV_chosenidx]>2.6)&&(SV_mass[SV_chosenidx]<3.4)")
        df_denominator = df_denominator.Define("Dimuon_system", "make_dimuon_system(SV_chosenidx, SV_mass, SV_dxy, SV_mu1idx, SV_mu2idx, SV_x, SV_y, SV_z, Muon_pt, Muon_eta, Muon_phi, Muon_dxy)")
        df_denominator = df_denominator.Define("subleading_pt", "Dimuon_system.subleading_pt").Define("min_dxy", "Dimuon_system.min_dxy").Define("dimuon_lxy", "Dimuon_system.dimuon_lxy").Define("dimuon_dr", "Dimuon_system.dimuon_dr").Define("dimuon_mass", "Dimuon_system.dimuon_mass").Define("dimuon_pt", "Dimuon_system.dimuon_pt").Define("dimuon_eta", "Dimuon_system.dimuon_eta")

        df_pass = df_denominator.Filter("MuTrigger==true")

        lxy_bins = [
            "abs(dimuon_lxy) > 0.0 && abs(dimuon_lxy) < 0.2",
            "abs(dimuon_lxy) > 0.2 && abs(dimuon_lxy) < 1.0",
            "abs(dimuon_lxy) > 1.0 && abs(dimuon_lxy) < 2.4",
            "abs(dimuon_lxy) > 2.4 && abs(dimuon_lxy) < 3.1",
            "abs(dimuon_lxy) > 3.1 && abs(dimuon_lxy) < 7.0",
            "abs(dimuon_lxy) > 7.0 && abs(dimuon_lxy) < 11.0",
            "abs(dimuon_lxy) > 11.0 && abs(dimuon_lxy) < 16.0",
            "abs(dimuon_lxy) > 16.0 && abs(dimuon_lxy) < 70.0",
        ]

        pt_bins = [
            "subleading_pt > 3.0 && subleading_pt < 5.0",
            "subleading_pt > 5.0 && subleading_pt < 10.0",
            "subleading_pt > 10.0 && subleading_pt < 15.0",
            "subleading_pt > 15.0 && subleading_pt < 30.0",
        ]

        dr_bins = [
            "dimuon_dr > 0.0 && dimuon_dr < 0.2",
            "dimuon_dr > 0.2 && dimuon_dr < 0.4",
            "dimuon_dr > 0.4 && dimuon_dr < 0.6",
            "dimuon_dr > 0.6 && dimuon_dr < 1.2",
            "dimuon_dr > 1.2 && dimuon_dr < 2.0",
        ]

        histos = {}

        #2D efficiencies
        #pT-Lxy
        for lxy_bin, lxy_filter in enumerate(lxy_bins):
            for pt_bin, pt_filter in enumerate(pt_bins):
                h_lxy_pT_total = df_denominator.Filter(lxy_filter).Filter(pt_filter).Histo1D(("h_lxy_%s_pT_%s_total"%(lxy_bin, pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")
                h_lxy_pT_pass = df_pass.Filter(lxy_filter).Filter(pt_filter).Histo1D(("h_lxy_%s_pT_%s_pass"%(lxy_bin, pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")

                h_lxy_pT_fail = h_lxy_pT_total.GetPtr() - h_lxy_pT_pass.GetPtr()
                h_lxy_pT_fail.SetName("h_lxy_%s_pT_%s_fail"% (lxy_bin, pt_bin)) 

                histos["h_lxy_%s_pT_%s_pass"%(lxy_bin, pt_bin)] = h_lxy_pT_pass
                histos["h_lxy_%s_pT_%s_fail"%(lxy_bin, pt_bin)] = h_lxy_pT_fail

        #pT-dR
        for dr_bin, dr_filter in enumerate(dr_bins):
            for pt_bin, pt_filter in enumerate(pt_bins):
                h_dr_pT_total = df_denominator.Filter(dr_filter).Filter(pt_filter).Histo1D(("h_dr_%s_pT_%s_total"%(dr_bin, pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")
                h_dr_pT_pass = df_pass.Filter(dr_filter).Filter(pt_filter).Histo1D(("h_dr_%s_pT_%s_pass"%(dr_bin, pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")

                h_dr_pT_fail = h_dr_pT_total.GetPtr() - h_dr_pT_pass.GetPtr()
                h_dr_pT_fail.SetName("h_dr_%s_pT_%s_fail"% (dr_bin, pt_bin))

                histos["h_dr_%s_pT_%s_pass"%(dr_bin, pt_bin)] = h_dr_pT_pass
                histos["h_dr_%s_pT_%s_fail"%(dr_bin, pt_bin)] = h_dr_pT_fail
        
        #1D efficiencies
        #pT
        for pt_bin, pt_filter in enumerate(pt_bins):
            h_lxy_pT_total = df_denominator.Filter(pt_filter).Histo1D(("h_lxy_100_pT_%s_total"%(pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")
            h_lxy_pT_pass = df_pass.Filter(pt_filter).Histo1D(("h_lxy_100_pT_%s_pass"%(pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")

            h_lxy_pT_fail = h_lxy_pT_total.GetPtr() - h_lxy_pT_pass.GetPtr()
            h_lxy_pT_fail.SetName("h_lxy_100_pT_%s_fail"% (pt_bin))

            histos["h_lxy_100_pT_%s_pass"%(pt_bin)] = h_lxy_pT_pass
            histos["h_lxy_100_pT_%s_fail"%(pt_bin)] = h_lxy_pT_fail
        
        #Lxy
        for lxy_bin, lxy_filter in enumerate(lxy_bins):
            h_lxy_pT_total = df_denominator.Filter(lxy_filter).Histo1D(("h_lxy_%s_pT_100_total"%(lxy_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")
            h_lxy_pT_pass = df_pass.Filter(lxy_filter).Histo1D(("h_lxy_%s_pT_100_pass"%(lxy_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")

            h_lxy_pT_fail = h_lxy_pT_total.GetPtr() - h_lxy_pT_pass.GetPtr()
            h_lxy_pT_fail.SetName("h_lxy_%s_pT_100_fail"% (lxy_bin))

            histos["h_lxy_%s_pT_100_pass"%(lxy_bin)] = h_lxy_pT_pass
            histos["h_lxy_%s_pT_100_fail"%(lxy_bin)] = h_lxy_pT_fail

        #Overall kinematics
        histos["h_sublead_pT"] = df_denominator.Histo1D(("h_sublead_pT", "; Subleading muon pT (GeV); Events/0.33 GeV", 100, 0, 30), "subleading_pt")
        histos["h_dimuon_pT"] = df_denominator.Histo1D(("h_dimuon_pT", "; Dimuon pT (GeV); Events/0.3 GeV", 100, 0, 200), "dimuon_pt")
        histos["h_dimuon_eta"] = df_denominator.Histo1D(("h_dimuon_eta", "; Dimuon eta; Events/0.1", 50, -2.5, 2.5), "dimuon_eta")
        #Check kinematics of J/Psi (pT and eta) in the lxy [3.1, 7.0] region
        histos["h_dimuon_pT_eta_lxy_4"] = df_denominator.Filter("abs(dimuon_lxy) > 3.1 && abs(dimuon_lxy) < 7.0").Histo2D(("h_dimuon_pT_eta_lxy_4", "; Dimuon pT (GeV); Dimuon eta; Events", 100, 0, 30, 50, -2.5, 2.5), "dimuon_pt", "dimuon_eta")
        #Draw 2D histogram of J/Psi pT vs Lxy using non uniform bins
        pt_bins_2d = array('d',[3.0, 5.0, 10.0, 15.0, 30.0])
        lxy_bins_2d = array('d',[0.0, 0.2, 1.0, 2.4, 3.1, 7.0, 11.0, 16.0, 70.0])
        histos["h_dimuon_pT_lxy"] = df_denominator.Histo2D(("h_dimuon_pT_lxy", "; Dimuon pT (GeV); Dimuon Lxy (cm); Events", len(pt_bins_2d)-1, pt_bins_2d, len(lxy_bins_2d)-1, lxy_bins_2d), "dimuon_pt", "dimuon_lxy")


        histo_file = ROOT.TFile.Open(create_file_dir(self.output().path), "RECREATE")
        for histo in histos.values():
            histo.Write()
        histo_file.Close()

#MC task which doesn't require the orthogonal trigger and inherits from TriggerSFScouting
class TriggerSFScoutingMC(TriggerSFScouting):
    def run(self):
        ROOT = import_root()
        self.add_to_root(ROOT)

        df = ROOT.RDataFrame("Events", self.input()["data"][0].path)

        eg_trigs = "(L1_DoubleEG_LooseIso16_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso18_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso20_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso22_LooseIso12_er1p5==true)||(L1_SingleLooseIsoEG28er2p1==true)||(L1_SingleLooseIsoEG28er1p5==true)||(L1_SingleLooseIsoEG30er1p5==true)||(L1_SingleIsoEG28er2p1==true)||(L1_SingleIsoEG30er2p1==true)||(L1_SingleIsoEG32er2p1==true)"
        jet_trigs = "(L1_HTT200er==true)||(L1_HTT255er==true)||(L1_HTT280er==true)||(L1_HTT320er==true)||(L1_HTT360er==true)||(L1_HTT400er==true)||(L1_HTT450er==true)||(L1_ETT2000==true)||(L1_SingleJet180==true)||(L1_SingleJet200==true)||(L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5==true)"
        ortho_trigs = eg_trigs + "||" + jet_trigs
        muscouting_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"

        #Define orthogonal dataframe
        df_orthogonal = df.Define("OrthoTrigger", ortho_trigs).Define("MuTrigger", muscouting_trigs)

        #df_orthogonal = df_orthogonal.Filter("OrthoTrigger==true")

        #Now get the denominator
        #Muon and SV selection
        df_denominator = df_orthogonal.Filter("nMuon==2").Filter("(Muon_pt[0]>3.0)&&(Muon_pt[1]>3.0)").Filter("(std::abs(Muon_eta[0])<2.4)&&(std::abs(Muon_eta[1])<2.4)").Filter("Muon_charge[0]!=Muon_charge[1]")
        df_denominator = df_denominator.Filter("ROOT::VecOps::Sum(SV_nMuonMatched==2)==1")
        #Get the index of the SV that has nMuonMatched==2
        df_denominator = df_denominator.Define("SV_chosenidx", "get_sv_idx(SV_nMuonMatched)")
        df_denominator = df_denominator.Filter("std::abs(SV_xError[SV_chosenidx])<0.05").Filter("std::abs(SV_yError[SV_chosenidx])<0.05").Filter("std::abs(SV_zError[SV_chosenidx])<0.1").Filter("(SV_chi2[SV_chosenidx]/SV_ndof[SV_chosenidx])<5")
        #Set up dimuon quantities
        df_denominator = df_denominator.Filter("(SV_mass[SV_chosenidx]>2.6)&&(SV_mass[SV_chosenidx]<3.4)")
        df_denominator = df_denominator.Define("Dimuon_system", "make_dimuon_system(SV_chosenidx, SV_mass, SV_dxy, SV_mu1idx, SV_mu2idx, SV_x, SV_y, SV_z, Muon_pt, Muon_eta, Muon_phi, Muon_dxy)")
        df_denominator = df_denominator.Define("subleading_pt", "Dimuon_system.subleading_pt").Define("min_dxy", "Dimuon_system.min_dxy").Define("dimuon_lxy", "Dimuon_system.dimuon_lxy").Define("dimuon_dr", "Dimuon_system.dimuon_dr").Define("dimuon_mass", "Dimuon_system.dimuon_mass").Define("dimuon_pt", "Dimuon_system.dimuon_pt").Define("dimuon_eta", "Dimuon_system.dimuon_eta")

        df_pass = df_denominator.Filter("MuTrigger==true")

        lxy_bins = [
            "abs(dimuon_lxy) > 0.0 && abs(dimuon_lxy) < 0.2",
            "abs(dimuon_lxy) > 0.2 && abs(dimuon_lxy) < 1.0",
            "abs(dimuon_lxy) > 1.0 && abs(dimuon_lxy) < 2.4",
            "abs(dimuon_lxy) > 2.4 && abs(dimuon_lxy) < 3.1",
            "abs(dimuon_lxy) > 3.1 && abs(dimuon_lxy) < 7.0",
            "abs(dimuon_lxy) > 7.0 && abs(dimuon_lxy) < 11.0",
            "abs(dimuon_lxy) > 11.0 && abs(dimuon_lxy) < 16.0",
            "abs(dimuon_lxy) > 16.0 && abs(dimuon_lxy) < 70.0",
        ]

        pt_bins = [
            "subleading_pt > 3.0 && subleading_pt < 5.0",
            "subleading_pt > 5.0 && subleading_pt < 10.0",
            "subleading_pt > 10.0 && subleading_pt < 15.0",
            "subleading_pt > 15.0 && subleading_pt < 30.0",
        ]

        dr_bins = [
            "dimuon_dr > 0.0 && dimuon_dr < 0.2",
            "dimuon_dr > 0.2 && dimuon_dr < 0.4",
            "dimuon_dr > 0.4 && dimuon_dr < 0.6",
            "dimuon_dr > 0.6 && dimuon_dr < 1.2",
            "dimuon_dr > 1.2 && dimuon_dr < 2.0",
        ]

        histos = {}

        #2D efficiencies
        #pT-Lxy
        for lxy_bin, lxy_filter in enumerate(lxy_bins):
            for pt_bin, pt_filter in enumerate(pt_bins):
                h_lxy_pT_total = df_denominator.Filter(lxy_filter).Filter(pt_filter).Histo1D(("h_lxy_%s_pT_%s_total"%(lxy_bin, pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")
                h_lxy_pT_pass = df_pass.Filter(lxy_filter).Filter(pt_filter).Histo1D(("h_lxy_%s_pT_%s_pass"%(lxy_bin, pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")

                h_lxy_pT_fail = h_lxy_pT_total.GetPtr() - h_lxy_pT_pass.GetPtr()
                h_lxy_pT_fail.SetName("h_lxy_%s_pT_%s_fail"% (lxy_bin, pt_bin)) 

                histos["h_lxy_%s_pT_%s_pass"%(lxy_bin, pt_bin)] = h_lxy_pT_pass
                histos["h_lxy_%s_pT_%s_fail"%(lxy_bin, pt_bin)] = h_lxy_pT_fail

        #pT-dR
        for dr_bin, dr_filter in enumerate(dr_bins):
            for pt_bin, pt_filter in enumerate(pt_bins):
                h_dr_pT_total = df_denominator.Filter(dr_filter).Filter(pt_filter).Histo1D(("h_dr_%s_pT_%s_total"%(dr_bin, pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")
                h_dr_pT_pass = df_pass.Filter(dr_filter).Filter(pt_filter).Histo1D(("h_dr_%s_pT_%s_pass"%(dr_bin, pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")

                h_dr_pT_fail = h_dr_pT_total.GetPtr() - h_dr_pT_pass.GetPtr()
                h_dr_pT_fail.SetName("h_dr_%s_pT_%s_fail"% (dr_bin, pt_bin))

                histos["h_dr_%s_pT_%s_pass"%(dr_bin, pt_bin)] = h_dr_pT_pass
                histos["h_dr_%s_pT_%s_fail"%(dr_bin, pt_bin)] = h_dr_pT_fail
        
        #1D efficiencies
        #pT
        for pt_bin, pt_filter in enumerate(pt_bins):
            h_lxy_pT_total = df_denominator.Filter(pt_filter).Histo1D(("h_lxy_100_pT_%s_total"%(pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")
            h_lxy_pT_pass = df_pass.Filter(pt_filter).Histo1D(("h_lxy_100_pT_%s_pass"%(pt_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")

            h_lxy_pT_fail = h_lxy_pT_total.GetPtr() - h_lxy_pT_pass.GetPtr()
            h_lxy_pT_fail.SetName("h_lxy_100_pT_%s_fail"% (pt_bin))

            histos["h_lxy_100_pT_%s_pass"%(pt_bin)] = h_lxy_pT_pass
            histos["h_lxy_100_pT_%s_fail"%(pt_bin)] = h_lxy_pT_fail
        
        #Lxy
        for lxy_bin, lxy_filter in enumerate(lxy_bins):
            h_lxy_pT_total = df_denominator.Filter(lxy_filter).Histo1D(("h_lxy_%s_pT_100_total"%(lxy_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")
            h_lxy_pT_pass = df_pass.Filter(lxy_filter).Histo1D(("h_lxy_%s_pT_100_pass"%(lxy_bin), "; Dimuon mass (GeV); Events/0.027 GeV", 30, 2.6, 3.4), "dimuon_mass")

            h_lxy_pT_fail = h_lxy_pT_total.GetPtr() - h_lxy_pT_pass.GetPtr()
            h_lxy_pT_fail.SetName("h_lxy_%s_pT_100_fail"% (lxy_bin))

            histos["h_lxy_%s_pT_100_pass"%(lxy_bin)] = h_lxy_pT_pass
            histos["h_lxy_%s_pT_100_fail"%(lxy_bin)] = h_lxy_pT_fail

        #Overall kinematics
        histos["h_sublead_pT"] = df_denominator.Histo1D(("h_sublead_pT", "; Subleading muon pT (GeV); Events/0.3 GeV", 100, 0, 30), "subleading_pt")
        histos["h_dimuon_pT"] = df_denominator.Histo1D(("h_dimuon_pT", "; Dimuon pT (GeV); Events/0.3 GeV", 100, 0, 200), "dimuon_pt")
        histos["h_dimuon_eta"] = df_denominator.Histo1D(("h_dimuon_eta", "; Dimuon eta; Events/0.1", 50, -2.5, 2.5), "dimuon_eta")
        #Check kinematics of J/Psi (pT and eta) in the lxy [3.1, 7.0] region
        histos["h_dimuon_pT_eta_lxy_4"] = df_denominator.Filter("abs(dimuon_lxy) > 3.1 && abs(dimuon_lxy) < 7.0").Histo2D(("h_dimuon_pT_eta_lxy_4", "; Dimuon pT (GeV); Dimuon eta; Events", 100, 0, 30, 50, -2.5, 2.5), "dimuon_pt", "dimuon_eta")
        #Draw 2D histogram of J/Psi pT vs Lxy using non uniform bins
        pt_bins_2d = array('d',[3.0, 5.0, 10.0, 15.0, 30.0])
        lxy_bins_2d = array('d',[0.0, 0.2, 1.0, 2.4, 3.1, 7.0, 11.0, 16.0, 70.0])
        histos["h_dimuon_pT_lxy"] = df_denominator.Histo2D(("h_dimuon_pT_lxy", "; Dimuon pT (GeV); Dimuon Lxy (cm); Events", len(pt_bins_2d)-1, pt_bins_2d, len(lxy_bins_2d)-1, lxy_bins_2d), "dimuon_pt", "dimuon_lxy")

        histo_file = ROOT.TFile.Open(create_file_dir(self.output().path), "RECREATE")
        for histo in histos.values():
            histo.Write()
        histo_file.Close()
