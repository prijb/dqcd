import law
import os
from array import array

from cmt.base_tasks.base import DatasetTaskWithCategory, HTCondorWorkflow, SGEWorkflow, InputData

from analysis_tools.utils import create_file_dir, import_root

class TriggerSFScoutingSnT(DatasetTaskWithCategory, law.LocalWorkflow, HTCondorWorkflow, SGEWorkflow):
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
            using Vbool = ROOT::RVec<bool>;
            using Vint = ROOT::RVec<int>;
            using Vfloat = ROOT::RVec<float>;
                                  
            float getdR(float eta1, float phi1, float eta2, float phi2) {
                float deta = eta1 - eta2;
                float dphi = phi1 - phi2;
                float pi = float(TMath::Pi());
                if(dphi > pi) dphi -= 2*pi;
                return std::sqrt(deta*deta + dphi*dphi);
            }

            struct Dimuon_trigger_system {
                float subleading_pt;
                float min_dxy;
                float dimuon_lxy;
                float dimuon_dr;
                float dimuon_mass;
                float dimuon_pt;
                float dimuon_eta;
            };

            //Get the dimuon trigger system for the event 
            Dimuon_trigger_system make_dimuon_system(Vfloat SV_lxy, Vfloat SV_x, Vfloat SV_y, Vfloat SV_z,
            Vfloat Muon_pt, Vfloat Muon_eta, Vfloat Muon_phiCorr, Vfloat Muon_dxy, unsigned int CommonVertexIdx){
                Dimuon_trigger_system event_dimuon;
                TLorentzVector muon1, muon2;
                TLorentzVector dimuon;
                TVector3 SV_vec(SV_x[0], SV_y[0], SV_z[0]);
                muon1.SetPtEtaPhiM(Muon_pt[0], Muon_eta[0], Muon_phiCorr[0], 0.1055);
                muon2.SetPtEtaPhiM(Muon_pt[1], Muon_eta[1], Muon_phiCorr[1], 0.1055);
                dimuon = muon1 + muon2;

                //Assign the dimuon trigger system variables
                event_dimuon.subleading_pt = (Muon_pt[0] > Muon_pt[1]) ? Muon_pt[1] : Muon_pt[0];
                event_dimuon.min_dxy = (std::abs(Muon_dxy[0]) < std::abs(Muon_dxy[1])) ? std::abs(Muon_dxy[0]) : std::abs(Muon_dxy[1]);
                event_dimuon.dimuon_lxy = SV_lxy[0];
                event_dimuon.dimuon_dr = getdR(Muon_eta[0], Muon_phiCorr[0], Muon_eta[1], Muon_phiCorr[1]);
                event_dimuon.dimuon_mass = dimuon.M();
                event_dimuon.dimuon_pt = dimuon.Pt();
                event_dimuon.dimuon_eta = dimuon.Eta();
                return event_dimuon;                    
            }
        """)

    def run(self):
        ROOT = import_root()
        self.add_to_root(ROOT)

        df = ROOT.RDataFrame("tout", self.input()["data"][0].path)

        eg_trigs = "(L1_DoubleEG_LooseIso16_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso18_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso20_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso22_LooseIso12_er1p5==true)||(L1_SingleLooseIsoEG28er2p1==true)||(L1_SingleLooseIsoEG28er1p5==true)||(L1_SingleLooseIsoEG30er1p5==true)||(L1_SingleIsoEG28er2p1==true)||(L1_SingleIsoEG30er2p1==true)||(L1_SingleIsoEG32er2p1==true)"
        jet_trigs = "(L1_HTT200er==true)||(L1_HTT255er==true)||(L1_HTT280er==true)||(L1_HTT320er==true)||(L1_HTT360er==true)||(L1_HTT400er==true)||(L1_HTT450er==true)||(L1_ETT2000==true)||(L1_SingleJet180==true)||(L1_SingleJet200==true)||(L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5==true)"
        ortho_trigs = eg_trigs + "||" + jet_trigs
        muscouting_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"

        #df_orthogonal = df.Define("MuTrigger", muscouting_trigs)
        
        df_orthogonal = df.Define("OrthoTrigger", ortho_trigs).Define("MuTrigger", muscouting_trigs)
        df_orthogonal = df_orthogonal.Filter("OrthoTrigger==true")

        # Denominator definition
        # Kinematic cuts on muons
        df_denominator = df_orthogonal.Define("nMuon", "Muon_pt.size()").Filter("nMuon==2").Filter("(Muon_pt[0]>3.0)&&(Muon_pt[1]>3.0)").Filter("(std::abs(Muon_eta[0])<2.4)&&(std::abs(Muon_eta[1])<2.4)").Filter("Muon_ch[0]!=Muon_ch[1]")
        df_denominator = df_denominator.Define("nSV", "SV_x.size()").Filter("nSV==1")
        # Vertex matching requirement
        df_denominator = df_denominator.Filter("Muon_bestAssocSVIdx[0]==Muon_bestAssocSVIdx[1]").Define("CommonVertexIdx", "Muon_bestAssocSVIdx[0]").Filter("CommonVertexIdx==SV_index[0]")
        # Cut on first (only) SV
        df_denominator = df_denominator.Filter("SV_chi2Ndof[0] < 5.0").Filter("(std::abs(SV_xe[0]) < 0.05) && (std::abs(SV_ye[0]) < 0.05) && (std::abs(SV_ze[0]) < 0.1)")
        df_denominator = df_denominator.Filter("SV_onModuleWithinUnc[0]==false")
        # Construct dimuon system
        df_denominator = df_denominator.Define("dimuon_system", "make_dimuon_system(SV_lxy, SV_x, SV_y, SV_z, Muon_pt, Muon_eta, Muon_phiCorr, Muon_dxy, CommonVertexIdx)")
        df_denominator = df_denominator.Define("subleading_pt", "dimuon_system.subleading_pt").Define("min_dxy", "dimuon_system.min_dxy").Define("dimuon_lxy", "dimuon_system.dimuon_lxy").Define("dimuon_dr", "dimuon_system.dimuon_dr").Define("dimuon_mass", "dimuon_system.dimuon_mass").Define("dimuon_pt", "dimuon_system.dimuon_pt").Define("dimuon_eta", "dimuon_system.dimuon_eta")

        # Choose resonance (J/Psi)
        df_denominator = df_denominator.Filter("(dimuon_mass>2.6) && (dimuon_mass<3.4)")
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
        pt_bins_2d_kine = array('d',[0.0, 3.0, 5.0, 10.0, 15.0, 30.0, 50.0])
        jpsi_pt_bins_2d_kine = array('d',[0.0, 3.0, 5.0, 10.0, 15.0, 30.0, 50.0, 80.0, 200.0])
        eta_bins_2d_kine = array('d',[-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
        lxy_bins_2d = array('d',[0.0, 0.2, 1.0, 2.4, 3.1, 7.0, 11.0, 16.0, 70.0])

        #These are kinematics for events passing the Mu trigger and SV selections (for reweighting)
        histos["h_sublead_pT_dimuon_eta"] = df_pass.Histo2D(("h_sublead_pT_dimuon_eta", "; Subleading muon pT (GeV); Dimuon eta; Events", len(pt_bins_2d_kine)-1, pt_bins_2d_kine, len(eta_bins_2d_kine)-1, eta_bins_2d_kine), "subleading_pt", "dimuon_eta")

        # J/Psi kinematics
        histos["h_sublead_pT"] = df_denominator.Histo1D(("h_sublead_pT", "; Subleading muon pT (GeV); Events/0.3 GeV", 100, 0, 30), "subleading_pt")
        histos["h_dimuon_pT"] = df_denominator.Histo1D(("h_dimuon_pT", "; Dimuon pT (GeV); Events/0.3 GeV", 100, 0, 200), "dimuon_pt")
        histos["h_dimuon_eta"] = df_denominator.Histo1D(("h_dimuon_eta", "; Dimuon eta; Events/0.1", 50, -2.5, 2.5), "dimuon_eta")
        #Check kinematics of J/Psi (pT and eta) in the lxy [3.1, 7.0] region
        histos["h_dimuon_pT_eta_lxy_4"] = df_denominator.Filter("abs(dimuon_lxy) > 3.1 && abs(dimuon_lxy) < 7.0").Histo2D(("h_dimuon_pT_eta_lxy_4", "; Dimuon pT (GeV); Dimuon eta; Events", len(jpsi_pt_bins_2d_kine)-1, jpsi_pt_bins_2d_kine, len(eta_bins_2d_kine)-1, eta_bins_2d_kine), "dimuon_pt", "dimuon_eta")
        #Draw 2D histogram of J/Psi pT vs Lxy using non uniform bins
        histos["h_dimuon_pT_lxy"] = df_denominator.Histo2D(("h_dimuon_pT_lxy", "; Dimuon pT (GeV); Dimuon Lxy (cm); Events", len(jpsi_pt_bins_2d_kine)-1, jpsi_pt_bins_2d_kine, len(lxy_bins_2d)-1, lxy_bins_2d), "dimuon_pt", "dimuon_lxy")

        histo_file = ROOT.TFile.Open(create_file_dir(self.output().path), "RECREATE")
        for histo in histos.values():
            histo.Write()
        histo_file.Close()


class TriggerSFScoutingSnTMC(TriggerSFScoutingSnT):
    def run(self):
        ROOT = import_root()
        self.add_to_root(ROOT)

        df = ROOT.RDataFrame("tout", self.input()["data"][0].path)

        eg_trigs = "(L1_DoubleEG_LooseIso16_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso18_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso20_LooseIso12_er1p5==true)||(L1_DoubleEG_LooseIso22_LooseIso12_er1p5==true)||(L1_SingleLooseIsoEG28er2p1==true)||(L1_SingleLooseIsoEG28er1p5==true)||(L1_SingleLooseIsoEG30er1p5==true)||(L1_SingleIsoEG28er2p1==true)||(L1_SingleIsoEG30er2p1==true)||(L1_SingleIsoEG32er2p1==true)"
        jet_trigs = "(L1_HTT200er==true)||(L1_HTT255er==true)||(L1_HTT280er==true)||(L1_HTT320er==true)||(L1_HTT360er==true)||(L1_HTT400er==true)||(L1_HTT450er==true)||(L1_ETT2000==true)||(L1_SingleJet180==true)||(L1_SingleJet200==true)||(L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5==true)||(L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5==true)"
        ortho_trigs = eg_trigs + "||" + jet_trigs
        muscouting_trigs = "(L1_DoubleMu_15_7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7==true)||(L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18==true)||(L1_DoubleMu4_SQ_OS_dR_Max1p2==true)||(L1_DoubleMu4p5_SQ_OS_dR_Max1p2==true)"

        df_orthogonal = df.Define("OrthoTrigger", ortho_trigs).Define("MuTrigger", muscouting_trigs)

        # Normally commented for stats
        #df_orthogonal = df_orthogonal.Filter("OrthoTrigger==true")

        # Denominator definition
        # Kinematic cuts on muons
        df_denominator = df_orthogonal.Define("nMuon", "Muon_pt.size()").Filter("nMuon==2").Filter("(Muon_pt[0]>3.0)&&(Muon_pt[1]>3.0)").Filter("(std::abs(Muon_eta[0])<2.4)&&(std::abs(Muon_eta[1])<2.4)").Filter("Muon_ch[0]!=Muon_ch[1]")
        df_denominator = df_denominator.Define("nSV", "SV_x.size()").Filter("nSV==1")
        # Vertex matching requirement
        df_denominator = df_denominator.Filter("Muon_bestAssocSVIdx[0]==Muon_bestAssocSVIdx[1]").Define("CommonVertexIdx", "Muon_bestAssocSVIdx[0]").Filter("CommonVertexIdx==SV_index[0]")
        # Cut on first (only) SV
        df_denominator = df_denominator.Filter("SV_chi2Ndof[0] < 5.0").Filter("(std::abs(SV_xe[0]) < 0.05) && (std::abs(SV_ye[0]) < 0.05) && (std::abs(SV_ze[0]) < 0.1)")
        df_denominator = df_denominator.Filter("SV_onModuleWithinUnc[0]==false")
        # Construct dimuon system
        df_denominator = df_denominator.Define("dimuon_system", "make_dimuon_system(SV_lxy, SV_x, SV_y, SV_z, Muon_pt, Muon_eta, Muon_phiCorr, Muon_dxy, CommonVertexIdx)")
        df_denominator = df_denominator.Define("subleading_pt", "dimuon_system.subleading_pt").Define("min_dxy", "dimuon_system.min_dxy").Define("dimuon_lxy", "dimuon_system.dimuon_lxy").Define("dimuon_dr", "dimuon_system.dimuon_dr").Define("dimuon_mass", "dimuon_system.dimuon_mass").Define("dimuon_pt", "dimuon_system.dimuon_pt").Define("dimuon_eta", "dimuon_system.dimuon_eta")

        # Choose resonance (J/Psi)
        df_denominator = df_denominator.Filter("(dimuon_mass>2.6) && (dimuon_mass<3.4)")
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
        pt_bins_2d_kine = array('d',[0.0, 3.0, 5.0, 10.0, 15.0, 30.0, 50.0])
        jpsi_pt_bins_2d_kine = array('d',[0.0, 3.0, 5.0, 10.0, 15.0, 30.0, 50.0, 80.0, 200.0])
        eta_bins_2d_kine = array('d',[-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
        lxy_bins_2d = array('d',[0.0, 0.2, 1.0, 2.4, 3.1, 7.0, 11.0, 16.0, 70.0])

        #These are kinematics for events passing the Mu trigger and SV selections (for reweighting)
        histos["h_sublead_pT_dimuon_eta"] = df_pass.Histo2D(("h_sublead_pT_dimuon_eta", "; Subleading muon pT (GeV); Dimuon eta; Events", len(pt_bins_2d_kine)-1, pt_bins_2d_kine, len(eta_bins_2d_kine)-1, eta_bins_2d_kine), "subleading_pt", "dimuon_eta")

        # J/Psi kinematics
        histos["h_sublead_pT"] = df_denominator.Histo1D(("h_sublead_pT", "; Subleading muon pT (GeV); Events/0.3 GeV", 100, 0, 30), "subleading_pt")
        histos["h_dimuon_pT"] = df_denominator.Histo1D(("h_dimuon_pT", "; Dimuon pT (GeV); Events/0.3 GeV", 100, 0, 200), "dimuon_pt")
        histos["h_dimuon_eta"] = df_denominator.Histo1D(("h_dimuon_eta", "; Dimuon eta; Events/0.1", 50, -2.5, 2.5), "dimuon_eta")
        #Check kinematics of J/Psi (pT and eta) in the lxy [3.1, 7.0] region
        histos["h_dimuon_pT_eta_lxy_4"] = df_denominator.Filter("abs(dimuon_lxy) > 3.1 && abs(dimuon_lxy) < 7.0").Histo2D(("h_dimuon_pT_eta_lxy_4", "; Dimuon pT (GeV); Dimuon eta; Events", len(jpsi_pt_bins_2d_kine)-1, jpsi_pt_bins_2d_kine, len(eta_bins_2d_kine)-1, eta_bins_2d_kine), "dimuon_pt", "dimuon_eta")
        #Draw 2D histogram of J/Psi pT vs Lxy using non uniform bins
        histos["h_dimuon_pT_lxy"] = df_denominator.Histo2D(("h_dimuon_pT_lxy", "; Dimuon pT (GeV); Dimuon Lxy (cm); Events", len(jpsi_pt_bins_2d_kine)-1, jpsi_pt_bins_2d_kine, len(lxy_bins_2d)-1, lxy_bins_2d), "dimuon_pt", "dimuon_lxy")

        histo_file = ROOT.TFile.Open(create_file_dir(self.output().path), "RECREATE")
        for histo in histos.values():
            histo.Write()
        histo_file.Close()




