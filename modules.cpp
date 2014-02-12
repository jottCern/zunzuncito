#include "modules.hpp"
#include "utils.hpp"
#include "event.hpp"

#include "JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"

#include <sstream>

using namespace std;

// ----------------------------------------------------------------- //
bool dumper::process(event & evt){
   std::cout << "dumper(" << prefix << "): ";
   dump(evt, nmax, evt.is_mc);
   return true;
}

// ----------------------------------------------------------------- //
jec::jec(const std::string & path_, int njetsmax_, const std::string & globaltag_mc, const std::string & globaltag_data): njetsmax(njetsmax_), path(path_), gt_data(globaltag_data), gt_mc(globaltag_mc){
   // make sure path ends with a '/':
   if(path.size() > 0){
      if(path[path.size()-1]!='/'){
         path += '/';
      }
   }
}

void jec::start_dataset(const dataset & d, TFile & outfile){
   std::vector<JetCorrectorParameters> pars;
   const std::string & globaltag = d.mc ? gt_mc : gt_data;
   pars.push_back(path + globaltag + "_L1FastJet_AK5PFchs.txt");
   pars.push_back(path + globaltag + "_L2Relative_AK5PFchs.txt");
   pars.push_back(path + globaltag + "_L3Absolute_AK5PFchs.txt");
  //  pars.push_back(path + globaltag + "_L1FastJet_AK5Calo.txt");
//    pars.push_back(path + globaltag + "_L2Relative_AK5Calo.txt");
//    pars.push_back(path + globaltag + "_L3Absolute_AK5Calo.txt");
   if(!d.mc){
      pars.push_back(path + globaltag + "_L2L3Residual_AK5PFchs.txt");
      //  pars.push_back(path + globaltag + "_L2L3Residual_AK5Calo.txt");
   }
   fjc.reset(new FactorizedJetCorrector(pars));
}

bool jec::process(event & evt){
   // correct all jets:
   int nmax = evt.NobjJet;
   if(njetsmax > 0 && njetsmax < nmax){
      nmax = njetsmax;
   }
   for(int i=0; i<nmax; ++i){
      fjc->setJetPt(evt.JetPt[i]);
      fjc->setJetEta(evt.JetEta[i]);
      fjc->setJetE(evt.JetE[i]);
      fjc->setJetA(evt.JetArea[i]);
      fjc->setRho(evt.Rho);
      float corr = fjc->getCorrection();
      evt.JetPt[i] *= corr;
      evt.JetE[i] *= corr;
   }
   evt.sort_jets_pt(nmax);
   return true;
}

jec::~jec(){}

// ----------------------------------------------------------------- //
jec_uncertainty::jec_uncertainty(const string & path_txtfile, syst_dir::e_syst_dir sdir_, int njetsmax_): sdir(sdir_), njetsmax(njetsmax_){
   jec_unc.reset(new JetCorrectionUncertainty(path_txtfile));
}

bool jec_uncertainty::process(event & evt){
   int nmax = evt.NobjJet;
   if(njetsmax > 0 && njetsmax < nmax){
      nmax = njetsmax;
   }
   for(int i=0; i < nmax; ++i){
      jec_unc->setJetEta(evt.JetEta[i]);
      jec_unc->setJetPt(evt.JetPt[i]);
      float correctionfactor;
      if(sdir == syst_dir::up){
         float unc = jec_unc->getUncertainty(1);
         correctionfactor = (1 + fabs(unc));
      } else {
         float unc = jec_unc->getUncertainty(-1);
         correctionfactor = (1 - fabs(unc));
      }
      evt.JetPt[i] *= correctionfactor;
      evt.JetE[i] *= correctionfactor;
   }
   evt.sort_jets_pt(nmax);
   return true;
}

jec_uncertainty::~jec_uncertainty(){}

// ----------------------------------------------------------------- //
bool checktrigger::process(event & evt){
   if(evt.NobjJet < 2) return false;

   float ptave = 0.5 * (evt.JetPt[0] + evt.JetPt[1]);
   if(ptave < 62) return false;  // min pt threshold in data and mc

   // keep all mc events above min pt
   if(ptave >= 62 && evt.is_mc) return true;

   // do trigger selection in data depending on dataset and ptave bin
   if(!evt.is_mc) {
      if(evt.is_jetdata || evt.is_jetmondata) {
         if(ptave >= 62   && ptave < 107 && evt.HltDiPFJetAve40) return true;
         if(ptave >= 107  && ptave < 175 && evt.HltDiPFJetAve80) return true;
         if(ptave >= 175  && ptave < 242 && evt.HltDiPFJetAve140) return true;
         if(ptave >= 242  && ptave < 310 && evt.HltDiPFJetAve200) return true;
         if(ptave >= 310  && ptave < 379 && evt.HltDiPFJetAve260) return true;
      }

      if(evt.is_jetdata || evt.is_jethtdata) {
         if(ptave >= 379 && ptave < 467 && evt.HltDiPFJetAve320) return true;
         if(ptave >= 467 && evt.HltDiPFJetAve400) return true;
      }
   }
  
   // reject all other events
   return false;
}

checktrigger::~checktrigger(){}

// ----------------------------------------------------------------- //
pureweighting::pureweighting(const std::string & path_, const std::string & true_mc_, const std::string & pu_data_): path(path_), true_mc(true_mc_), pu_data(pu_data_){
   // make sure path ends with a '/':
   if(path.size() > 0){
      if(path[path.size()-1]!='/'){
         path += '/';
      }
   }
}  

void pureweighting::start_dataset(const dataset & d, TFile & outfile){
   if(d.mc) {
      TString fileNamePU_data = path + pu_data;
      TString fileNamePU_truemc = path + true_mc + "Summer12S10CMSSW53_TrueMCPUDistributions.root";

      TFile file_true(fileNamePU_truemc, "READ");
      TH1 *h_truemc = 0;
      file_true.GetObject("pileup", h_truemc);
      if (!h_truemc) {
         std::cout << "ERROR in module PU Reweighting: Histogram 'pileup' does not exist in file for true MC distributions '" << fileNamePU_truemc << "'\n." << std::endl;
         exit(1);
      }
   
      TFile file_data_DiPFJetAve40(fileNamePU_data + "MyDataPileupHistogramHLT_DiPFJetAve40.root", "READ");
      TFile file_data_DiPFJetAve80(fileNamePU_data + "MyDataPileupHistogramHLT_DiPFJetAve80.root", "READ");
      TFile file_data_DiPFJetAve140(fileNamePU_data + "MyDataPileupHistogramHLT_DiPFJetAve140.root", "READ");
      TFile file_data_DiPFJetAve200(fileNamePU_data + "MyDataPileupHistogramHLT_DiPFJetAve200.root", "READ");
      TFile file_data_DiPFJetAve260(fileNamePU_data + "MyDataPileupHistogramHLT_DiPFJetAve260.root", "READ");
      TFile file_data_DiPFJetAve320(fileNamePU_data + "MyDataPileupHistogramHLT_DiPFJetAve320.root", "READ");
      TFile file_data_DiPFJetAve400(fileNamePU_data + "MyDataPileupHistogramHLT_DiPFJetAve400.root", "READ");

      h_DiPFJetAve40 = 0;
      h_DiPFJetAve80 = 0;
      h_DiPFJetAve140 = 0;
      h_DiPFJetAve200 = 0;
      h_DiPFJetAve260 = 0;
      h_DiPFJetAve320 = 0;
      h_DiPFJetAve400 = 0;

      h_DiPFJetAve40 = (TH1F*)file_data_DiPFJetAve40.Get("pileup"); 
      h_DiPFJetAve40->SetDirectory(0);
      h_DiPFJetAve80 = (TH1F*)file_data_DiPFJetAve80.FindObjectAny("pileup"); 
      h_DiPFJetAve80->SetDirectory(0);
      h_DiPFJetAve140 = (TH1F*)file_data_DiPFJetAve140.FindObjectAny("pileup");
      h_DiPFJetAve140->SetDirectory(0);
      h_DiPFJetAve200 = (TH1F*)file_data_DiPFJetAve200.FindObjectAny("pileup");
      h_DiPFJetAve200->SetDirectory(0);
      h_DiPFJetAve260 = (TH1F*)file_data_DiPFJetAve260.FindObjectAny("pileup");
      h_DiPFJetAve260->SetDirectory(0);
      h_DiPFJetAve320 = (TH1F*)file_data_DiPFJetAve320.FindObjectAny("pileup");
      h_DiPFJetAve320->SetDirectory(0);
      h_DiPFJetAve400 = (TH1F*)file_data_DiPFJetAve400.FindObjectAny("pileup");
      h_DiPFJetAve400->SetDirectory(0);
   
      if (!(h_DiPFJetAve40 || h_DiPFJetAve80 || h_DiPFJetAve140 ||  h_DiPFJetAve200 || h_DiPFJetAve260  || h_DiPFJetAve320  || h_DiPFJetAve400)) {
         std::cout << "ERROR in module PU Reweighting: Histogram 'pileup' does not exist in file for data distributions '" << fileNamePU_data << "'\n." << std::endl;
         exit(1);
      }

      std::cout << "Computing pu reweighting weights" << std::endl;

      // normalize histos first
      h_truemc->Scale(1.0/h_truemc->Integral() );
      h_DiPFJetAve40->Scale(1.0/h_DiPFJetAve40->Integral() );
      h_DiPFJetAve80->Scale(1.0/h_DiPFJetAve80->Integral() );
      h_DiPFJetAve140->Scale(1.0/h_DiPFJetAve140->Integral() );
      h_DiPFJetAve200->Scale(1.0/h_DiPFJetAve200->Integral() );
      h_DiPFJetAve260->Scale(1.0/h_DiPFJetAve260->Integral() );
      h_DiPFJetAve320->Scale(1.0/h_DiPFJetAve320->Integral() );
      h_DiPFJetAve400->Scale(1.0/h_DiPFJetAve400->Integral() );
    
      //  MC * data/MC = data, so the weights are data/MC:
      h_DiPFJetAve40->Divide(h_truemc);
      h_DiPFJetAve80->Divide(h_truemc);
      h_DiPFJetAve140->Divide(h_truemc);
      h_DiPFJetAve200->Divide(h_truemc);
      h_DiPFJetAve260->Divide(h_truemc);
      h_DiPFJetAve320->Divide(h_truemc);
      h_DiPFJetAve400->Divide(h_truemc);
   }
}

  
bool pureweighting::process(event & evt){
   if(!evt.is_mc) return true;
   if(evt.NobjJet < 2) return false;

   //  cout << "Evt weight before pu: " << evt.Weight << endl;
 
   float ptave = 0.5 * (evt.JetPt[0] + evt.JetPt[1]);
   if(ptave < 62) return false;
   if(ptave >= 62   && ptave < 107) evt.Weight = evt.Weight* h_DiPFJetAve40->GetBinContent(h_DiPFJetAve40->FindBin(evt.PUMCNumTruth));
   if(ptave >= 107  && ptave < 175) evt.Weight = evt.Weight* h_DiPFJetAve80->GetBinContent(h_DiPFJetAve80->FindBin(evt.PUMCNumTruth));
   if(ptave >= 175  && ptave < 242) evt.Weight = evt.Weight* h_DiPFJetAve140->GetBinContent(h_DiPFJetAve140->FindBin(evt.PUMCNumTruth));
   if(ptave >= 242  && ptave < 310) evt.Weight = evt.Weight* h_DiPFJetAve200->GetBinContent(h_DiPFJetAve200->FindBin(evt.PUMCNumTruth));
   if(ptave >= 310  && ptave < 379) evt.Weight = evt.Weight* h_DiPFJetAve260->GetBinContent(h_DiPFJetAve260->FindBin(evt.PUMCNumTruth));
   if(ptave >= 379  && ptave < 467) evt.Weight = evt.Weight* h_DiPFJetAve320->GetBinContent(h_DiPFJetAve320->FindBin(evt.PUMCNumTruth));
   if(ptave >= 467) evt.Weight = evt.Weight* h_DiPFJetAve400->GetBinContent(h_DiPFJetAve400->FindBin(evt.PUMCNumTruth));

   //  cout << "Evt weight after pu: " << evt.Weight << endl;
   
   return true;
}
   

pureweighting::~pureweighting(){}

// ----------------------------------------------------------------- //
bool eventfilter::process(event & evt){
   // 1. we must have a vertex with rho < 2cm, abs(z) < 24cm and ndof >= 4, so reject the event otherwise:
   if(evt.VtxNDof < 4.0f or std::abs(evt.VtxPosZ) > 24.0f or hypotf(evt.VtxPosX, evt.VtxPosY) > 2.0f) return false;
      
   // if event passed this: keep it:
   return true;
}

eventfilter::~eventfilter(){}

// ----------------------------------------------------------------- //
bool eventcuts::process(event & evt){
   // 1. there must be at least 2 jets
   if(evt.NobjJet < 2) return false;

   // 2. (corrected) jets relevant for the analysis must fulfill loose JetID
   if(!(evt.JetIDLoose[0] && evt.JetIDLoose[1])) return false;
   if(evt.NobjJet > 2 && !evt.JetIDLoose[2]) return false;

   // 3. select DiJet-like structure
   float deltaPhi =  fabs( evt.JetPhi[0] - evt.JetPhi[1] );
   if( deltaPhi > TMath::Pi() ) deltaPhi = 2* TMath::Pi() - deltaPhi;
   if( deltaPhi <= 2.7 ) return false;
   
   // if event passed this: keep it:
   return true;
}

eventcuts::~eventcuts(){}

// ----------------------------------------------------------------- //
bool reasonablemceventchecker::process(event & evt){
   if(!evt.is_mc) return true;

   double res1 = evt.JetPt[0]/evt.GenJetPt[0]; 
   double dr1 = evt.JetGenJetDeltaR[0]; 
   double res2 = evt.JetPt[1]/evt.GenJetPt[1];
   double dr2 = evt.JetGenJetDeltaR[1]; 

   //cut on pthat introduced according to JES-mail by Mikko 11 Aug 2011
   //NOW: Introduce additional deltaR-matching cut...
   if(!( res1 > 0.2 && res1 < 2.0 && res2 > 0.2 && res2 < 2.0
         && (evt.JetPt[0] < 2.0 * evt.GenEvtScale)  
         && (evt.JetPt[1] < 2.0 * evt.GenEvtScale)
         && (dr1 < 0.25) && (dr2 < 0.25)) ) return false;
     
   // if event passed this: keep it:
   return true;
}

reasonablemceventchecker::~reasonablemceventchecker(){}

// ----------------------------------------------------------------- //
void alphareweighting::start_dataset(const dataset & d, TFile & outfile)
{
   f = new TF1("func", "0.5*[0]*(TMath::Erf([1]*x-[2])+1)", 0, 1.0);
 //   f->SetParameter(0, 1.15);
//    f->SetParameter(1, 15);
//    f->SetParameter(2, 0.1);

   f->SetParameter(0, 1.09);
   f->SetParameter(1, 13.5);
   f->SetParameter(2, 0.02);
}

bool alphareweighting::process(event & evt){
   if(!evt.is_mc) return true;

   float alpha = 0.0;
   // if there is a third jet --> third jet pt should not be too small
   if(evt.NobjJet > 2 && evt.JetPt[2] > 10.0f) alpha = evt.JetPt[2] / (0.5f * (evt.JetPt[0] + evt.JetPt[1]));

   if( alpha > 0. && alpha < 0.3) evt.Weight = evt.Weight * (f->Eval(alpha));  
  
   return true;
}

alphareweighting::~alphareweighting(){}

// ----------------------------------------------------------------- //
smearmc::smearmc(int njetsmax_): njetsmax(njetsmax_){}

bool smearmc::process(event & evt){
   if(!evt.is_mc) return true;

   float diff = 0; //difference between corr reco pt and genpt
   float smearfactor = 1; 

   // smear nmax jets:
   int nmax = evt.NobjJet;
   if(njetsmax > 0 && njetsmax < nmax){
      nmax = njetsmax;
   }
   for(int i = 0; i < nmax; ++i){
      if( evt.JetGenJetDeltaR[i] > 0.25 || evt.GenJetPt[i] < 10.) continue;
      diff = evt.JetPt[i] - evt.GenJetPt[i];
      smearfactor = getsmearfactor(abs(evt.JetEta[i]));

    //   cout << "jet pt: " << evt.JetPt[i] << endl;
//       cout << "jet eta: " << evt.JetEta[i] << endl;
//       cout << "gen jet pt: " << evt.GenJetPt[i] << endl;
//       cout << "diff: " << diff << endl;
//       cout << "smearfactor: " << smearfactor << endl;

      evt.JetPt[i] = evt.GenJetPt[i] + smearfactor * diff;

      //   cout << "corrected jet pt: " << evt.JetPt[i] << endl;
     
   }
   evt.sort_jets_pt(nmax);
   return true;
}

float smearmc::getsmearfactor(float eta)
{
   const float eta_high[] = {0.5f, 1.1f, 1.7f, 2.3f, 5.2f}; // eta bins defining smearfactors
   const int neta = sizeof(eta_high) / sizeof(float);
   // const float smearval[] = {1.080, 1.103, 1.124, 1.222, 1.206};
   const float smearval[] = {1.1, 1.1, 1.1, 1.1, 1.1};
   //const float smearval[] = {1, 1, 1, 1, 1};

   for(int ieta = 0; ieta < neta; ++ieta){
      float eta_hi = eta_high[ieta];
      float eta_lo = ieta == 0 ? 0.0f : eta_high[ieta-1];
      if(eta >= eta_lo && eta < eta_hi) return smearval[ieta];
   }

   return 1;
}

smearmc::~smearmc(){}

// ----------------------------------------------------------------- //
controlplots::controlplots(const std::string & suffix_, const std::string & dir_): suffix(suffix_), dir(dir_){
}  

void controlplots::start_dataset(const dataset & d, TFile & outfile){
   // create all histograms in the output file:
   outfile.cd();
   if(!dir.empty()){
      outfile.mkdir(dir.c_str());
      outfile.cd(dir.c_str());
   }
   // create the control histograms
   // HltDiPFJetAve40
   NVtx_HltDiPFJetAve40 = new TH1F(("NVtx_HltDiPFJetAve40_" + suffix).c_str(), ("NVtx_" + suffix).c_str(), 60, 0, 60);
   NVtx_HltDiPFJetAve40->Sumw2();
   Jet1Pt_HltDiPFJetAve40 = new TH1F(("Jet1Pt_HltDiPFJetAve40_" + suffix).c_str(), ("Jet1Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet1Pt_HltDiPFJetAve40->Sumw2();
   Jet2Pt_HltDiPFJetAve40 = new TH1F(("Jet2Pt_HltDiPFJetAve40_" + suffix).c_str(), ("Jet2Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet2Pt_HltDiPFJetAve40->Sumw2();
   Jet3Pt_HltDiPFJetAve40 = new TH1F(("Jet3Pt_HltDiPFJetAve40_" + suffix).c_str(), ("Jet3Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet3Pt_HltDiPFJetAve40->Sumw2();
   GenJet1Pt_HltDiPFJetAve40 = new TH1F(("GenJet1Pt_HltDiPFJetAve40_" + suffix).c_str(), ("GenJet1Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet1Pt_HltDiPFJetAve40->Sumw2();
   GenJet2Pt_HltDiPFJetAve40 = new TH1F(("GenJet2Pt_HltDiPFJetAve40_" + suffix).c_str(), ("GenJet2Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet2Pt_HltDiPFJetAve40->Sumw2();
   GenJet3Pt_HltDiPFJetAve40 = new TH1F(("GenJet3Pt_HltDiPFJetAve40_" + suffix).c_str(), ("GenJet3Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet3Pt_HltDiPFJetAve40->Sumw2();
   PtAve_HltDiPFJetAve40 = new TH1F(("PtAve_HltDiPFJetAve40_" + suffix).c_str(), ("PtAve_" + suffix).c_str(), 2500, 0, 2500);
   PtAve_HltDiPFJetAve40->Sumw2();
   GenPtAve_HltDiPFJetAve40 = new TH1F(("GenPtAve_HltDiPFJetAve40_" + suffix).c_str(), ("GenPtAve_" + suffix).c_str(), 2500, 0, 2500);
   GenPtAve_HltDiPFJetAve40->Sumw2();
   Alpha_HltDiPFJetAve40 = new TH1F(("Alpha_HltDiPFJetAve40_" + suffix).c_str(), ("Alpha_" + suffix).c_str(), 100, 0, 1);
   Alpha_HltDiPFJetAve40->Sumw2();
   GenAlpha_HltDiPFJetAve40 = new TH1F(("GenAlpha_HltDiPFJetAve40_" + suffix).c_str(), ("GenAlpha_" + suffix).c_str(), 100, 0, 1);
   GenAlpha_HltDiPFJetAve40->Sumw2();
   DeltaPhi_HltDiPFJetAve40 = new TH1F(("DeltaPhi_HltDiPFJetAve40_" + suffix).c_str(), ("DeltaPhi_" + suffix).c_str(), 320, 0, 3.2);
   DeltaPhi_HltDiPFJetAve40->Sumw2();

   // HltDiPFJetAve80
   NVtx_HltDiPFJetAve80 = new TH1F(("NVtx_HltDiPFJetAve80_" + suffix).c_str(), ("NVtx_" + suffix).c_str(), 60, 0, 60);
   NVtx_HltDiPFJetAve80->Sumw2();
   Jet1Pt_HltDiPFJetAve80 = new TH1F(("Jet1Pt_HltDiPFJetAve80_" + suffix).c_str(), ("Jet1Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet1Pt_HltDiPFJetAve80->Sumw2();
   Jet2Pt_HltDiPFJetAve80 = new TH1F(("Jet2Pt_HltDiPFJetAve80_" + suffix).c_str(), ("Jet2Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet2Pt_HltDiPFJetAve80->Sumw2();
   Jet3Pt_HltDiPFJetAve80 = new TH1F(("Jet3Pt_HltDiPFJetAve80_" + suffix).c_str(), ("Jet3Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet3Pt_HltDiPFJetAve80->Sumw2();
   GenJet1Pt_HltDiPFJetAve80 = new TH1F(("GenJet1Pt_HltDiPFJetAve80_" + suffix).c_str(), ("GenJet1Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet1Pt_HltDiPFJetAve80->Sumw2();
   GenJet2Pt_HltDiPFJetAve80 = new TH1F(("GenJet2Pt_HltDiPFJetAve80_" + suffix).c_str(), ("GenJet2Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet2Pt_HltDiPFJetAve80->Sumw2();
   GenJet3Pt_HltDiPFJetAve80 = new TH1F(("GenJet3Pt_HltDiPFJetAve80_" + suffix).c_str(), ("GenJet3Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet3Pt_HltDiPFJetAve80->Sumw2();
   PtAve_HltDiPFJetAve80 = new TH1F(("PtAve_HltDiPFJetAve80_" + suffix).c_str(), ("PtAve_" + suffix).c_str(), 2500, 0, 2500);
   PtAve_HltDiPFJetAve80->Sumw2();
   GenPtAve_HltDiPFJetAve80 = new TH1F(("GenPtAve_HltDiPFJetAve80_" + suffix).c_str(), ("GenPtAve_" + suffix).c_str(), 2500, 0, 2500);
   GenPtAve_HltDiPFJetAve80->Sumw2();
   Alpha_HltDiPFJetAve80 = new TH1F(("Alpha_HltDiPFJetAve80_" + suffix).c_str(), ("Alpha_" + suffix).c_str(), 100, 0, 1);
   Alpha_HltDiPFJetAve80->Sumw2();
   GenAlpha_HltDiPFJetAve80 = new TH1F(("GenAlpha_HltDiPFJetAve80_" + suffix).c_str(), ("GenAlpha_" + suffix).c_str(), 100, 0, 1);
   GenAlpha_HltDiPFJetAve80->Sumw2();
   DeltaPhi_HltDiPFJetAve80 = new TH1F(("DeltaPhi_HltDiPFJetAve80_" + suffix).c_str(), ("DeltaPhi_" + suffix).c_str(), 320, 0, 3.2);
   DeltaPhi_HltDiPFJetAve80->Sumw2();

   // HltDiPFJetAve140
   NVtx_HltDiPFJetAve140 = new TH1F(("NVtx_HltDiPFJetAve140_" + suffix).c_str(), ("NVtx_" + suffix).c_str(), 60, 0, 60);
   NVtx_HltDiPFJetAve140->Sumw2();
   Jet1Pt_HltDiPFJetAve140 = new TH1F(("Jet1Pt_HltDiPFJetAve140_" + suffix).c_str(), ("Jet1Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet1Pt_HltDiPFJetAve140->Sumw2();
   Jet2Pt_HltDiPFJetAve140 = new TH1F(("Jet2Pt_HltDiPFJetAve140_" + suffix).c_str(), ("Jet2Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet2Pt_HltDiPFJetAve140->Sumw2();
   Jet3Pt_HltDiPFJetAve140 = new TH1F(("Jet3Pt_HltDiPFJetAve140_" + suffix).c_str(), ("Jet3Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet3Pt_HltDiPFJetAve140->Sumw2();
   GenJet1Pt_HltDiPFJetAve140 = new TH1F(("GenJet1Pt_HltDiPFJetAve140_" + suffix).c_str(), ("GenJet1Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet1Pt_HltDiPFJetAve140->Sumw2();
   GenJet2Pt_HltDiPFJetAve140 = new TH1F(("GenJet2Pt_HltDiPFJetAve140_" + suffix).c_str(), ("GenJet2Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet2Pt_HltDiPFJetAve140->Sumw2();
   GenJet3Pt_HltDiPFJetAve140 = new TH1F(("GenJet3Pt_HltDiPFJetAve140_" + suffix).c_str(), ("GenJet3Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet3Pt_HltDiPFJetAve140->Sumw2();
   PtAve_HltDiPFJetAve140 = new TH1F(("PtAve_HltDiPFJetAve140_" + suffix).c_str(), ("PtAve_" + suffix).c_str(), 2500, 0, 2500);
   PtAve_HltDiPFJetAve140->Sumw2();
   GenPtAve_HltDiPFJetAve140 = new TH1F(("GenPtAve_HltDiPFJetAve140_" + suffix).c_str(), ("GenPtAve_" + suffix).c_str(), 2500, 0, 2500);
   GenPtAve_HltDiPFJetAve140->Sumw2();
   Alpha_HltDiPFJetAve140 = new TH1F(("Alpha_HltDiPFJetAve140_" + suffix).c_str(), ("Alpha_" + suffix).c_str(), 100, 0, 1);
   Alpha_HltDiPFJetAve140->Sumw2();
   GenAlpha_HltDiPFJetAve140 = new TH1F(("GenAlpha_HltDiPFJetAve140_" + suffix).c_str(), ("GenAlpha_" + suffix).c_str(), 100, 0, 1);
   GenAlpha_HltDiPFJetAve140->Sumw2();
   DeltaPhi_HltDiPFJetAve140 = new TH1F(("DeltaPhi_HltDiPFJetAve140_" + suffix).c_str(), ("DeltaPhi_" + suffix).c_str(), 320, 0, 3.2);
   DeltaPhi_HltDiPFJetAve140->Sumw2();

   // HltDiPFJetAve200
   NVtx_HltDiPFJetAve200 = new TH1F(("NVtx_HltDiPFJetAve200_" + suffix).c_str(), ("NVtx_" + suffix).c_str(), 60, 0, 60);
   NVtx_HltDiPFJetAve200->Sumw2();
   Jet1Pt_HltDiPFJetAve200 = new TH1F(("Jet1Pt_HltDiPFJetAve200_" + suffix).c_str(), ("Jet1Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet1Pt_HltDiPFJetAve200->Sumw2();
   Jet2Pt_HltDiPFJetAve200 = new TH1F(("Jet2Pt_HltDiPFJetAve200_" + suffix).c_str(), ("Jet2Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet2Pt_HltDiPFJetAve200->Sumw2();
   Jet3Pt_HltDiPFJetAve200 = new TH1F(("Jet3Pt_HltDiPFJetAve200_" + suffix).c_str(), ("Jet3Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet3Pt_HltDiPFJetAve200->Sumw2();
   GenJet1Pt_HltDiPFJetAve200 = new TH1F(("GenJet1Pt_HltDiPFJetAve200_" + suffix).c_str(), ("GenJet1Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet1Pt_HltDiPFJetAve200->Sumw2();
   GenJet2Pt_HltDiPFJetAve200 = new TH1F(("GenJet2Pt_HltDiPFJetAve200_" + suffix).c_str(), ("GenJet2Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet2Pt_HltDiPFJetAve200->Sumw2();
   GenJet3Pt_HltDiPFJetAve200 = new TH1F(("GenJet3Pt_HltDiPFJetAve200_" + suffix).c_str(), ("GenJet3Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet3Pt_HltDiPFJetAve200->Sumw2();
   PtAve_HltDiPFJetAve200 = new TH1F(("PtAve_HltDiPFJetAve200_" + suffix).c_str(), ("PtAve_" + suffix).c_str(), 2500, 0, 2500);
   PtAve_HltDiPFJetAve200->Sumw2();
   GenPtAve_HltDiPFJetAve200 = new TH1F(("GenPtAve_HltDiPFJetAve200_" + suffix).c_str(), ("GenPtAve_" + suffix).c_str(), 2500, 0, 2500);
   GenPtAve_HltDiPFJetAve200->Sumw2();
   Alpha_HltDiPFJetAve200 = new TH1F(("Alpha_HltDiPFJetAve200_" + suffix).c_str(), ("Alpha_" + suffix).c_str(), 100, 0, 1);
   Alpha_HltDiPFJetAve200->Sumw2();
   GenAlpha_HltDiPFJetAve200 = new TH1F(("GenAlpha_HltDiPFJetAve200_" + suffix).c_str(), ("GenAlpha_" + suffix).c_str(), 100, 0, 1);
   GenAlpha_HltDiPFJetAve200->Sumw2();
   DeltaPhi_HltDiPFJetAve200 = new TH1F(("DeltaPhi_HltDiPFJetAve200_" + suffix).c_str(), ("DeltaPhi_" + suffix).c_str(), 320, 0, 3.2);
   DeltaPhi_HltDiPFJetAve200->Sumw2();

   // HltDiPFJetAve260
   NVtx_HltDiPFJetAve260 = new TH1F(("NVtx_HltDiPFJetAve260_" + suffix).c_str(), ("NVtx_" + suffix).c_str(), 60, 0, 60);
   NVtx_HltDiPFJetAve260->Sumw2();
   Jet1Pt_HltDiPFJetAve260 = new TH1F(("Jet1Pt_HltDiPFJetAve260_" + suffix).c_str(), ("Jet1Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet1Pt_HltDiPFJetAve260->Sumw2();
   Jet2Pt_HltDiPFJetAve260 = new TH1F(("Jet2Pt_HltDiPFJetAve260_" + suffix).c_str(), ("Jet2Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet2Pt_HltDiPFJetAve260->Sumw2();
   Jet3Pt_HltDiPFJetAve260 = new TH1F(("Jet3Pt_HltDiPFJetAve260_" + suffix).c_str(), ("Jet3Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet3Pt_HltDiPFJetAve260->Sumw2();
   GenJet1Pt_HltDiPFJetAve260 = new TH1F(("GenJet1Pt_HltDiPFJetAve260_" + suffix).c_str(), ("GenJet1Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet1Pt_HltDiPFJetAve260->Sumw2();
   GenJet2Pt_HltDiPFJetAve260 = new TH1F(("GenJet2Pt_HltDiPFJetAve260_" + suffix).c_str(), ("GenJet2Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet2Pt_HltDiPFJetAve260->Sumw2();
   GenJet3Pt_HltDiPFJetAve260 = new TH1F(("GenJet3Pt_HltDiPFJetAve260_" + suffix).c_str(), ("GenJet3Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet3Pt_HltDiPFJetAve260->Sumw2();
   PtAve_HltDiPFJetAve260 = new TH1F(("PtAve_HltDiPFJetAve260_" + suffix).c_str(), ("PtAve_" + suffix).c_str(), 2500, 0, 2500);
   PtAve_HltDiPFJetAve260->Sumw2();
   GenPtAve_HltDiPFJetAve260 = new TH1F(("GenPtAve_HltDiPFJetAve260_" + suffix).c_str(), ("GenPtAve_" + suffix).c_str(), 2500, 0, 2500);
   GenPtAve_HltDiPFJetAve260->Sumw2();
   Alpha_HltDiPFJetAve260 = new TH1F(("Alpha_HltDiPFJetAve260_" + suffix).c_str(), ("Alpha_" + suffix).c_str(), 100, 0, 1);
   Alpha_HltDiPFJetAve260->Sumw2();
   GenAlpha_HltDiPFJetAve260 = new TH1F(("GenAlpha_HltDiPFJetAve260_" + suffix).c_str(), ("GenAlpha_" + suffix).c_str(), 100, 0, 1);
   GenAlpha_HltDiPFJetAve260->Sumw2();
   DeltaPhi_HltDiPFJetAve260 = new TH1F(("DeltaPhi_HltDiPFJetAve260_" + suffix).c_str(), ("DeltaPhi_" + suffix).c_str(), 320, 0, 3.2);
   DeltaPhi_HltDiPFJetAve260->Sumw2();

   // HltDiPFJetAve320
   NVtx_HltDiPFJetAve320 = new TH1F(("NVtx_HltDiPFJetAve320_" + suffix).c_str(), ("NVtx_" + suffix).c_str(), 60, 0, 60);
   NVtx_HltDiPFJetAve320->Sumw2();
   Jet1Pt_HltDiPFJetAve320 = new TH1F(("Jet1Pt_HltDiPFJetAve320_" + suffix).c_str(), ("Jet1Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet1Pt_HltDiPFJetAve320->Sumw2();
   Jet2Pt_HltDiPFJetAve320 = new TH1F(("Jet2Pt_HltDiPFJetAve320_" + suffix).c_str(), ("Jet2Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet2Pt_HltDiPFJetAve320->Sumw2();
   Jet3Pt_HltDiPFJetAve320 = new TH1F(("Jet3Pt_HltDiPFJetAve320_" + suffix).c_str(), ("Jet3Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet3Pt_HltDiPFJetAve320->Sumw2();
   GenJet1Pt_HltDiPFJetAve320 = new TH1F(("GenJet1Pt_HltDiPFJetAve320_" + suffix).c_str(), ("GenJet1Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet1Pt_HltDiPFJetAve320->Sumw2();
   GenJet2Pt_HltDiPFJetAve320 = new TH1F(("GenJet2Pt_HltDiPFJetAve320_" + suffix).c_str(), ("GenJet2Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet2Pt_HltDiPFJetAve320->Sumw2();
   GenJet3Pt_HltDiPFJetAve320 = new TH1F(("GenJet3Pt_HltDiPFJetAve320_" + suffix).c_str(), ("GenJet3Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet3Pt_HltDiPFJetAve320->Sumw2();
   PtAve_HltDiPFJetAve320 = new TH1F(("PtAve_HltDiPFJetAve320_" + suffix).c_str(), ("PtAve_" + suffix).c_str(), 2500, 0, 2500);
   PtAve_HltDiPFJetAve320->Sumw2();
   GenPtAve_HltDiPFJetAve320 = new TH1F(("GenPtAve_HltDiPFJetAve320_" + suffix).c_str(), ("GenPtAve_" + suffix).c_str(), 2500, 0, 2500);
   GenPtAve_HltDiPFJetAve320->Sumw2();
   Alpha_HltDiPFJetAve320 = new TH1F(("Alpha_HltDiPFJetAve320_" + suffix).c_str(), ("Alpha_" + suffix).c_str(), 100, 0, 1);
   Alpha_HltDiPFJetAve320->Sumw2();
   GenAlpha_HltDiPFJetAve320 = new TH1F(("GenAlpha_HltDiPFJetAve320_" + suffix).c_str(), ("GenAlpha_" + suffix).c_str(), 100, 0, 1);
   GenAlpha_HltDiPFJetAve320->Sumw2();
   DeltaPhi_HltDiPFJetAve320 = new TH1F(("DeltaPhi_HltDiPFJetAve320_" + suffix).c_str(), ("DeltaPhi_" + suffix).c_str(), 320, 0, 3.2);
   DeltaPhi_HltDiPFJetAve320->Sumw2();

   // HltDiPFJetAve400
   NVtx_HltDiPFJetAve400 = new TH1F(("NVtx_HltDiPFJetAve400_" + suffix).c_str(), ("NVtx_" + suffix).c_str(), 60, 0, 60);
   NVtx_HltDiPFJetAve400->Sumw2();
   Jet1Pt_HltDiPFJetAve400 = new TH1F(("Jet1Pt_HltDiPFJetAve400_" + suffix).c_str(), ("Jet1Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet1Pt_HltDiPFJetAve400->Sumw2();
   Jet2Pt_HltDiPFJetAve400 = new TH1F(("Jet2Pt_HltDiPFJetAve400_" + suffix).c_str(), ("Jet2Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet2Pt_HltDiPFJetAve400->Sumw2();
   Jet3Pt_HltDiPFJetAve400 = new TH1F(("Jet3Pt_HltDiPFJetAve400_" + suffix).c_str(), ("Jet3Pt_" + suffix).c_str(), 250, 0, 2500);
   Jet3Pt_HltDiPFJetAve400->Sumw2();
   GenJet1Pt_HltDiPFJetAve400 = new TH1F(("GenJet1Pt_HltDiPFJetAve400_" + suffix).c_str(), ("GenJet1Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet1Pt_HltDiPFJetAve400->Sumw2();
   GenJet2Pt_HltDiPFJetAve400 = new TH1F(("GenJet2Pt_HltDiPFJetAve400_" + suffix).c_str(), ("GenJet2Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet2Pt_HltDiPFJetAve400->Sumw2();
   GenJet3Pt_HltDiPFJetAve400 = new TH1F(("GenJet3Pt_HltDiPFJetAve400_" + suffix).c_str(), ("GenJet3Pt_" + suffix).c_str(), 250, 0, 2500);
   GenJet3Pt_HltDiPFJetAve400->Sumw2();
   PtAve_HltDiPFJetAve400 = new TH1F(("PtAve_HltDiPFJetAve400_" + suffix).c_str(), ("PtAve_" + suffix).c_str(), 2500, 0, 2500);
   PtAve_HltDiPFJetAve400->Sumw2();
   GenPtAve_HltDiPFJetAve400 = new TH1F(("GenPtAve_HltDiPFJetAve400_" + suffix).c_str(), ("GenPtAve_" + suffix).c_str(), 2500, 0, 2500);
   GenPtAve_HltDiPFJetAve400->Sumw2();
   Alpha_HltDiPFJetAve400 = new TH1F(("Alpha_HltDiPFJetAve400_" + suffix).c_str(), ("Alpha_" + suffix).c_str(), 100, 0, 1);
   Alpha_HltDiPFJetAve400->Sumw2();
   GenAlpha_HltDiPFJetAve400 = new TH1F(("GenAlpha_HltDiPFJetAve400_" + suffix).c_str(), ("GenAlpha_" + suffix).c_str(), 100, 0, 1);
   GenAlpha_HltDiPFJetAve400->Sumw2();
   DeltaPhi_HltDiPFJetAve400 = new TH1F(("DeltaPhi_HltDiPFJetAve400_" + suffix).c_str(), ("DeltaPhi_" + suffix).c_str(), 320, 0, 3.2);
   DeltaPhi_HltDiPFJetAve400->Sumw2();
}

bool controlplots::process(event & evt){
   float ptave = 0.5 * (evt.JetPt[0] + evt.JetPt[1]);
   float genptave = 0.5 * (evt.GenJetPt[0] + evt.GenJetPt[1]);
    
   // HltDiPFJetAve40
   if(ptave >= 62 && ptave < 107) {
      NVtx_HltDiPFJetAve40->Fill(evt.VtxN, evt.Weight);
      Jet1Pt_HltDiPFJetAve40->Fill(evt.JetPt[0], evt.Weight);
      Jet2Pt_HltDiPFJetAve40->Fill(evt.JetPt[1], evt.Weight);
      GenJet1Pt_HltDiPFJetAve40->Fill(evt.GenJetPt[0], evt.Weight);
      GenJet2Pt_HltDiPFJetAve40->Fill(evt.GenJetPt[1], evt.Weight);
      PtAve_HltDiPFJetAve40->Fill(ptave, evt.Weight);
      GenPtAve_HltDiPFJetAve40->Fill(genptave, evt.Weight);
      if(evt.NobjJet > 2) {
         Jet3Pt_HltDiPFJetAve40->Fill(evt.JetPt[2], evt.Weight);
         GenJet3Pt_HltDiPFJetAve40->Fill(evt.GenJetPt[2], evt.Weight);
         Alpha_HltDiPFJetAve40->Fill(evt.JetPt[2]/ptave, evt.Weight);
         GenAlpha_HltDiPFJetAve40->Fill(evt.GenJetPt[2]/genptave, evt.Weight);
      }
      float deltaPhi =  fabs( evt.JetPhi[0] - evt.JetPhi[1] );
      if( deltaPhi > TMath::Pi() ) deltaPhi = 2* TMath::Pi() - deltaPhi;
      DeltaPhi_HltDiPFJetAve40->Fill(deltaPhi, evt.Weight);
   }
   // HltDiPFJetAve80
   if(ptave >= 107 && ptave < 175) {
      NVtx_HltDiPFJetAve80->Fill(evt.VtxN, evt.Weight);
      Jet1Pt_HltDiPFJetAve80->Fill(evt.JetPt[0], evt.Weight);
      Jet2Pt_HltDiPFJetAve80->Fill(evt.JetPt[1], evt.Weight);
      GenJet1Pt_HltDiPFJetAve80->Fill(evt.GenJetPt[0], evt.Weight);
      GenJet2Pt_HltDiPFJetAve80->Fill(evt.GenJetPt[1], evt.Weight);
      PtAve_HltDiPFJetAve80->Fill(ptave, evt.Weight);
      GenPtAve_HltDiPFJetAve80->Fill(genptave, evt.Weight);
      if(evt.NobjJet > 2) {
         Jet3Pt_HltDiPFJetAve80->Fill(evt.JetPt[2], evt.Weight);
         GenJet3Pt_HltDiPFJetAve80->Fill(evt.GenJetPt[2], evt.Weight);
         Alpha_HltDiPFJetAve80->Fill(evt.JetPt[2]/ptave, evt.Weight);
         GenAlpha_HltDiPFJetAve80->Fill(evt.GenJetPt[2]/genptave, evt.Weight);
      }
      float deltaPhi =  fabs( evt.JetPhi[0] - evt.JetPhi[1] );
      if( deltaPhi > TMath::Pi() ) deltaPhi = 2* TMath::Pi() - deltaPhi;
      DeltaPhi_HltDiPFJetAve80->Fill(deltaPhi, evt.Weight);
   }
   // HltDiPFJetAve140
   if(ptave >= 175 && ptave < 242) {
      NVtx_HltDiPFJetAve140->Fill(evt.VtxN, evt.Weight);
      Jet1Pt_HltDiPFJetAve140->Fill(evt.JetPt[0], evt.Weight);
      Jet2Pt_HltDiPFJetAve140->Fill(evt.JetPt[1], evt.Weight);
      GenJet1Pt_HltDiPFJetAve140->Fill(evt.GenJetPt[0], evt.Weight);
      GenJet2Pt_HltDiPFJetAve140->Fill(evt.GenJetPt[1], evt.Weight);
      PtAve_HltDiPFJetAve140->Fill(ptave, evt.Weight);
      GenPtAve_HltDiPFJetAve140->Fill(genptave, evt.Weight);
      if(evt.NobjJet > 2) {
         Jet3Pt_HltDiPFJetAve140->Fill(evt.JetPt[2], evt.Weight);
         GenJet3Pt_HltDiPFJetAve140->Fill(evt.GenJetPt[2], evt.Weight);
         Alpha_HltDiPFJetAve140->Fill(evt.JetPt[2]/ptave, evt.Weight);
         GenAlpha_HltDiPFJetAve140->Fill(evt.GenJetPt[2]/genptave, evt.Weight);
      }
      float deltaPhi =  fabs( evt.JetPhi[0] - evt.JetPhi[1] );
      if( deltaPhi > TMath::Pi() ) deltaPhi = 2* TMath::Pi() - deltaPhi;
      DeltaPhi_HltDiPFJetAve140->Fill(deltaPhi, evt.Weight);
   }
   // HltDiPFJetAve200
   if(ptave >= 242 && ptave < 310) {
      NVtx_HltDiPFJetAve200->Fill(evt.VtxN, evt.Weight);
      Jet1Pt_HltDiPFJetAve200->Fill(evt.JetPt[0], evt.Weight);
      Jet2Pt_HltDiPFJetAve200->Fill(evt.JetPt[1], evt.Weight);
      GenJet1Pt_HltDiPFJetAve200->Fill(evt.GenJetPt[0], evt.Weight);
      GenJet2Pt_HltDiPFJetAve200->Fill(evt.GenJetPt[1], evt.Weight);
      PtAve_HltDiPFJetAve200->Fill(ptave, evt.Weight);
      GenPtAve_HltDiPFJetAve200->Fill(genptave, evt.Weight);
      if(evt.NobjJet > 2) {
         Jet3Pt_HltDiPFJetAve200->Fill(evt.JetPt[2], evt.Weight);
         GenJet3Pt_HltDiPFJetAve200->Fill(evt.GenJetPt[2], evt.Weight);
         Alpha_HltDiPFJetAve200->Fill(evt.JetPt[2]/ptave, evt.Weight);
         GenAlpha_HltDiPFJetAve200->Fill(evt.GenJetPt[2]/genptave, evt.Weight);
      }
      float deltaPhi =  fabs( evt.JetPhi[0] - evt.JetPhi[1] );
      if( deltaPhi > TMath::Pi() ) deltaPhi = 2* TMath::Pi() - deltaPhi;
      DeltaPhi_HltDiPFJetAve200->Fill(deltaPhi, evt.Weight);
   }
   // HltDiPFJetAve260
   if(ptave >= 310  && ptave < 379) {
      NVtx_HltDiPFJetAve260->Fill(evt.VtxN, evt.Weight);
      Jet1Pt_HltDiPFJetAve260->Fill(evt.JetPt[0], evt.Weight);
      Jet2Pt_HltDiPFJetAve260->Fill(evt.JetPt[1], evt.Weight);
      GenJet1Pt_HltDiPFJetAve260->Fill(evt.GenJetPt[0], evt.Weight);
      GenJet2Pt_HltDiPFJetAve260->Fill(evt.GenJetPt[1], evt.Weight);
      PtAve_HltDiPFJetAve260->Fill(ptave, evt.Weight);
      GenPtAve_HltDiPFJetAve260->Fill(genptave, evt.Weight);
      if(evt.NobjJet > 2) {
         Jet3Pt_HltDiPFJetAve260->Fill(evt.JetPt[2], evt.Weight);
         GenJet3Pt_HltDiPFJetAve260->Fill(evt.GenJetPt[2], evt.Weight);
         Alpha_HltDiPFJetAve260->Fill(evt.JetPt[2]/ptave, evt.Weight);
         GenAlpha_HltDiPFJetAve260->Fill(evt.GenJetPt[2]/genptave, evt.Weight);
      }
      float deltaPhi =  fabs( evt.JetPhi[0] - evt.JetPhi[1] );
      if( deltaPhi > TMath::Pi() ) deltaPhi = 2* TMath::Pi() - deltaPhi;
      DeltaPhi_HltDiPFJetAve260->Fill(deltaPhi, evt.Weight);
   }
   // HltDiPFJetAve320
   if(ptave >= 379 && ptave < 467) {
      NVtx_HltDiPFJetAve320->Fill(evt.VtxN, evt.Weight);
      Jet1Pt_HltDiPFJetAve320->Fill(evt.JetPt[0], evt.Weight);
      Jet2Pt_HltDiPFJetAve320->Fill(evt.JetPt[1], evt.Weight);
      GenJet1Pt_HltDiPFJetAve320->Fill(evt.GenJetPt[0], evt.Weight);
      GenJet2Pt_HltDiPFJetAve320->Fill(evt.GenJetPt[1], evt.Weight);
      PtAve_HltDiPFJetAve320->Fill(ptave, evt.Weight);
      GenPtAve_HltDiPFJetAve320->Fill(genptave, evt.Weight);
      if(evt.NobjJet > 2) {
         Jet3Pt_HltDiPFJetAve320->Fill(evt.JetPt[2], evt.Weight);
         GenJet3Pt_HltDiPFJetAve320->Fill(evt.GenJetPt[2], evt.Weight);
         Alpha_HltDiPFJetAve320->Fill(evt.JetPt[2]/ptave, evt.Weight);
         GenAlpha_HltDiPFJetAve320->Fill(evt.GenJetPt[2]/genptave, evt.Weight);
      }
      float deltaPhi =  fabs( evt.JetPhi[0] - evt.JetPhi[1] );
      if( deltaPhi > TMath::Pi() ) deltaPhi = 2* TMath::Pi() - deltaPhi;
      DeltaPhi_HltDiPFJetAve320->Fill(deltaPhi, evt.Weight);
   }
   // HltDiPFJetAve400
   if(ptave >= 467) {
      NVtx_HltDiPFJetAve400->Fill(evt.VtxN, evt.Weight);
      Jet1Pt_HltDiPFJetAve400->Fill(evt.JetPt[0], evt.Weight);
      Jet2Pt_HltDiPFJetAve400->Fill(evt.JetPt[1], evt.Weight);
      GenJet1Pt_HltDiPFJetAve400->Fill(evt.GenJetPt[0], evt.Weight);
      GenJet2Pt_HltDiPFJetAve400->Fill(evt.GenJetPt[1], evt.Weight);
      PtAve_HltDiPFJetAve400->Fill(ptave, evt.Weight);
      GenPtAve_HltDiPFJetAve400->Fill(genptave, evt.Weight);
      if(evt.NobjJet > 2) {
         Jet3Pt_HltDiPFJetAve400->Fill(evt.JetPt[2], evt.Weight);
         GenJet3Pt_HltDiPFJetAve400->Fill(evt.GenJetPt[2], evt.Weight);
         Alpha_HltDiPFJetAve400->Fill(evt.JetPt[2]/ptave, evt.Weight);
         GenAlpha_HltDiPFJetAve400->Fill(evt.GenJetPt[2]/genptave, evt.Weight);
      }
      float deltaPhi =  fabs( evt.JetPhi[0] - evt.JetPhi[1] );
      if( deltaPhi > TMath::Pi() ) deltaPhi = 2* TMath::Pi() - deltaPhi;
      DeltaPhi_HltDiPFJetAve400->Fill(deltaPhi, evt.Weight);
   }

   return true;
}

controlplots::~controlplots(){}


// ----------------------------------------------------------------- //
int standard_ptave_binning::nbins() const {
   return 13;
}

int standard_ptave_binning::ibin(const event & evt) const{
   if(evt.NobjJet < 2) return -1;
   float ptave = 0.5 * (evt.JetPt[0] + evt.JetPt[1]);
  
   // define same pt bins for data and mc
   if(ptave < 62) return -1;
   if(ptave >= 62    && ptave < 107)  return 0;  
   if(ptave >= 107   && ptave < 175)  return 1; 
   if(ptave >= 175   && ptave < 205)  return 2; 
   if(ptave >= 205   && ptave < 242)  return 3; 
   if(ptave >= 242   && ptave < 270)  return 4; 
   if(ptave >= 270   && ptave < 310)  return 5; 
   if(ptave >= 310   && ptave < 335)  return 6; 
   if(ptave >= 335   && ptave < 379)  return 7; 
   if(ptave >= 379   && ptave < 410)  return 8;  
   if(ptave >= 410   && ptave < 467)  return 9; 
   if(ptave >= 467   && ptave < 600)  return 10;  
   if(ptave >= 600   && ptave < 1000) return 11;  
   if(ptave >= 1000  ) return 12;   

   return -1;
    
}

standard_ptave_binning::~standard_ptave_binning(){}

// ----------------------------------------------------------------- //
asymm_histos::asymm_histos(const boost::shared_ptr<pt_binning> & binning, const std::string & dir_): dir(dir_){
   ptbinning = binning;
}

int asymm_histos::nbins_eta() const{
   return 5;
}

int asymm_histos::ibin_eta(const event & evt) const{
   if(evt.NobjJet < 2) return -1; 
   float eta0 = abs(evt.JetEta[0]);
   float eta1 = abs(evt.JetEta[1]);
   const float eta_high[] = {0.5f, 1.1f, 1.7f, 2.3f, 5.2f};
   const int neta = sizeof(eta_high) / sizeof(float);
   for(int ieta = 0; ieta < neta; ++ieta){
      // both the leading and the sub-leading jet must be within the same eta bin extending from eta_lo to eta_hi:
      float eta_hi = eta_high[ieta];
      float eta_lo = ieta == 0 ? 0.0f : eta_high[ieta-1];
      if(eta0 >= eta_lo and eta0 < eta_hi and eta1 >= eta_lo and eta1 < eta_hi) return ieta;
   }
   return -1;
}

int asymm_histos::nbins_geneta() const{
   return 5;
}

int asymm_histos::ibin_geneta(const event & evt) const{
   if(evt.NobjJet < 2) return -1; 
   float eta0 = abs(evt.GenJetEta[0]);
   float eta1 = abs(evt.GenJetEta[1]);
   const float eta_high[] = {0.5f, 1.1f, 1.7f, 2.3f, 5.2f};
   const int neta = sizeof(eta_high) / sizeof(float);
   for(int ieta = 0; ieta < neta; ++ieta){
      // both the leading and the sub-leading jet must be within the same eta bin extending from eta_lo to eta_hi:
      float eta_hi = eta_high[ieta];
      float eta_lo = ieta == 0 ? 0.0f : eta_high[ieta-1];
      if(eta0 >= eta_lo and eta0 < eta_hi and eta1 >= eta_lo and eta1 < eta_hi) return ieta;
   }
   return -1;
}
   
int asymm_histos::nbins_alpha() const{
   return 8;
}

int asymm_histos::ibin_alpha(const event & evt) const{
   float alpha = 0.0;
   // if there is a third jet --> third jet pt should not be too small
   if(evt.NobjJet > 2) {
      if(evt.JetPt[2] < 10.0f) return -1;
      alpha = evt.JetPt[2] / (0.5f * (evt.JetPt[0] + evt.JetPt[1]));
   }
   const float alpha_min[] = {0.0f, 0.1f, 0.125f, 0.15f, 0.175f, 0.20f, 0.225f, 0.25f}; // lower bin borders of alpha binning
   const int nalpha = sizeof(alpha_min) / sizeof(float);
   // find the bin index in alpha for this event. The last bin extends to +infinity implicitly
   for(int ialpha=0; ialpha < nalpha; ++ialpha){
      if(alpha >= alpha_min[ialpha] and (ialpha + 1 == nalpha or alpha < alpha_min[ialpha + 1])){
         //   cout << "alpha: " << alpha << endl;
         //   cout << "ialpha: " << ialpha << endl;
         return ialpha;
      }
   }
   return -1;
}

int asymm_histos::nbins_genalpha() const{
   return 8;
}

int asymm_histos::ibin_genalpha(const event & evt) const{
   // leading reco jets should match leading gen jets
   if( evt.JetGenJetDeltaR[0] > 0.25 || evt.JetGenJetDeltaR[1] > 0.25 ) {
      //   cout << "leading gen jets do not match leading reco jets" << endl;
      return -1;
   }

   // find 3rd gen-jet
   float pt1gen = evt.GenJetPt[0];
   float pt2gen = evt.GenJetPt[1];
   if(evt.GenJetPt[0] < pt2gen) {
      pt1gen = evt.GenJetPt[1];
      pt2gen = evt.GenJetPt[0];
   }
   float pt3gen = 0.;
   for (int i = 2; i < 50; i++) {
      if(evt.GenJetPt[i] > pt2gen && !(evt.GenJetPt[i] == pt1gen)) {
         return -1;
      }
      if(evt.GenJetPt[i] < pt2gen && evt.GenJetPt[i] > pt3gen) {
         pt3gen = evt.GenJetPt[i];
      }
   }

   if(pt3gen < 10.) return -1;
   float genalpha = pt3gen / (0.5f * (evt.GenJetPt[0] + evt.GenJetPt[1]));   
 
   const float alpha_min[] = {0.0f, 0.1f, 0.125f, 0.15f, 0.175f, 0.20f, 0.225f, 0.25f}; // lower bin borders of alpha binning
   const int nalpha = sizeof(alpha_min) / sizeof(float);
   // find the bin index in alpha for this event. The last bin extends to +infinity implicitly
   for(int ialpha=0; ialpha < nalpha; ++ialpha){
      if(genalpha >= alpha_min[ialpha] and (ialpha + 1 == nalpha or genalpha < alpha_min[ialpha + 1])){
      //    cout << "genalpha: " << genalpha << endl;
//          cout << "ialpha: " << ialpha << endl;
         return ialpha;
      }
   }
   return -1;
}

int asymm_histos::nbins_genpt() const{
   return 13;
}

int asymm_histos::ibin_genpt(const event & evt) const{
   if(evt.NobjJet < 2) return -1;

   // leading reco jets should match leading gen jets
   if( evt.JetGenJetDeltaR[0] > 0.25 || evt.JetGenJetDeltaR[1] > 0.25 ) {
      //   cout << "leading gen jets do not match leading reco jets" << endl;
      return -1;
   }

   float genptave = 0.5 * (evt.GenJetPt[0] + evt.GenJetPt[1]);
  
   // define same pt bins for data and mc
   if(genptave < 62) return -1;
   if(genptave >= 62    && genptave < 107)  return 0;  
   if(genptave >= 107   && genptave < 175)  return 1; 
   if(genptave >= 175   && genptave < 205)  return 2; 
   if(genptave >= 205   && genptave < 242)  return 3; 
   if(genptave >= 242   && genptave < 270)  return 4; 
   if(genptave >= 270   && genptave < 310)  return 5; 
   if(genptave >= 310   && genptave < 335)  return 6; 
   if(genptave >= 335   && genptave < 379)  return 7; 
   if(genptave >= 379   && genptave < 410)  return 8;  
   if(genptave >= 410   && genptave < 467)  return 9; 
   if(genptave >= 467   && genptave < 600)  return 10;  
   if(genptave >= 600   && genptave < 1000) return 11;  
   if(genptave >= 1000  ) return 12;   

   return -1;
}

void asymm_histos::start_dataset(const dataset & d, TFile & outfile){
   // clear the histograms from previous dataset:
   histos_asymm.clear();
   histos_genasymm.clear();
   histos_alphaspectrum.clear();
   // create all histograms in the output file:
   outfile.cd();
   if(!dir.empty()){
      outfile.mkdir(dir.c_str());
      outfile.cd(dir.c_str());
   }
   // create the (gen)asymmetry histograms: histograms from 0 to 1 with 1000 bins:
   const int n_pt = ptbinning->nbins();
   const int n_alpha = nbins_alpha();
   const int n_genalpha = nbins_genalpha();
   const int n_eta = nbins_eta();
   for(int ipt=0; ipt<n_pt; ++ipt){
      for(int ieta=0; ieta < n_eta; ++ieta){
         stringstream ss2;
         ss2 << "AlphaSpectrum_Pt" << ipt << "_eta" << ieta;
         TH1F *tmp_alpha = new TH1F(ss2.str().c_str(), ss2.str().c_str(), 1000, 0.0, 1.0);
         tmp_alpha->Sumw2();
         histos_alphaspectrum.push_back(tmp_alpha);

         for(int ialpha=0; ialpha < n_alpha; ++ialpha){
            stringstream ss;
            ss << "Pt" << ipt << "_eta" << ieta << "_alpha" << ialpha;
            TH1F *tmp = new TH1F(ss.str().c_str(), ss.str().c_str(), 1000, 0.0, 1.0);
            tmp->Sumw2();
            histos_asymm.push_back(tmp);
         }

         if(d.mc){
            for(int igenalpha=0; igenalpha < n_genalpha; ++igenalpha){
               stringstream ss;
               ss << "GenPt" << ipt << "_geneta" << ieta << "_genalpha" << igenalpha;
               TH1F *tmp = new TH1F(ss.str().c_str(), ss.str().c_str(), 1000, 0.0, 1.0);
               tmp->Sumw2();
               histos_genasymm.push_back(tmp);
            }
         }
      }
   }
}

bool asymm_histos::process(event & evt){
   // if(evt.NobjJet < 3) return false;
   int ptbin = ptbinning->ibin(evt);
   if(ptbin < 0) return false;
   assert(ptbin < ptbinning->nbins());
    
   int etabin = ibin_eta(evt);
   if(etabin < 0) return false;
   assert(etabin < nbins_eta());
    
   int alphabin = ibin_alpha(evt);
   if(alphabin < 0) return false;
   assert(alphabin < nbins_alpha());

   int genalphabin = ibin_genalpha(evt);
   assert(genalphabin < nbins_genalpha());

   int genptbin = ibin_genpt(evt);
   assert(genptbin < nbins_genpt());

   int genetabin = ibin_geneta(evt);
   assert(genetabin < nbins_geneta());
    
   const int n_eta = nbins_eta();
   const int n_alpha = nbins_alpha();
   const int n_genalpha = nbins_genalpha();
   const int n_geneta = nbins_geneta();

   // fill histos with alpha spectrum
   float alpha = 0.0;
   // if there is a third jet --> third jet pt should not be too small
   if(evt.NobjJet > 2 && evt.JetPt[2] > 10.0f) {
      alpha = evt.JetPt[2] / (0.5f * (evt.JetPt[0] + evt.JetPt[1]));
   }
   int histo_index_alpha = (ptbin * n_eta) + etabin;
   histos_alphaspectrum[histo_index_alpha]->Fill(alpha, evt.Weight);
 
   // fill asymmetry histos inclusive in alpha
   float asymmetry = (evt.JetPt[0] - evt.JetPt[1]) / (evt.JetPt[0] + evt.JetPt[1]);
   // asymmetry defined that way should always be positive, if the pt-sorting was right:
   assert(asymmetry >= 0.0f);
   for(int i = 0; i < n_alpha-alphabin; i++) {
      int histo_index = ptbin * (n_eta * n_alpha) + etabin * n_alpha + (i+alphabin);
      assert(histo_index < int(histos_asymm.size()));
      histos_asymm[histo_index]->Fill(asymmetry, evt.Weight);
   }

   // fill genasymmetry histos inclusive in genalpha for PLI (only for MC)
   if(evt.is_mc && !(genalphabin < 0 || genptbin < 0 || genetabin < 0) ) {
      float genasymmetry = TMath::Abs((evt.GenJetPt[0] - evt.GenJetPt[1]) / (evt.GenJetPt[0] + evt.GenJetPt[1]));
      // genasymmetry defined that way should always be positive
      assert(genasymmetry >= 0.0f);
      for(int i = 0; i < n_genalpha-genalphabin; i++) {
         int histo_index = genptbin * (n_geneta * n_genalpha) + genetabin * n_genalpha + (i+genalphabin);
         assert(histo_index < int(histos_genasymm.size()));
         histos_genasymm[histo_index]->Fill(genasymmetry, evt.Weight);
      }
   }

   // fill asymmetry histos exclusive in alpha
   //  int histo_index = ptbin * (n_eta * n_alpha) + etabin * n_alpha + alphabin;
   //  assert(histo_index < int(histos.size()));
   //  float asymmetry = (evt.JetPt[0] - evt.JetPt[1]) / (evt.JetPt[0] + evt.JetPt[1]);
   //  // asymmetry defined that way should always be positive, if the pt-sorting was right:
   //  assert(asymmetry >= 0.0f);
   //  histos[histo_index]->Fill(asymmetry, evt.Weight);
   return true;
}

// ----------------------------------------------------------------- //
response_histos::response_histos(const std::string & dir_): dir(dir_){
}

int response_histos::nbins_etagen() const{
   return 5;
}

int response_histos::ibin_etagen(const float & eta) const{
   const float eta_high[] = {0.5f, 1.1f, 1.7f, 2.3f, 5.2f};
   const int neta = sizeof(eta_high) / sizeof(float);
   for(int ieta = 0; ieta < neta; ++ieta){
      float eta_hi = eta_high[ieta];
      float eta_lo = ieta == 0 ? 0.0f : eta_high[ieta-1];
      if(eta >= eta_lo and eta < eta_hi) return ieta;
   }
   return -1;
}

int response_histos::nbins_ptgen() const{
   return 13;
}

int response_histos::ibin_ptgen(const float & ptgen) const{
   // define ptgen bins for mc
   if(ptgen < 62) return -1;
   if(ptgen >= 62    && ptgen < 107)  return 0;  
   if(ptgen >= 107   && ptgen < 175)  return 1; 
   if(ptgen >= 175   && ptgen < 205)  return 2; 
   if(ptgen >= 205   && ptgen < 242)  return 3; 
   if(ptgen >= 242   && ptgen < 270)  return 4; 
   if(ptgen >= 270   && ptgen < 310)  return 5; 
   if(ptgen >= 310   && ptgen < 335)  return 6; 
   if(ptgen >= 335   && ptgen < 379)  return 7; 
   if(ptgen >= 379   && ptgen < 410)  return 8;  
   if(ptgen >= 410   && ptgen < 467)  return 9; 
   if(ptgen >= 467   && ptgen < 600)  return 10;  
   if(ptgen >= 600   && ptgen < 1000) return 11;  
   if(ptgen >= 1000  ) return 12;  
   
   return -1; 
}
   
void response_histos::start_dataset(const dataset & d, TFile & outfile){
   // clear the histograms from previous dataset:
   histos_response.clear();
   // create all histograms in the output file:
   outfile.cd();
   if(!dir.empty()){
      outfile.mkdir(dir.c_str());
      outfile.cd(dir.c_str());
   }
   // create the response histograms: histograms from 0 to 2 with 2000 bins:
   const int n_pt = nbins_ptgen();
   const int n_eta = nbins_etagen();
   for(int ipt=0; ipt<n_pt; ++ipt){
      for(int ieta=0; ieta < n_eta; ++ieta){
         stringstream ss;
         ss << "Response_Pt" << ipt << "_eta" << ieta;
         TH1F *tmp = new TH1F(ss.str().c_str(), ss.str().c_str(), 2000, 0.0, 2.0);
         tmp->Sumw2();
         histos_response.push_back(tmp);
      }
   }
   ResponseCorr_eta0 = new TH2F("ResponseCorr_eta0", "Response", 2000, 0.0, 2.0, 2000, 0.0, 2.0 );
   ResponseCorr_eta1 = new TH2F("ResponseCorr_eta1", "Response", 2000, 0.0, 2.0, 2000, 0.0, 2.0 );
   ResponseCorr_eta2 = new TH2F("ResponseCorr_eta2", "Response", 2000, 0.0, 2.0, 2000, 0.0, 2.0 );
   ResponseCorr_eta3 = new TH2F("ResponseCorr_eta3", "Response", 2000, 0.0, 2.0, 2000, 0.0, 2.0 );
   ResponseCorr_eta4 = new TH2F("ResponseCorr_eta4", "Response", 2000, 0.0, 2.0, 2000, 0.0, 2.0 );
}

bool response_histos::process(event & evt){
   if(!evt.is_mc) return true;

   int ptbin1 = ibin_ptgen(evt.GenJetPt[0]);
   int ptbin2 = ibin_ptgen(evt.GenJetPt[1]);
   assert(ptbin1 < nbins_ptgen());
   assert(ptbin2 < nbins_ptgen());
    
   int etabin1 = ibin_etagen(abs(evt.GenJetEta[0]));
   int etabin2 = ibin_etagen(abs(evt.GenJetEta[1]));
   assert(etabin1 < nbins_etagen());
   assert(etabin2 < nbins_etagen());
        
   const int n_eta = nbins_etagen();

   // fill response histos
   // cout << "GenEta1: " << evt.GenJetEta[0] << endl;
   //cout << "ptbin1: " << ptbin1 << "  n_eta: " << n_eta << "  etabin1: " << etabin1 << endl;
   int histo_index1 = (ptbin1 * n_eta) + etabin1;
   int histo_index2 = (ptbin2 * n_eta) + etabin2;
   //cout << "histos_index1 : " << histo_index1 << endl;
   assert(histo_index1 < int(histos_response.size()));
   assert(histo_index2 < int(histos_response.size()));
   float res1 = evt.JetPt[0] / evt.GenJetPt[0];
   float res2 = evt.JetPt[1] / evt.GenJetPt[1];
   if((evt.JetGenJetDeltaR[0] < 0.25 && !(ptbin1 < 0) && !(etabin1 < 0)) &&
      (evt.JetGenJetDeltaR[1] < 0.25 && !(ptbin2 < 0) && !(etabin2 < 0))) {
      histos_response[histo_index1]->Fill(res1, evt.Weight);
      histos_response[histo_index2]->Fill(res2, evt.Weight);

      if(etabin1 == 0) ResponseCorr_eta0->Fill(res1, res2, evt.Weight);
      if(etabin1 == 1) ResponseCorr_eta1->Fill(res1, res2, evt.Weight);
      if(etabin1 == 2) ResponseCorr_eta2->Fill(res1, res2, evt.Weight);
      if(etabin1 == 3) ResponseCorr_eta3->Fill(res1, res2, evt.Weight);
      if(etabin1 == 4) ResponseCorr_eta4->Fill(res1, res2, evt.Weight);
   }
   return true;
}
// ----------------------------------------------------------------- //
