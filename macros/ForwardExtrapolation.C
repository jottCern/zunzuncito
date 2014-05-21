#include "/afs/desy.de/user/k/kheine/xxl-af-cms/Kalibri/scripts/tdrstyle_mod.C"
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TFitResult.h>
#include <TVirtualFitter.h>

#include "util/utils.h"
#include "util/HistOps.h"
#include "util/FileOps.h"
#include "util/LabelFactory.h"
#include "util/StyleSettings.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>


// --------------------------------- //
std::pair <float,float> determineMinMax(TGraphErrors* graph)
{
   std::pair <float,float> minMaxPair(graph->GetY()[0],graph->GetY()[0]);
   for(Int_t i=0;i<graph->GetN();i++){
      if(graph->GetY()[i]<minMaxPair.first)minMaxPair.first=graph->GetY()[i];
      if(graph->GetY()[i]>minMaxPair.second)minMaxPair.second=graph->GetY()[i];
   }
   return minMaxPair;
}

// --------------------------------- //
// any data to be used by chi2_linear for evaluating the function:
struct fit_data{
   // x, y values:
   std::vector<double> x_val, y_val;
   // variance-covariance matrix for y values:
   TMatrixD y_cov;
   // inverted cov matrix; calculated by chi2_linear "on demand".
   TMatrixD y_cov_inv;
   
   void reset(){
      x_val.clear();
      y_val.clear();
      y_cov.ResizeTo(0,0);
      y_cov_inv.ResizeTo(0,0);
   }

   void CheckPoints(){
      std::vector<int> RemovedPoints;
      TMatrixD y_cov_new;
      int j = 0;
    
      for(unsigned int i = 0; i < y_val.size(); i++) {
         //         std::cout << "i: " << i << "   j: " << j << std::endl;
         if( y_val.at(i) == 0) {
            x_val.erase(x_val.begin()+i);
            y_val.erase(y_val.begin()+i);
            RemovedPoints.push_back(j);
            i = i-1;
         }
         j++;
      }    
    //   for(int i = 0; i < x_val.size(); i++) {
//           std::cout << "x: " << x_val.at(i) << std::endl;

//       }

      //  std::cout<< "Removed Points: " << RemovedPoints.size() << std::endl;
      //  std::cout<< "Remaining Points: " << x_val.size() << std::endl;

      y_cov_new.ResizeTo(x_val.size(),x_val.size());
      for(unsigned int i=0; i < x_val.size(); i++) {
         for(unsigned int k= 0; k < x_val.size(); k++) {
            y_cov_new(i,k) = y_cov(i+RemovedPoints.size(),k+RemovedPoints.size());
         }
      }
      y_cov.ResizeTo(0,0);
      y_cov.ResizeTo(x_val.size(),x_val.size());
      y_cov = y_cov_new;  
   }
};

fit_data data;

// --------------------------------- //
// the chi^2 to minimize for fitting a linear function
//   y = p[0]*x + p[1]
// with fit parameters p[0], p[1] to data with known x and y and covariance
// matrix for y.
void chi2_linear(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* p, Int_t status){
    if(data.y_cov_inv.GetNcols()==0){
        double dummy;
        int ncols = data.y_cov.GetNcols();
        data.y_cov_inv.ResizeTo(ncols, ncols);
        data.y_cov_inv = data.y_cov.Invert(&dummy);
    }
    const size_t ndata = data.x_val.size(); // number of data points in x,y graph to fit to
    std::vector<double> delta_y(ndata);
    for(size_t i=0; i<ndata; ++i){
        delta_y[i] = data.x_val[i]*p[0] + p[1] - data.y_val[i];
    }
    // now calculate the chi2, i.e.
    //  dy^T * C^{-1} * dy
    // where C is the variance--covariance matrix and dy = (y_data - y_pred)
    // This could probably be implemented in ROOT, but it's so simple, we just do it here:
    fval = 0.0;
    for(size_t i=0; i<ndata; ++i){
        for(size_t j=0; j<ndata; ++j){
            fval += delta_y[i] * delta_y[j] * data.y_cov_inv(i,j);
        }
    }
}

// --------------------------------- //
double make_lin_fit(double & slope, double & d_slope, double & offset, double & d_offset){
    TMinuit min;
    min.SetPrintLevel(-1);
    // min.SetPrintLevel(0);
    int err = min.DefineParameter(0, "slope", slope, d_slope, 0.05, 1.0);
    assert(err==0);
    err = min.DefineParameter(1, "offset", offset, d_offset, 0.001, 0.2);
    assert(err==0);
    min.SetFCN(chi2_linear);
    min.mnmigr();
    min.GetParameter(0, slope, d_slope);
    min.GetParameter(1, offset, d_offset);

    int n = 2;
    double fval;
    double pars[2] = {slope, offset};
    chi2_linear(n, 0, fval, pars, 0);
    return fval;
}

// --------------------------------- //
float GetAsymmWidth(TH1F* htemp, double * xq_IQW, double * yq_IQW) 
{
   const int nq = 2;

   float width = 0.;
   
   if( htemp->GetEntries() > 100 ) {
      htemp->GetXaxis()->SetRange(0,-1);
      htemp->ComputeIntegral();
      htemp->GetQuantiles(nq,yq_IQW,xq_IQW);
      Int_t IQW_high_bin_i = htemp->FindBin(yq_IQW[1]);
 
      // cout << "Truncated Integral in %: " << htemp->Integral(1, IQW_high_bin_i)/htemp->Integral() << endl;

      for(int i=1; i <= IQW_high_bin_i; i++) {
         width += htemp->GetBinContent(i)* std::pow(htemp->GetBinCenter(i), 2);
      }

      width = TMath::Sqrt(1/(htemp->Integral(1, IQW_high_bin_i))*width);
   }

   return width;  
}

// --------------------------------- //
float GetTruthRes(TH1F* htemp, double trunc_val) 
{
   float width = 0.;

   float integral_tot = htemp->Integral();
   
   if( htemp->GetEntries() > 100 ) {

      Int_t MeanBin =  htemp->FindBin(htemp->GetMean());
      float integral = htemp->GetBinContent(MeanBin);

      cout << "Mean: " << htemp->GetMean() << endl;

      Int_t IQW_bin_i_low = MeanBin;
      Int_t IQW_bin_i_high = MeanBin;

      for(int i = 1; i <= htemp->GetNbinsX(); i++) {
         integral += htemp->GetBinContent(MeanBin+i);
         integral += htemp->GetBinContent(MeanBin-i);
         if(integral/integral_tot > trunc_val) {
            IQW_bin_i_low = MeanBin-i;
            IQW_bin_i_high = MeanBin+i;
            break;
         }
      }

      cout << "Truncated Integral in %: " << integral/htemp->Integral() << endl;

      htemp->GetXaxis()->SetRange(IQW_bin_i_low, IQW_bin_i_high); 
      width = htemp->GetRMS();
   }

   return width;  
}

// --------------------------------- //
float GetTruthResErr(TH1F* htemp, double trunc_val) 
{
   float width_err = 0.;

   float integral_tot = htemp->Integral();
   
   if( htemp->GetEntries() > 100 ) {

      Int_t MeanBin =  htemp->FindBin(htemp->GetMean());
      float integral = htemp->GetBinContent(MeanBin);

      // cout << "Mean: " << htemp->GetMean() << endl;

      Int_t IQW_bin_i_low = MeanBin;
      Int_t IQW_bin_i_high = MeanBin;

      for(int i = 1; i <= htemp->GetNbinsX(); i++) {
         integral += htemp->GetBinContent(MeanBin+i);
         integral += htemp->GetBinContent(MeanBin-i);
         if(integral/integral_tot > trunc_val) {
            IQW_bin_i_low = MeanBin-i;
            IQW_bin_i_high = MeanBin+i;
            break;
         }
      }

      // cout << "Truncated Integral in %: " << integral/htemp->Integral() << endl;

      htemp->GetXaxis()->SetRange(IQW_bin_i_low, IQW_bin_i_high); 
      width_err = htemp->GetRMSError();
   }

   return width_err;  
}

// --------------------------------- //
void ForwardExtrapolation()
{
   setTDRStyle();
   gROOT->ForceStyle();

   // input files
   // --------------------- //
   TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_ForwardExtension_final_v4.root", "READ");
   // TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v4.root", "READ");
   // TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_ForwardExtensionNextToCentral_final_v3.root", "READ");
   // TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v3.root", "READ");
   // TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_ForwardExtensionSecondNextToCentral_final_v1.root", "READ");
   // TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_ForwardExtensionSecondNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v1.root", "READ");
  
   // --------------------- //
   TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD_ForwardExtension_final_v4.root", "READ");
   // TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v4.root", "READ");
   // TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD_ForwardExtensionNextToCentral_final_v3.root", "READ");
   // TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v3.root", "READ");
   // TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD_ForwardExtensionSecondNextToCentral_final_v1.root", "READ");
   // TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD_ForwardExtensionSecondNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v1.root", "READ");
  
   // --------------------- //
   TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD_ForwardExtension_final_v4.root", "READ");
   // TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v4.root", "READ");
   // TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD_ForwardExtensionNextToCentral_final_v3.root", "READ");
   // TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v3.root", "READ");
   // TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD_ForwardExtensionSecondNextToCentral_final_v1.root", "READ");
   // TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD_ForwardExtensionSecondNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v1.root", "READ");
  
   // --------------------- //
   TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_v4.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_MBXS73500_v4.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_JECup_v4.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_JECdown_v4.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_v4.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_ReweightAlphaSpectrum_v4.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_ReweightGluonSplitting_v4.root", "READ");

   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionNextToCentral_v3.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionNextToCentral_MBXS73500_v3.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionNextToCentral_JECup_v3.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionNextToCentral_JECdown_v3.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_v3.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionNextToCentral_ReweightAlphaSpectrum_v3.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionNextToCentral_ReweightGluonSplitting_v3.root", "READ");

   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionSecondNextToCentral_v1.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionSecondNextToCentral_MBXS73500_v1.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionSecondNextToCentral_JECup_v1.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionSecondNextToCentral_JECdown_v1.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionSecondNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_v1.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionSecondNextToCentral_ReweightAlphaSpectrum_v1.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionSecondNextToCentral_ReweightGluonSplitting_v1.root", "READ");

   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_ForwardExtension_v4.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_ForwardExtensionNextToCentral_v3.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_ForwardExtensionSecondNextToCentral_v1.root", "READ");
  
   // --------------------- //
 
   // --------------------- //
   // TString suffix = "_test";

   TString suffix = "_ForwardExtension_v4";
   // TString suffix = "_ForwardExtension_MBXS73500_v4";
   // TString suffix = "_ForwardExtension_JECup_v4";
   // TString suffix = "_ForwardExtension_JECdown_v4";
   // TString suffix = "_ForwardExtension_PLIup_v4";
   // TString suffix = "_ForwardExtension_PLIdown_v4";
   // TString suffix = "_ForwardExtension_95Truncation_v4";
   // TString suffix = "_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_v4";
   // TString suffix = "_ForwardExtension_ReweightAlphaSpectrum_v4";
   // TString suffix = "_ForwardExtension_ReweightGluonSplitting_v4";

   // TString suffix = "_ForwardExtensionNextToCentral_v3";
   // TString suffix = "_ForwardExtensionNextToCentral_MBXS73500_v3";
   // TString suffix = "_ForwardExtensionNextToCentral_JECup_v3";
   // TString suffix = "_ForwardExtensionNextToCentral_JECdown_v3";
   // TString suffix = "_ForwardExtensionNextToCentral_PLIup_v3";
   // TString suffix = "_ForwardExtensionNextToCentral_PLIdown_v3";
   // TString suffix = "_ForwardExtensionNextToCentral_95Truncation_v3";
   // TString suffix = "_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_v3";
   // TString suffix = "_ForwardExtensionNextToCentral_ReweightAlphaSpectrum_v3";
   // TString suffix = "_ForwardExtensionNextToCentral_ReweightGluonSplitting_v3";

   // TString suffix = "_ForwardExtensionSecondNextToCentral_v1";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_MBXS73500_v1";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_JECup_v1";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_JECdown_v1";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_PLIup_v1";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_PLIdown_v1";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_95Truncation_v1";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_v1";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_ReweightAlphaSpectrum_v1";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_ReweightGluonSplitting_v1";

   // TString suffix = "_herwigpp_ForwardExtension_v4";
   // TString suffix = "_herwigpp_ForwardExtensionNextToCentral_v3";
   // TString suffix = "_herwigpp_ForwardExtensionSecondNextToCentral_v1";
  
   // --------------------- //

   bool MakeAllControlPlots = false;
   bool useTruth = false;
    
   // define scale factor for PLI correction --> can be used for systematic variation
   double PLI_scale = 1.0; // 0.75, 1.0 or 1.25

   // define helper histos
   TH1F *tmp_data1 = new TH1F();
   TH1F *tmp_data2 = new TH1F();
   TH1F *tmp_mc = new TH1F();
   TH1F *tmp_gen = new TH1F();
   TH1F *tmp_res = new TH1F();
   tmp_data1->Sumw2();
   tmp_data2->Sumw2();
   tmp_mc->Sumw2();
   tmp_gen->Sumw2();
   tmp_res->Sumw2();

   std::vector<float> alpha;
   // alpha.push_back(0.05); 
   alpha.push_back(0.1); 
   alpha.push_back(0.125); 
   alpha.push_back(0.15); 
   alpha.push_back(0.175); 
   alpha.push_back(0.20); 
   alpha.push_back(0.225);
   alpha.push_back(0.25); 

   const int Neta = 7;
   const int Npt = 13;
   int Nalpha = alpha.size();
      
   float pt_bins[Npt+1] = {62, 107, 175, 205, 242, 270, 310, 335, 379, 410, 467, 600, 1000, 2000};
   // float eta_bins[6] = {0, 0.5, 1.1, 1.7, 2.3, 5.0};
   // float eta_bins[7] = {0, 0.5, 1.1, 1.7, 2.3, 3.2, 5.0};
   float eta_bins[Neta+1] = {0, 0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
   TH1F *extrapolated_data = new TH1F("extrapolated_data", "extrapolated_data", Npt, pt_bins);
   TH1F *extrapolated_data_slope = new TH1F("extrapolated_data_slope", "extrapolated_data", Npt, pt_bins);
   TH1F *extrapolated_mc = new TH1F("extrapolated_mc", "extrapolated_mc", Npt, pt_bins);
   TH1F *extrapolated_mc_slope = new TH1F("extrapolated_mc_slope", "extrapolated_mc", Npt, pt_bins);
   TH1F *extrapolated_gen = new TH1F("extrapolated_gen", "extrapolated_gen", Npt, pt_bins);
   TH1F *extrapolated_gen_slope = new TH1F("extrapolated_gen_slope", "extrapolated_gen", Npt, pt_bins);
   TH1F *extrapolated_data_with_pli = new TH1F("extrapolated_data_with_pli", "extrapolated_data", Npt, pt_bins);
   TH1F *extrapolated_mc_with_pli = new TH1F("extrapolated_mc_with_pli", "extrapolated_mc", Npt, pt_bins);
   TH1F *truth_resolution = new TH1F("truth_resolution", "truth_resolution", Npt, pt_bins);
   TH1F *chi2_data = new TH1F("chi2_data", "chi2/ndf", Npt, pt_bins);
   TH1F *chi2_mc = new TH1F("chi2_mc", "chi2/ndf", Npt, pt_bins);
   TH1F *chi2_gen = new TH1F("chi2_gen", "chi2/ndf", Npt, pt_bins);
   // TH1F *RatioVsEta = new TH1F("RatioVsEta", "", Neta, eta_bins);
   // TH1F *RatioVsEta_with_pli = new TH1F("RatioVsEta_with_pli", "", Neta, eta_bins);
   // TH1F *RatioVsEta = new TH1F("RatioVsEta", "", Neta, eta_bins);
   // TH1F *RatioVsEta_with_pli = new TH1F("RatioVsEta_with_pli", "", Neta, eta_bins);
   TH1F *RatioVsEta = new TH1F("RatioVsEta", "", Neta, eta_bins);
   TH1F *RatioVsEta_with_pli = new TH1F("RatioVsEta_with_pli", "", Neta, eta_bins);
   extrapolated_data->Sumw2();
   extrapolated_data_slope->Sumw2();
   extrapolated_mc->Sumw2();
   extrapolated_mc_slope->Sumw2();
   extrapolated_gen->Sumw2();
   extrapolated_gen_slope->Sumw2();
   extrapolated_data_with_pli->Sumw2();
   extrapolated_mc_with_pli->Sumw2();
   truth_resolution->Sumw2();
   chi2_data->Sumw2();
   chi2_mc->Sumw2();
   chi2_gen->Sumw2();
   RatioVsEta->Sumw2();
   RatioVsEta_with_pli->Sumw2();

   // how much should be truncated?
   double xq_IQW = 0.985;
   
   //// get asymmetry histos
   for(int ieta=0; ieta < Neta; ++ieta){
      //  cout << "eta Bin: " << ieta << endl;

      extrapolated_data->Reset();
      extrapolated_data_slope->Reset();
      extrapolated_mc->Reset();
      extrapolated_mc_slope->Reset();
      extrapolated_gen->Reset();
      extrapolated_gen_slope->Reset();
      extrapolated_data_with_pli->Reset();
      extrapolated_mc_with_pli->Reset();
      truth_resolution->Reset();
      chi2_data->Reset();
      chi2_mc->Reset();
      chi2_gen->Reset();

      for(int ipt=0; ipt < Npt; ++ipt){  
 
         // hack !!!
         // ------ //
         if(suffix.Contains("herwigpp_ForwardExtensionSecondNextToCentral") && ipt == 4 && ieta == 4) continue;
         // ------ //

         //  cout << "pt Bin: " << ipt << endl;
         std::vector<double> x,x_e,MCy,MCy_e,Datay,Datay_e,Geny,Geny_e;
         for(int ialpha=0; ialpha < Nalpha; ++ialpha){        // nominal
            //  cout << "alpha Bin: " << ialpha << endl;
            TString hname = Form("Forward_Pt%i_eta%i_alpha%i", ipt, ieta, ialpha);
            TString hname_gen = Form("Forward_GenAsymm_Pt%i_eta%i_alpha%i", ipt, ieta, ialpha);

            //  cout << "hname: " << hname << endl;

            // make sure to use histos from correct dataset
            if( ipt < 9) {
               jetmon_file->cd();
               tmp_data1 = 0;
               tmp_data1 = (TH1F*) gDirectory->FindObjectAny(hname);
               jet_file->cd();
               tmp_data2 = 0;
               tmp_data2 = (TH1F*) gDirectory->FindObjectAny(hname);
            }
            else {
               jetht_file->cd();
               tmp_data1 = 0;
               tmp_data1 = (TH1F*) gDirectory->FindObjectAny(hname);
               jet_file->cd();
               tmp_data2 = 0;
               tmp_data2 = (TH1F*) gDirectory->FindObjectAny(hname);
            }

            // add histos from different datasets
            tmp_data1->Add(tmp_data2);

            mc_file->cd();
            tmp_mc = 0;
            tmp_gen = 0;
            tmp_mc = (TH1F*) gDirectory->FindObjectAny(hname);
            tmp_gen = (TH1F*) gDirectory->FindObjectAny(hname_gen);
                               
            x.push_back(alpha.at(ialpha));
            x_e.push_back(0.);

            cout << "Pt: " << ipt << " eta: " << ieta << " alpha: " << ialpha << endl;
              
            tmp_mc->Scale(tmp_data1->Integral()/tmp_mc->Integral());
         
            double mc_width = GetTruthRes(tmp_mc, xq_IQW);
            double mc_width_err = GetTruthResErr(tmp_mc, xq_IQW);
            double data_width = GetTruthRes(tmp_data1, xq_IQW);
            double data_width_err = GetTruthResErr(tmp_data1, xq_IQW);
            double gen_width = GetTruthRes(tmp_gen, xq_IQW);
            double gen_width_err = GetTruthResErr(tmp_gen, xq_IQW);

            tmp_mc->Rebin(10);
            tmp_gen->Rebin(10);
            tmp_gen->Scale(tmp_data1->Integral()/tmp_gen->Integral());
            tmp_data1->Rebin(10);

            // define Gaussian functions using width of truncated RMS
            TF1 *gauss_mc = new TF1("gauss_mc", "gaus(0)", 0, 1);
            TF1 *gauss_gen = new TF1("gauss_gen", "gaus(0)", 0, 1);
            TF1 *gauss_data = new TF1("gauss_data", "gaus(0)", 0, 1);
            gauss_mc->SetParameter(0, (1/(mc_width*TMath::Sqrt(2*TMath::Pi())))*tmp_mc->Integral()*tmp_mc->GetBinWidth(5)*2);
            gauss_mc->SetParameter(1, 0);
            gauss_mc->SetParameter(2, mc_width);
            gauss_gen->SetParameter(0, (1/(gen_width*TMath::Sqrt(2*TMath::Pi())))*tmp_gen->Integral()*tmp_gen->GetBinWidth(5)*2);
            gauss_gen->SetParameter(1, 0);
            gauss_gen->SetParameter(2, gen_width);
            gauss_data->SetParameter(0, (1/(data_width*TMath::Sqrt(2*TMath::Pi())))*tmp_data1->Integral()*tmp_data1->GetBinWidth(5)*2);
            gauss_data->SetParameter(1, 0);
            gauss_data->SetParameter(2, data_width);
          
            // use truncated RMS
            MCy.push_back( mc_width );
            MCy_e.push_back( mc_width_err );
            Datay.push_back( data_width );
            Datay_e.push_back( data_width_err );
            Geny.push_back( gen_width );
            Geny_e.push_back( gen_width_err );

            // --------------------- //
            // draw asymmetry histos
            TCanvas *c5 = new TCanvas("c5", "", 600, 600);
            c5->SetLogy();
            tmp_mc->GetYaxis()->SetRangeUser(0.1, 1000.*tmp_mc->GetMaximum());
            tmp_mc->GetXaxis()->SetTitle("|A|");
            tmp_mc->GetYaxis()->SetTitle("Events");
            tmp_mc->SetLineColor(30);
            tmp_mc->SetFillColor(30);
            tmp_mc->Draw("hist");
            tmp_data1->SetMarkerStyle(20);
            tmp_data1->Draw("same");

            TPaveText *label = util::LabelFactory::createPaveTextWithOffset(5,1.0,0.01);
            label->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
            label->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq  p_{T}^{ave} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label->AddText( Form("#alpha #leq %0.3f", alpha.at(ialpha)) );
            label->AddText( Form("mean data: %0.3f", tmp_data1->GetMean()) );
            label->AddText( Form("mean mc: %0.5f", tmp_mc->GetMean()) );
            label->Draw("same");
            
            TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(2,1.0,0.25);
            leg1->AddEntry(tmp_data1,"Data","LP");
            leg1->AddEntry(tmp_mc,"Simulation","LF");
            leg1->Draw();

            if(ieta == 0 && ipt == 0 && ialpha == 0 ) c5->Print("ForwardExtrapolation/AsymmHistos" + suffix + ".eps(");
            else if(ieta == Neta-1 && ipt == Npt-1 && ialpha == Nalpha-1) c5->Print("ForwardExtrapolation/AsymmHistos" + suffix + ".eps)");
            else c5->Print("ForwardExtrapolation/AsymmHistos" + suffix + ".eps"); 

            if(MakeAllControlPlots) {
               TString asymm_name;
               asymm_name = Form("ForwardExtrapolation/AsymmHistos_Eta%i_pt%i_alpha%i" + suffix + ".eps", ieta, ipt, ialpha);
               c5->Print(asymm_name);
            }
            // --------------------- //

            // --------------------- //
            // draw gen asymmetry histos
            TCanvas *c5b = new TCanvas("c5b", "", 600, 600);
            c5b->SetLogy();
            tmp_gen->GetYaxis()->SetRangeUser(0.1, 1000.*tmp_gen->GetMaximum());
            tmp_gen->GetXaxis()->SetTitle("|A_{gen}|");
            tmp_gen->GetYaxis()->SetTitle("Events");
            tmp_gen->SetLineColor(kOrange-3);
            tmp_gen->SetFillColor(kOrange-3);
            tmp_gen->Draw("hist");
                                
            TPaveText *label2 = util::LabelFactory::createPaveTextWithOffset(4,1.0,0.01);
            label2->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
            label2->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq  p_{T}^{ave} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label2->AddText( Form("#alpha #leq %0.3f", alpha.at(ialpha)) );
            label->AddText( Form("mean gen: %0.5f", tmp_gen->GetMean()) );
            label2->Draw("same");
            
            TLegend* leg2 = util::LabelFactory::createLegendColWithOffset(1,1.0,0.25);
            leg2->AddEntry(tmp_gen,"Simulation","LF");
            leg2->Draw();

            if(ieta == 0 && ipt == 0 && ialpha == 0 ) c5b->Print("ForwardExtrapolation/GenAsymmHistos" + suffix + ".eps(");
            else if(ieta == Neta-1 && ipt == Npt-1 && ialpha == Nalpha-1) c5b->Print("ForwardExtrapolation/GenAsymmHistos" + suffix + ".eps)");
            else c5b->Print("ForwardExtrapolation/GenAsymmHistos" + suffix + ".eps"); 

            if(MakeAllControlPlots) {
               TString genasymm_name;
               genasymm_name = Form("ForwardExtrapolation/GenAsymmHistos_Eta%i_pt%i_alpha%i" + suffix + ".eps", ieta, ipt, ialpha);
               c5b->Print(genasymm_name);
            }
            // --------------------- //

            // --------------------- //
            TCanvas *bb = new TCanvas("bb", "", 600, 600);
            bb->SetLogy();
            bb->SetBottomMargin(0.25 + 0.75*bb->GetBottomMargin()-0.25*bb->GetTopMargin());
            bb->cd();
            tmp_data1->GetYaxis()->SetRangeUser(0.1, 100.*tmp_data1->GetMaximum());
            tmp_data1->GetXaxis()->SetLabelSize(0);
            gauss_data->GetXaxis()->SetLabelSize(0);
            tmp_data1->GetXaxis()->SetTitle("");
            tmp_data1->GetYaxis()->SetTitle("Events");
            tmp_data1->SetLineColor(kOrange-3);
            tmp_data1->SetFillColor(kOrange-3);
            tmp_data1->Draw("hist");
            gauss_data->SetLineWidth(2);
            gauss_data->SetLineColor(kBlue+1);
            gauss_data->Draw("same");

            TPaveText *label4 = util::LabelFactory::createPaveTextWithOffset(3,0.8,0.01);
            label4->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
            label4->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq  p_{T}^{ave} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label4->AddText( Form("#alpha #leq %0.3f", alpha.at(ialpha)) );
            label4->Draw("same");
            
            TLegend* leg4 = util::LabelFactory::createLegendColWithOffset(2,0.65,0.2);
            leg4->AddEntry(tmp_data1,"Data","LF");
            leg4->AddEntry(gauss_data,"Gaussian Function","LF");
            leg4->Draw();

            TPad *pad3 = new TPad("pad3a", "pad3a", 0, 0, 1, 1);
            pad3->SetTopMargin(0.75 - 0.75*pad3->GetBottomMargin()+0.25*pad3->GetTopMargin());
            pad3->SetFillStyle(0);
            pad3->SetFrameFillColor(10);
            pad3->SetFrameBorderMode(0);
            pad3->Draw();
            pad3->cd();

            TH1F* r1 = new TH1F(*tmp_data1);
            r1->Reset();
            r1->SetTitle("");
            r1->GetXaxis()->SetTitle("|A|");
            r1->GetYaxis()->SetTitle("Data / Gauss.");
            r1->SetLineColor(kBlack);
            r1->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
            r1->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y"));
            r1->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.3);
            r1->GetYaxis()->CenterTitle();
            r1->GetYaxis()->SetNdivisions(505);
            r1->SetStats(0);
            r1->SetMarkerStyle(20);
            r1->SetMarkerSize(1.12);
            r1->SetMarkerColor(kBlack);
            r1->Add(tmp_data1, 1);
            r1->Divide(gauss_data);
            r1->GetYaxis()->SetRangeUser(0.7, 1.3);
            r1->Draw("ep");
            TLine l2;
            l2.DrawLine(0., 1., 1., 1.);
            bb->cd();
           
            if(ieta == 0 && ipt == 0 && ialpha == 0 ) bb->Print("ForwardExtrapolation/AsymmHistosDataWithRatio" + suffix + ".eps(");
            else if(ieta == Neta-1 && ipt == Npt-1 && ialpha == Nalpha-1) bb->Print("ForwardExtrapolation/AsymmHistosDataWithRatio" + suffix + ".eps)");
            else bb->Print("ForwardExtrapolation/AsymmHistosDataWithRatio" + suffix + ".eps"); 

            if(MakeAllControlPlots) {
               TString ratiodata_name;
               ratiodata_name = Form("ForwardExtrapolation/AsymmHistosDataWithRatio_Eta%i_pt%i_alpha%i" + suffix + ".eps", ieta, ipt, ialpha);
               bb->Print(ratiodata_name);
            }
            // --------------------- //

            // --------------------- //
            TCanvas *cc = new TCanvas("cc", "", 600, 600);
            cc->SetLogy();
            cc->SetBottomMargin(0.25 + 0.75*cc->GetBottomMargin()-0.25*cc->GetTopMargin());
            cc->cd();
            tmp_mc->GetYaxis()->SetRangeUser(0.1, 100.*tmp_mc->GetMaximum());
            tmp_mc->GetXaxis()->SetLabelSize(0);
            gauss_mc->GetXaxis()->SetLabelSize(0);
            tmp_mc->GetXaxis()->SetTitle("");
            tmp_mc->GetYaxis()->SetTitle("Events");
            tmp_mc->SetLineColor(30);
            tmp_mc->SetFillColor(30);
            tmp_mc->Draw("hist");
            gauss_mc->SetLineWidth(2);
            gauss_mc->SetLineColor(kRed);
            gauss_mc->Draw("same");

            TPaveText *label3 = util::LabelFactory::createPaveTextWithOffset(3,0.8,0.01);
            label3->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
            label3->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq  p_{T}^{ave} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label3->AddText( Form("#alpha #leq %0.3f", alpha.at(ialpha)) );
            label3->Draw("same");
            
            TLegend* leg3 = util::LabelFactory::createLegendColWithOffset(2,0.65,0.2);
            leg3->AddEntry(tmp_mc,"Simulation","LF");
            leg3->AddEntry(gauss_mc,"Gaussian Function","LF");
            leg3->Draw();

            TPad *pad2 = new TPad("pad2a", "pad2a", 0, 0, 1, 1);
            pad2->SetTopMargin(0.75 - 0.75*pad2->GetBottomMargin()+0.25*pad2->GetTopMargin());
            pad2->SetFillStyle(0);
            pad2->SetFrameFillColor(10);
            pad2->SetFrameBorderMode(0);
            pad2->Draw();
            pad2->cd();

            TH1F* r = new TH1F(*tmp_mc);
            r->Reset();
            r->SetTitle("");
            r->GetXaxis()->SetTitle("|A|");
            r->GetYaxis()->SetTitle("Sim. / Gauss.");
            r->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
            r->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y"));
            r->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.3);
            r->GetYaxis()->CenterTitle();
            r->GetYaxis()->SetNdivisions(505);
            r->SetStats(0);
            r->SetMarkerStyle(20);
            r->SetMarkerSize(1.12);
            r->SetMarkerColor(kBlack);
            r->SetLineColor(kBlack);
            r->Add(tmp_mc, 1);
            r->Divide(gauss_mc);
            r->GetYaxis()->SetRangeUser(0.7, 1.3);
            r->Draw("ep");
            TLine l;
            l.DrawLine(0., 1., 1., 1.);
            cc->cd();
           
            if(ieta == 0 && ipt == 0 && ialpha == 0 ) cc->Print("ForwardExtrapolation/AsymmHistosSimWithRatio" + suffix + ".eps(");
            else if(ieta == Neta-1 && ipt == Npt-1 && ialpha == Nalpha-1) cc->Print("ForwardExtrapolation/AsymmHistosSimWithRatio" + suffix + ".eps)");
            else cc->Print("ForwardExtrapolation/AsymmHistosSimWithRatio" + suffix + ".eps"); 

            if(MakeAllControlPlots) {
               TString ratiomc_name;
               ratiomc_name = Form("ForwardExtrapolation/AsymmHistosSimWithRatio_Eta%i_pt%i_alpha%i" + suffix + ".eps", ieta, ipt, ialpha);
               cc->Print(ratiomc_name);
            }
            // --------------------- //

            // --------------------- //
            TCanvas *dd = new TCanvas("dd", "", 600, 600);
            dd->SetLogy();
            dd->SetBottomMargin(0.25 + 0.75*dd->GetBottomMargin()-0.25*dd->GetTopMargin());
            dd->cd();
            tmp_gen->GetYaxis()->SetRangeUser(0.1, 100.*tmp_gen->GetMaximum());
            tmp_gen->GetXaxis()->SetLabelSize(0);
            gauss_gen->GetXaxis()->SetLabelSize(0);
            tmp_gen->GetXaxis()->SetTitle("");
            tmp_gen->GetYaxis()->SetTitle("Events");
            tmp_gen->SetLineColor(kBlue-10);
            tmp_gen->SetFillColor(kBlue-10);
            tmp_gen->Draw("hist");
            gauss_gen->SetLineWidth(2);
            gauss_gen->SetLineColor(kRed);
            gauss_gen->Draw("same");

            TPaveText *label3d = util::LabelFactory::createPaveTextWithOffset(3,0.8,0.01);
            label3d->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
            label3d->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq  p_{T}^{ave} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label3d->AddText( Form("#alpha #leq %0.3f", alpha.at(ialpha)) );
            label3d->Draw("same");
            
            TLegend* leg3d = util::LabelFactory::createLegendColWithOffset(2,0.65,0.2);
            leg3d->AddEntry(tmp_gen,"Simulation","LF");
            leg3d->AddEntry(gauss_gen,"Gaussian Function","LF");
            leg3d->Draw();

            TPad *pad2d = new TPad("pad2da", "pad2da", 0, 0, 1, 1);
            pad2d->SetTopMargin(0.75 - 0.75*pad2d->GetBottomMargin()+0.25*pad2d->GetTopMargin());
            pad2d->SetFillStyle(0);
            pad2d->SetFrameFillColor(10);
            pad2d->SetFrameBorderMode(0);
            pad2d->Draw();
            pad2d->cd();

            TH1F* rd = new TH1F(*tmp_gen);
            rd->Reset();
            rd->SetTitle("");
            rd->GetXaxis()->SetTitle("|A_{Gen}|");
            rd->GetYaxis()->SetTitle("Gen. / Gauss.");
            rd->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
            rd->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y"));
            rd->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.3);
            rd->GetYaxis()->CenterTitle();
            rd->GetYaxis()->SetNdivisions(505);
            rd->SetStats(0);
            rd->SetMarkerStyle(20);
            rd->SetMarkerSize(1.12);
            rd->SetMarkerColor(kBlack);
            rd->SetLineColor(kBlack);
            rd->Add(tmp_gen, 1);
            rd->Divide(gauss_gen);
            rd->GetYaxis()->SetRangeUser(0.7, 1.3);
            rd->Draw("ep");
            TLine ld;
            ld.DrawLine(0., 1., 1., 1.);
            dd->cd();
           
            if(ieta == 0 && ipt == 0 && ialpha == 0 ) dd->Print("ForwardExtrapolation/AsymmHistosGenWithRatio" + suffix + ".eps(");
            else if(ieta == Neta-1 && ipt == Npt-1 && ialpha == Nalpha-1) dd->Print("ForwardExtrapolation/AsymmHistosGenWithRatio" + suffix + ".eps)");
            else dd->Print("ForwardExtrapolation/AsymmHistosGenWithRatio" + suffix + ".eps"); 

            if(MakeAllControlPlots) {
               TString ratiogen_name;
               ratiogen_name = Form("ForwardExtrapolation/AsymmHistosGenWithRatio_Eta%i_pt%i_alpha%i" + suffix + ".eps", ieta, ipt, ialpha);
               dd->Print(ratiogen_name);
            }
            // --------------------- //
         }

         // Covariance matrices needed for fitting 
         TMatrixD y_cov_mc;
         TMatrixD y_cov_data;
         TMatrixD y_cov_gen;
         y_cov_mc.ResizeTo(alpha.size(), alpha.size());
         y_cov_data.ResizeTo(alpha.size(), alpha.size());
         y_cov_gen.ResizeTo(alpha.size(), alpha.size());

         // fill covariance matrix for data and mc
         for(unsigned int ialpha=0; ialpha < alpha.size(); ++ialpha){
            for(unsigned int jalpha =0; jalpha < alpha.size(); jalpha++){
               if( ialpha <= jalpha ) {
                  double n1_mc = pow(MCy.at(ialpha),2)/(2*pow(MCy_e.at(ialpha),2));
                  double n2_mc = pow(MCy.at(jalpha),2)/(2*pow(MCy_e.at(jalpha),2));
          
                  double n1_data = pow(Datay.at(ialpha),2)/(2*pow(Datay_e.at(ialpha),2));
                  double n2_data = pow(Datay.at(jalpha),2)/(2*pow(Datay_e.at(jalpha),2));

                  double n1_gen = pow(Geny.at(ialpha),2)/(2*pow(Geny_e.at(ialpha),2));
                  double n2_gen = pow(Geny.at(jalpha),2)/(2*pow(Geny_e.at(jalpha),2));

                  y_cov_mc(ialpha, jalpha) = pow(MCy_e.at(ialpha),2) * pow((n1_mc/n2_mc),2)*
                     (MCy.at(ialpha)/MCy.at(jalpha));
                  y_cov_data(ialpha, jalpha) = pow(Datay_e.at(ialpha),2) * pow((n1_data/n2_data),2)*
                     (Datay.at(ialpha)/Datay.at(jalpha));  
                  y_cov_gen(ialpha, jalpha) = pow(Geny_e.at(ialpha),2) * pow((n1_gen/n2_gen),2)*
                     (Geny.at(ialpha)/Geny.at(jalpha));  
               }
               else {
                  double n1_mc = pow(MCy.at(jalpha),2)/(2*pow(MCy_e.at(jalpha),2));
                  double n2_mc = pow(MCy.at(ialpha),2)/(2*pow(MCy_e.at(ialpha),2));

                  double n1_data = pow(Datay.at(jalpha),2)/(2*pow(Datay_e.at(jalpha),2));
                  double n2_data = pow(Datay.at(ialpha),2)/(2*pow(Datay_e.at(ialpha),2));

                  double n1_gen = pow(Geny.at(jalpha),2)/(2*pow(Geny_e.at(jalpha),2));
                  double n2_gen = pow(Geny.at(ialpha),2)/(2*pow(Geny_e.at(ialpha),2));

                  y_cov_mc(ialpha, jalpha) = pow(MCy_e.at(jalpha),2) * pow((n1_mc/n2_mc),2)*
                     (MCy.at(jalpha)/MCy.at(ialpha));
                  y_cov_data(ialpha, jalpha) = pow(Datay_e.at(jalpha),2) * pow((n1_data/n2_data),2)*
                     (Datay.at(jalpha)/Datay.at(ialpha));
                  y_cov_gen(ialpha, jalpha) = pow(Geny_e.at(jalpha),2) * pow((n1_gen/n2_gen),2)*
                     (Geny.at(jalpha)/Geny.at(ialpha));
               }
            }
         }        
  
         //create TGraphErrors from previously defined vectors
         TGraphErrors *extrapol_MC = new TGraphErrors(alpha.size(),&x[0],&MCy[0],&x_e[0],&MCy_e[0]);
         TGraphErrors *extrapol_Data = new TGraphErrors(alpha.size(),&x[0],&Datay[0],&x_e[0],&Datay_e[0]);
         TGraphErrors *extrapol_Gen = new TGraphErrors(alpha.size(),&x[0],&Geny[0],&x_e[0],&Geny_e[0]);
 
         // fit linear extrapolation function
         TF1 *lin_extrapol_mc = new TF1("lin_extrapol_mc","[0]+[1]*x",0,alpha.back()+0.05); 
         TF1 *lin_extrapol_data = new TF1("lin_extrapol_data","[0]+[1]*x",0,alpha.back()+0.05);
         TF1 *lin_extrapol_gen = new TF1("lin_extrapol_gen","[0]+[1]*x",0,alpha.back()+0.05);
                  
         //fit extrapolation function to the TGraphErrors for data and MC  
         // extrapol_MC->Fit("lin_extrapol_mc","Q","same",0,alpha.back()+0.05);
         // extrapol_Data->Fit("lin_extrapol_data","Q","same",0,alpha.back()+0.05);
         // extrapol_Gen->Fit("lin_extrapol_gen","Q","same",0,alpha.back()+0.05);
           
         // --------------------- //
         // fit mc
         data.reset();
         data.x_val = x;
         data.y_val = MCy;
         data.y_cov.ResizeTo(alpha.size(), alpha.size());
         data.y_cov = y_cov_mc;
         data.CheckPoints();
         
         // choose start values for MC fit
         double slope = (MCy.at(MCy.size()-1) - MCy.at(MCy.size()-4))/(x.at(x.size()-1) - x.at(x.size()-4));
         double d_slope = slope;
         double offset = MCy.at(MCy.size()-1) - (slope*x.at(x.size()-1));
         double d_offset = offset;
             
         std::cout << "MC start values: " << "slope: " << slope << " offset: " << offset << std::endl; 
         double chi2MC = make_lin_fit(slope, d_slope, offset, d_offset);
         std::cout << "MC fit values: " << "slope: " << slope << " offset: " << offset << std::endl; 

         chi2_mc->SetBinContent(ipt+1, chi2MC/(data.x_val.size()-2));

         lin_extrapol_mc->SetParameter(0, offset);
         lin_extrapol_mc->SetParError(0, d_offset);
         lin_extrapol_mc->SetParameter(1, slope);
         lin_extrapol_mc->SetParError(1, d_slope);
         extrapol_MC->GetListOfFunctions()->Add(lin_extrapol_mc);
         
         data.reset();
         // --------------------- //
         
         // --------------------- //
         // fit data
         data.x_val = x;
         data.y_val = Datay;
         data.y_cov.ResizeTo(alpha.size(), alpha.size());
         data.y_cov = y_cov_data;
         data.CheckPoints();
         
         // choose start values for data fit
         slope = (Datay.at(Datay.size()-1) - Datay.at(Datay.size()-4))/(x.at(x.size()-1) - x.at(x.size()-4));
         d_slope = slope;
         offset = Datay.at(Datay.size()-1) - (slope*x.at(x.size()-1));
         d_offset = offset;
            
         std::cout << "Data start values: " << "slope: " << slope << " offset: " << offset << std::endl; 
         double chi2Data = make_lin_fit(slope, d_slope, offset, d_offset);
         std::cout << "Data fit values: " << "slope: " << slope << " offset: " << offset << std::endl; 

         chi2_data->SetBinContent(ipt+1, chi2Data/(data.x_val.size()-2));
         
         lin_extrapol_data->SetParameter(0, offset);
         lin_extrapol_data->SetParError(0, d_offset);
         lin_extrapol_data->SetParameter(1, slope);
         lin_extrapol_data->SetParError(1, d_slope);
         extrapol_Data->GetListOfFunctions()->Add(lin_extrapol_data);
         
         data.reset();
         // --------------------- //

         // --------------------- //
         // fit gen
         data.x_val = x;
         data.y_val = Geny;
         data.y_cov.ResizeTo(alpha.size(), alpha.size());
         data.y_cov = y_cov_gen;
         data.CheckPoints();
         
         // choose start values for gen fit
         slope = (Geny.at(Geny.size()-1) - Geny.at(Geny.size()-4))/(x.at(x.size()-1) - x.at(x.size()-4));
         d_slope = slope;
         offset = Geny.at(Geny.size()-1) - (slope*x.at(x.size()-1));
         d_offset = offset;
            
         std::cout << "Gen start values: " << "slope: " << slope << " offset: " << offset << std::endl; 
         double chi2Gen = make_lin_fit(slope, d_slope, offset, d_offset);
         std::cout << "Gen fit values: " << "slope: " << slope << " offset: " << offset << std::endl; 

         chi2_gen->SetBinContent(ipt+1, chi2Gen/(data.x_val.size()-2));
         
         lin_extrapol_gen->SetParameter(0, offset);
         lin_extrapol_gen->SetParError(0, d_offset);
         lin_extrapol_gen->SetParameter(1, slope);
         lin_extrapol_gen->SetParError(1, d_slope);
         extrapol_Gen->GetListOfFunctions()->Add(lin_extrapol_gen);
         
         data.reset();
         // --------------------- //
         
         // --------------------- //
         // draw extrapolations data + mc
         TCanvas *c = new TCanvas("c","",600,600);
         std::pair <float,float> minMaxPair = determineMinMax(extrapol_Data);
         c->DrawFrame(0,minMaxPair.first*0.5-0.05,alpha.back()+0.05,minMaxPair.second*1.47,(";#alpha_{max};#sigma_{A}"));
         extrapol_MC->SetMarkerStyle(20);
         extrapol_MC->SetMarkerColor(kRed+1);
         extrapol_MC->SetLineColor(kRed+1);
         extrapol_MC->Draw("P");
         extrapol_Data->SetMarkerStyle(20);
         extrapol_Data->SetMarkerColor(kBlack);
         extrapol_Data->SetLineColor(kBlack);
         extrapol_Data->Draw("Psame");
         TF1* MCTemp = new TF1();
         TF1* DataTemp = new TF1();
         extrapol_MC->GetFunction("lin_extrapol_mc")->SetLineColor(kRed+1);
         extrapol_MC->GetFunction("lin_extrapol_mc")->SetLineStyle(2);
         extrapol_Data->GetFunction("lin_extrapol_data")->SetLineColor(kBlack);
         extrapol_Data->GetFunction("lin_extrapol_data")->SetLineStyle(2);
         MCTemp=(TF1*) extrapol_MC->GetFunction("lin_extrapol_mc")->Clone();
         DataTemp=(TF1*) extrapol_Data->GetFunction("lin_extrapol_data")->Clone();
         MCTemp->SetRange(0.1,1);
         MCTemp->SetLineStyle(1);
         MCTemp->Draw("same");
         DataTemp->SetRange(0.1,1);
         DataTemp->SetLineStyle(1);
         DataTemp->Draw("same");
       
         TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.5);
         label->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
         label->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq p_{T}^{ave} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
         label->Draw("same");
  
         TLegend* leg1 = util::LabelFactory::createLegendWithOffset(2,0.6);
         leg1->AddEntry(extrapol_Data,"Extrapolation (Data)","LP");
         leg1->AddEntry(extrapol_MC,"Extrapolation (MC)","LP");

         leg1->Draw();
         //   cmsPrel(-1, false , 8);

         TString name;
         name = Form("ForwardExtrapolation/Extrapol_Eta%i_pt%i" + suffix + ".eps", ieta, ipt);
         c->Print(name);
         // --------------------- //

         // --------------------- //
         // draw extrapolations gen (PLI)
         TCanvas *cb = new TCanvas("c","",600,600);
         std::pair <float,float> minMaxPair2 = determineMinMax(extrapol_Data);
         cb->DrawFrame(0,minMaxPair2.first*0.5-0.05,alpha.back()+0.05,minMaxPair2.second*1.47,(";#alpha_{max};#sigma_{A, gen}"));
         extrapol_Gen->SetMarkerStyle(20);
         extrapol_Gen->SetMarkerColor(kBlue+1);
         extrapol_Gen->SetLineColor(kBlue+1);
         extrapol_Gen->Draw("P");
         TF1* GenTemp = new TF1();
         extrapol_Gen->GetFunction("lin_extrapol_gen")->SetLineColor(kBlue+1);
         extrapol_Gen->GetFunction("lin_extrapol_gen")->SetLineStyle(2);
         GenTemp=(TF1*) extrapol_Gen->GetFunction("lin_extrapol_gen")->Clone();
         GenTemp->SetRange(0.1,1);
         GenTemp->SetLineStyle(1);
         GenTemp->Draw("same");

         TPaveText *label2 = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.05);
         label2->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
         label2->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq p_{T}^{ave} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
         label2->Draw("same");
  
         TLegend* leg2 = util::LabelFactory::createLegendWithOffset(2,0.15);
         leg2->AddEntry(extrapol_Gen,"Extrapolation (PLI)","LP");

         leg2->Draw();
         //  cmsPrel(-1, false , 8);

         TString name2;
         name2 = Form("ForwardExtrapolation/Extrapol_Eta%i_pt%i_gen" + suffix + ".eps", ieta, ipt);
         cb->Print(name2);
         // --------------------- //

         // get values from extrapolation fit
         float par_data = lin_extrapol_data->GetParameter(0);
         float par_data_err = lin_extrapol_data->GetParError(0);
         float par_mc = lin_extrapol_mc->GetParameter(0);
         float par_mc_err = lin_extrapol_mc->GetParError(0);
         float par_gen = lin_extrapol_gen->GetParameter(0);
         float par_gen_err = lin_extrapol_gen->GetParError(0);

         // scale pli
         par_gen = par_gen * PLI_scale;
             
         cout << "ieta : " << ieta << "  ipt : " << ipt << endl;
         cout << "Parameter data: " << par_data << endl;
         cout << "Parameter error data: " << par_data_err << endl;
         cout << "Parameter mc: " << par_mc << endl;
         cout << "Parameter error mc: " << par_mc_err << endl;
         cout << "Parameter gen: " << par_gen << endl;
         cout << "Parameter error gen: " << par_gen_err << endl;

         // fill extrapolated histos
         extrapolated_data->SetBinContent(ipt+1, par_data);
         extrapolated_data->SetBinError(ipt+1, par_data_err);
         extrapolated_data_slope->SetBinContent(ipt+1, lin_extrapol_data->GetParameter(1));
         extrapolated_data_slope->SetBinError(ipt+1, lin_extrapol_data->GetParError(1));
         extrapolated_mc->SetBinContent(ipt+1, par_mc);
         extrapolated_mc->SetBinError(ipt+1, par_mc_err);  
         extrapolated_mc_slope->SetBinContent(ipt+1, lin_extrapol_mc->GetParameter(1));
         extrapolated_mc_slope->SetBinError(ipt+1, lin_extrapol_mc->GetParError(1));
         extrapolated_gen->SetBinContent(ipt+1, par_gen);
         extrapolated_gen->SetBinError(ipt+1, par_gen_err);  
         extrapolated_gen_slope->SetBinContent(ipt+1, lin_extrapol_gen->GetParameter(1));
         extrapolated_gen_slope->SetBinError(ipt+1, lin_extrapol_gen->GetParError(1));

         float par_data_pli_corr = 0;
         float par_data_pli_corr_err = 0;
         float par_mc_pli_corr = 0;
         float par_mc_pli_corr_err = 0;

         // use only bins where extrapolation worked + calc values with pli correction
         if(par_gen > 0 && par_data > 0 && par_mc > 0 && par_data > par_gen && par_mc > par_gen) {
            par_data_pli_corr = TMath::Sqrt(pow(par_data,2) - pow(par_gen,2));
            par_data_pli_corr_err = TMath::Sqrt( pow(par_data,2)/(pow(par_data,2) - pow(par_gen,2)) * pow(par_data_err,2) +  pow(par_gen,2)/(pow(par_data,2) - pow(par_gen,2)) * pow(par_gen_err,2));
            par_mc_pli_corr = TMath::Sqrt(pow(par_mc,2) - pow(par_gen,2));
            par_mc_pli_corr_err = TMath::Sqrt( pow(par_mc,2)/(pow(par_mc,2) - pow(par_gen,2)) * pow(par_mc_err,2) +  pow(par_gen,2)/(pow(par_mc,2) - pow(par_gen,2)) * pow(par_gen_err,2));
         }
            
         // fill extrapolated histos after pli correction
         extrapolated_data_with_pli->SetBinContent(ipt+1, par_data_pli_corr);
         extrapolated_data_with_pli->SetBinError(ipt+1, par_data_pli_corr_err);
         extrapolated_mc_with_pli->SetBinContent(ipt+1, par_mc_pli_corr);
         extrapolated_mc_with_pli->SetBinError(ipt+1, par_mc_pli_corr_err);

         cout << "Parameter data after pli: " << par_data_pli_corr << endl;
         cout << "Parameter error data after pli: " << par_data_pli_corr_err << endl;
         cout << "Parameter mc after pli: " << par_mc_pli_corr << endl;
         cout << "Parameter error mc after pli: " << par_mc_pli_corr_err << endl;

         // calc MC truth resolution
         if(useTruth) {
            // get mc truth resolution
            TString hname_res = Form("Response_Pt%i_eta%i", ipt, ieta);
    
            mc_file->cd();
            tmp_res = 0;
            tmp_res = (TH1F*) gDirectory->FindObjectAny(hname_res);

            double truth_res = GetTruthRes(tmp_res, 0.985);
            double truth_res_err = truth_res/(TMath::Sqrt(2*tmp_res->GetEffectiveEntries()));
            //double truth_res = 0;
            //double truth_res_err = 0;

            double chi2 = 0;
            double ndf = 0;
            double mean = tmp_res->GetMean();
            double mean_err = 0;
            // fit gauss function to core of response
            if( tmp_res->GetEntries() > 100 ) {
               mean = tmp_res->GetMean(); 
               truth_res = tmp_res->GetRMS();
               tmp_res->Fit("gaus","QNI","", mean - 2.5 * truth_res,mean + 2.5 * truth_res);
            
               TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
               f->SetLineColor(kRed);
               f->SetLineWidth(3);
               mean = f->GetParameter(1);
               truth_res = f->GetParameter(2);
            
               if( (tmp_res->Fit(f,"QI","same",mean - 2.0 * truth_res, mean + 2.0 * truth_res) == 0) ) { 
                  mean = f->GetParameter(1);
                  mean_err = f->GetParError(1);
                  truth_res = f->GetParameter(2);
                  truth_res_err = f->GetParError(2);
                  chi2 = f->GetChisquare();
                  ndf = f->GetNDF();
               }
            }
                  
            truth_resolution->SetBinContent(ipt+1, truth_res);
            truth_resolution->SetBinError(ipt+1, truth_res_err);

            tmp_res->Rebin(10);

            TF1 *gauss_res = new TF1("gauss_res", "gaus(0)", 0, 2);
            gauss_res->SetParameter(0, (1/(truth_res*TMath::Sqrt(2*TMath::Pi())))*tmp_res->Integral()*tmp_res->GetBinWidth(5));
            gauss_res->SetParameter(1, mean );
            gauss_res->SetParameter(2, truth_res);  
         
            // draw truth response histos
            // --------------------- //
            TCanvas *res = new TCanvas("res", "", 600, 600);
            res->SetLogy();
            res->SetBottomMargin(0.25 + 0.75*res->GetBottomMargin()-0.25*res->GetTopMargin());
            res->cd();
            tmp_res->GetYaxis()->SetRangeUser(0.000001, 500 * tmp_res->GetMaximum());
            tmp_res->SetLabelSize(0);
            tmp_res->GetXaxis()->SetTitle("");
            tmp_res->Draw("hist");
            gauss_res->Draw("same");

            TPaveText *label3d = util::LabelFactory::createPaveTextWithOffset(3,1.0,0.01);
            label3d->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
            label3d->AddText( Form("%0.1f #leq |#eta_{gen}| #leq %0.1f, %3.0f #leq p_{T, gen} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label3d->AddText("");
            label3d->Draw("same");

            TLegend* legres = util::LabelFactory::createLegendColWithOffset(1,1.0,0.12);
            legres->AddEntry(tmp_res,"Simulation","LF");
            legres->Draw("same");

            TPad *pad2b = new TPad("pad2ba", "pad2ba", 0, 0, 1, 1);
            pad2b->SetTopMargin(0.75 - 0.75*pad2b->GetBottomMargin()+0.25*pad2b->GetTopMargin());
            pad2b->SetFillStyle(0);
            pad2b->SetFrameFillColor(10);
            pad2b->SetFrameBorderMode(0);
            pad2b->Draw();
            pad2b->cd();

            TH1F* r2b = new TH1F(*tmp_res);
            r2b->Reset();
            r2b->SetTitle("");
            r2b->GetXaxis()->SetTitle("p_{T, reco}/p_{T, gen}");
            r2b->GetYaxis()->SetTitle("Gauss. / Sim.");
            r2b->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
            r2b->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y"));
            r2b->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.3);
            r2b->GetYaxis()->CenterTitle();
            r2b->GetYaxis()->SetNdivisions(505);
            r2b->SetStats(0);
            r2b->SetMarkerStyle(20);
            r2b->SetMarkerSize(1.12);
            r2b->SetMarkerColor(kBlack);
            r2b->SetLineColor(kBlack);
            r2b->Add(gauss_res, 1);
            r2b->Divide(tmp_res);
            r2b->GetYaxis()->SetRangeUser(0.7, 1.3);
            r2b->Draw("ep");
            TLine l2b;
            l2b.DrawLine(0., 1., 2., 1.);
            res->cd();

            if(ieta == 0 && ipt == 0 ) res->Print("ForwardExtrapolation/TruthResponse" + suffix + ".eps(");
            else if(ieta == Neta-1 && ipt == Npt-1) res->Print("ForwardExtrapolation/TruthResponse" + suffix + ".eps)");
            else res->Print("ForwardExtrapolation/TruthResponse" + suffix + ".eps");

            if(ieta == 0 && ipt == 8) res->Print("ForwardExtrapolation/TruthResponse_example" + suffix + ".eps"); 
            // --------------------- //
         }
      }

      // --------------------- //
      // write pli to root file for each eta-bin
      TString name_pli;
      name_pli = Form("ForwardExtrapolation/PLI_Eta%i" + suffix + ".root", ieta);
      TFile* output_pli = new TFile(name_pli, "RECREATE");
      extrapolated_gen->Write();
      output_pli->Write();
      // --------------------- //

      // --------------------- //
      // write truth resolution to root file for each eta-bin
      TString name_truth;
      name_truth = Form("ForwardExtrapolation/TruthResolution_Eta%i" + suffix + "_2SigmaGaussFit.root", ieta);
      // name_truth = Form("ForwardExtrapolation/TruthResolution_Eta%i" + suffix + "_1p5SigmaGaussFit.root", ieta);
      // name_truth = Form("ForwardExtrapolation/TruthResolution_Eta%i" + suffix + "_98p5TruncatedRMS.root", ieta);
      // name_truth = Form("ForwardExtrapolation/TruthResolution_Eta%i" + suffix + "_97TruncatedRMS.root", ieta);
      // name_truth = Form("ForwardExtrapolation/TruthResolution_Eta%i" + suffix + "_95TruncatedRMS.root", ieta);
      TFile* output_truth = new TFile(name_truth, "RECREATE");
      truth_resolution->Write();
      output_truth->Write();
      // --------------------- //

      // --------------------- //
      // write absolute resolution after extrapolation to root file for each eta-bin
      TString name_resolution;
      name_resolution = Form("ForwardExtrapolation/AbsoluteResolution_Eta%i" + suffix + ".root", ieta);
      TFile* output_resolution = new TFile(name_resolution, "RECREATE");
      extrapolated_mc_with_pli->Write();
      extrapolated_data_with_pli->Write();
      output_resolution->Write();
      // --------------------- //

      // --------------------------------------- //    
      // draw res after extrapolations (--> closure test for absolute resolution)
      TCanvas *c2 = new TCanvas("c2","",600,600);
      c2->SetLogx();
      c2->SetBottomMargin(0.25 + 0.75*c2->GetBottomMargin()-0.25*c2->GetTopMargin());
      c2->cd();
      truth_resolution->GetXaxis()->SetLabelSize(0);
      truth_resolution->GetYaxis()->SetRangeUser(0, 0.3);
      truth_resolution->GetYaxis()->SetTitle("#sqrt{2}#sigma_{A}");
      truth_resolution->SetMarkerStyle(26);
      truth_resolution->SetMarkerColor(kRed+1);
      truth_resolution->Draw();
      extrapolated_mc_with_pli->Scale(TMath::Sqrt(2));
      extrapolated_mc_with_pli->SetMarkerStyle(25);
      extrapolated_mc_with_pli->SetMarkerColor(kBlue+2);
      extrapolated_mc_with_pli->Draw("same");

      // scale data accordingly to get ratio correct
      extrapolated_data_with_pli->Scale(TMath::Sqrt(2));
    
      TLegend* leg3 = util::LabelFactory::createLegendWithOffset(2,0.175);
      leg3->AddEntry(truth_resolution,"Truth resolution","P");
      leg3->AddEntry(extrapolated_mc_with_pli,"Measured resolution (MC)","P");      
      leg3->Draw();
      // cmsPrel(-1, false , 8);

      TPad *pad2 = new TPad("pad2a", "pad2a", 0, 0, 1, 1);
      pad2->SetLogx();
      pad2->SetTopMargin(0.75 - 0.75*pad2->GetBottomMargin()+0.25*pad2->GetTopMargin());
      pad2->SetFillStyle(0);
      pad2->SetFrameFillColor(10);
      pad2->SetFrameBorderMode(0);
      pad2->Draw();
      pad2->cd();

      Double_t xMin1 = extrapolated_mc_with_pli->GetXaxis()->GetXmin();
      Double_t xMax1 = extrapolated_mc_with_pli->GetXaxis()->GetXmax();

      TH1F* r = new TH1F(*extrapolated_mc_with_pli);
      r->Sumw2();
      r->SetXTitle("p_{T}^{ref} [GeV]");
      r->SetYTitle("#sigma - #sigma_{truth} ) / #sigma");
      r->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
      r->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y"));
      r->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.3);
      r->GetYaxis()->SetNdivisions(505);
      r->SetStats(0);
      r->SetMarkerStyle(20);
      r->SetMarkerColor(kBlack);
      r->Reset();
      r->Add(extrapolated_mc_with_pli, 1);
      r->Add(truth_resolution, -1);
      r->Divide(extrapolated_mc_with_pli);
      r->SetMaximum(0.35);
      r->SetMinimum(-0.35);
      r->Draw("ep");
      TLine l;
      l.DrawLine(xMin1, 0., xMax1, 0.);
      c2->cd();

      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.05);
      label->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
      label->AddText( Form("%0.1f #leq |#eta| #leq %0.1f", eta_bins[ieta], eta_bins[ieta+1]) );
      label->Draw("same");

      TString name2;
      name2 = Form("ForwardExtrapolation/AbsoluteResolutionClosure_Eta%i" + suffix + ".eps", ieta);
      c2->Print(name2);
      // --------------------- //

      // --------------------- //
      // plot results of extrapolation
      TCanvas *c2b = new TCanvas("c2b","",600,600);
      c2b->SetLogx();
      extrapolated_mc->GetYaxis()->SetTitle("#sigma_{A}");
      extrapolated_mc->GetYaxis()->SetRangeUser(0., 0.2);
      extrapolated_mc->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
      extrapolated_mc->SetLineColor(kRed+1);
      extrapolated_mc->Draw("hist");
      extrapolated_gen->SetLineColor(kBlue+1);
      extrapolated_gen->SetMarkerColor(kBlue+1);
      extrapolated_gen->Draw("samehist");
      extrapolated_data->Draw("same");

      TLegend* leg3b = util::LabelFactory::createLegendWithOffset(3,0.175);
      leg3b->AddEntry(extrapolated_data,"Data","PL");
      leg3b->AddEntry(extrapolated_mc,"MC","L");
      leg3b->AddEntry(extrapolated_gen,"PLI","L");
      
      leg3b->Draw();
      // cmsPrel(-1, false , 8);

      TPaveText *labelb = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.05);
      labelb->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
      labelb->AddText( Form("%0.1f #leq |#eta| #leq %0.1f", eta_bins[ieta], eta_bins[ieta+1]) );
      labelb->Draw("same");

      TString name2b;
      name2b = Form("ForwardExtrapolation/Extrapol_Eta%i" + suffix + ".eps", ieta);
      c2b->Print(name2b);
      // --------------------- //

      // --------------------- //
      // plot slopes of extrapolation
      TCanvas *c4b = new TCanvas("c4b","",600,600);
      c4b->SetLogx();
      extrapolated_mc_slope->GetYaxis()->SetTitle("#Delta#sigma_{A}");
      extrapolated_mc_slope->GetYaxis()->SetRangeUser(0.8, 1.2);
      extrapolated_mc_slope->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
      extrapolated_mc_slope->SetMarkerColor(kRed+1);
      extrapolated_mc_slope->SetLineColor(kRed+1);
      extrapolated_gen_slope->SetLineColor(kBlue+1);
      extrapolated_gen_slope->SetMarkerColor(kBlue+1);
      extrapolated_data_slope->GetYaxis()->SetRangeUser(0.8, 1.4);
      extrapolated_data_slope->Divide(extrapolated_mc_slope);
      extrapolated_data_slope->Draw();

      TLegend* leg4c = util::LabelFactory::createLegendWithOffset(3,0.175);
      leg4c->AddEntry(extrapolated_data_slope,"Data","PL");
      leg4c->AddEntry(extrapolated_mc_slope,"MC","L");
      leg4c->AddEntry(extrapolated_gen_slope,"PLI","L");
      
      leg4c->Draw();
      // cmsPrel(-1, false , 8);

      labelb->Draw("same");

      TString name3b;
      name3b = Form("ForwardExtrapolation/ExtrapolSlope_Eta%i" + suffix + ".eps", ieta);
      c4b->Print(name3b);
      // --------------------- //

      // reset color
      extrapolated_mc->SetLineColor(kBlack);

      // --------------------- //
      // plot chi2/ndf for extrapolations
      TCanvas *c2c = new TCanvas("c2c","",600,600);
      c2c->SetLogx();
      chi2_mc->GetYaxis()->SetTitle("#chi^{2}/ndf");
      chi2_mc->GetYaxis()->SetRangeUser(0., 8.0);
      chi2_mc->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
      chi2_mc->SetLineColor(kRed+1);
      chi2_mc->Draw();
      chi2_gen->SetLineColor(kBlue+1);
      chi2_gen->SetMarkerColor(kBlue+1);
      chi2_data->Draw("same");

      TLegend* leg3c = util::LabelFactory::createLegendWithOffset(2,0.175);
      leg3c->AddEntry(chi2_data,"Data","L");
      leg3c->AddEntry(chi2_mc,"MC","L");
      
      leg3c->Draw();
      // cmsPrel(-1, false , 8);

      TPaveText *labelc = util::LabelFactory::createPaveTextWithOffset(2,1.0,0.05);
      labelc->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
      labelc->AddText( Form("%0.1f #leq |#eta| #leq %0.1f", eta_bins[ieta], eta_bins[ieta+1]) );
      labelc->Draw("same");

      TString name2c;
      name2c = Form("ForwardExtrapolation/GoodnessOfFit_Eta%i" + suffix + ".eps", ieta);
      c2c->Print(name2c);
      // --------------------- //

      // --------------------------------------- //
      // calc data/mc ratio and fit with constant
      if(extrapolated_mc->FindFirstBinAbove(0) > 0 && extrapolated_mc_with_pli->FindFirstBinAbove(0) > 0) {
         TH1F* ratio = new TH1F(*extrapolated_mc);
         TH1F* ratio_with_pli = new TH1F(*extrapolated_mc);
         ratio->Divide(extrapolated_data, extrapolated_mc, 1, 1);
         ratio_with_pli->Divide(extrapolated_data_with_pli, extrapolated_mc_with_pli, 1, 1);
    
         TF1 *fit_const = new TF1("fit_const", "[0]", ratio->GetXaxis()->GetXmin(), ratio->GetXaxis()->GetXmax());
         double chi2_const = fit_const->GetChisquare();
         double ndf_const = fit_const->GetNDF();
         fit_const->SetFillColor(kGray);
         fit_const->SetParameters(0, 1.1);
         fit_const->SetParName(0, "const");
         // fit ratio without pli corrrection
         ratio->Fit("fit_const", "", "same");
         ratio->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
         ratio->GetYaxis()->SetRangeUser(0.8, 2.0);
         ratio->GetYaxis()->SetTitle("Data/MC ratio (const fit)");

         RatioVsEta->SetBinContent(ieta+1, ratio->GetFunction("fit_const")->GetParameter(0));
         RatioVsEta->SetBinError(ieta+1, ratio->GetFunction("fit_const")->GetParError(0));

         double fit_qual = chi2_const/ndf_const;
      
         //Create a histogram to hold the confidence intervals
         TH1D *hint = new TH1D("hint","Fitted function with .95 conf.band", 100, ratio->GetXaxis()->GetXmin(), ratio->GetXaxis()->GetXmax());
         hint->Sumw2();
         (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint, 0.68);
         //Now the "hint" histogram has the fitted function values as the 
         //bin contents and the confidence intervals as bin errors
         hint->SetMarkerStyle(1);
         hint->SetLineColor(kRed);
         hint->SetFillColor(kGray);
    
         // --------------------- //
         // plot ratio without pli correction
         TCanvas *c3 = new TCanvas("c3","",600,600);
         c3->SetLogx();
         ratio->Draw();
         hint->DrawClone("e4 same");
         ratio->Draw("same");
         label->Draw("same");

         TPaveText *labeld = util::LabelFactory::createPaveTextWithOffset(1,1.0,0.35);
         labeld->AddText( Form("#chi^{2}/ndf : %0.2f", fit_qual) );
         labeld->Draw("same");

         TLegend* leg4 = util::LabelFactory::createLegendWithOffset(2,0.175);
         leg4->AddEntry(ratio,"Measured Ratio","P");
         leg4->AddEntry(hint,"Constant Fit","LF");
         leg4->Draw("same");

         TString name3;
         name3 = Form("ForwardExtrapolation/ExtrapolRatio_Eta%i" + suffix + ".eps", ieta);
         c3->Print(name3);
         // --------------------- //

         // fit ratio with pli correction
         ratio_with_pli->Fit("fit_const", "", "same");
         ratio_with_pli->GetXaxis()->SetTitle("p_{T}^{ave} [GeV]");
         ratio_with_pli->GetYaxis()->SetRangeUser(0.8, 2.0);
         ratio_with_pli->GetYaxis()->SetTitle("Data/MC ratio (const fit)");

         RatioVsEta_with_pli->SetBinContent(ieta+1, ratio_with_pli->GetFunction("fit_const")->GetParameter(0));
         RatioVsEta_with_pli->SetBinError(ieta+1, ratio_with_pli->GetFunction("fit_const")->GetParError(0));

         chi2_const = fit_const->GetChisquare();
         ndf_const = fit_const->GetNDF();
         fit_qual = chi2_const/ndf_const;

         (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint, 0.68);
    
         // --------------------- //
         // plot ratio with pli correction --> final result
         TCanvas *c3b = new TCanvas("c3","",600,600);
         c3b->SetLogx();
         ratio_with_pli->Draw();
         hint->DrawClone("e4 same");
         ratio_with_pli->Draw("same");
         leg4->Draw("same");

         label->Draw("same");
     
         TPaveText *labele = util::LabelFactory::createPaveTextWithOffset(1,0.3,0.7);
         labele->AddText( Form("#chi^{2}/ndf : %0.2f", fit_qual) );
         labele->Draw("same");

         ratio_with_pli->Draw("same");

         TString name4;
         name4 = Form("ForwardExtrapolation/ExtrapolRatio_Eta%i_with_pli" + suffix + ".eps", ieta);
         c3b->Print(name4);
         // --------------------- //
      }
   }

   // --------------------- //
   // draw data/mc scaling factors vs. eta
   TCanvas *c4 = new TCanvas();
   RatioVsEta->GetXaxis()->SetRangeUser(0., 5.0);
   RatioVsEta->GetYaxis()->SetRangeUser(0.7, 1.3);
   RatioVsEta->GetXaxis()->SetTitle("|#eta|");
   RatioVsEta->GetYaxis()->SetTitle("Data /MC ratio (const fit)");
   RatioVsEta->Draw();
   c4->Print("ForwardExtrapolation/ScalingFactorsVsEta" + suffix + ".eps");
   // --------------------- //

   // --------------------- //
   TCanvas *c4b = new TCanvas();
   RatioVsEta_with_pli->GetXaxis()->SetRangeUser(0., 5.0);
   RatioVsEta_with_pli->GetYaxis()->SetRangeUser(0.7, 1.4);
   RatioVsEta_with_pli->GetXaxis()->SetTitle("|#eta|");
   RatioVsEta_with_pli->GetYaxis()->SetTitle("Data /MC ratio (const fit)");
   RatioVsEta_with_pli->Draw();
   c4b->Print("ForwardExtrapolation/ScalingFactorsVsEta_with_pli" + suffix + ".eps");
   // --------------------- //

   cout << "//----------------------------------------------//" << endl;
   cout << "Scaling factors without PLI: " << endl;
   cout << "Ratio eta1: " << RatioVsEta->GetBinContent(1) << " +- " << RatioVsEta->GetBinError(1) << endl;
   cout << "Ratio eta2: " << RatioVsEta->GetBinContent(2) << " +- " << RatioVsEta->GetBinError(2) << endl;
   cout << "Ratio eta3: " << RatioVsEta->GetBinContent(3) << " +- " << RatioVsEta->GetBinError(3) << endl;
   cout << "Ratio eta4: " << RatioVsEta->GetBinContent(4) << " +- " << RatioVsEta->GetBinError(4) << endl;
   cout << "Ratio eta5: " << RatioVsEta->GetBinContent(5) << " +- " << RatioVsEta->GetBinError(5) << endl;
   cout << "Ratio eta6: " << RatioVsEta->GetBinContent(6) << " +- " << RatioVsEta->GetBinError(6) << endl;
   cout << "Ratio eta7: " << RatioVsEta->GetBinContent(7) << " +- " << RatioVsEta->GetBinError(7) << endl;
   cout << "//----------------------------------------------//" << endl;

   cout << "//----------------------------------------------//" << endl;
   cout << "Scaling factors with PLI: " << endl;
   cout << "Ratio eta1: " << RatioVsEta_with_pli->GetBinContent(1) << " +- " << RatioVsEta_with_pli->GetBinError(1) << endl;
   cout << "Ratio eta2: " << RatioVsEta_with_pli->GetBinContent(2) << " +- " << RatioVsEta_with_pli->GetBinError(2) << endl;
   cout << "Ratio eta3: " << RatioVsEta_with_pli->GetBinContent(3) << " +- " << RatioVsEta_with_pli->GetBinError(3) << endl;
   cout << "Ratio eta4: " << RatioVsEta_with_pli->GetBinContent(4) << " +- " << RatioVsEta_with_pli->GetBinError(4) << endl;
   cout << "Ratio eta5: " << RatioVsEta_with_pli->GetBinContent(5) << " +- " << RatioVsEta_with_pli->GetBinError(5) << endl;
   cout << "Ratio eta6: " << RatioVsEta_with_pli->GetBinContent(6) << " +- " << RatioVsEta_with_pli->GetBinError(6) << endl;
   cout << "Ratio eta7: " << RatioVsEta_with_pli->GetBinContent(7) << " +- " << RatioVsEta_with_pli->GetBinError(7) << endl;
   cout << "//----------------------------------------------//" << endl;

   // --------------------- //
   // write final scaling factors to root file
   TFile* output = new TFile("ForwardExtrapolation/JER_RatioVsEta" + suffix + ".root", "RECREATE");
   RatioVsEta_with_pli->Write();

   output->Write();
   // --------------------- //
}



