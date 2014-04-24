#include "/afs/desy.de/user/k/kheine/xxl-af-cms/Kalibri/scripts/tdrstyle_mod.C"
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TString.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TLegend.h>

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

void CalcSysUncJER()
{
   setTDRStyle();
   // gStyle->SetHatchesSpacing(1.0);
   //  gStyle->SetHatchesLineWidth(10);

   // gStyle->SetErrorX();
   //  gROOT->ForceStyle();

   // ------------------------------------------------------------------ //
   // get files with nominal value and variations
   // nominal values
   TFile* file_nominal = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_v2.root", "READ");

   // jec down
   TFile* file_jecDOWN = new TFile("Extrapolation/JER_RatioVsEta_JECdown_final_v2.root", "READ");

   // jec up
   TFile* file_jecUP = new TFile("Extrapolation/JER_RatioVsEta_JECup_final_v2.root", "READ");

   // pli down
   TFile* file_pliDOWN = new TFile("Extrapolation/JER_RatioVsEta_PLIdown_final_v2.root", "READ");

   // pli up
   TFile* file_pliUP = new TFile("Extrapolation/JER_RatioVsEta_PLIup_final_v2.root", "READ");

   // alpha spectrum reweighting
   TFile* file_alpha = new TFile("Extrapolation/JER_RatioVsEta_ReweightAlphaSpectrum_final_v2.root", "READ");

   // PU reweighting
   TFile* file_pu = new TFile("Extrapolation/JER_RatioVsEta_MBXS73500_final_v2.root", "READ");

   // extrapolation range (low alpha)
   //  TFile* file_alpha_low = new TFile("Extrapolation/JER_RatioVsEta_FirstThreeAlphaPoints_final_v2.root", "READ");

   // extrapolation range (high alpha)
   //  TFile* file_alpha_high = new TFile("Extrapolation/JER_RatioVsEta_LastThreeAlphaPoints_final_v2.root", "READ");

   // non-gaussian tails
   TFile* file_tails = new TFile("Extrapolation/JER_RatioVsEta_95Truncation_final_v2.root", "READ");

   // gluon splitting reweighting
   TFile* file_gluon = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_GluonSplittingReweighting_v2.root", "READ");

   // additional alpha-bin 0.0-0.05
   TFile* file_addalpha = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_v2_NoMinPtCutForThirdJet_AddNewAlphaBin.root", "READ");

   // ------------------------------------------------------------------ //
   // define histos needed for uncertainties
   TH1F *RatioNominal = new TH1F();
   TH1F *RatioJECDown = new TH1F();
   TH1F *RatioJECUp = new TH1F();
   TH1F *RatioPLIDown = new TH1F();
   TH1F *RatioPLIUp = new TH1F();
   TH1F *RatioAlphaReweightUp = new TH1F();
   TH1F *RatioPUUp = new TH1F();
   //  TH1F *RatioAlphaLow = new TH1F();
   //  TH1F *RatioAlphaHigh = new TH1F();
   TH1F *RatioTailsDown = new TH1F();
   TH1F *RatioGluonUp = new TH1F();
   TH1F *RatioAddAlphaDown = new TH1F();

   // ------------------------------------------------------------------ //
   // get histos from files
   file_nominal->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioNominal);

   file_jecDOWN->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioJECDown);

   file_jecUP->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioJECUp);

   file_pliDOWN->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioPLIDown);

   file_pliUP->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioPLIUp);

   file_alpha->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioAlphaReweightUp);
  
   file_pu->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioPUUp);

   //  file_alpha_low->cd();
   //  gDirectory->GetObject("RatioVsEta_with_pli;1", RatioAlphaLow);

   //  file_alpha_high->cd();
   //  gDirectory->GetObject("RatioVsEta_with_pli;1", RatioAlphaHigh);

   file_tails->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioTailsDown);

   file_gluon->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioGluonUp);

   file_addalpha->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioAddAlphaDown);

   // ------------------------------------------------------------------ //
   // calc symmetric PU errors
   TH1F *RatioPUDown = new TH1F(*RatioPUUp);
   RatioPUDown->Reset();

   for (int i = 1; i < RatioPUUp->GetNbinsX()+1; i++) {
      double bin_content_pu = RatioPUUp->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_pu - bin_content_nominal;
      if(i == 2) {
         diff = 0.003;
         bin_content_pu = bin_content_nominal + 0.003;
      }

      //  cout << "PU i diff: " << i << "  "  << diff << endl;

      if( diff > 0 ) {
         RatioPUUp->SetBinContent(i, bin_content_pu);
         RatioPUDown->SetBinContent(i, bin_content_nominal - diff);
      }
      else if ( diff < 0 ) {
         RatioPUDown->SetBinContent(i, bin_content_pu);
         RatioPUUp->SetBinContent(i, bin_content_nominal - diff);
      }
   }

   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "PU Uncertainty : " << endl;
   cout << "Eta0 : " << (RatioPUUp->GetBinContent(1) - RatioNominal->GetBinContent(1))/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << (RatioPUUp->GetBinContent(2) - RatioNominal->GetBinContent(2))/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << (RatioPUUp->GetBinContent(3) - RatioNominal->GetBinContent(3))/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << (RatioPUUp->GetBinContent(4) - RatioNominal->GetBinContent(4))/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << (RatioPUUp->GetBinContent(5) - RatioNominal->GetBinContent(5))/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;

   // calc symmetric AlphaReweight errors
   TH1F *RatioAlphaReweightDown = new TH1F(*RatioAlphaReweightUp);
   RatioAlphaReweightDown->Reset();

   for (int i = 1; i < RatioAlphaReweightUp->GetNbinsX()+1; i++) {
      double bin_content_alpha = RatioAlphaReweightUp->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_alpha - bin_content_nominal;

      //  cout << "Alpha reweight i diff: " << i << "  "  << diff << endl;

      if( diff > 0 ) {
         RatioAlphaReweightUp->SetBinContent(i, bin_content_alpha);
         RatioAlphaReweightDown->SetBinContent(i, bin_content_nominal - diff);
      }
      else if ( diff < 0 ) {
         RatioAlphaReweightDown->SetBinContent(i, bin_content_alpha);
         RatioAlphaReweightUp->SetBinContent(i, bin_content_nominal - diff);
      }
   }

   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "Alpha Spectrum : " << endl;
   cout << "Eta0 : " << (RatioAlphaReweightUp->GetBinContent(1) - RatioNominal->GetBinContent(1))/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << (RatioAlphaReweightUp->GetBinContent(2) - RatioNominal->GetBinContent(2))/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << (RatioAlphaReweightUp->GetBinContent(3) - RatioNominal->GetBinContent(3))/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << (RatioAlphaReweightUp->GetBinContent(4) - RatioNominal->GetBinContent(4))/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << (RatioAlphaReweightUp->GetBinContent(5) - RatioNominal->GetBinContent(5))/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;

   // calc symmetric JEC errors
   for (int i = 1; i < RatioJECDown->GetNbinsX()+1; i++) {
      double bin_content_jec_low = RatioJECDown->GetBinContent(i);
      double bin_content_jec_high = RatioJECUp->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff_low = bin_content_jec_low - bin_content_nominal;
      double diff_high = bin_content_jec_high - bin_content_nominal;
      double diff = 0;

      //   cout << "JEC i diff low: " << i << "  "  << diff_low << endl;
      //  cout << "JEC i diff high: " << i << "  "  << diff_high << endl;

      if( diff_low > 0 && diff_high > 0 || diff_low < 0 && diff_high < 0) {
         if(diff_low < 0 && diff_high < 0) {
            if( diff_low > diff_high) diff = TMath::Abs(diff_high);
            else diff = TMath::Abs(diff_low);
         }
         else {
            if( diff_low < diff_high) diff = diff_high;
            else diff = diff_low;
         }
      }
      else {
         diff = 0.5 * (bin_content_jec_high - bin_content_jec_low);
      }

      //  cout << "JEC i diff : " << i << "  "  << diff << endl;

      if( diff > 0 ) {
         RatioJECDown->SetBinContent(i, bin_content_nominal - diff);
         RatioJECUp->SetBinContent(i, bin_content_nominal + diff);
      }
      else if( diff < 0) {
         RatioJECDown->SetBinContent(i, bin_content_nominal - TMath::Abs(diff));
         RatioJECUp->SetBinContent(i, bin_content_nominal + TMath::Abs(diff));
      }
   }

   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "JEC : " << endl;
   cout << "Eta0 : " << (RatioJECUp->GetBinContent(1) - RatioNominal->GetBinContent(1))/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << (RatioJECUp->GetBinContent(2) - RatioNominal->GetBinContent(2))/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << (RatioJECUp->GetBinContent(3) - RatioNominal->GetBinContent(3))/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << (RatioJECUp->GetBinContent(4) - RatioNominal->GetBinContent(4))/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << (RatioJECUp->GetBinContent(5) - RatioNominal->GetBinContent(5))/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;

   // calc symmetric PLI errors
   for (int i = 1; i < RatioPLIDown->GetNbinsX()+1; i++) {
      double bin_content_pli_low = RatioPLIDown->GetBinContent(i);
      double bin_content_pli_high = RatioPLIUp->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff_low = bin_content_pli_low - bin_content_nominal;
      double diff_high = bin_content_pli_high - bin_content_nominal;
      double diff = 0;

      //   cout << "PLI i diff low: " << i << "  "  << diff_low << endl;
      //   cout << "PLI i diff high: " << i << "  "  << diff_high << endl;

      if( diff_low > 0 && diff_high > 0 || diff_low < 0 && diff_high < 0) {
         if(diff_low < 0 && diff_high < 0) {
            if( diff_low > diff_high) diff = TMath::Abs(diff_high);
            else diff = TMath::Abs(diff_low);
         }
         else {
            if( diff_low < diff_high) diff = diff_high;
            else diff = diff_low;
         }
      }
      else {
         diff = 0.5 * (bin_content_pli_high - bin_content_pli_low);
      }

      if( diff > 0 ) {
         RatioPLIDown->SetBinContent(i, bin_content_nominal - diff);
         RatioPLIUp->SetBinContent(i, bin_content_nominal + diff);
      }
      else if( diff < 0) {
         RatioPLIDown->SetBinContent(i, bin_content_nominal - TMath::Abs(diff));
         RatioPLIUp->SetBinContent(i, bin_content_nominal + TMath::Abs(diff));
      }
   }

   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "PLI : " << endl;
   cout << "Eta0 : " << (RatioPLIUp->GetBinContent(1) - RatioNominal->GetBinContent(1))/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << (RatioPLIUp->GetBinContent(2) - RatioNominal->GetBinContent(2))/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << (RatioPLIUp->GetBinContent(3) - RatioNominal->GetBinContent(3))/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << (RatioPLIUp->GetBinContent(4) - RatioNominal->GetBinContent(4))/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << (RatioPLIUp->GetBinContent(5) - RatioNominal->GetBinContent(5))/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;


   // calc symmetric AlphaRange errors
   /*  for (int i = 1; i < RatioAlphaLow->GetNbinsX()+1; i++) {
      double bin_content_alpha_low = RatioAlphaLow->GetBinContent(i);
      double bin_content_alpha_high = RatioAlphaHigh->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      cout << "Alpha Range Low Eta: " << i << "  " << bin_content_alpha_low << " +- " << RatioAlphaLow->GetBinError(i) << endl;
      cout << "Alpha Range High Eta: " << i << "  " << bin_content_alpha_high << " +- " << RatioAlphaHigh->GetBinError(i) << endl;
 
      double diff_low = bin_content_alpha_low - bin_content_nominal;
      double diff_high = bin_content_alpha_high - bin_content_nominal;
      double diff = 0;

      //  cout << "Alpha Range i diff low: " << i << "  "  << diff_low << endl;
      //   cout << "Alpha Range i diff high: " << i << "  "  << diff_high << endl;

      if( diff_low > 0 && diff_high > 0 || diff_low < 0 && diff_high < 0) {
         if(diff_low < 0 && diff_high < 0) {
            if( diff_low > diff_high) diff = TMath::Abs(diff_high);
            else diff = TMath::Abs(diff_low);
         }
         else {
            if( diff_low < diff_high) diff = diff_high;
            else diff = diff_low;
         }
      }
      else {
         diff = 0.5 * (bin_content_alpha_high - bin_content_alpha_low);
      }

      if( diff > 0 ) {
         RatioAlphaLow->SetBinContent(i, bin_content_nominal - diff);
         RatioAlphaHigh->SetBinContent(i, bin_content_nominal + diff);
      }
      else if( diff < 0 ) {
         RatioAlphaLow->SetBinContent(i, bin_content_nominal - TMath::Abs(diff));
         RatioAlphaHigh->SetBinContent(i, bin_content_nominal + TMath::Abs(diff));
      }
   }

   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "Alpha Range: " << endl;
   cout << "Eta0 : " << (RatioAlphaHigh->GetBinContent(1) - RatioNominal->GetBinContent(1))/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << (RatioAlphaHigh->GetBinContent(2) - RatioNominal->GetBinContent(2))/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << (RatioAlphaHigh->GetBinContent(3) - RatioNominal->GetBinContent(3))/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << (RatioAlphaHigh->GetBinContent(4) - RatioNominal->GetBinContent(4))/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << (RatioAlphaHigh->GetBinContent(5) - RatioNominal->GetBinContent(5))/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;*/


   // calc symmetric Tail errors
   TH1F *RatioTailsUp = new TH1F(*RatioTailsDown);
   RatioTailsUp->Reset();

   for (int i = 1; i < RatioTailsUp->GetNbinsX()+1; i++) {
      double bin_content_tails = RatioTailsDown->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_tails - bin_content_nominal;
      if(i == 2) {
         diff = 0.003;
         bin_content_tails = bin_content_nominal + 0.003;
      }

      //  cout << "PU i diff: " << i << "  "  << diff << endl;

      if( diff > 0 ) {
         RatioTailsUp->SetBinContent(i, bin_content_tails);
         RatioTailsDown->SetBinContent(i, bin_content_nominal - diff);
      }
      else if ( diff < 0 ) {
         RatioTailsDown->SetBinContent(i, bin_content_tails);
         RatioTailsUp->SetBinContent(i, bin_content_nominal - diff);
      }
   }

   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "Tail Uncertainty : " << endl;
   cout << "Eta0 : " << (RatioTailsUp->GetBinContent(1) - RatioNominal->GetBinContent(1))/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << (RatioTailsUp->GetBinContent(2) - RatioNominal->GetBinContent(2))/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << (RatioTailsUp->GetBinContent(3) - RatioNominal->GetBinContent(3))/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << (RatioTailsUp->GetBinContent(4) - RatioNominal->GetBinContent(4))/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << (RatioTailsUp->GetBinContent(5) - RatioNominal->GetBinContent(5))/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;

   // calc symmetric GluonReweighting errors
   TH1F *RatioGluonDown = new TH1F(*RatioGluonUp);
   RatioGluonDown->Reset();

   for (int i = 1; i < RatioGluonUp->GetNbinsX()+1; i++) {
      double bin_content_alpha = RatioGluonUp->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_alpha - bin_content_nominal;

      //  cout << "Alpha reweight i diff: " << i << "  "  << diff << endl;

      if( diff > 0 ) {
         RatioGluonUp->SetBinContent(i, bin_content_alpha);
         RatioGluonDown->SetBinContent(i, bin_content_nominal - diff);
      }
      else if ( diff < 0 ) {
         RatioGluonDown->SetBinContent(i, bin_content_alpha);
         RatioGluonUp->SetBinContent(i, bin_content_nominal - diff);
      }
   }

   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "Gluon Splitting Reweighting : " << endl;
   cout << "Eta0 : " << (RatioGluonUp->GetBinContent(1) - RatioNominal->GetBinContent(1))/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << (RatioGluonUp->GetBinContent(2) - RatioNominal->GetBinContent(2))/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << (RatioGluonUp->GetBinContent(3) - RatioNominal->GetBinContent(3))/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << (RatioGluonUp->GetBinContent(4) - RatioNominal->GetBinContent(4))/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << (RatioGluonUp->GetBinContent(5) - RatioNominal->GetBinContent(5))/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;

   // calc symmetric errors for additional alpha bin --> alpha range
   TH1F *RatioAddAlphaUp = new TH1F(*RatioAddAlphaDown);
   RatioAddAlphaUp->Reset();

   for (int i = 1; i < RatioAddAlphaUp->GetNbinsX()+1; i++) {
      double bin_content_alpha = RatioAddAlphaDown->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_alpha - bin_content_nominal;

      //  cout << "Alpha reweight i diff: " << i << "  "  << diff << endl;

      if( diff > 0 ) {
         RatioAddAlphaUp->SetBinContent(i, bin_content_alpha);
         RatioAddAlphaDown->SetBinContent(i, bin_content_nominal - diff);
      }
      else if ( diff < 0 ) {
         RatioAddAlphaDown->SetBinContent(i, bin_content_alpha);
         RatioAddAlphaUp->SetBinContent(i, bin_content_nominal - diff);
      }
   }

   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "Alpha Range (Add low Alpha-Bin) : " << endl;
   cout << "Eta0 : " << (RatioAddAlphaUp->GetBinContent(1) - RatioNominal->GetBinContent(1))/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << (RatioAddAlphaUp->GetBinContent(2) - RatioNominal->GetBinContent(2))/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << (RatioAddAlphaUp->GetBinContent(3) - RatioNominal->GetBinContent(3))/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << (RatioAddAlphaUp->GetBinContent(4) - RatioNominal->GetBinContent(4))/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << (RatioAddAlphaUp->GetBinContent(5) - RatioNominal->GetBinContent(5))/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;

   // calc symmetric errors for possible pt-trend
   TH1F *RatioPtTrendDown = new TH1F(*RatioAddAlphaDown);
   TH1F *RatioPtTrendUp = new TH1F(*RatioAddAlphaDown);
   RatioPtTrendDown->Reset();
   RatioPtTrendUp->Reset();

   for (int i = 1; i < RatioPtTrendDown->GetNbinsX()+1; i++) {
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      RatioPtTrendUp->SetBinContent(i, bin_content_nominal + 0.02*bin_content_nominal);
      RatioPtTrendDown->SetBinContent(i, bin_content_nominal - 0.02*bin_content_nominal);
    
   }

   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "Pt Trend : " << endl;
   cout << "Eta0 : " << (RatioPtTrendUp->GetBinContent(1) - RatioNominal->GetBinContent(1))/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << (RatioPtTrendUp->GetBinContent(2) - RatioNominal->GetBinContent(2))/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << (RatioPtTrendUp->GetBinContent(3) - RatioNominal->GetBinContent(3))/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << (RatioPtTrendUp->GetBinContent(4) - RatioNominal->GetBinContent(4))/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << (RatioPtTrendUp->GetBinContent(5) - RatioNominal->GetBinContent(5))/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;

   // ------------------------------------------------------------------ //
   // remove hist errors for unc. variations
   for (int i = 1; i < RatioNominal->GetNbinsX()+1; i++) {
      RatioJECDown->SetBinError(i, 0.0001);
      RatioJECUp->SetBinError(i, 0.0001);
      RatioPLIDown->SetBinError(i, 0.0001);
      RatioPLIUp->SetBinError(i, 0.0001);
      RatioAlphaReweightDown->SetBinError(i, 0.0001);
      RatioAlphaReweightUp->SetBinError(i, 0.0001);
      RatioPUDown->SetBinError(i, 0.0001);
      RatioPUUp->SetBinError(i, 0.0001);
      //    RatioAlphaLow->SetBinError(i, 0.0001);
      //   RatioAlphaHigh->SetBinError(i, 0.0001);
      RatioTailsUp->SetBinError(i, 0.0001);
      RatioTailsDown->SetBinError(i, 0.0001);
   }

   // ------------------------------------------------------------------ //
   // calc lower bounds
   double eta0_sysUncDown = TMath::Sqrt(pow(RatioNominal->GetBinContent(1) - 
                                            RatioJECDown->GetBinContent(1),2) + 
                                        pow(RatioNominal->GetBinContent(1) - 
                                            RatioPLIDown->GetBinContent(1),2) + 
                                        pow(RatioNominal->GetBinContent(1) - 
                                            RatioAlphaReweightDown->GetBinContent(1),2) +
                                        //    pow(RatioNominal->GetBinContent(1) - 
                                        //     RatioAlphaLow->GetBinContent(1),2) + 
                                        pow(RatioNominal->GetBinContent(1) - 
                                            RatioPUDown->GetBinContent(1),2) +
                                        pow(RatioNominal->GetBinContent(1) - 
                                            RatioGluonDown->GetBinContent(1),2) +
                                        pow(RatioNominal->GetBinContent(1) - 
                                            RatioAddAlphaDown->GetBinContent(1),2) +
                                        pow(RatioNominal->GetBinContent(1) - 
                                            RatioPtTrendDown->GetBinContent(1),2) +
                                        pow(RatioNominal->GetBinContent(1) - 
                                            RatioTailsDown->GetBinContent(1),2));
   double eta1_sysUncDown = TMath::Sqrt(pow(RatioNominal->GetBinContent(2) - 
                                            RatioJECDown->GetBinContent(2),2) + 
                                        pow(RatioNominal->GetBinContent(2) - 
                                            RatioPLIDown->GetBinContent(2),2) +
                                        pow(RatioNominal->GetBinContent(2) - 
                                            RatioAlphaReweightDown->GetBinContent(2),2) +
                                        //   pow(RatioNominal->GetBinContent(2) - 
                                        //       RatioAlphaLow->GetBinContent(2),2) + 
                                        pow(RatioNominal->GetBinContent(2) - 
                                            RatioPUDown->GetBinContent(2),2) +
                                        pow(RatioNominal->GetBinContent(2) - 
                                            RatioGluonDown->GetBinContent(2),2) +
                                        pow(RatioNominal->GetBinContent(2) - 
                                            RatioAddAlphaDown->GetBinContent(2),2) +
                                        pow(RatioNominal->GetBinContent(2) - 
                                            RatioPtTrendDown->GetBinContent(2),2) +
                                        pow(RatioNominal->GetBinContent(2) - 
                                            RatioTailsDown->GetBinContent(2),2));
   double eta2_sysUncDown = TMath::Sqrt(pow(RatioNominal->GetBinContent(3) - 
                                            RatioJECDown->GetBinContent(3),2) + 
                                        pow(RatioNominal->GetBinContent(3) - 
                                            RatioPLIDown->GetBinContent(3),2) + 
                                        pow(RatioNominal->GetBinContent(3) - 
                                            RatioAlphaReweightDown->GetBinContent(3),2) +
                                        //   pow(RatioNominal->GetBinContent(3) - 
                                        //       RatioAlphaLow->GetBinContent(3),2) + 
                                        pow(RatioNominal->GetBinContent(3) - 
                                            RatioPUDown->GetBinContent(3),2) +
                                        pow(RatioNominal->GetBinContent(3) - 
                                            RatioGluonDown->GetBinContent(3),2) +
                                        pow(RatioNominal->GetBinContent(3) - 
                                            RatioAddAlphaDown->GetBinContent(3),2) +
                                        pow(RatioNominal->GetBinContent(3) - 
                                            RatioPtTrendDown->GetBinContent(3),2) +
                                        pow(RatioNominal->GetBinContent(3) - 
                                            RatioTailsDown->GetBinContent(3),2));
   double eta3_sysUncDown = TMath::Sqrt(pow(RatioNominal->GetBinContent(4) - 
                                            RatioJECDown->GetBinContent(4),2) + 
                                        pow(RatioNominal->GetBinContent(4) - 
                                            RatioPLIDown->GetBinContent(4),2) + 
                                        pow(RatioNominal->GetBinContent(4) - 
                                            RatioAlphaReweightDown->GetBinContent(4),2) +
                                        //   pow(RatioNominal->GetBinContent(4) - 
                                        //        RatioAlphaLow->GetBinContent(4),2) + 
                                        pow(RatioNominal->GetBinContent(4) - 
                                            RatioPUDown->GetBinContent(4),2) +
                                        pow(RatioNominal->GetBinContent(4) - 
                                            RatioGluonDown->GetBinContent(4),2) +
                                        pow(RatioNominal->GetBinContent(4) - 
                                            RatioAddAlphaDown->GetBinContent(4),2) +
                                        pow(RatioNominal->GetBinContent(4) - 
                                            RatioPtTrendDown->GetBinContent(4),2) +
                                        pow(RatioNominal->GetBinContent(4) - 
                                            RatioTailsDown->GetBinContent(4),2));
   double eta4_sysUncDown = TMath::Sqrt(pow(RatioNominal->GetBinContent(5) - 
                                            RatioJECDown->GetBinContent(5),2) + 
                                        pow(RatioNominal->GetBinContent(5) - 
                                            RatioPLIDown->GetBinContent(5),2) + 
                                        pow(RatioNominal->GetBinContent(5) - 
                                            RatioAlphaReweightDown->GetBinContent(5),2) +
                                        //     pow(RatioNominal->GetBinContent(5) - 
                                        //   RatioAlphaLow->GetBinContent(5),2) + 
                                        pow(RatioNominal->GetBinContent(5) - 
                                            RatioPUDown->GetBinContent(5),2) +
                                        pow(RatioNominal->GetBinContent(5) - 
                                            RatioGluonDown->GetBinContent(5),2) +
                                        pow(RatioNominal->GetBinContent(5) - 
                                            RatioAddAlphaDown->GetBinContent(5),2) +
                                        pow(RatioNominal->GetBinContent(5) - 
                                            RatioPtTrendDown->GetBinContent(5),2) +
                                        pow(RatioNominal->GetBinContent(5) - 
                                            RatioTailsDown->GetBinContent(5),2));

   // ------------------------------------------------------------------ //
   // calc upper bounds
   double eta0_sysUncUp = TMath::Sqrt(pow(RatioJECUp->GetBinContent(1) - 
                                          RatioNominal->GetBinContent(1),2) + 
                                      pow(RatioPLIUp->GetBinContent(1) - 
                                          RatioNominal->GetBinContent(1),2) + 
                                      pow(RatioAlphaReweightUp->GetBinContent(1) - 
                                          RatioNominal->GetBinContent(1),2) +
                                      //    pow(RatioAlphaHigh->GetBinContent(1) - 
                                      //        RatioNominal->GetBinContent(1),2) + 
                                      pow(RatioPUUp->GetBinContent(1) - 
                                          RatioNominal->GetBinContent(1),2) +
                                      pow(RatioGluonUp->GetBinContent(1) - 
                                          RatioNominal->GetBinContent(1),2) +
                                      pow(RatioAddAlphaUp->GetBinContent(1) - 
                                          RatioNominal->GetBinContent(1),2) +
                                      pow(RatioPtTrendUp->GetBinContent(1) - 
                                          RatioNominal->GetBinContent(1),2) +
                                      pow(RatioTailsUp->GetBinContent(1) - 
                                          RatioNominal->GetBinContent(1),2));
   double eta1_sysUncUp = TMath::Sqrt(pow(RatioJECUp->GetBinContent(2) - 
                                          RatioNominal->GetBinContent(2),2) + 
                                      pow(RatioPLIUp->GetBinContent(2) - 
                                          RatioNominal->GetBinContent(2),2) + 
                                      pow(RatioAlphaReweightUp->GetBinContent(2) - 
                                          RatioNominal->GetBinContent(2),2) +
                                      //    pow(RatioAlphaHigh->GetBinContent(2) - 
                                      //        RatioNominal->GetBinContent(2),2) + 
                                      pow(RatioPUUp->GetBinContent(2) - 
                                          RatioNominal->GetBinContent(2),2) + 
                                      pow(RatioGluonUp->GetBinContent(2) - 
                                          RatioNominal->GetBinContent(2),2) +
                                      pow(RatioAddAlphaUp->GetBinContent(2) - 
                                          RatioNominal->GetBinContent(2),2) +
                                      pow(RatioPtTrendUp->GetBinContent(2) - 
                                          RatioNominal->GetBinContent(2),2) +
                                      pow(RatioTailsUp->GetBinContent(2) - 
                                          RatioNominal->GetBinContent(2),2)); 
   double eta2_sysUncUp = TMath::Sqrt(pow(RatioJECUp->GetBinContent(3) - 
                                          RatioNominal->GetBinContent(3),2) + 
                                      pow(RatioPLIUp->GetBinContent(3) - 
                                          RatioNominal->GetBinContent(3),2) + 
                                      pow(RatioAlphaReweightUp->GetBinContent(3) - 
                                          RatioNominal->GetBinContent(3),2) +
                                      //   pow(RatioAlphaHigh->GetBinContent(3) - 
                                      //       RatioNominal->GetBinContent(3),2) + 
                                      pow(RatioPUUp->GetBinContent(3) - 
                                          RatioNominal->GetBinContent(3),2) +
                                      pow(RatioGluonUp->GetBinContent(3) - 
                                          RatioNominal->GetBinContent(3),2) +
                                      pow(RatioAddAlphaUp->GetBinContent(3) - 
                                          RatioNominal->GetBinContent(3),2) +
                                      pow(RatioPtTrendUp->GetBinContent(3) - 
                                          RatioNominal->GetBinContent(3),2) +
                                      pow(RatioTailsUp->GetBinContent(3) - 
                                          RatioNominal->GetBinContent(3),2));
   double eta3_sysUncUp = TMath::Sqrt(pow(RatioJECUp->GetBinContent(4) - 
                                          RatioNominal->GetBinContent(4),2) + 
                                      pow(RatioPLIUp->GetBinContent(4) - 
                                          RatioNominal->GetBinContent(4),2) + 
                                      pow(RatioAlphaReweightUp->GetBinContent(4) - 
                                          RatioNominal->GetBinContent(4),2) +
                                      //   pow(RatioAlphaHigh->GetBinContent(4) - 
                                      //       RatioNominal->GetBinContent(4),2) + 
                                      pow(RatioPUUp->GetBinContent(4) - 
                                          RatioNominal->GetBinContent(4),2) +
                                      pow(RatioGluonUp->GetBinContent(4) - 
                                          RatioNominal->GetBinContent(4),2) +
                                      pow(RatioAddAlphaUp->GetBinContent(4) - 
                                          RatioNominal->GetBinContent(4),2) +
                                      pow(RatioPtTrendUp->GetBinContent(4) - 
                                          RatioNominal->GetBinContent(4),2) +
                                      pow(RatioTailsUp->GetBinContent(4) - 
                                          RatioNominal->GetBinContent(4),2));
   double eta4_sysUncUp = TMath::Sqrt(pow(RatioJECUp->GetBinContent(5) - 
                                          RatioNominal->GetBinContent(5),2) + 
                                      pow(RatioPLIUp->GetBinContent(5) - 
                                          RatioNominal->GetBinContent(5),2) + 
                                      pow(RatioAlphaReweightUp->GetBinContent(5) - 
                                          RatioNominal->GetBinContent(5),2) +
                                      //     pow(RatioAlphaHigh->GetBinContent(5) - 
                                      //         RatioNominal->GetBinContent(5),2) + 
                                      pow(RatioPUUp->GetBinContent(5) - 
                                          RatioNominal->GetBinContent(5),2) +
                                      pow(RatioGluonUp->GetBinContent(5) - 
                                          RatioNominal->GetBinContent(5),2) +
                                      pow(RatioAddAlphaUp->GetBinContent(5) - 
                                          RatioNominal->GetBinContent(5),2) +
                                      pow(RatioPtTrendUp->GetBinContent(5) - 
                                          RatioNominal->GetBinContent(5),2) +
                                      pow(RatioTailsUp->GetBinContent(5) - 
                                          RatioNominal->GetBinContent(5),2));


   // ------------------------------------------------------------------ //
   // show values on screen total syst. uncertainty in %
   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "Total syst. uncertainty (in '%') : " << endl;
   cout << "Eta0 : " << eta0_sysUncUp/RatioNominal->GetBinContent(1) << endl;
   cout << "Eta1 : " << eta1_sysUncUp/RatioNominal->GetBinContent(2) << endl;
   cout << "Eta2 : " << eta2_sysUncUp/RatioNominal->GetBinContent(3) << endl;
   cout << "Eta3 : " << eta3_sysUncUp/RatioNominal->GetBinContent(4) << endl;
   cout << "Eta4 : " << eta4_sysUncUp/RatioNominal->GetBinContent(5) << endl;
   cout << " // ------------------------------------------------------------------ // " << endl;

   // ------------------------------------------------------------------ //
   // show values on screen with stat. uncertainty
   cout << "eta 0 stat.: " << RatioNominal->GetBinContent(1) << " +- " << RatioNominal->GetBinError(1) << endl;
   cout << "eta 1 stat.: " << RatioNominal->GetBinContent(2) << " +- " << RatioNominal->GetBinError(2) << endl;
   cout << "eta 2 stat.: " << RatioNominal->GetBinContent(3) << " +- " << RatioNominal->GetBinError(3) << endl;
   cout << "eta 3 stat.: " << RatioNominal->GetBinContent(4) << " +- " << RatioNominal->GetBinError(4) << endl;
   cout << "eta 4 stat.: " << RatioNominal->GetBinContent(5) << " +- " << RatioNominal->GetBinError(5) << endl;

   // ------------------------------------------------------------------ //

   // ------------------------------------------------------------------ //
   // show values on screen with sys. uncertainty
   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "eta 0 sys.: " << RatioNominal->GetBinContent(1) << " + " << eta0_sysUncUp << " - " << eta0_sysUncDown << endl;
   cout << "eta 1 sys.: " << RatioNominal->GetBinContent(2) << " + " << eta1_sysUncUp << " - " << eta1_sysUncDown << endl;
   cout << "eta 2 sys.: " << RatioNominal->GetBinContent(3) << " + " << eta2_sysUncUp << " - " << eta2_sysUncDown << endl;
   cout << "eta 3 sys.: " << RatioNominal->GetBinContent(4) << " + " << eta3_sysUncUp << " - " << eta3_sysUncDown << endl;
   cout << "eta 4 sys.: " << RatioNominal->GetBinContent(5) << " + " << eta4_sysUncUp << " - " << eta4_sysUncDown << endl;

   // ------------------------------------------------------------------ //
   // calc total errors
   double eta0_TotalUncDown = TMath::Sqrt(pow(RatioNominal->GetBinError(1),2) + pow(eta0_sysUncDown,2));
   double eta1_TotalUncDown = TMath::Sqrt(pow(RatioNominal->GetBinError(2),2) + pow(eta1_sysUncDown,2));
   double eta2_TotalUncDown = TMath::Sqrt(pow(RatioNominal->GetBinError(3),2) + pow(eta2_sysUncDown,2));
   double eta3_TotalUncDown = TMath::Sqrt(pow(RatioNominal->GetBinError(4),2) + pow(eta3_sysUncDown,2));
   double eta4_TotalUncDown = TMath::Sqrt(pow(RatioNominal->GetBinError(5),2) + pow(eta4_sysUncDown,2));  

   double eta0_TotalUncUp = TMath::Sqrt(pow(RatioNominal->GetBinError(1),2) + pow(eta0_sysUncUp,2));
   double eta1_TotalUncUp = TMath::Sqrt(pow(RatioNominal->GetBinError(2),2) + pow(eta1_sysUncUp,2));
   double eta2_TotalUncUp = TMath::Sqrt(pow(RatioNominal->GetBinError(3),2) + pow(eta2_sysUncUp,2));
   double eta3_TotalUncUp = TMath::Sqrt(pow(RatioNominal->GetBinError(4),2) + pow(eta3_sysUncUp,2));
   double eta4_TotalUncUp = TMath::Sqrt(pow(RatioNominal->GetBinError(5),2) + pow(eta4_sysUncUp,2));  

   // ------------------------------------------------------------------ //
   // fill total uncertainties to TGraph
   TGraphAsymmErrors *Res_2012_total = new TGraphAsymmErrors(RatioNominal);
   Res_2012_total->SetPointError(0, 0.25, 0.25, eta0_TotalUncDown, eta0_TotalUncUp);
   Res_2012_total->SetPointError(1, 0.3, 0.3, eta1_TotalUncDown, eta1_TotalUncUp);
   Res_2012_total->SetPointError(2, 0.3, 0.3, eta2_TotalUncDown, eta2_TotalUncUp);
   Res_2012_total->SetPointError(3, 0.3, 0.3, eta3_TotalUncDown, eta3_TotalUncUp);
   Res_2012_total->SetPointError(4, 1.45, 1.45, eta4_TotalUncDown, eta4_TotalUncUp);

   // ------------------------------------------------------------------ //
   // show values on screen with tot. uncertainty
   cout << " // ------------------------------------------------------------------ // " << endl;
   cout << "eta 0 tot.: " << RatioNominal->GetBinContent(1) << " + " << eta0_TotalUncUp << " - " << eta0_TotalUncDown << endl;
   cout << "eta 1 tot.: " << RatioNominal->GetBinContent(2) << " + " << eta1_TotalUncUp << " - " << eta1_TotalUncDown << endl;
   cout << "eta 2 tot.: " << RatioNominal->GetBinContent(3) << " + " << eta2_TotalUncUp << " - " << eta2_TotalUncDown << endl;
   cout << "eta 3 tot.: " << RatioNominal->GetBinContent(4) << " + " << eta3_TotalUncUp << " - " << eta3_TotalUncDown << endl;
   cout << "eta 4 tot.: " << RatioNominal->GetBinContent(5) << " + " << eta4_TotalUncUp << " - " << eta4_TotalUncDown << endl;

   // ------------------------------------------------------------------ //
   // plot nominal value + separate uncertainties 
   TCanvas *c = new TCanvas();
   RatioNominal->GetXaxis()->SetTitle("|#eta|");
   RatioNominal->GetYaxis()->SetTitle("Data/MC ratio (const fit)");
   RatioNominal->GetXaxis()->SetRangeUser(0., 5.);
   RatioNominal->GetYaxis()->SetRangeUser(0.8, 1.5);
   RatioNominal->SetMarkerSize(1.);
   RatioNominal->Draw();
   RatioJECDown->GetXaxis()->SetRangeUser(0., 5.);
   RatioJECDown->SetMarkerSize(1.);
   RatioJECDown->SetMarkerStyle(22);
   RatioJECDown->SetMarkerColor(kRed+1);
   RatioJECDown->Draw("same");
   RatioJECUp->GetXaxis()->SetRangeUser(0., 5.);
   RatioJECUp->SetMarkerSize(1.);
   RatioJECUp->SetMarkerStyle(22);
   RatioJECUp->SetMarkerColor(kRed+1);
   RatioJECUp->Draw("same");
   RatioPLIDown->GetXaxis()->SetRangeUser(0., 5.);
   RatioPLIDown->SetMarkerSize(1.);
   RatioPLIDown->SetMarkerStyle(24);
   RatioPLIDown->SetMarkerColor(kGreen+1);
   RatioPLIDown->Draw("same");
   RatioPLIUp->GetXaxis()->SetRangeUser(0., 5.);
   RatioPLIUp->SetMarkerSize(1.);
   RatioPLIUp->SetMarkerStyle(24);
   RatioPLIUp->SetMarkerColor(kGreen+1);
   RatioPLIUp->Draw("same");
   RatioAlphaReweightDown->GetXaxis()->SetRangeUser(0., 5.);
   RatioAlphaReweightDown->SetMarkerSize(1.);
   RatioAlphaReweightDown->SetMarkerStyle(21);
   RatioAlphaReweightDown->SetMarkerColor(kSpring+4);
   RatioAlphaReweightDown->Draw("same");
   RatioAlphaReweightUp->GetXaxis()->SetRangeUser(0., 5.);
   RatioAlphaReweightUp->SetMarkerSize(1.);
   RatioAlphaReweightUp->SetMarkerStyle(21);
   RatioAlphaReweightUp->SetMarkerColor(kSpring+4);
   RatioAlphaReweightUp->Draw("same");
   RatioPUDown->GetXaxis()->SetRangeUser(0., 5.);
   RatioPUDown->SetMarkerSize(1.);
   RatioPUDown->SetMarkerStyle(23);
   RatioPUDown->SetMarkerColor(kCyan+1);
   RatioPUDown->Draw("same");
   RatioPUUp->GetXaxis()->SetRangeUser(0., 5.);
   RatioPUUp->SetMarkerSize(1.);
   RatioPUUp->SetMarkerStyle(23);
   RatioPUUp->SetMarkerColor(kCyan+1);
   RatioPUUp->Draw("same");
   //  RatioAlphaLow->GetXaxis()->SetRangeUser(0., 5.);
   //RatioAlphaLow->SetMarkerSize(1.);
   //RatioAlphaLow->SetMarkerStyle(24);
   //RatioAlphaLow->SetMarkerColor(kMagenta+1);
   //RatioAlphaLow->Draw("same");
   //RatioAlphaHigh->GetXaxis()->SetRangeUser(0., 5.);
   //RatioAlphaHigh->SetMarkerSize(1.);
   //RatioAlphaHigh->SetMarkerStyle(24);
   //RatioAlphaHigh->SetMarkerColor(kMagenta+1);
   //RatioAlphaHigh->Draw("same");
   RatioTailsDown->GetXaxis()->SetRangeUser(0., 5.);
   RatioTailsDown->SetMarkerSize(1.);
   RatioTailsDown->SetMarkerStyle(25);
   RatioTailsDown->SetMarkerColor(kPink+1);
   RatioTailsDown->Draw("same");
   RatioTailsUp->GetXaxis()->SetRangeUser(0., 5.);
   RatioTailsUp->SetMarkerSize(1.);
   RatioTailsUp->SetMarkerStyle(25);
   RatioTailsUp->SetMarkerColor(kPink+1);
   RatioTailsUp->Draw("same");
   RatioNominal->Draw("same");

   TLegend *leg = new TLegend(0.17, 0.17, 0.54, 0.4);
   leg->SetBorderSize(0);
   // leg->SetBorderMode(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.035);

   leg->AddEntry(RatioNominal,"Nominal JER + stat. unc.", "pfl");
   leg->AddEntry(RatioJECDown,"JEC", "P");
   leg->AddEntry(RatioPLIDown,"PLI", "P");
   leg->AddEntry(RatioAlphaReweightDown,"Alpha Spectrum", "P");
   //  leg->AddEntry(RatioAlphaLow,"Alpha Range", "P");
   leg->AddEntry(RatioPUUp, "PU reweighting", "P");
   leg->AddEntry(RatioTailsUp, "Non-Gaussian tails", "P");
  
   leg->Draw("same");

   c->Print("Results/JER_2012_uncertainties_final_v2.eps");
   c->Print("Results/JER_2012_uncertainties_final_v2.png");


   // ------------------------------------------------------------------ //
   // plot nominal values with total unc.
   TH1F *dummy = new TH1F(*RatioNominal);
   dummy->Reset();

   TCanvas *c2 = new TCanvas();
   dummy->GetXaxis()->SetTitle("|#eta|");
   dummy->GetYaxis()->SetTitle("Data/MC ratio (const fit)");
   dummy->GetXaxis()->SetRangeUser(0., 5.0);
   dummy->GetYaxis()->SetRangeUser(0.9, 1.7);
   dummy->Draw();
  
   Res_2012_total->SetMarkerStyle(20);
   Res_2012_total->SetMarkerSize(1.4);
   //   Res_2012_total->SetFillColor(kYellow-3);
   Res_2012_total->SetFillColor(kCyan-2);
   Res_2012_total->SetLineColor(kCyan-2);
   Res_2012_total->SetFillStyle(1001);
   //  Res_2012_total->SetFillStyle(3001);
   //  Res_2012_total->SetLineColor(kYellow-3);
   Res_2012_total->DrawClone("2psame");

   Res_2012_total->SetPointError(0, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(1, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(2, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(3, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(4, 0., 0., 0., 0.);

   TPaveText *label = util::LabelFactory::createPaveTextWithOffset(1,1.05,0.01);
   label->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
   label->Draw("same");

   TLegend* leg2 = util::LabelFactory::createLegendColWithOffset(1, 1.05, 0.07);
   leg2->AddEntry(Res_2012_total,"JER 2012 Dijets (with tot. unc.)","PF");
    
   leg2->Draw("same");

   Res_2012_total->GetXaxis()->SetRangeUser(0., 5);
   Res_2012_total->Draw("psame");
   gPad->RedrawAxis();

   c2->Print("Results/JER_2012_final_v2.eps");
   c2->Print("Results/JER_2012_final_v2.png");
   c2->Print("Results/JER_2012_final_v2.pdf");

   // ------------------------------------------------------------------ //
   // plot nominal values with total unc. + 2011 comparison
   TH1F *Res_2011Final = new TH1F(*RatioNominal);
   Res_2011Final->Reset();
   Res_2011Final->SetBinContent(1, 1.052);
   Res_2011Final->SetBinContent(2, 1.057);
   Res_2011Final->SetBinContent(3, 1.096);
   Res_2011Final->SetBinContent(4, 1.134);
   Res_2011Final->SetBinContent(5, 1.288);

   TGraphAsymmErrors *Res_2011 = new TGraphAsymmErrors(Res_2011Final);
   Res_2011->SetPointError(0, 025., 0.25, 0.063, 0.062);
   Res_2011->SetPointError(1, 0.3, 0.3, 0.057, 0.056);
   Res_2011->SetPointError(2, 0.3, 0.3, 0.065, 0.064);
   Res_2011->SetPointError(3, 0.3, 0.3, 0.094, 0.092);
   Res_2011->SetPointError(4, 1.35, 1.35, 0.200, 0.199);
 
   TCanvas *c3 = new TCanvas();
   dummy->GetXaxis()->SetTitle("|#eta|");
   dummy->GetYaxis()->SetTitle("Data/MC ratio (const fit)");
   dummy->GetXaxis()->SetRangeUser(0., 5.);
   dummy->GetYaxis()->SetRangeUser(0.9, 2.0);
   dummy->Draw();

   Res_2012_total->SetPointError(0, 0.25, 0.25, eta0_TotalUncDown, eta0_TotalUncUp);
   Res_2012_total->SetPointError(1, 0.3, 0.3, eta1_TotalUncDown, eta1_TotalUncUp);
   Res_2012_total->SetPointError(2, 0.3, 0.3, eta2_TotalUncDown, eta2_TotalUncUp);
   Res_2012_total->SetPointError(3, 0.3, 0.3, eta3_TotalUncDown, eta3_TotalUncUp);
   Res_2012_total->SetPointError(4, 1.35, 1.35, eta4_TotalUncDown, eta4_TotalUncUp);
  
   Res_2012_total->SetMarkerStyle(20);
   Res_2012_total->SetMarkerSize(1.4);
   // Res_2012_total->SetFillColor(kYellow-3);
   Res_2012_total->SetFillColor(kCyan-2);
   Res_2012_total->SetFillStyle(1001);
   //Res_2012_total->SetLineColor(kYellow-3);
   Res_2012_total->SetLineColor(kCyan-2);
   Res_2012_total->DrawClone("2Psame");

   Res_2011->SetMarkerStyle(24);
   Res_2011->SetMarkerSize(1.4);
   Res_2011->SetFillColor(kGray+2);
   Res_2011->SetFillStyle(3744);
   Res_2011->SetLineColor(kGray+2);
   Res_2011->DrawClone("2psame");

   // cmsPrel();

   Res_2012_total->SetPointError(0, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(1, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(2, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(3, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(4, 0., 0., 0., 0.);

   Res_2012_total->GetXaxis()->SetRangeUser(0., 5.);
   Res_2012_total->Draw("psame");

   Res_2011->SetPointError(0, 0., 0., 0., 0.);
   Res_2011->SetPointError(1, 0., 0., 0., 0.);
   Res_2011->SetPointError(2, 0., 0., 0., 0.);
   Res_2011->SetPointError(3, 0., 0., 0., 0.);
   Res_2011->SetPointError(4, 0., 0., 0., 0.);

   Res_2011->Draw("psame");

   label->Draw("same");

   TLegend* leg3 = util::LabelFactory::createLegendColWithOffset(1, 1.05, 0.07);
   leg3->AddEntry(Res_2012_total,"JER 2012 Dijets (with tot. unc.)","PF");
     
   leg3->Draw("same");

   TPaveText *label2 = util::LabelFactory::createPaveTextWithOffset(1,1.05,0.15);
   label2->AddText("Anti-k_{T} (R=0.5) PF Jets");
   label2->Draw("same");

   TLegend* leg6 = util::LabelFactory::createLegendColWithOffset(1, 1.05, 0.21);
   leg6->AddEntry(Res_2011,"JER 2011 Dijets (with tot. unc.)","PF");
    
   leg6->Draw("same");

   c3->Print("Results/JER_2012_comp2011_final_v2.eps");
   c3->Print("Results/JER_2012_comp2011_final_v2.png");
   c3->Print("Results/JER_2012_comp2011_final_v2.pdf");


   // ------------------------------------------------------------------ //
   // plot nominal values with total unc. + photon comparison
   double photon_bins[5] = {0, 0.5, 1.1, 1.7, 2.3};
   TH1F *Res_2012Photon_hist = new TH1F("PhotonJER", "", 4, photon_bins);
   Res_2012Photon_hist->Reset();
   Res_2012Photon_hist->SetBinContent(1, 1.067);
   Res_2012Photon_hist->SetBinContent(2, 1.087);
   Res_2012Photon_hist->SetBinContent(3, 1.104);
   Res_2012Photon_hist->SetBinContent(4, 1.199);

   TGraphAsymmErrors *Res_2012Photon = new TGraphAsymmErrors(Res_2012Photon_hist);
   Res_2012Photon->SetPointError(0, 0.25, 0.25, 0.027, 0.028);
   Res_2012Photon->SetPointError(1, 0.3, 0.3, 0.041, 0.041);
   Res_2012Photon->SetPointError(2, 0.3, 0.3, 0.052, 0.052);
   Res_2012Photon->SetPointError(3, 0.3, 0.3, 0.084, 0.084);

   Res_2012_total->SetPointError(0, 0.25, 0.25, eta0_TotalUncDown, eta0_TotalUncUp);
   Res_2012_total->SetPointError(1, 0.3, 0.3, eta1_TotalUncDown, eta1_TotalUncUp);
   Res_2012_total->SetPointError(2, 0.3, 0.3, eta2_TotalUncDown, eta2_TotalUncUp);
   Res_2012_total->SetPointError(3, 0.3, 0.3, eta3_TotalUncDown, eta3_TotalUncUp);
   Res_2012_total->SetPointError(4, 1.35, 1.35, eta4_TotalUncDown, eta4_TotalUncUp);

   TCanvas *c4 = new TCanvas();
   dummy->GetXaxis()->SetTitle("|#eta|");
   dummy->GetYaxis()->SetTitle("Data/MC ratio (const fit)");
   dummy->GetXaxis()->SetRangeUser(0., 5.);
   dummy->GetYaxis()->SetRangeUser(0.9, 1.7);
   dummy->Draw();

   Res_2012_total->SetMarkerStyle(20);
   Res_2012_total->SetMarkerSize(1.4);
   // Res_2012_total->SetFillColor(kYellow-3);
   Res_2012_total->SetFillColor(kCyan-2);
   Res_2012_total->SetFillStyle(1001);
   //Res_2012_total->SetLineColor(kYellow-3);
   Res_2012_total->SetLineColor(kCyan-2);
   Res_2012_total->DrawClone("2Psame");

   Res_2012Photon->SetMarkerStyle(22);
   Res_2012Photon->SetMarkerSize(1.4);
   // Res_2012Photon->SetLineColor(kYellow-3);
   // Res_2012Photon->SetFillColor(kYellow-3);
   Res_2012Photon->SetLineColor(kPink-8);
   Res_2012Photon->SetFillColor(kPink-8);
   Res_2012Photon->SetFillStyle(3744);
   Res_2012Photon->DrawClone("2Psame");

   Res_2012_total->SetPointError(0, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(1, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(2, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(3, 0., 0., 0., 0.);
   Res_2012_total->SetPointError(4, 0., 0., 0., 0.);

   Res_2012_total->GetXaxis()->SetRangeUser(0., 5.);
   Res_2012_total->Draw("psame");
   gPad->RedrawAxis();

   label->Draw("same");

   TLegend* leg4 = util::LabelFactory::createLegendColWithOffset(2, 1.05, 0.07);
   leg4->AddEntry(Res_2012_total,"JER 2012 Dijets (with tot. unc.)","PF");
   leg4->AddEntry(Res_2012Photon,"JER 2012 Photon + Jet (with tot. unc.)","PF");
    
   leg4->Draw("same");

   // cmsPrel();

   c4->Print("Results/JER_2012_compPhoton_final_v2.eps");
   c4->Print("Results/JER_2012_compPhoton_final_v2.png");
   c4->Print("Results/JER_2012_compPhoton_final_v2.pdf");
}
