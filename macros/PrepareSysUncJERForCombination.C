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

void PrepareSysUncJERForCombination()
{
   setTDRStyle();

   TString suffix = "_final_SameEta_v4";
   // TString suffix = "_ForwardExtension_v4";
   // TString suffix = "_ForwardExtensionNextToCentral_v3";
   // TString suffix = "_ForwardExtensionSecondNextToCentral_v1";

   // ------------------------------------------------------------------ //
   // get files with nominal value and variations
   TFile* file_nominal = new TFile();
   TFile* file_jecDOWN = new TFile();
   TFile* file_jecUP = new TFile();
   TFile* file_pliDOWN = new TFile();
   TFile* file_pliUP = new TFile();
   TFile* file_alpha = new TFile();
   TFile* file_pu = new TFile();
   TFile* file_tails = new TFile();
   TFile* file_gluon = new TFile();
   TFile* file_addalpha = new TFile();

   if(suffix.Contains("final_SameEta_v4")) {
      file_nominal = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_v4.root", "READ");
      file_jecDOWN = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_JECdown_v4.root", "READ");
      file_jecUP = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_JECup_v4.root", "READ");
      file_pliDOWN = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_PLIdown_v4.root", "READ");
      file_pliUP = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_PLIup_v4.root", "READ");
      file_alpha = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_ReweightAlphaSpectrum_v4.root", "READ");
      file_pu = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_MBXS73500_v4.root", "READ");
      file_tails = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_95Truncation_v4.root", "READ");
      file_gluon = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_ReweightGluonSplitting_v4.root", "READ");
      file_addalpha = new TFile("Extrapolation/JER_RatioVsEta_final_nominal_NoMinPtCutForThirdJet_AddNewAlphaBin_v4.root", "READ");
   }

   if(suffix.Contains("ForwardExtension_")) {
      file_nominal = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_v4.root", "READ");
      file_jecDOWN = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_JECdown_v4.root", "READ");
      file_jecUP = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_JECup_v4.root", "READ");
      file_pliDOWN = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_PLIdown_v4.root", "READ");
      file_pliUP = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_PLIup_v4.root", "READ");
      file_alpha = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_ReweightAlphaSpectrum_v4.root", "READ");
      file_pu = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_MBXS73500_v4.root", "READ");
      file_tails = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_95Truncation_v4.root", "READ");
      file_gluon = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_ReweightGluonSplitting_v4.root", "READ");
      file_addalpha = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_v4.root", "READ");
   }

   if(suffix.Contains("ForwardExtensionNextToCentral")) {
      file_nominal = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_v3.root", "READ");
      file_jecDOWN = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_JECdown_v3.root", "READ");
      file_jecUP = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_JECup_v3.root", "READ");
      file_pliDOWN = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_PLIdown_v3.root", "READ");
      file_pliUP = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_PLIup_v3.root", "READ");
      file_alpha = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_ReweightAlphaSpectrum_v3.root", "READ");
      file_pu = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_MBXS73500_v3.root", "READ");
      file_tails = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_95Truncation_v3.root", "READ");
      file_gluon = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_ReweightGluonSplitting_v3.root", "READ");
      file_addalpha = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_v3.root", "READ");
   }

   if(suffix.Contains("ForwardExtensionSecondNextToCentral")) {
      file_nominal = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_v1.root", "READ");
      file_jecDOWN = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_JECdown_v1.root", "READ");
      file_jecUP = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_JECup_v1.root", "READ");
      file_pliDOWN = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_PLIdown_v1.root", "READ");
      file_pliUP = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_PLIup_v1.root", "READ");
      file_alpha = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_ReweightAlphaSpectrum_v1.root", "READ");
      file_pu = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_MBXS73500_v1.root", "READ");
      file_tails = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_95Truncation_v1.root", "READ");
      file_gluon = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_ReweightGluonSplitting_v1.root", "READ");
      file_addalpha = new TFile("ForwardExtension/JER_RatioVsEta_ForwardExtensionSecondNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_v1.root", "READ");
   }
  

   // ------------------------------------------------------------------ //
   // define histos needed for uncertainties
   TH1F *RatioNominal = new TH1F();
   TH1F *RatioJECDown = new TH1F();
   TH1F *RatioJECUp = new TH1F();
   TH1F *RatioPLIDown = new TH1F();
   TH1F *RatioPLIUp = new TH1F();
   TH1F *RatioAlphaReweightUp = new TH1F();
   TH1F *RatioPUUp = new TH1F();
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

   file_tails->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioTailsDown);

   file_gluon->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioGluonUp);

   file_addalpha->cd();
   gDirectory->GetObject("RatioVsEta_with_pli;1", RatioAddAlphaDown);

   // ------------------------------------------------------------------ //
   // calc PU errors
   TH1F *RatioPUDown = new TH1F(*RatioPUUp);
   RatioPUDown->Reset();

   for (int i = 1; i < RatioPUUp->GetNbinsX()+1; i++) {
      double bin_content_pu = RatioPUUp->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_pu - bin_content_nominal;
        
      //  cout << "PU i diff: " << i << "  "  << diff << endl;

      RatioPUDown->SetBinContent(i, 0);
      RatioPUDown->SetBinError(i, 0);

      /*  if( diff > 0 ) {
         RatioPUUp->SetBinContent(i, bin_content_pu);
         RatioPUUp->SetBinError(i, RatioPUUp->GetBinError(i));
         RatioPUDown->SetBinContent(i, bin_content_nominal - diff);
         RatioPUDown->SetBinError(i, RatioPUUp->GetBinError(i));
      }
      else if ( diff < 0 ) {
         RatioPUDown->SetBinContent(i, bin_content_pu);
         RatioPUDown->SetBinError(i, RatioPUUp->GetBinError(i));
         RatioPUUp->SetBinContent(i, bin_content_nominal - diff);
         RatioPUUp->SetBinError(i, RatioPUUp->GetBinError(i));
         }*/
   }  

   // calc AlphaReweight errors
   TH1F *RatioAlphaReweightDown = new TH1F(*RatioAlphaReweightUp);
   RatioAlphaReweightDown->Reset();

   for (int i = 1; i < RatioAlphaReweightUp->GetNbinsX()+1; i++) {
      double bin_content_alpha = RatioAlphaReweightUp->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_alpha - bin_content_nominal;
    
      //  cout << "Alpha reweight i diff: " << i << "  "  << diff << endl;

      RatioAlphaReweightDown->SetBinContent(i, 0);
      RatioAlphaReweightDown->SetBinError(i, 0);

      /* if( diff > 0 ) {
         RatioAlphaReweightUp->SetBinContent(i, bin_content_alpha);
         RatioAlphaReweightUp->SetBinError(i, RatioAlphaReweightUp->GetBinError(i));
         RatioAlphaReweightDown->SetBinContent(i, bin_content_nominal - diff);
         RatioAlphaReweightDown->SetBinError(i, RatioAlphaReweightUp->GetBinError(i));
      }
      else if ( diff < 0 ) {
         RatioAlphaReweightDown->SetBinContent(i, bin_content_alpha);
         RatioAlphaReweightDown->SetBinError(i, RatioAlphaReweightUp->GetBinError(i));
         RatioAlphaReweightUp->SetBinContent(i, bin_content_nominal - diff);
         RatioAlphaReweightUp->SetBinError(i, RatioAlphaReweightUp->GetBinError(i));
         }*/
   }

   // calc symmetric JEC errors
   /*  for (int i = 1; i < RatioJECDown->GetNbinsX()+1; i++) {
      double bin_content_jec_low = RatioJECDown->GetBinContent(i);
      double bin_error_jec_low = RatioJECDown->GetBinError(i);
      double bin_content_jec_high = RatioJECUp->GetBinContent(i);
      double bin_error_jec_high = RatioJECUp->GetBinError(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff_low = bin_content_jec_low - bin_content_nominal;
      double diff_high = bin_content_jec_high - bin_content_nominal;

     //  cout << "JEC Up in: " << bin_content_jec_high << endl;
//       cout << "JEC Down in: " << bin_content_jec_low << endl;
//       cout << "diff up: " << diff_high << endl;
//       cout << "diff down: " << diff_low << endl;
     
      if( diff_low > 0 && diff_high > 0 || diff_low < 0 && diff_high < 0) {
         if(diff_low < 0 && diff_high < 0) {
            if( diff_low > diff_high) {
               bin_content_jec_low = bin_content_nominal - TMath::Abs(diff_high);
               bin_content_jec_high = bin_content_nominal + TMath::Abs(diff_high);
               bin_error_jec_low = bin_error_jec_high;
            }
            else {
               bin_content_jec_low = bin_content_nominal - TMath::Abs(diff_low);
               bin_content_jec_high = bin_content_nominal + TMath::Abs(diff_low);
               bin_error_jec_high = bin_error_jec_low;
            }
         }
         else {
            if( diff_low < diff_high) {
               bin_content_jec_low = bin_content_nominal - TMath::Abs(diff_high);
               bin_content_jec_high = bin_content_nominal + TMath::Abs(diff_high);
               bin_error_jec_low = bin_error_jec_high;
            }
            else {
               bin_content_jec_low = bin_content_nominal - TMath::Abs(diff_low);
               bin_content_jec_high = bin_content_nominal + TMath::Abs(diff_low);
               bin_error_jec_high = bin_error_jec_low;
            }
         }
      }

     //  cout << "JEC Up outn: " << bin_content_jec_high << endl;
//       cout << "JEC Down out: " << bin_content_jec_low << endl;

      if( diff_low > 0 && diff_high < 0) {
         bin_content_jec_high = bin_content_nominal + diff_low;
         bin_error_jec_high = bin_error_jec_low;
         bin_content_jec_low = bin_content_nominal + diff_high;
         bin_error_jec_low = bin_error_jec_high;
      }
        
      RatioJECDown->SetBinContent(i, bin_content_jec_low);
      RatioJECDown->SetBinError(i, bin_error_jec_low);
      RatioJECUp->SetBinContent(i, bin_content_jec_high);
      RatioJECUp->SetBinError(i, bin_error_jec_high); 
      }*/

   // calc symmetric PLI errors
   /*for (int i = 1; i < RatioPLIDown->GetNbinsX()+1; i++) {
      double bin_content_pli_low = RatioPLIDown->GetBinContent(i);
      double bin_error_pli_low = RatioPLIDown->GetBinError(i);
      double bin_content_pli_high = RatioPLIUp->GetBinContent(i);
      double bin_error_pli_high = RatioPLIUp->GetBinError(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff_low = bin_content_pli_low - bin_content_nominal;
      double diff_high = bin_content_pli_high - bin_content_nominal;
     
      if( diff_low > 0 && diff_high > 0 || diff_low < 0 && diff_high < 0) {
         if(diff_low < 0 && diff_high < 0) {
            if( diff_low > diff_high) {
               bin_content_pli_low = bin_content_nominal - TMath::Abs(diff_high);
               bin_content_pli_high = bin_content_nominal + TMath::Abs(diff_high);
               bin_error_pli_low = bin_error_pli_high;
            }
            else {
               bin_content_pli_low = bin_content_nominal - TMath::Abs(diff_low);
               bin_content_pli_high = bin_content_nominal + TMath::Abs(diff_low);
               bin_error_pli_high = bin_error_pli_low;
            }
         }
         else {
            if( diff_low < diff_high) {
               bin_content_pli_low = bin_content_nominal - TMath::Abs(diff_high);
               bin_content_pli_high = bin_content_nominal + TMath::Abs(diff_high);
               bin_error_pli_low = bin_error_pli_high;
            }
            else {
               bin_content_pli_low = bin_content_nominal - TMath::Abs(diff_low);
               bin_content_pli_high = bin_content_nominal + TMath::Abs(diff_low);
               bin_error_pli_high = bin_error_pli_low;
            }
         }
      }

      if( diff_low > 0 && diff_high < 0) {
         bin_content_pli_high = bin_content_nominal + diff_low;
         bin_error_pli_high = bin_error_pli_low;
         bin_content_pli_low = bin_content_nominal + diff_high;
         bin_error_pli_low = bin_error_pli_high;
      }
     
      RatioPLIDown->SetBinContent(i, bin_content_pli_low);
      RatioPLIDown->SetBinError(i, bin_error_pli_low);
      RatioPLIUp->SetBinContent(i, bin_content_pli_high);
      RatioPLIUp->SetBinError(i, bin_error_pli_high); 
      }*/

   // calc symmetric Tail errors
   TH1F *RatioTailsUp = new TH1F(*RatioTailsDown);
   RatioTailsUp->Reset();

   for (int i = 1; i < RatioTailsUp->GetNbinsX()+1; i++) {
      double bin_content_tails = RatioTailsDown->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_tails - bin_content_nominal;

      RatioTailsUp->SetBinContent(i, 0);
      RatioTailsUp->SetBinError(i, 0);
     
      /* if( diff > 0 ) {
         RatioTailsUp->SetBinContent(i, bin_content_tails);
         RatioTailsUp->SetBinError(i, RatioTailsDown->GetBinError(i));
         RatioTailsDown->SetBinContent(i, bin_content_nominal - diff);
         RatioTailsDown->SetBinError(i, RatioTailsDown->GetBinError(i));
      }
      else if ( diff < 0 ) {
         RatioTailsDown->SetBinContent(i, bin_content_tails);
         RatioTailsDown->SetBinError(i, RatioTailsDown->GetBinError(i));
         RatioTailsUp->SetBinContent(i, bin_content_nominal - diff);
         RatioTailsUp->SetBinError(i, RatioTailsDown->GetBinError(i));
         }*/
   }

   // calc symmetric GluonReweighting errors
   TH1F *RatioGluonDown = new TH1F(*RatioGluonUp);
   RatioGluonDown->Reset();

   for (int i = 1; i < RatioGluonUp->GetNbinsX()+1; i++) {
      double bin_content_alpha = RatioGluonUp->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_alpha - bin_content_nominal;

      RatioGluonDown->SetBinContent(i, 0);
      RatioGluonDown->SetBinError(i, 0);

      /* if( diff > 0 ) {
         RatioGluonUp->SetBinContent(i, bin_content_alpha);
         RatioGluonUp->SetBinError(i, RatioGluonUp->GetBinError(i));
         RatioGluonDown->SetBinContent(i, bin_content_nominal - diff);
         RatioGluonDown->SetBinError(i, RatioGluonUp->GetBinError(i));
      }
      else if ( diff < 0 ) {
         RatioGluonDown->SetBinContent(i, bin_content_alpha);
         RatioGluonDown->SetBinError(i, RatioGluonUp->GetBinError(i));
         RatioGluonUp->SetBinContent(i, bin_content_nominal - diff);
         RatioGluonUp->SetBinError(i, RatioGluonUp->GetBinError(i));
         }*/
   }

   // calc symmetric errors for additional alpha bin --> alpha range
   TH1F *RatioAddAlphaUp = new TH1F(*RatioAddAlphaDown);
   RatioAddAlphaUp->Reset();

   for (int i = 1; i < RatioAddAlphaUp->GetNbinsX()+1; i++) {
      double bin_content_alpha = RatioAddAlphaDown->GetBinContent(i);
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      double diff = bin_content_alpha - bin_content_nominal;

      RatioAddAlphaUp->SetBinContent(i, 0);
      RatioAddAlphaUp->SetBinError(i, 0);

      /* if( diff > 0 ) {
         RatioAddAlphaUp->SetBinContent(i, bin_content_alpha);
         RatioAddAlphaUp->SetBinError(i, RatioAddAlphaDown->GetBinError(i));
         RatioAddAlphaDown->SetBinContent(i, bin_content_nominal - diff);
         RatioAddAlphaDown->SetBinError(i, RatioAddAlphaDown->GetBinError(i));
      }
      else if ( diff < 0 ) {
         RatioAddAlphaDown->SetBinContent(i, bin_content_alpha);
         RatioAddAlphaDown->SetBinError(i, RatioAddAlphaDown->GetBinError(i));
         RatioAddAlphaUp->SetBinContent(i, bin_content_nominal - diff);
         RatioAddAlphaUp->SetBinError(i, RatioAddAlphaDown->GetBinError(i));
         }*/
   }

   // calc symmetric errors for possible pt-trend
   TH1F *RatioPtTrendDown = new TH1F(*RatioAddAlphaDown);
   TH1F *RatioPtTrendUp = new TH1F(*RatioAddAlphaDown);
   RatioPtTrendDown->Reset();
   RatioPtTrendUp->Reset();

   for (int i = 1; i < RatioPtTrendDown->GetNbinsX()+1; i++) {
      double bin_content_nominal = RatioNominal->GetBinContent(i);

      RatioPtTrendUp->SetBinContent(i, 0);
      RatioPtTrendUp->SetBinError(i, 0);
      RatioPtTrendDown->SetBinContent(i, 0);
      RatioPtTrendDown->SetBinError(i, 0);

      /* RatioPtTrendUp->SetBinContent(i, bin_content_nominal + 0.02*bin_content_nominal);
      RatioPtTrendUp->SetBinError(i, RatioNominal->GetBinError(i));
      RatioPtTrendDown->SetBinContent(i, bin_content_nominal - 0.02*bin_content_nominal);
      RatioPtTrendDown->SetBinError(i, RatioNominal->GetBinError(i));*/
   }

   // write results to output root-file
   TFile* output_file = new TFile("Results/Results_With_Uncertainties" + suffix + ".root", "RECREATE");

   RatioNominal->SetName("RatioNominal");
   RatioNominal->Write(); 
   RatioJECDown->SetName("RatioJECDown");
   RatioJECDown->Write(); 
   RatioJECUp->SetName("RatioJECUp");
   RatioJECUp->Write(); 
   RatioPLIDown->SetName("RatioPLIDown");
   RatioPLIDown->Write(); 
   RatioPLIUp->SetName("RatioPLIUp");
   RatioPLIUp->Write(); 
   RatioAlphaReweightDown->SetName("RatioAlphaReweightDown");
   RatioAlphaReweightDown->Write(); 
   RatioAlphaReweightUp->SetName("RatioAlphaReweightUp");
   RatioAlphaReweightUp->Write(); 
   RatioPUDown->SetName("RatioPUDown");
   RatioPUDown->Write(); 
   RatioPUUp->SetName("RatioPUUp");
   RatioPUUp->Write(); 
   RatioTailsDown->SetName("RatioTailsDown");
   RatioTailsDown->Write();
   RatioTailsUp->SetName("RatioTailsUp");
   RatioTailsUp->Write();
   RatioGluonDown->SetName("RatioGluonDown");
   RatioGluonDown->Write();
   RatioGluonUp->SetName("RatioGluonUp");
   RatioGluonUp->Write();
   RatioAddAlphaDown->SetName("RatioAddAlphaDown");
   RatioAddAlphaDown->Write(); 
   RatioAddAlphaUp->SetName("RatioAddAlphaUp");
   RatioAddAlphaUp->Write(); 
   RatioPtTrendDown->SetName("RatioPtTrendDown");
   RatioPtTrendDown->Write();
   RatioPtTrendUp->SetName("RatioPtTrendUp");
   RatioPtTrendUp->Write();

   output_file->Write();
}
