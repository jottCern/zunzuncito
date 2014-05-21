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

double CalcWeightedResult(TString histo_name, int bin, bool isMean)
{
   // ------------------------------------------------------------------ //
   // get files with nominal value and variations
   TFile* file_SameEta = new TFile("Results/Results_With_Uncertainties_final_SameEta_v4.root", "READ");
   TFile* file_ForwardExtension = new TFile("Results/Results_With_Uncertainties_ForwardExtension_v4.root", "READ");
   TFile* file_ForwardExtensionNextToCentral = new TFile("Results/Results_With_Uncertainties_ForwardExtensionNextToCentral_v3.root", "READ");
   TFile* file_ForwardExtensionSecondNextToCentral = new TFile("Results/Results_With_Uncertainties_ForwardExtensionSecondNextToCentral_v1.root", "READ");

   TH1F *histo_SameEta = new TH1F();
   TH1F *histo_ForwardExtension = new TH1F();
   TH1F *histo_ForwardExtensionNextToCentral = new TH1F();
   TH1F *histo_ForwardExtensionSecondNextToCentral = new TH1F();

   file_SameEta->cd();
   gDirectory->GetObject(histo_name, histo_SameEta);

   file_ForwardExtension->cd();
   gDirectory->GetObject(histo_name, histo_ForwardExtension);

   file_ForwardExtensionNextToCentral->cd();
   gDirectory->GetObject(histo_name, histo_ForwardExtensionNextToCentral);

   file_ForwardExtensionSecondNextToCentral->cd();
   gDirectory->GetObject(histo_name, histo_ForwardExtensionSecondNextToCentral);

   cout << "//---------------------------------//" << endl;
   cout << "Histo Name" << histo_name << endl;
   cout << "Same Eta: " << histo_SameEta->GetBinContent(bin) << " +- " << histo_SameEta->GetBinError(bin) << endl;
   cout << "Forward extension: " << histo_ForwardExtension->GetBinContent(bin) << " +- " << histo_ForwardExtension->GetBinError(bin) << endl;
   cout << "Forward extension next to central: " << histo_ForwardExtensionNextToCentral->GetBinContent(bin) << " +- " << histo_ForwardExtensionNextToCentral->GetBinError(bin) << endl;
   cout << "Forward extension second next to central: " << histo_ForwardExtensionSecondNextToCentral->GetBinContent(bin) << " +- " << histo_ForwardExtensionSecondNextToCentral->GetBinError(bin) << endl;
   cout << "//---------------------------------//" << endl;

   double content_SameEta = 0.;
   double unc_SameEta = 0.;
   double wi_SameEta = 0.;
   if(histo_SameEta->GetBinContent(bin) > 0 && bin < 6){
      content_SameEta = histo_SameEta->GetBinContent(bin);
      unc_SameEta = histo_SameEta->GetBinError(bin);
      wi_SameEta = 1/(unc_SameEta*unc_SameEta);
   }

   double content_ForwardExtension = 0.;
   double unc_ForwardExtension = 0.;
   double wi_ForwardExtension = 0.;
   if(histo_ForwardExtension->GetBinContent(bin) > 0 && bin != 1) {
      content_ForwardExtension = histo_ForwardExtension->GetBinContent(bin);
      unc_ForwardExtension = histo_ForwardExtension->GetBinError(bin);
      wi_ForwardExtension = 1/(unc_ForwardExtension*unc_ForwardExtension);
   }

   double content_ForwardExtensionNextToCentral = 0.;
   double unc_ForwardExtensionNextToCentral = 0.;
   double wi_ForwardExtensionNextToCentral = 0.;
   if(histo_ForwardExtensionNextToCentral->GetBinContent(bin) > 0 && bin != 2) {
      content_ForwardExtensionNextToCentral = histo_ForwardExtensionNextToCentral->GetBinContent(bin);
      unc_ForwardExtensionNextToCentral = histo_ForwardExtensionNextToCentral->GetBinError(bin);
      wi_ForwardExtensionNextToCentral = 1/(unc_ForwardExtensionNextToCentral*unc_ForwardExtensionNextToCentral);
   }

   double content_ForwardExtensionSecondNextToCentral = 0.;
   double unc_ForwardExtensionSecondNextToCentral = 0.;
   double wi_ForwardExtensionSecondNextToCentral = 0.;
   if(histo_ForwardExtensionSecondNextToCentral->GetBinContent(bin) > 0 && bin != 3) {
      content_ForwardExtensionSecondNextToCentral = histo_ForwardExtensionSecondNextToCentral->GetBinContent(bin);
      unc_ForwardExtensionSecondNextToCentral = histo_ForwardExtensionSecondNextToCentral->GetBinError(bin);
      wi_ForwardExtensionSecondNextToCentral = 1/(unc_ForwardExtensionSecondNextToCentral*unc_ForwardExtensionSecondNextToCentral);
   }
   
   // double variance = 1/(wi_SameEta + wi_ForwardExtension + wi_ForwardExtensionNextToCentral);
   double variance = 1/(wi_SameEta + wi_ForwardExtension + wi_ForwardExtensionNextToCentral + wi_ForwardExtensionSecondNextToCentral);

   // double mean = (content_SameEta*wi_SameEta + content_ForwardExtension*wi_ForwardExtension + content_ForwardExtensionNextToCentral*wi_ForwardExtensionNextToCentral) * variance;
   double mean = (content_SameEta*wi_SameEta + content_ForwardExtension*wi_ForwardExtension + content_ForwardExtensionNextToCentral*wi_ForwardExtensionNextToCentral + content_ForwardExtensionSecondNextToCentral*wi_ForwardExtensionSecondNextToCentral) * variance;

   if(isMean) return mean;
   else return TMath::Sqrt(variance);

   return -1;

}

void CalcSysUncJER_weighted()
{
   setTDRStyle();

   // ------------------------------------------------------------------ //
   double mean = 0.;
   double stat_unc = 0.;
   double JEC_unc_up = 0.;
   double JEC_unc_down = 0.;
   double PLI_unc_up = 0.;
   double PLI_unc_down = 0.;
   double AlphaReweight_unc_up = 0.;
   double AlphaReweight_unc_down = 0.;
   double PU_unc_up = 0.;
   double PU_unc_down = 0.;
   double Tail_unc_up = 0.;
   double Tail_unc_down = 0.;
   double Gluon_unc_up = 0.;
   double Gluon_unc_down = 0.;
   double AlphaRange_unc_up = 0.;
   double AlphaRange_unc_down = 0.;
   double PtTrend_unc_up = 0.;
   double PtTrend_unc_down = 0.;
   double tot_sys_up = 0.;
   double tot_sys_down = 0.;

   for(int i = 1; i < 8; i++) {

      mean = CalcWeightedResult("RatioNominal", i, true);
      stat_unc = CalcWeightedResult("RatioNominal", i, false);

      JEC_unc_up = CalcWeightedResult("RatioJECUp", i, true) - mean;
      JEC_unc_down = mean - CalcWeightedResult("RatioJECDown", i, true);
      PLI_unc_up = CalcWeightedResult("RatioPLIUp", i, true) - mean;
      PLI_unc_down = mean - CalcWeightedResult("RatioPLIDown", i, true);
      AlphaReweight_unc_up = CalcWeightedResult("RatioAlphaReweightUp", i, true) - mean;
      AlphaReweight_unc_down = mean - CalcWeightedResult("RatioAlphaReweightDown", i, true);
      PU_unc_up = CalcWeightedResult("RatioPUUp", i, true) - mean;
      PU_unc_down = mean - CalcWeightedResult("RatioPUDown", i, true);
      Tail_unc_up = CalcWeightedResult("RatioTailsUp", i, true) - mean;
      Tail_unc_down = mean - CalcWeightedResult("RatioTailsDown", i, true);
      Gluon_unc_up = CalcWeightedResult("RatioGluonUp", i, true) - mean;
      Gluon_unc_down = mean - CalcWeightedResult("RatioGluonDown", i, true);
      AlphaRange_unc_up = CalcWeightedResult("RatioAddAlphaUp", i, true) - mean;
      AlphaRange_unc_down = mean - CalcWeightedResult("RatioAddAlphaDown", i, true);
      PtTrend_unc_up = CalcWeightedResult("RatioPtTrendUp", i, true) - mean;
      PtTrend_unc_down = mean - CalcWeightedResult("RatioPtTrendDown", i, true);

      tot_sys_up = TMath::Sqrt(JEC_unc_up*JEC_unc_up + PLI_unc_up*PLI_unc_up + AlphaReweight_unc_up*AlphaReweight_unc_up + PU_unc_up*PU_unc_up + Tail_unc_up*Tail_unc_up + Gluon_unc_up*Gluon_unc_up + AlphaRange_unc_up*AlphaRange_unc_up + PtTrend_unc_up*PtTrend_unc_up);

      tot_sys_down = TMath::Sqrt(JEC_unc_down*JEC_unc_down + PLI_unc_down*PLI_unc_down + AlphaReweight_unc_down*AlphaReweight_unc_down + PU_unc_down*PU_unc_down + Tail_unc_down*Tail_unc_down + Gluon_unc_down*Gluon_unc_down + AlphaRange_unc_down*AlphaRange_unc_down + PtTrend_unc_down*PtTrend_unc_down);

      /* cout << "//---------------------------------//" << endl;
      cout << "//----- RAW --------- //" << endl;
      cout << "Eta Bin : " << i << endl;
      cout << "//---------------------------------//" << endl;
      cout << "mean: " << mean << endl;
      cout << "stat unc = " << stat_unc << endl;
      cout << "//---------------------------------//" << endl;
      cout << "JEC: " << JEC_unc_up << " - " << JEC_unc_down << endl;
      cout << "PLI: " << PLI_unc_up << " - " << PLI_unc_down << endl;
      cout << "Alpha Reweight: " << AlphaReweight_unc_up << " - " << AlphaReweight_unc_down << endl;
      cout << "PU_unc: " << PU_unc_up << " - " << PU_unc_down << endl;
      cout << "Tail_unc: " << Tail_unc_up << " - " << Tail_unc_down << endl;
      cout << "Gluon_unc: " << Gluon_unc_up << " - " << Gluon_unc_down << endl;
      cout << "Alpha Range: " << AlphaRange_unc_up << " - " << AlphaRange_unc_down << endl;
      cout << "Pt Trend: " << PtTrend_unc_up << " - " << PtTrend_unc_down << endl;
      cout << "//---------------------------------//" << endl;
      cout << "Syst. unc = " << tot_sys_up << " - " << tot_sys_down << endl;
      cout << "//---------------------------------//" << endl;
      cout << "Total unc: " << TMath::Sqrt(tot_sys_up*tot_sys_up + stat_unc*stat_unc) << endl;*/

      // JEC
      if(JEC_unc_up < 0 || JEC_unc_down < 0) {
         if(TMath::Abs(JEC_unc_up) > TMath::Abs(JEC_unc_down) ) {
            JEC_unc_down = TMath::Abs(JEC_unc_up);
            JEC_unc_up = TMath::Abs(JEC_unc_up);
         }
         else {
            JEC_unc_up = TMath::Abs(JEC_unc_down);
            JEC_unc_down = TMath::Abs(JEC_unc_down);
         }
      }
      if(JEC_unc_up > 0 && JEC_unc_down > 0) {
         JEC_unc_up = 0.5*(JEC_unc_up + JEC_unc_down);
         JEC_unc_down = 0.5*(JEC_unc_up + JEC_unc_down);
      }

      // PLI
      if(PLI_unc_up < 0 || PLI_unc_down < 0) {
         if(TMath::Abs(PLI_unc_up) > TMath::Abs(PLI_unc_down) ) {
            PLI_unc_down = TMath::Abs(PLI_unc_up);
            PLI_unc_up = TMath::Abs(PLI_unc_up);
         }
         else {
            PLI_unc_up = TMath::Abs(PLI_unc_down);
            PLI_unc_down = TMath::Abs(PLI_unc_down);
         }
      }
      if(PLI_unc_up > 0 && PLI_unc_down > 0) {
         PLI_unc_up = 0.5*(PLI_unc_up + PLI_unc_down);
         PLI_unc_down = 0.5*(PLI_unc_up + PLI_unc_down);
      }

      // Alpha Reweight
      if(TMath::Abs(AlphaReweight_unc_up) > 0) {
         AlphaReweight_unc_up = TMath::Abs(AlphaReweight_unc_up);
         AlphaReweight_unc_down = TMath::Abs(AlphaReweight_unc_up);
      }
      else if(TMath::Abs(AlphaReweight_unc_down) > 0) {
         AlphaReweight_unc_up = TMath::Abs(AlphaReweight_unc_down);
         AlphaReweight_unc_down = TMath::Abs(AlphaReweight_unc_down);
      }

      // PU
      if(TMath::Abs(PU_unc_up) > 0) {
         PU_unc_up = TMath::Abs(PU_unc_up);
         PU_unc_down = TMath::Abs(PU_unc_up);
      }
      else if(TMath::Abs(PU_unc_down) > 0) {
         PU_unc_up = TMath::Abs(PU_unc_down);
         PU_unc_down = TMath::Abs(PU_unc_down);
      }

      // Tail
      if(TMath::Abs(Tail_unc_up) > 0) {
         Tail_unc_up = TMath::Abs(Tail_unc_up);
         Tail_unc_down = TMath::Abs(Tail_unc_up);
      }
      else if(TMath::Abs(Tail_unc_down) > 0) {
         Tail_unc_up = TMath::Abs(Tail_unc_down);
         Tail_unc_down = TMath::Abs(Tail_unc_down);
      }

      // Gluon
      if(TMath::Abs(Gluon_unc_up) > 0) {
         Gluon_unc_up = TMath::Abs(Gluon_unc_up);
         Gluon_unc_down = TMath::Abs(Gluon_unc_up);
      }
      else if(TMath::Abs(Gluon_unc_down) > 0) {
         Gluon_unc_up = TMath::Abs(Gluon_unc_down);
         Gluon_unc_down = TMath::Abs(Gluon_unc_down);
      }

      // Alpha Range
      if(TMath::Abs(AlphaRange_unc_up) > 0) {
         AlphaRange_unc_up = TMath::Abs(AlphaRange_unc_up);
         AlphaRange_unc_down = TMath::Abs(AlphaRange_unc_up);
      }
      else if(TMath::Abs(AlphaRange_unc_down) > 0) {
         AlphaRange_unc_up = TMath::Abs(AlphaRange_unc_down);
         AlphaRange_unc_down = TMath::Abs(AlphaRange_unc_down);
      }

      // Pt Trend
      PtTrend_unc_up = 0.02 * mean;
      PtTrend_unc_down = 0.02 * mean;

      tot_sys_up = TMath::Sqrt(JEC_unc_up*JEC_unc_up + PLI_unc_up*PLI_unc_up + AlphaReweight_unc_up*AlphaReweight_unc_up + PU_unc_up*PU_unc_up + Tail_unc_up*Tail_unc_up + Gluon_unc_up*Gluon_unc_up + AlphaRange_unc_up*AlphaRange_unc_up + PtTrend_unc_up*PtTrend_unc_up);

      tot_sys_down = TMath::Sqrt(JEC_unc_down*JEC_unc_down + PLI_unc_down*PLI_unc_down + AlphaReweight_unc_down*AlphaReweight_unc_down + PU_unc_down*PU_unc_down + Tail_unc_down*Tail_unc_down + Gluon_unc_down*Gluon_unc_down + AlphaRange_unc_down*AlphaRange_unc_down + PtTrend_unc_down*PtTrend_unc_down);

      cout << "//---------------------------------//" << endl;
      cout << "//----- CORRECTED --------- //" << endl;
      cout << "Eta Bin : " << i << endl;
      cout << "//---------------------------------//" << endl;
      cout << "mean: " << mean << endl;
      cout << "stat unc = " << stat_unc << endl;
      cout << "//---------------------------------//" << endl;
      cout << "JEC: " << JEC_unc_up << " - " << JEC_unc_down << endl;
      cout << "PLI: " << PLI_unc_up << " - " << PLI_unc_down << endl;
      cout << "Alpha Reweight: " << AlphaReweight_unc_up << " - " << AlphaReweight_unc_down << endl;
      cout << "PU_unc: " << PU_unc_up << " - " << PU_unc_down << endl;
      cout << "Tail_unc: " << Tail_unc_up << " - " << Tail_unc_down << endl;
      cout << "Gluon_unc: " << Gluon_unc_up << " - " << Gluon_unc_down << endl;
      cout << "Alpha Range: " << AlphaRange_unc_up << " - " << AlphaRange_unc_down << endl;
      cout << "Pt Trend: " << PtTrend_unc_up << " - " << PtTrend_unc_down << endl;
      cout << "//---------------------------------//" << endl;
      cout << "Syst. unc = " << tot_sys_up << " - " << tot_sys_down << endl;
      cout << "//---------------------------------//" << endl;
      cout << "Total unc: " << TMath::Sqrt(tot_sys_up*tot_sys_up + stat_unc*stat_unc) << endl;
   }
}
