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
   TFile* file_SameEta = new TFile("Extrapolation/JER_RatioVsEta_herwigpp_final_v4.root", "READ");
   TFile* file_ForwardExtension = new TFile("ForwardExtension/JER_RatioVsEta_herwigpp_ForwardExtension_v4.root", "READ");
   TFile* file_ForwardExtensionNextToCentral = new TFile("ForwardExtension/JER_RatioVsEta_herwigpp_ForwardExtensionNextToCentral_v3.root", "READ");
   TFile* file_ForwardExtensionSecondNextToCentral = new TFile("ForwardExtension/JER_RatioVsEta_herwigpp_ForwardExtensionSecondNextToCentral_v1.root", "READ");

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

void CalcWeightedMeanHerwig()
{
   setTDRStyle();

   // ------------------------------------------------------------------ //
   double mean = 0.;
   double stat_unc = 0.;
  
   for(int i = 1; i < 8; i++) {

      mean = CalcWeightedResult("RatioVsEta_with_pli", i, true);
      stat_unc = CalcWeightedResult("RatioVsEta_with_pli", i, false);

      cout << "//---------------------------------//" << endl;
      cout << "Eta Bin : " << i << endl;
      cout << "//---------------------------------//" << endl;
      cout << "mean: " << mean << endl;
      cout << "stat unc = " << stat_unc << endl;
   }
}
