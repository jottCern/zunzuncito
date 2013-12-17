#include "/nfs/dust/test/cms/user/kheine/Kalibri/scripts/tdrstyle_mod.C"
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
    
      for(int i = 0; i < y_val.size(); i++) {
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
      for( int i=0; i < x_val.size(); i++) {
         for(int k= 0; k < x_val.size(); k++) {
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
void make_lin_fit(double & slope, double & d_slope, double & offset, double & d_offset){
    TMinuit min;
    min.SetPrintLevel(-1);
    //min.SetPrintLevel(0);
    int err = min.DefineParameter(0, "slope", slope, d_slope, 0.05, 1.0);
    assert(err==0);
    err = min.DefineParameter(1, "offset", offset, d_offset, 0.001, 0.2);
    assert(err==0);
    min.SetFCN(chi2_linear);
    min.mnmigr();
    min.GetParameter(0, slope, d_slope);
    min.GetParameter(1, offset, d_offset);
}

// --------------------------------- //
float GetAsymmWidth(TH1F* htemp, double * xq_IQW, double * yq_IQW) 
{
   const int nq = 2;

   float width = 0.;
   
   if( htemp->GetEntries() > 150 ) {
      htemp->GetXaxis()->SetRange(0,-1);
      htemp->ComputeIntegral();
      htemp->GetQuantiles(nq,yq_IQW,xq_IQW);
      Int_t IQW_high_bin_i = htemp->FindBin(yq_IQW[1]);
      cout << "Bin where to cut: " << IQW_high_bin_i << endl;
      cout << "Bins before truncation: " << htemp->GetNbinsX() << endl;
      cout << "Integral before truncation: " << htemp->Integral() << endl;
      //  htemp->GetXaxis()->SetRange(0,IQW_high_bin_i);
      // htemp->SetBins(IQW_high_bin_i, 0, IQW_high_bin_i);
      // cout << "Bins after truncation: " << htemp->GetNbinsX() << endl;
      cout << "Integral after truncation: " << htemp->Integral(1, IQW_high_bin_i) << endl;
      cout << "Truncated Integral in %: " << htemp->Integral(1, IQW_high_bin_i)/htemp->Integral() << endl;
      //  for(int i=1; i <= htemp->GetNbinsX(); i++) {
      for(int i=1; i <= IQW_high_bin_i; i++) {
         width += htemp->GetBinContent(i)* std::pow(htemp->GetBinCenter(i), 2);
      }

      // width = TMath::Sqrt(1/(htemp->Integral())*width);
      width = TMath::Sqrt(1/(htemp->Integral(1, IQW_high_bin_i))*width);
   }

   // cout << "width: " << width << endl;
   //  htemp->GetXaxis()->SetRange(0,-1);

   return width;
   
}

// --------------------------------- //
void Extrapolation()
{
   setTDRStyle();
   gROOT->ForceStyle();

   // input files
   TFile* jet_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/Jet_ReRecoA.root", "READ");
   TFile* jetht_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD.root", "READ");
   TFile* jetmon_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD.root", "READ");
   // TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat.root", "READ");
   //TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_SmearedWithMeasuredValues.root", "READ");
   //TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_ReweightAlphaSpectrum.root", "READ");
   //TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_MBXS73500.root", "READ");
   //TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_JECup.root", "READ");
   //TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_JECdown.root", "READ");
   TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_ChangedPLIDefinitionToGenPtAve.root", "READ");

   //TString suffix = "_nominal";
   //TString suffix = "_MCSmearedWithMeasuredValues";
   //TString suffix = "_ReweightAlphaSpectrum";
   //TString suffix = "_MBXS73500";
   //TString suffix = "_JECup";
   //TString suffix = "_JECdown";
   // TString suffix = "_PLIup";
   //TString suffix = "_PLIdown";
   // TString suffix = "_FirstThreeAlphaPoints";
   TString suffix = "_ChangedPLIDefinitionToGenPtAve";
  
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
   alpha.push_back(0.1); 
   alpha.push_back(0.125); 
   alpha.push_back(0.15); 
   alpha.push_back(0.175); 
   alpha.push_back(0.20); 
   alpha.push_back(0.225);
   alpha.push_back(0.25); 
      
   float pt_bins[14] = {62, 107, 175, 205, 242, 270, 310, 335, 379, 410, 467, 600, 1000, 2000};
   float eta_bins[6] = {0, 0.5, 1.1, 1.7, 2.3, 5.2};
   TH1F *extrapolated_data = new TH1F("extrapolated_data", "extrapolated_data", 13, pt_bins);
   TH1F *extrapolated_mc = new TH1F("extrapolated_mc", "extrapolated_mc", 13, pt_bins);
   TH1F *extrapolated_gen = new TH1F("extrapolated_gen", "extrapolated_gen", 13, pt_bins);
   TH1F *extrapolated_data_with_pli = new TH1F("extrapolated_data_with_pli", "extrapolated_data", 13, pt_bins);
   TH1F *extrapolated_mc_with_pli = new TH1F("extrapolated_mc_with_pli", "extrapolated_mc", 13, pt_bins);
   TH1F *truth_resolution = new TH1F("truth_resolution", "truth_resolution", 13, pt_bins);
   TH1F* RatioVsEta = new TH1F("RatioVsEta", "", 5, eta_bins);
   TH1F* RatioVsEta_with_pli = new TH1F("RatioVsEta_with_pli", "", 5, eta_bins);
   extrapolated_data->Sumw2();
   extrapolated_mc->Sumw2();
   extrapolated_gen->Sumw2();
   extrapolated_data_with_pli->Sumw2();
   extrapolated_mc_with_pli->Sumw2();
   truth_resolution->Sumw2();
   RatioVsEta->Sumw2();
   RatioVsEta_with_pli->Sumw2();

   // how much should be truncated?
   double yq_IQW[2],xq_IQW[2];
   xq_IQW[0] = 0.0;
   xq_IQW[1] = 0.97;
   
   //// get asymmetry histos
   for(int ieta=0; ieta < 5; ++ieta){
      //  cout << "eta Bin: " << ieta << endl;

      extrapolated_data->Reset();
      extrapolated_mc->Reset();
      extrapolated_gen->Reset();
      extrapolated_data_with_pli->Reset();
      extrapolated_mc_with_pli->Reset();
      truth_resolution->Reset();

      for(int ipt=0; ipt < 13; ++ipt){     
         //  cout << "pt Bin: " << ipt << endl;
         std::vector<double> x,x_e,MCy,MCy_e,Datay,Datay_e,Geny,Geny_e;

         for(int ialpha=0; ialpha < 7; ++ialpha){
            //  cout << "alpha Bin: " << ialpha << endl;
            TString hname = Form("Pt%i_eta%i_alpha%i", ipt, ieta, ialpha);
            TString hname_gen = Form("Pt%i_eta%i_genalpha%i", ipt, ieta, ialpha);

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

            double mc_width = GetAsymmWidth(tmp_mc, xq_IQW, yq_IQW);
            double mc_width_err = GetAsymmWidth(tmp_mc, xq_IQW, yq_IQW)/(TMath::Sqrt(2*tmp_mc->GetEffectiveEntries()));
            double data_width = GetAsymmWidth(tmp_data1, xq_IQW, yq_IQW);
            double data_width_err = GetAsymmWidth(tmp_data1, xq_IQW, yq_IQW)/(TMath::Sqrt(2*tmp_data1->GetEffectiveEntries()));
            double gen_width = GetAsymmWidth(tmp_gen, xq_IQW, yq_IQW);
            double gen_width_err = GetAsymmWidth(tmp_gen, xq_IQW, yq_IQW)/(TMath::Sqrt(2*tmp_gen->GetEffectiveEntries()));

          
            MCy.push_back( mc_width );
            MCy_e.push_back( mc_width_err );
            Datay.push_back( data_width );
            Datay_e.push_back( data_width_err );
            Geny.push_back( gen_width );
            Geny_e.push_back( gen_width_err );

            TCanvas *c5 = new TCanvas("c5", "", 600, 600);
            c5->SetLogy();
            tmp_mc->Rebin(10);
            tmp_mc->Scale(tmp_data1->Integral()/tmp_mc->Integral());
            //  tmp_mc->GetYaxis()->SetRangeUser(tmp_mc->GetMinimum(), 100.*tmp_mc->GetMaximum());
            tmp_mc->GetXaxis()->SetTitle("(p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");
            tmp_mc->SetLineColor(38);
            tmp_mc->SetFillColor(38);
            tmp_mc->Draw("hist");
            tmp_data1->Rebin(10);
            tmp_data1->SetMarkerStyle(20);
            tmp_data1->Draw("same");

            TPaveText *label = util::LabelFactory::createPaveTextWithOffset(3,0.8,0.01);
            label->AddText("Anti-k_{T} (R=0.5) PFchs Jets");
            label->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq #bar{ p}_{T} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label->AddText( Form("#alpha #leq %0.3f", alpha.at(ialpha)) );
            label->Draw("same");
            
            TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(2,0.65,0.2);
            leg1->AddEntry(tmp_data1,"Asymmetry (Data)","LP");
            leg1->AddEntry(tmp_mc,"Asymmetry (MC)","LF");
            leg1->Draw();

            if(ieta == 0 && ipt == 0 && ialpha == 0 ) c5->Print("Extrapolation/AsymmHistos" + suffix + ".eps(");
            else if(ieta == 4 && ipt == 12 && ialpha == 6) c5->Print("Extrapolation/AsymmHistos" + suffix + ".eps)");
            else c5->Print("Extrapolation/AsymmHistos" + suffix + ".eps"); 

            TCanvas *c5b = new TCanvas("c5b", "", 600, 600);
            c5b->SetLogy();
            tmp_gen->Rebin(10);
            //  tmp_gen->Scale(tmp_data->Integral()/tmp_gen->Integral());
            //  tmp_gen->GetYaxis()->SetRangeUser(tmp_gen->GetMinimum(), 100.*tmp_gen->GetMaximum());
            tmp_gen->GetXaxis()->SetTitle("(p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");
            tmp_gen->SetLineColor(39);
            tmp_gen->SetFillColor(39);
            tmp_gen->Draw("hist");
           
            TPaveText *label2 = util::LabelFactory::createPaveTextWithOffset(3,0.8,0.01);
            label2->AddText("Anti-k_{T} (R=0.5) PFchs Jets");
            label2->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq #bar{ p, gen}_{T} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label2->AddText( Form("#alpha_{gen} #leq %0.3f", alpha.at(ialpha)) );
            label2->Draw("same");
            
            TLegend* leg2 = util::LabelFactory::createLegendColWithOffset(2,0.65,0.2);
            leg2->AddEntry(tmp_gen,"Asymmetry (Gen)","LF");
            leg2->Draw();

            if(ieta == 0 && ipt == 0 && ialpha == 0 ) c5b->Print("Extrapolation/GenAsymmHistos" + suffix + ".eps(");
            else if(ieta == 4 && ipt == 12 && ialpha == 6) c5b->Print("Extrapolation/GenAsymmHistos" + suffix + ".eps)");
            else c5b->Print("Extrapolation/GenAsymmHistos" + suffix + ".eps"); 
         }

         // Covariance matrices needed for fitting 
         TMatrixD y_cov_mc;
         TMatrixD y_cov_data;
         TMatrixD y_cov_gen;
         y_cov_mc.ResizeTo(alpha.size(), alpha.size());
         y_cov_data.ResizeTo(alpha.size(), alpha.size());
         y_cov_gen.ResizeTo(alpha.size(), alpha.size());

         // fill covariance matrix for data and mc
         for(int ialpha=0; ialpha < alpha.size(); ++ialpha){
            for (Int_t jalpha =0; jalpha < alpha.size(); jalpha++){
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
         //  extrapol_Data->Fit("lin_extrapol_data","Q","same",0,alpha.back()+0.05);
                 
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
         make_lin_fit(slope, d_slope, offset, d_offset);
         std::cout << "MC fit values: " << "slope: " << slope << " offset: " << offset << std::endl; 

         lin_extrapol_mc->SetParameter(0, offset);
         lin_extrapol_mc->SetParError(0, d_offset);
         lin_extrapol_mc->SetParameter(1, slope);
         lin_extrapol_mc->SetParError(1, d_slope);
         extrapol_MC->GetListOfFunctions()->Add(lin_extrapol_mc);
         
         data.reset();
         
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
         make_lin_fit(slope, d_slope, offset, d_offset);
         std::cout << "Data fit values: " << "slope: " << slope << " offset: " << offset << std::endl; 
         
         lin_extrapol_data->SetParameter(0, offset);
         lin_extrapol_data->SetParError(0, d_offset);
         lin_extrapol_data->SetParameter(1, slope);
         lin_extrapol_data->SetParError(1, d_slope);
         extrapol_Data->GetListOfFunctions()->Add(lin_extrapol_data);
         
         data.reset();

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
         make_lin_fit(slope, d_slope, offset, d_offset);
         std::cout << "Gen fit values: " << "slope: " << slope << " offset: " << offset << std::endl; 
         
         lin_extrapol_gen->SetParameter(0, offset);
         lin_extrapol_gen->SetParError(0, d_offset);
         lin_extrapol_gen->SetParameter(1, slope);
         lin_extrapol_gen->SetParError(1, d_slope);
         extrapol_Gen->GetListOfFunctions()->Add(lin_extrapol_gen);
         
         data.reset();

         // draw extrapolations data + mc
         TCanvas *c = new TCanvas("c","",600,600);
         std::pair <float,float> minMaxPair = determineMinMax(extrapol_Data);
         c->DrawFrame(0,minMaxPair.first*0.5-0.05,alpha.back()+0.05,minMaxPair.second*1.47,(";Threshold #alpha_{max};#sigma_{A}"));
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
         label->AddText("Anti-k_{T} (R=0.5) PFchs Jets");
         label->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq #bar{ p}_{T} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
         label->Draw("same");
  
         TLegend* leg1 = util::LabelFactory::createLegendWithOffset(2,0.6);
         leg1->AddEntry(extrapol_Data,"Extrapolation (Data)","LP");
         leg1->AddEntry(extrapol_MC,"Extrapolation (MC)","LP");

         leg1->Draw();
         cmsPrel(-1, false , 8);

         TString name;
         name = Form("Extrapolation/Extrapol_Eta%i_pt%i" + suffix + ".eps", ieta, ipt);
         c->Print(name);

         TCanvas *cb = new TCanvas("c","",600,600);
         std::pair <float,float> minMaxPair2 = determineMinMax(extrapol_Data);
         cb->DrawFrame(0,minMaxPair2.first*0.5-0.05,alpha.back()+0.05,minMaxPair2.second*1.47,(";Threshold #alpha_{max, gen};#sigma_{A, gen}"));
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
         label2->AddText("Anti-k_{T} (R=0.5) PFchs Jets");
         label2->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq #bar{ p}_{T, gen} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
         label2->Draw("same");
  
         TLegend* leg2 = util::LabelFactory::createLegendWithOffset(2,0.15);
         leg2->AddEntry(extrapol_Gen,"Extrapolation (PLI)","LP");

         leg2->Draw();
         cmsPrel(-1, false , 8);

         TString name2;
         name2 = Form("Extrapolation/Extrapol_Eta%i_pt%i_gen" + suffix + ".eps", ieta, ipt);
         cb->Print(name2);

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

         extrapolated_data->SetBinContent(ipt+1, par_data);
         extrapolated_data->SetBinError(ipt+1, par_data_err);
         extrapolated_mc->SetBinContent(ipt+1, par_mc);
         extrapolated_mc->SetBinError(ipt+1, par_mc_err);   
         extrapolated_gen->SetBinContent(ipt+1, par_gen);
         extrapolated_gen->SetBinError(ipt+1, par_gen_err);    

         float par_data_pli_corr = 0;
         float par_data_pli_corr_err = 0;
         float par_mc_pli_corr = 0;
         float par_mc_pli_corr_err = 0;

         if(par_gen > 0 && par_data > 0 && par_mc > 0 && par_data > par_gen && par_mc > par_gen) {
            par_data_pli_corr = TMath::Sqrt(pow(par_data,2) - pow(par_gen,2));
            par_data_pli_corr_err = TMath::Sqrt( pow(par_data,2)/(pow(par_data,2) - pow(par_gen,2)) * pow(par_data_err,2) +  pow(par_gen,2)/(pow(par_data,2) - pow(par_gen,2)) * pow(par_gen_err,2));
            par_mc_pli_corr = TMath::Sqrt(pow(par_mc,2) - pow(par_gen,2));
            par_mc_pli_corr_err = TMath::Sqrt( pow(par_mc,2)/(pow(par_mc,2) - pow(par_gen,2)) * pow(par_mc_err,2) +  pow(par_gen,2)/(pow(par_mc,2) - pow(par_gen,2)) * pow(par_gen_err,2));
         }
            
         extrapolated_data_with_pli->SetBinContent(ipt+1, par_data_pli_corr);
         extrapolated_data_with_pli->SetBinError(ipt+1, par_data_pli_corr_err);
         extrapolated_mc_with_pli->SetBinContent(ipt+1, par_mc_pli_corr);
         extrapolated_mc_with_pli->SetBinError(ipt+1, par_mc_pli_corr_err);

         cout << "Parameter data after pli: " << par_data_pli_corr << endl;
         cout << "Parameter error data after pli: " << par_data_pli_corr_err << endl;
         cout << "Parameter mc after pli: " << par_mc_pli_corr << endl;
         cout << "Parameter error mc after pli: " << par_mc_pli_corr_err << endl;

         // get mc truth resolution
         TString hname_res = Form("Response_Pt%i_eta%i", ipt, ieta);
         //TString hname_res = Form("ResponseRecoPtAve_Pt%i_eta%i", ipt, ieta);
    
         mc_file->cd();
         tmp_res = 0;
         tmp_res = (TH1F*) gDirectory->FindObjectAny(hname_res);
         tmp_res->Rebin(10);

         double truth_res = 0; 
         double truth_res_err = 0;
         double chi2 = 0;
         double ndf = 0;
         double mean = 0;
         double mean_err = 0;
         if( tmp_res->GetEntries() > 200 ) {
            mean = tmp_res->GetMean(); 
            truth_res = tmp_res->GetRMS();
            tmp_res->Fit("gaus","QNI","", mean - 3 * truth_res,mean + 3 * truth_res);
            
            TF1 *f = (TF1*)gROOT->GetFunction("gaus")->Clone();
            f->SetLineColor(kRed);
            f->SetLineWidth(3);
            mean = f->GetParameter(1);
            truth_res = f->GetParameter(2);
            
            if( (tmp_res->Fit(f,"QI","same",mean - 1.5 * truth_res, mean + 1.5 * truth_res) == 0) ) { 
               mean = f->GetParameter(1);
               mean_err = f->GetParError(1);
               truth_res = f->GetParameter(2);
               truth_res_err = f->GetParError(2);
               chi2 = f->GetChisquare();
               ndf = f->GetNDF();
            }
         }

         if(mean && truth_res) {
            truth_res = truth_res/mean;
            truth_res_err = TMath::Sqrt(pow(1/mean, 2) * pow(truth_res_err,2) + pow(truth_res/(mean*mean), 2) * pow(mean_err,2));
         }

         double fit_quality = chi2/ndf;
         
         truth_resolution->SetBinContent(ipt+1, truth_res);
         truth_resolution->SetBinError(ipt+1, truth_res_err);

         TCanvas *res = new TCanvas("res", "", 600, 600);
         res->SetLogy();
         tmp_res->GetXaxis()->SetTitle("p_{T, reco}/p_{T, gen}");
         // tmp_res->SetMarkerSize(0.7);
         tmp_res->Draw();

         TPaveText *label3 = util::LabelFactory::createPaveTextWithOffset(3,1.0,0.5);
         label3->AddText("Anti-k_{T} (R=0.5) PFchs Jets");
         //   label3->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq p_{T, gen} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
         label3->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq #bar{ p}_{T} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
         label3->AddText( Form("#chi^{2}/ndf : %0.2f", fit_quality) );
         label3->Draw("same");

         if(ieta == 0 && ipt == 0 ) res->Print("Extrapolation/TruthResponse" + suffix + ".eps(");
         else if(ieta == 4 && ipt == 12) res->Print("Extrapolation/TruthResponse" + suffix + ".eps)");
         else res->Print("Extrapolation/TruthResponse" + suffix + ".eps"); 
      }

      // fit function for truth resolution
      TF1 *res_func = new TF1("res_func", "TMath::Sqrt(sign([0])*pow([0]/x,2)+pow([1],2)*pow(x,[2]-1.)+pow([3],2))" );
      res_func->SetParameters(1.2, 0.5, 0.03, 0.1);
    
      // draw res after extrapolations
      TCanvas *c2 = new TCanvas("c2","",600,600);
      c2->SetLogx();
      c2->SetBottomMargin(0.25 + 0.75*c2->GetBottomMargin()-0.25*c2->GetTopMargin());
      c2->cd();
      //  truth_resolution->GetXaxis()->SetTitle("#bar{ p}_{T} [GeV]");
      truth_resolution->GetXaxis()->SetLabelSize(0);
      truth_resolution->GetYaxis()->SetRangeUser(0, 0.3);
      truth_resolution->GetYaxis()->SetTitle("#sqrt{2}#sigma_{A}");
      truth_resolution->SetMarkerStyle(26);
      truth_resolution->SetMarkerColor(kRed+1);
      truth_resolution->Draw();
      truth_resolution->Fit("res_func", "WLQ", "same", truth_resolution->GetBinCenter(truth_resolution->FindFirstBinAbove()), truth_resolution->GetBinCenter(truth_resolution->FindLastBinAbove()) );
      // extrapolated_mc->Draw("hist");
      // extrapolated_gen->Scale(TMath::Sqrt(2));
      //extrapolated_gen->SetLineColor(kBlue+1);
      //  extrapolated_gen->Draw("histsame");
      //  extrapolated_data->Draw("same");
      extrapolated_mc_with_pli->Scale(TMath::Sqrt(2));
      extrapolated_mc_with_pli->SetMarkerStyle(25);
      extrapolated_mc_with_pli->SetMarkerColor(kBlue+2);
      extrapolated_mc_with_pli->Draw("same");
      extrapolated_data_with_pli->Scale(TMath::Sqrt(2));
      extrapolated_data_with_pli->SetMarkerColor(41);
      extrapolated_data_with_pli->SetMarkerStyle(27);
      //  extrapolated_data_with_pli->Draw("same");

      TLegend* leg3 = util::LabelFactory::createLegendWithOffset(2,0.175);
      leg3->AddEntry(truth_resolution,"Truth resolution","P");
      leg3->AddEntry(extrapolated_mc_with_pli,"Measured resolution","P");
      // leg3->AddEntry(extrapolated_data_with_pli,"Data resolution","P");
      
      leg3->Draw();
      cmsPrel(-1, false , 8);

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
      r->SetXTitle("#bar{ p}_{T} [GeV]");
      r->SetYTitle("#sigma - #sigma_{truth} ) / #sigma");
      r->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
      r->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y"));
      r->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.3);
      r->GetYaxis()->SetNdivisions(505);
      r->SetStats(0);
      r->SetMarkerStyle(20);
      //      r->SetMarkerSize(1.12);
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
      label->AddText("Anti-k_{T} (R=0.5) PFchs Jets");
      label->AddText( Form("%0.1f #leq |#eta| #leq %0.1f", eta_bins[ieta], eta_bins[ieta+1]) );
      label->Draw("same");

      TString name2;
      name2 = Form("Extrapolation/Extrapol_Eta%i" + suffix + ".eps", ieta);
      c2->Print(name2);

      // --------------------------------------- //
      // calc data/mc ratio and fit with constant
      TH1F* ratio = new TH1F(*extrapolated_mc);
      TH1F* ratio_with_pli = new TH1F(*extrapolated_mc);
      ratio->Divide(extrapolated_data, extrapolated_mc, 1, 1);
      ratio_with_pli->Divide(extrapolated_data_with_pli, extrapolated_mc_with_pli, 1, 1);
    
      TF1 *fit_const = new TF1("fit_const", "[0]", ratio->GetXaxis()->GetXmin(), ratio->GetXaxis()->GetXmax());
      fit_const->SetFillColor(kGray);
      // fit_const->SetFillStyle(1001);
      fit_const->SetParameters(0, 1.1);
      fit_const->SetParName(0, "const");
      ratio->Fit("fit_const", "", "same");
      ratio->GetXaxis()->SetTitle("#bar{ p}_{T} [GeV]");
      ratio->GetYaxis()->SetRangeUser(0.8, 1.5);
      ratio->GetYaxis()->SetTitle("Data/MC ratio (const fit)");

      RatioVsEta->SetBinContent(ieta+1, ratio->GetFunction("fit_const")->GetParameter(0));
      RatioVsEta->SetBinError(ieta+1, ratio->GetFunction("fit_const")->GetParError(0));
      
      //Create a histogram to hold the confidence intervals
      TH1D *hint = new TH1D("hint","Fitted function with .95 conf.band", 100, ratio->GetXaxis()->GetXmin(), ratio->GetXaxis()->GetXmax());
      hint->Sumw2();
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint, 0.68);
      //Now the "hint" histogram has the fitted function values as the 
      //bin contents and the confidence intervals as bin errors
      hint->SetMarkerStyle(1);
      hint->SetLineColor(kRed);
      hint->SetFillColor(kGray);
    
      TCanvas *c3 = new TCanvas("c3","",600,600);
      c3->SetLogx();
      ratio->Draw();
      hint->DrawClone("e4 same");
      ratio->Draw("same");

      label->Draw("same");

      TLegend* leg4 = util::LabelFactory::createLegendWithOffset(2,0.175);
      leg4->AddEntry(ratio,"Measured Ratio","P");
      leg4->AddEntry(hint,"Constant Fit","LF");
      // leg4->AddEntry(extrapolated_data_with_pli,"Data resolution","P");

      leg4->Draw("same");

      TString name3;
      name3 = Form("Extrapolation/ExtrapolRatio_Eta%i" + suffix + ".eps", ieta);
      c3->Print(name3);

      ratio_with_pli->Fit("fit_const", "", "same");
      ratio_with_pli->GetXaxis()->SetTitle("#bar{ p}_{T} [GeV]");
      ratio_with_pli->GetYaxis()->SetRangeUser(0.8, 1.5);
      ratio_with_pli->GetYaxis()->SetTitle("Data/MC ratio (const fit)");

      RatioVsEta_with_pli->SetBinContent(ieta+1, ratio_with_pli->GetFunction("fit_const")->GetParameter(0));
      RatioVsEta_with_pli->SetBinError(ieta+1, ratio_with_pli->GetFunction("fit_const")->GetParError(0));

      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint, 0.68);
    
      TCanvas *c3b = new TCanvas("c3","",600,600);
      c3b->SetLogx();
      ratio_with_pli->Draw();
      hint->DrawClone("e4 same");
      ratio_with_pli->Draw("same");
      leg4->Draw("same");

      label->Draw("same");

      TString name4;
      name4 = Form("Extrapolation/ExtrapolRatio_Eta%i_with_pli" + suffix + ".eps", ieta);
      c3b->Print(name4);
   }

   // draw data/mc scaling factors vs. eta
   TCanvas *c4 = new TCanvas();
   RatioVsEta->GetYaxis()->SetRangeUser(0.7, 1.3);
   RatioVsEta->GetXaxis()->SetTitle("|#eta|");
   RatioVsEta->GetYaxis()->SetTitle("Data /MC ratio (const fit)");
   RatioVsEta->Draw();
   c4->Print("Extrapolation/ScalingFactorsVsEta" + suffix + ".eps");

   TCanvas *c4b = new TCanvas();
   RatioVsEta_with_pli->GetYaxis()->SetRangeUser(0.7, 1.4);
   RatioVsEta_with_pli->GetXaxis()->SetTitle("|#eta|");
   RatioVsEta_with_pli->GetYaxis()->SetTitle("Data /MC ratio (const fit)");
   RatioVsEta_with_pli->Draw();
   c4b->Print("Extrapolation/ScalingFactorsVsEta_with_pli" + suffix + ".eps");

   cout << "Ratio eta1: " << RatioVsEta_with_pli->GetBinContent(1) << " +- " << RatioVsEta_with_pli->GetBinError(1) << endl;
   cout << "Ratio eta2: " << RatioVsEta_with_pli->GetBinContent(2) << " +- " << RatioVsEta_with_pli->GetBinError(2) << endl;
   cout << "Ratio eta3: " << RatioVsEta_with_pli->GetBinContent(3) << " +- " << RatioVsEta_with_pli->GetBinError(3) << endl;
   cout << "Ratio eta4: " << RatioVsEta_with_pli->GetBinContent(4) << " +- " << RatioVsEta_with_pli->GetBinError(4) << endl;
   cout << "Ratio eta5: " << RatioVsEta_with_pli->GetBinContent(5) << " +- " << RatioVsEta_with_pli->GetBinError(5) << endl;

   TFile* output = new TFile("Extrapolation/JER_RatioVsEta" + suffix + ".root", "RECREATE");
   RatioVsEta_with_pli->Write();

   output->Write();

}



