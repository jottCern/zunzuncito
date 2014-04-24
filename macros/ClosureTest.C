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
   
   if( htemp->GetEntries() > 100 ) {
      htemp->GetXaxis()->SetRange(0,-1);
      htemp->ComputeIntegral();
      htemp->GetQuantiles(nq,yq_IQW,xq_IQW);
      Int_t IQW_high_bin_i = htemp->FindBin(yq_IQW[1]);
    //   cout << "x value: " << htemp->GetBinCenter(IQW_high_bin_i) << endl;
//       cout << "Bin where to cut: " << IQW_high_bin_i << endl;
//       cout << "Bins before truncation: " << htemp->GetNbinsX() << endl;
//       cout << "Integral before truncation: " << htemp->Integral() << endl;
//       cout << "Integral after truncation: " << htemp->Integral(1, IQW_high_bin_i) << endl;
      cout << "Truncated Integral in %: " << htemp->Integral(1, IQW_high_bin_i)/htemp->Integral() << endl;

      for(int i=1; i <= IQW_high_bin_i; i++) {
         width += htemp->GetBinContent(i)* std::pow(htemp->GetBinCenter(i), 2);
      }

      width = TMath::Sqrt(1/(htemp->Integral(1, IQW_high_bin_i))*width);
   }

   return width;  
}

// --------------------------------- //
void ClosureTest()
{
   setTDRStyle();
   gROOT->ForceStyle();

   // input files
   // TFile* mc_smeared_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_ClosureSecondHalf_v2.root", "READ");

   /* TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_v2.root", "READ");
      TFile* mc_smeared_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_ClosureFirstHalf_v2.root", "READ");*/

   TFile* mc_smeared_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_v2.root", "READ");
   TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_v2.root", "READ");

   //TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final.root", "READ");

   //TString suffix = "_PythiaSmearedVsPythia";
   TString suffix = "_HerwigVsPythia";
   // TString suffix = "_PythiaVsHerwig";

 
   // define helper histos
   TH1F *tmp_mcsmeared = new TH1F();
   TH1F *tmp_mc = new TH1F();
   TH1F *tmp_gensmeared = new TH1F();
   TH1F *tmp_gen = new TH1F();
   tmp_mcsmeared->Sumw2();
   tmp_mc->Sumw2();
   tmp_gensmeared->Sumw2();
   tmp_gen->Sumw2();
 
   std::vector<float> alpha;
   alpha.push_back(0.1); 
   alpha.push_back(0.125); 
   alpha.push_back(0.15); 
   alpha.push_back(0.175); 
   alpha.push_back(0.20); 
   alpha.push_back(0.225);
   alpha.push_back(0.25); 
      
   float pt_bins[14] = {62, 107, 175, 205, 242, 270, 310, 335, 379, 410, 467, 600, 1000, 2000};
   float eta_bins[6] = {0, 0.5, 1.1, 1.7, 2.3, 5.0};
   TH1F *extrapolated_mcsmeared = new TH1F("extrapolated_mcsmeared", "extrapolated_mcsmeared", 13, pt_bins);
   TH1F *extrapolated_mc = new TH1F("extrapolated_mc", "extrapolated_mc", 13, pt_bins);
   TH1F *extrapolated_gen = new TH1F("extrapolated_gen", "extrapolated_gen", 13, pt_bins);
   TH1F *extrapolated_mcsmeared_with_pli = new TH1F("extrapolated_mcsmeared_with_pli", "extrapolated_mcsmeared", 13, pt_bins);
   TH1F *extrapolated_mc_with_pli = new TH1F("extrapolated_mc_with_pli", "extrapolated_mc", 13, pt_bins);

   TH1F* RatioVsEta = new TH1F("RatioVsEta", "", 5, eta_bins);
   TH1F* RatioVsEta_with_pli = new TH1F("RatioVsEta_with_pli", "", 5, eta_bins);
   extrapolated_mcsmeared->Sumw2();
   extrapolated_mc->Sumw2();
   extrapolated_gen->Sumw2();
   extrapolated_mcsmeared_with_pli->Sumw2();
   extrapolated_mc_with_pli->Sumw2();
 
   RatioVsEta->Sumw2();
   RatioVsEta_with_pli->Sumw2();

   // how much should be truncated?
   double yq_IQW[2],xq_IQW[2];
   xq_IQW[0] = 0.0;
   xq_IQW[1] = 0.985;
   
   //// get asymmetry histos
   for(int ieta=0; ieta < 5; ++ieta){
      //  cout << "eta Bin: " << ieta << endl;

      extrapolated_mcsmeared->Reset();
      extrapolated_mc->Reset();
      extrapolated_gen->Reset();
      extrapolated_mcsmeared_with_pli->Reset();
      extrapolated_mc_with_pli->Reset();
    
      for(int ipt=0; ipt < 13; ++ipt){     
         //  cout << "pt Bin: " << ipt << endl;
         std::vector<double> x,x_e,MCy,MCy_e,Datay,Datay_e,Geny,Geny_e;

         for(int ialpha=0; ialpha < 7; ++ialpha){
         //for(int ialpha=0; ialpha < 3; ++ialpha){
            //  cout << "alpha Bin: " << ialpha << endl;
            TString hname = Form("Pt%i_eta%i_alpha%i", ipt, ieta, ialpha);
            TString hname_gen = Form("GenAsymm_Pt%i_eta%i_alpha%i", ipt, ieta, ialpha);

            //  cout << "hname: " << hname << endl;

            mc_file->cd();
            tmp_mc = 0;
            tmp_gen = 0;
            tmp_mc = (TH1F*) gDirectory->FindObjectAny(hname);
            tmp_gen = (TH1F*) gDirectory->FindObjectAny(hname_gen);

            mc_smeared_file->cd();
            tmp_mcsmeared = 0;
            tmp_gensmeared = 0;
            tmp_mcsmeared = (TH1F*) gDirectory->FindObjectAny(hname);
            tmp_gensmeared = (TH1F*) gDirectory->FindObjectAny(hname_gen);
                               
            x.push_back(alpha.at(ialpha));
            x_e.push_back(0.);

            double mc_width = GetAsymmWidth(tmp_mc, xq_IQW, yq_IQW);
            double mc_width_err = GetAsymmWidth(tmp_mc, xq_IQW, yq_IQW)/(TMath::Sqrt(2*tmp_mc->GetEffectiveEntries()));
            double data_width = GetAsymmWidth(tmp_mcsmeared, xq_IQW, yq_IQW);
            double data_width_err = GetAsymmWidth(tmp_mcsmeared, xq_IQW, yq_IQW)/(TMath::Sqrt(2*tmp_mcsmeared->GetEffectiveEntries()));
            double gen_width = GetAsymmWidth(tmp_gen, xq_IQW, yq_IQW);
            double gen_width_err = GetAsymmWidth(tmp_gen, xq_IQW, yq_IQW)/(TMath::Sqrt(2*tmp_gen->GetEffectiveEntries()));

            MCy.push_back( mc_width );
            MCy_e.push_back( mc_width_err );
            Datay.push_back( data_width );
            Datay_e.push_back( data_width_err );
            Geny.push_back( gen_width );
            Geny_e.push_back( gen_width_err );

            tmp_mc->Scale(tmp_mcsmeared->Integral()/tmp_mc->Integral());
            tmp_gen->Scale(tmp_gensmeared->Integral()/tmp_gen->Integral());

            // draw asymmetry histos
            TCanvas *c5 = new TCanvas("c5", "", 600, 600);
            c5->SetLogy();
            //  tmp_mc->GetYaxis()->SetRangeUser(0.1, 100.*tmp_mc->GetMaximum());
            // tmp_mc->GetXaxis()->SetTitle("(p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");
            tmp_mc->GetXaxis()->SetTitle("|A|");
            tmp_mc->GetYaxis()->SetTitle("Events");
            tmp_mc->SetLineColor(30);
            //  tmp_mc->SetFillColor(30);
            tmp_mc->Rebin(10);
            tmp_mc->Draw("hist");
            //   gauss_mc->Draw("same");
            tmp_mcsmeared->Rebin(10);
            tmp_mcsmeared->SetMarkerStyle(20);
            tmp_mcsmeared->Draw("histsame");
          
            TPaveText *label = util::LabelFactory::createPaveTextWithOffset(3,0.8,0.01);
            label->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
            label->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq  p_{T}^{ave} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label->AddText( Form("#alpha #leq %0.3f", alpha.at(ialpha)) );
            label->Draw("same");
            
            TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(2,0.65,0.2);
            leg1->AddEntry(tmp_mcsmeared,"Herwig","L");
            leg1->AddEntry(tmp_mc,"Pythia","L");
            leg1->Draw();

            if(ieta == 0 && ipt == 0 && ialpha == 0 ) c5->Print("ClosureTest/AsymmHistos" + suffix + ".eps(");
            else if(ieta == 4 && ipt == 12 && ialpha == 6) c5->Print("ClosureTest/AsymmHistos" + suffix + ".eps)");
            else c5->Print("ClosureTest/AsymmHistos" + suffix + ".eps"); 

            // draw gen-asymmetry histos
            TCanvas *c5b = new TCanvas("c5b", "", 600, 600);
            c5b->SetLogy();
            //  tmp_gen->GetYaxis()->SetRangeUser(0.1, 100.*tmp_gen->GetMaximum());
            // tmp_gen->GetXaxis()->SetTitle("(p_{T,1} - p_{T,2})/(p_{T,1} + p_{T,2})");
            tmp_gen->GetXaxis()->SetTitle("|A_{gen}|");
            tmp_gen->GetYaxis()->SetTitle("Events");
            tmp_gen->SetLineColor(30);
            //  tmp_gen->SetFillColor(30);
            tmp_gen->Rebin(10);
            tmp_gen->Draw("hist");
            //   gauss_gen->Draw("same");
            tmp_gensmeared->Rebin(10);
            tmp_gensmeared->SetMarkerStyle(20);
            tmp_gensmeared->Draw("histsame");
          
            TPaveText *label2 = util::LabelFactory::createPaveTextWithOffset(3,0.8,0.01);
            label2->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
            label2->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq  p_{T}^{ave} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
            label2->AddText( Form("#alpha #leq %0.3f", alpha.at(ialpha)) );
            label2->Draw("same");
            
            TLegend* leg2 = util::LabelFactory::createLegendColWithOffset(2,0.65,0.2);
            leg2->AddEntry(tmp_gensmeared,"Herwig","L");
            leg2->AddEntry(tmp_gen,"Pythia","L");
            leg2->Draw();

            if(ieta == 0 && ipt == 0 && ialpha == 0 ) c5b->Print("ClosureTest/GenAsymmHistos" + suffix + ".eps(");
            else if(ieta == 4 && ipt == 12 && ialpha == 6) c5b->Print("ClosureTest/GenAsymmHistos" + suffix + ".eps)");
            else c5b->Print("ClosureTest/GenAsymmHistos" + suffix + ".eps"); 
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
         leg1->AddEntry(extrapol_Data,"Extrapolation (MC smeared)","LP");
         leg1->AddEntry(extrapol_MC,"Extrapolation (MC)","LP");

         leg1->Draw();
         cmsPrel(-1, false , 8);

         TString name;
         name = Form("ClosureTest/Extrapol_Eta%i_pt%i" + suffix + ".eps", ieta, ipt);
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
         label2->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq #bar{ p}_{T} [GeV] #leq %3.0f", eta_bins[ieta], eta_bins[ieta+1], pt_bins[ipt], pt_bins[ipt+1]) );
         label2->Draw("same");
  
         TLegend* leg2 = util::LabelFactory::createLegendWithOffset(2,0.15);
         leg2->AddEntry(extrapol_Gen,"Extrapolation (PLI)","LP");

         leg2->Draw();
         cmsPrel(-1, false , 8);

         TString name2;
         name2 = Form("ClosureTest/Extrapol_Eta%i_pt%i_gen" + suffix + ".eps", ieta, ipt);
         cb->Print(name2);

         float par_data = lin_extrapol_data->GetParameter(0);
         float par_data_err = lin_extrapol_data->GetParError(0);
         float par_mc = lin_extrapol_mc->GetParameter(0);
         float par_mc_err = lin_extrapol_mc->GetParError(0);
         float par_gen = lin_extrapol_gen->GetParameter(0);
         float par_gen_err = lin_extrapol_gen->GetParError(0);
       
         cout << "ieta : " << ieta << "  ipt : " << ipt << endl;
         cout << "Parameter data: " << par_data << endl;
         cout << "Parameter error data: " << par_data_err << endl;
         cout << "Parameter mc: " << par_mc << endl;
         cout << "Parameter error mc: " << par_mc_err << endl;
         cout << "Parameter gen: " << par_gen << endl;
         cout << "Parameter error gen: " << par_gen_err << endl;

         extrapolated_mcsmeared->SetBinContent(ipt+1, par_data);
         extrapolated_mcsmeared->SetBinError(ipt+1, par_data_err);
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
            
         extrapolated_mcsmeared_with_pli->SetBinContent(ipt+1, par_data_pli_corr);
         extrapolated_mcsmeared_with_pli->SetBinError(ipt+1, par_data_pli_corr_err);
         extrapolated_mc_with_pli->SetBinContent(ipt+1, par_mc_pli_corr);
         extrapolated_mc_with_pli->SetBinError(ipt+1, par_mc_pli_corr_err);

         cout << "Parameter data after pli: " << par_data_pli_corr << endl;
         cout << "Parameter error data after pli: " << par_data_pli_corr_err << endl;
         cout << "Parameter mc after pli: " << par_mc_pli_corr << endl;
         cout << "Parameter error mc after pli: " << par_mc_pli_corr_err << endl;
      }

      // --------------------------------------- //
      // calc data/mc ratio and fit with constant
      TH1F* ratio = new TH1F(*extrapolated_mc);
      TH1F* ratio_with_pli = new TH1F(*extrapolated_mc);
      ratio->Divide(extrapolated_mcsmeared, extrapolated_mc, 1, 1);
      ratio_with_pli->Divide(extrapolated_mcsmeared_with_pli, extrapolated_mc_with_pli, 1, 1);
    
      TF1 *fit_const = new TF1("fit_const", "[0]", ratio->GetXaxis()->GetXmin(), ratio->GetXaxis()->GetXmax());
      fit_const->SetParameters(0, 1.1);
      fit_const->SetParName(0, "const");
      ratio->Fit("fit_const", "", "same");
      ratio->GetXaxis()->SetTitle("#bar{ p}_{T} [GeV]");
      ratio->GetYaxis()->SetRangeUser(0.7, 1.4);
      ratio->GetYaxis()->SetTitle("MC_{smeared}/MC ratio (const fit)");

      RatioVsEta->SetBinContent(ieta+1, ratio->GetFunction("fit_const")->GetParameter(0));
      RatioVsEta->SetBinError(ieta+1, ratio->GetFunction("fit_const")->GetParError(0));

      TCanvas *c3 = new TCanvas("c3","",600,600);
      c3->SetLogx();
      ratio->Draw();
      TString name3;
      name3 = Form("ClosureTest/ExtrapolRatio_Eta%i" + suffix + ".eps", ieta);
      c3->Print(name3);

      ratio_with_pli->Fit("fit_const", "", "same");
      ratio_with_pli->GetXaxis()->SetTitle("#bar{ p}_{T} [GeV]");
      ratio_with_pli->GetYaxis()->SetRangeUser(0.7, 1.4);
      ratio_with_pli->GetYaxis()->SetTitle("MC_{smeared}/MC ratio (const fit)");
    
      TCanvas *c3b = new TCanvas("c3","",600,600);
      c3b->SetLogx();
      ratio_with_pli->Draw();
      TString name4;
      name4 = Form("ClosureTest/ExtrapolRatio_Eta%i_with_pli" + suffix + ".eps", ieta);
      c3b->Print(name4);

      RatioVsEta_with_pli->SetBinContent(ieta+1, ratio_with_pli->GetFunction("fit_const")->GetParameter(0));
      RatioVsEta_with_pli->SetBinError(ieta+1, ratio_with_pli->GetFunction("fit_const")->GetParError(0));

   }

   // draw data/mc scaling factors vs. eta
   TCanvas *c4 = new TCanvas();
   RatioVsEta->GetYaxis()->SetRangeUser(0.7, 1.3);
   RatioVsEta->GetXaxis()->SetTitle("|#eta|");
   RatioVsEta->GetYaxis()->SetTitle("MC_{smeared} /MC ratio (const fit)");
   RatioVsEta->Draw();
   c4->Print("ClosureTest/ScalingFactorsVsEta" + suffix + ".eps");

   cout << "//----------------------------------------------//" << endl;
   cout << "Scaling factors without PLI: " << endl;
   cout << "Ratio eta1: " << RatioVsEta->GetBinContent(1) << " +- " << RatioVsEta->GetBinError(1) << endl;
   cout << "Ratio eta2: " << RatioVsEta->GetBinContent(2) << " +- " << RatioVsEta->GetBinError(2) << endl;
   cout << "Ratio eta3: " << RatioVsEta->GetBinContent(3) << " +- " << RatioVsEta->GetBinError(3) << endl;
   cout << "Ratio eta4: " << RatioVsEta->GetBinContent(4) << " +- " << RatioVsEta->GetBinError(4) << endl;
   cout << "Ratio eta5: " << RatioVsEta->GetBinContent(5) << " +- " << RatioVsEta->GetBinError(5) << endl;
   cout << "//----------------------------------------------//" << endl;

   TCanvas *c4b = new TCanvas();
   RatioVsEta_with_pli->GetYaxis()->SetRangeUser(0.7, 1.4);
   RatioVsEta_with_pli->GetXaxis()->SetTitle("|#eta|");
   RatioVsEta_with_pli->GetYaxis()->SetTitle("MC_{smeared} /MC ratio (const fit)");
   RatioVsEta_with_pli->Draw();
   c4b->Print("ClosureTest/ScalingFactorsVsEta_with_pli" + suffix + ".eps");

   cout << "//----------------------------------------------//" << endl;
   cout << "Scaling factors with PLI: " << endl;
   cout << "Ratio eta1: " << RatioVsEta_with_pli->GetBinContent(1) << " +- " << RatioVsEta_with_pli->GetBinError(1) << endl;
   cout << "Ratio eta2: " << RatioVsEta_with_pli->GetBinContent(2) << " +- " << RatioVsEta_with_pli->GetBinError(2) << endl;
   cout << "Ratio eta3: " << RatioVsEta_with_pli->GetBinContent(3) << " +- " << RatioVsEta_with_pli->GetBinError(3) << endl;
   cout << "Ratio eta4: " << RatioVsEta_with_pli->GetBinContent(4) << " +- " << RatioVsEta_with_pli->GetBinError(4) << endl;
   cout << "Ratio eta5: " << RatioVsEta_with_pli->GetBinContent(5) << " +- " << RatioVsEta_with_pli->GetBinError(5) << endl;
   cout << "//----------------------------------------------//" << endl;

}



