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

void MakePlot_Alpha(TFile *jet_file, TFile *jetht_file, TFile *jetmon_file, TFile *mc_file, TString xTitle, TString Label, TString histoName, int pt_bin, int eta_bin)
{
   float pt_bins[14] = {62, 107, 175, 205, 242, 270, 310, 335, 379, 410, 467, 600, 1000, 2000};
   float eta_bins[6] = {0, 0.5, 1.1, 1.7, 2.3, 5.2};

   TH1F *tmp_data = new TH1F();
   TH1F *tmp_mc = new TH1F();
   TH1F *tmp_help = new TH1F();
   
   jet_file->cd();  
   tmp_help = 0;
   tmp_help = (TH1F*) gDirectory->FindObjectAny(histoName);
   tmp_data = tmp_help;
      
   jetht_file->cd(); 
   tmp_help = 0;
   tmp_help = (TH1F*) gDirectory->FindObjectAny(histoName);
   tmp_data->Add(tmp_help);

   jetmon_file->cd();  
   tmp_help = 0;
   tmp_help = (TH1F*) gDirectory->FindObjectAny(histoName);
   tmp_data->Add(tmp_help);
      
   mc_file->cd();  
   tmp_help = 0;
   tmp_mc= (TH1F*) gDirectory->FindObjectAny(histoName);

   TCanvas *c = new TCanvas();
   c->SetLogy();
   c->SetBottomMargin(0.25 + 0.75*c->GetBottomMargin()-0.25*c->GetTopMargin());
   c->cd();
   tmp_data->GetXaxis()->SetLabelSize(0);
   tmp_data->GetXaxis()->SetTitle("");

   if(histoName.Contains("Alpha")) {
      tmp_data->Rebin(25);
      tmp_mc->Rebin(25);
   }

   tmp_data->GetYaxis()->SetTitle("Events");
   tmp_data->SetMarkerStyle(26);
   tmp_data->SetMarkerColor(kRed);
   tmp_data->Draw();
   tmp_mc->Scale(tmp_data->Integral()/tmp_mc->Integral());
   tmp_mc->SetLineColor(kBlue+2);
   tmp_mc->Draw("samehist");
      
   TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2,0.97,0.35);
   label->AddText(Label);
   label->AddText( Form("%0.1f #leq |#eta| #leq %0.1f, %3.0f #leq #bar{ p}_{T} [GeV] #leq %3.0f", eta_bins[eta_bin], eta_bins[eta_bin+1], pt_bins[pt_bin], pt_bins[pt_bin+1]) );
   label->Draw("same");
      
   TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(2,0.97,0.45);
   leg1->AddEntry(tmp_data,"Data","P");
   leg1->AddEntry(tmp_mc,"MC","LF");
   leg1->Draw();
      
   tmp_data->Draw("same");
   tmp_mc->Draw("samehist");

   TPad *pad = new TPad("pad", "pad", 0, 0, 1, 1);
   pad->SetTopMargin(0.75 - 0.75*pad->GetBottomMargin()+0.25*pad->GetTopMargin());
   pad->SetFillStyle(0);
   pad->SetFrameFillColor(10);
   pad->SetFrameBorderMode(0);
   pad->Draw();
   pad->cd();

   Double_t xMin1 = tmp_data->GetXaxis()->GetXmin();
   Double_t xMax1 = tmp_data->GetXaxis()->GetXmax();

   TH1F* r = new TH1F(*tmp_data);
   r->Sumw2();
   r->SetXTitle(xTitle);
   r->GetYaxis()->CenterTitle();
   // r->SetYTitle("(Data-MC)/Data");
   r->SetYTitle("Data/MC");
   r->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
   r->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y"));
   r->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.3);
   r->GetYaxis()->SetNdivisions(505);
   r->SetStats(0);
   r->SetMarkerStyle(20);
   //      r->SetMarkerSize(1.12);
   r->SetMarkerColor(kBlack);
   r->Reset();
   r->Add(tmp_data, 1);
   //  r->Add(tmp_mc, -1);
   // r->Divide(tmp_data);
   r->Divide(tmp_mc);
   r->SetMaximum(1.5);
   r->SetMinimum(0.5);
   r->Draw("ep");

   // define fit function
   TF1 *fit = new TF1("fit","0.5*[0]*(TMath::Erf([1]*x-[2])+1)", 0, 0.3);
   if( TFitResult *res = r->Fit( fit, "R0S" ).Get() ) {
      TF1 *f = new TF1( *r->GetFunction("fit") );
      f->SetLineColor( kBlue );
      f->SetLineWidth( 0.5 );
      f->Draw("same");
   }     

   TF1 *f2 = new TF1("func", "0.5*[0]*(TMath::Erf([1]*x-[2])+1)", 0, 0.3);
   f2->SetParameter(0, 1.15);
   f2->SetParameter(1, 15);
   f2->SetParameter(2, 0.1);

   f2->SetLineColor(kRed);
   f2->Draw("same");

   cout << "Par 0: " << fit->GetParameter(0) << endl; 
   cout << "Par 1: " << fit->GetParameter(1) << endl; 
   cout << "Par 2: " << fit->GetParameter(2) << endl; 

   TLine l;
   l.DrawLine(xMin1, 1., xMax1, 1.);
   c->cd();
      
   c->Print("ControlPlots/" + histoName + ".eps");
}

void MakePlot(TFile *jet_file, TFile *jetht_file, TFile *jetmon_file, TFile *mc_file, TString xTitle, TString Label, TString histoName)
{
   TH1F *tmp_data = new TH1F();
   TH1F *tmp_mc = new TH1F();
   TH1F *tmp_help = new TH1F();
   
   jet_file->cd();  
   tmp_help = 0;
   tmp_help = (TH1F*) gDirectory->FindObjectAny(histoName);
   tmp_data = tmp_help;
      
   jetht_file->cd(); 
   tmp_help = 0;
   tmp_help = (TH1F*) gDirectory->FindObjectAny(histoName);
   tmp_data->Add(tmp_help);

   jetmon_file->cd();  
   tmp_help = 0;
   tmp_help = (TH1F*) gDirectory->FindObjectAny(histoName);
   tmp_data->Add(tmp_help);
      
   mc_file->cd();  
   tmp_help = 0;
   tmp_mc= (TH1F*) gDirectory->FindObjectAny(histoName);

   TCanvas *c = new TCanvas();
   c->SetLogy();
   c->SetBottomMargin(0.25 + 0.75*c->GetBottomMargin()-0.25*c->GetTopMargin());
   c->cd();
   tmp_data->GetXaxis()->SetLabelSize(0);
   tmp_data->GetXaxis()->SetTitle("");
   tmp_data->GetYaxis()->SetTitle("Events");
   tmp_data->SetMarkerStyle(26);
   tmp_data->SetMarkerColor(kRed);
   tmp_data->Draw();
   tmp_mc->Scale(tmp_data->Integral()/tmp_mc->Integral());
   tmp_mc->SetLineColor(kBlue+2);
   tmp_mc->Draw("samehist");
      
   TPaveText *label = util::LabelFactory::createPaveTextWithOffset(1,0.925,0.3);
   label->AddText(Label);
   label->Draw("same");
      
   TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(2,0.925,0.38);
   leg1->AddEntry(tmp_data,"Data","P");
   leg1->AddEntry(tmp_mc,"MC","LF");
   leg1->Draw();
      
   tmp_data->Draw("same");
   tmp_mc->Draw("samehist");

   TPad *pad = new TPad("pad", "pad", 0, 0, 1, 1);
   pad->SetTopMargin(0.75 - 0.75*pad->GetBottomMargin()+0.25*pad->GetTopMargin());
   pad->SetFillStyle(0);
   pad->SetFrameFillColor(10);
   pad->SetFrameBorderMode(0);
   pad->Draw();
   pad->cd();

   Double_t xMin1 = tmp_data->GetXaxis()->GetXmin();
   Double_t xMax1 = tmp_data->GetXaxis()->GetXmax();

   TH1F* r = new TH1F(*tmp_data);
   r->Sumw2();
   r->SetXTitle(xTitle);
   r->GetYaxis()->CenterTitle();
   r->SetYTitle("(Data-MC)/Data");
   r->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X"));
   r->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y"));
   r->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/0.3);
   r->GetYaxis()->SetNdivisions(505);
   r->SetStats(0);
   r->SetMarkerStyle(20);
   //      r->SetMarkerSize(1.12);
   r->SetMarkerColor(kBlack);
   r->Reset();
   r->Add(tmp_data, 1);
   r->Add(tmp_mc, -1);
   r->Divide(tmp_data);
   r->SetMaximum(0.5);
   r->SetMinimum(-0.5);
   r->Draw("ep");
   TLine l;
   l.DrawLine(xMin1, 0., xMax1, 0.);
   c->cd();
      
   c->Print("ControlPlots/" + histoName + ".eps");
}

void ControlPlots()
{
   setTDRStyle();
   gROOT->ForceStyle();
   
   // input files
   TFile* jet_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/Jet_ReRecoA.root", "READ");
   TFile* jetht_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD.root", "READ");
   TFile* jetmon_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD.root", "READ");
   //TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat.root", "READ");
   //  TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_SmearedWithMeasuredValues.root", "READ");
   TFile* mc_file = new TFile("/nfs/dust/test/cms/user/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_ReweightAlphaSpectrum.root", "READ");
   
   MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "N_{Vtx}", "Before PU Reweighting", "NVtx_AfterTriggerSelection");
   MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "N_{Vtx}", "After PU Reweighting", "NVtx_AfterPUReweighting");
   MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "N_{Vtx}", "After Final Selection", "NVtx_AfterAsymmHistos");
 
   for(int ipt = 0; ipt < 13; ++ipt){
      for(int ieta = 0; ieta < 5; ++ieta){
         TString ss;
         ss.Form("AlphaSpectrum_Pt%i_eta%i", ipt, ieta);
         MakePlot_Alpha(jet_file, jetht_file, jetmon_file, mc_file, "#alpha", "After Final Selection", ss, ipt, ieta);
      }
   }
        
   //////////////////////////////////////////////////////////////

   TH1F *ptave_data = new TH1F();
   TH1F *ptave_mc = new TH1F();
   TH1F *ptave_mc_new = new TH1F();
   TH1F *tmp = new TH1F();

   ptave_data->Sumw2();
   ptave_mc->Sumw2();
   ptave_mc_new->Sumw2();
   tmp->Sumw2();

   jet_file->cd();  
   ptave_data = 0;
   tmp = 0;
   tmp = (TH1F*) gDirectory->FindObjectAny("PtAve_AfterAsymmHistos");
   ptave_data = tmp;

   jetht_file->cd();  
   tmp = 0;
   tmp = (TH1F*) gDirectory->FindObjectAny("PtAve_AfterAsymmHistos");
   ptave_data->Add(tmp);

   jetmon_file->cd();  
   tmp = 0;
   tmp = (TH1F*) gDirectory->FindObjectAny("PtAve_AfterAsymmHistos");
   ptave_data->Add(tmp);

   mc_file->cd();  
   ptave_mc = 0;
   ptave_mc = (TH1F*) gDirectory->FindObjectAny("PtAve_AfterAsymmHistos");

   TCanvas *c3 = new TCanvas();
   c3->SetLogx();
   c3->SetLogy();

   float int_ptbin1_data = ptave_data->Integral(ptave_data->FindBin(62), ptave_data->FindBin(107));
   float int_ptbin2_data = ptave_data->Integral(ptave_data->FindBin(107), ptave_data->FindBin(175));
   float int_ptbin3_data = ptave_data->Integral(ptave_data->FindBin(175), ptave_data->FindBin(242));
   float int_ptbin4_data = ptave_data->Integral(ptave_data->FindBin(242), ptave_data->FindBin(310));
   float int_ptbin5_data = ptave_data->Integral(ptave_data->FindBin(310), ptave_data->FindBin(379));
   float int_ptbin6_data = ptave_data->Integral(ptave_data->FindBin(379), ptave_data->FindBin(467));
   float int_ptbin7_data = ptave_data->Integral(ptave_data->FindBin(467), ptave_data->GetNbinsX());

   float int_ptbin1_mc = ptave_mc->Integral(ptave_mc->FindBin(62), ptave_mc->FindBin(107));
   float int_ptbin2_mc = ptave_mc->Integral(ptave_mc->FindBin(107), ptave_mc->FindBin(175));
   float int_ptbin3_mc = ptave_mc->Integral(ptave_mc->FindBin(175), ptave_mc->FindBin(242));
   float int_ptbin4_mc = ptave_mc->Integral(ptave_mc->FindBin(242), ptave_mc->FindBin(310));
   float int_ptbin5_mc = ptave_mc->Integral(ptave_mc->FindBin(310), ptave_mc->FindBin(379));
   float int_ptbin6_mc = ptave_mc->Integral(ptave_mc->FindBin(379), ptave_mc->FindBin(467));
   float int_ptbin7_mc = ptave_mc->Integral(ptave_mc->FindBin(467), ptave_mc->GetNbinsX());
  
   ptave_data->GetXaxis()->SetRangeUser(62, 2000);
   ptave_data->GetXaxis()->SetTitle("#bar{ p}_{T} [GeV]");
   ptave_data->Draw();

   for(int i = 1; i <= ptave_mc->GetNbinsX(); i++) {
      if(i >= ptave_mc->FindBin(62) && i < ptave_mc->FindBin(107) ) {
         ptave_mc->SetBinContent(i, ptave_mc->GetBinContent(i)*(int_ptbin1_data/int_ptbin1_mc));
      }
      if(i >= ptave_mc->FindBin(107) && i < ptave_mc->FindBin(175) ) {
         ptave_mc->SetBinContent(i, ptave_mc->GetBinContent(i)*int_ptbin2_data/int_ptbin2_mc);
      }
      if(i >= ptave_mc->FindBin(175) && i < ptave_mc->FindBin(242) ) {
         ptave_mc->SetBinContent(i, ptave_mc->GetBinContent(i)*int_ptbin3_data/int_ptbin3_mc);
      }
      if(i >= ptave_mc->FindBin(242) && i < ptave_mc->FindBin(310) ) {
         ptave_mc->SetBinContent(i, ptave_mc->GetBinContent(i)*int_ptbin4_data/int_ptbin4_mc);
      }
      if(i >= ptave_mc->FindBin(310) && i < ptave_mc->FindBin(379) ) {
         ptave_mc->SetBinContent(i, ptave_mc->GetBinContent(i)*int_ptbin5_data/int_ptbin5_mc);
      }
      if(i >= ptave_mc->FindBin(379) && i < ptave_mc->FindBin(467) ) {
         ptave_mc->SetBinContent(i, ptave_mc->GetBinContent(i)*int_ptbin6_data/int_ptbin6_mc);
      }
      if(i >= ptave_mc->FindBin(467) ) {
         ptave_mc->SetBinContent(i, ptave_mc->GetBinContent(i)*(int_ptbin7_data/int_ptbin7_mc));
      }
   }

   //  ptave_mc->Scale(int_ptbin1_data/ptave_mc_new->Integral());
   ptave_mc->SetLineColor(kRed);
   ptave_mc->Draw("samehist");

   TPaveText *label3 = util::LabelFactory::createPaveTextWithOffset(1,0.98,0.5);
   label3->AddText("After Final Selection");
   label3->Draw("same");

   //  leg1->Draw();

   ptave_data->Draw("same");
   ptave_mc->Draw("samehist");

   c3->Print("ControlPlots/PtAveSpectrum.eps");

}
