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

void MakePlot_Alpha(TFile *jet_file, TFile *jetht_file, TFile *jetmon_file, TFile *mc_file, TString xTitle, TString Label, TString histoName, int pt_bin, int eta_bin)
{
   float pt_bins[14] = {62, 107, 175, 205, 242, 270, 310, 335, 379, 410, 467, 600, 1000, 2000};
   float eta_bins[6] = {0, 0.5, 1.1, 1.7, 2.3, 5.0};

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
      //f->Draw("same");
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
      
   c->Print("ControlPlots/Nominal_" + histoName + ".eps");
}

void MakePlot_Trig(TFile *jet_file, TFile *jetht_file, TFile *jetmon_file, TFile *mc_file, TString xTitle, TString Label, TString histoName, TString suffix)
{
   TString trigger[7] = {"HltDiPFJetAve40", "HltDiPFJetAve80", "HltDiPFJetAve140", "HltDiPFJetAve200", "HltDiPFJetAve260", "HltDiPFJetAve320", "HltDiPFJetAve400"};

   TString pt_bins[8] = {"62", "107", "175", "242", "310", "379", "467", "#infty"};

   for(int i=0; i < 7; i++) {

      TString histname = histoName;
      histname += trigger[i];
      histname += suffix;

      TH1F *tmp_data = new TH1F();
      TH1F *tmp_mc = new TH1F();
      TH1F *tmp_help = new TH1F();
  
      jet_file->cd();  
      tmp_help = 0;
      tmp_help = (TH1F*) gDirectory->FindObjectAny(histname);
      tmp_data = tmp_help;
      
      jetht_file->cd(); 
      tmp_help = 0;
      tmp_help = (TH1F*) gDirectory->FindObjectAny(histname);
      tmp_data->Add(tmp_help);

      jetmon_file->cd();  
      tmp_help = 0;
      tmp_help = (TH1F*) gDirectory->FindObjectAny(histname);
      tmp_data->Add(tmp_help);
      
      mc_file->cd();  
      tmp_help = 0;
      tmp_mc= (TH1F*) gDirectory->FindObjectAny(histname);

      TCanvas *c = new TCanvas();
      c->SetLogy();
      // c->SetBottomMargin(0.25 + 0.75*c->GetBottomMargin()-0.25*c->GetTopMargin());
      c->cd();
      //  tmp_data->Rebin(100);
      tmp_data->GetXaxis()->SetTitle(xTitle);
      //tmp_data->GetXaxis()->SetLabelSize(0);
      tmp_data->GetYaxis()->SetTitle("Events");
      tmp_data->GetYaxis()->SetRangeUser(0.5, 10000* tmp_data->GetMaximum() );
      tmp_data->SetMarkerStyle(20);
      tmp_data->SetMarkerColor(kBlack);
      tmp_data->Draw();
      // tmp_mc->Rebin(100);
      tmp_mc->Scale(tmp_data->Integral()/tmp_mc->Integral());
      tmp_mc->SetLineColor(38);
      tmp_mc->SetFillColor(38);
      tmp_mc->Draw("samehist");
      tmp_data->Draw("same");
      gPad->RedrawAxis();
      
      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(1, 1.0, 0.01);
      if(i == 6) label->AddText(trigger[i] + ", "  + pt_bins[i] + " #leq p_{T}^{ave} [GeV] ");
      else label->AddText(trigger[i] + ", "  + pt_bins[i] + " #leq p_{T}^{ave} [GeV] #leq " + pt_bins[i+1]);
      label->Draw("same");
      
      TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(2,0.4,0.11);
      leg1->AddEntry(tmp_data,"Data","PL");
      leg1->AddEntry(tmp_mc,"Simulation","LF");
      leg1->Draw();
      
      //  tmp_data->Draw("same");
      //tmp_mc->Draw("samehist");

      /* TPad *pad = new TPad("pad", "pad", 0, 0, 1, 1);
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
      c->cd(); */// 
      
      c->Print("ControlPlots/" + histname + ".eps");
      // c->Print("ControlPlots/ForwardExtension" + histname + ".eps");
      //c->Print("ControlPlots/Herwig" + histname + ".eps");

   }
}

void MakePlot(TFile *jet_file, TFile *jetht_file, TFile *jetmon_file, TFile *mc_file, TString xTitle, TString Label, TString histoName, TString suffix)
{
   TString trigger[7] = {"HltDiPFJetAve40", "HltDiPFJetAve80", "HltDiPFJetAve140", "HltDiPFJetAve200", "HltDiPFJetAve260", "HltDiPFJetAve320", "HltDiPFJetAve400"};

   TString dummy_name = histoName;
   dummy_name += trigger[0];
   dummy_name += suffix;
   jet_file->cd();  
   TH1F *dummy_data = new TH1F();
   dummy_data = (TH1F*) gDirectory->FindObjectAny(dummy_name);
   mc_file->cd(); 
   TH1F *dummy_mc = new TH1F();
   dummy_mc = (TH1F*) gDirectory->FindObjectAny(dummy_name);
   // TH1F *dummy_gen = new TH1F();
   //dummy_gen = (TH1F*) gDirectory->FindObjectAny("Gen" + dummy_name);

   TH1F *tmp_data = new TH1F(*dummy_data);
   tmp_data->Clear();
   TH1F *tmp_mc = new TH1F(*dummy_mc);
   tmp_mc->Clear();
   // TH1F *tmp_gen = new TH1F(*dummy_gen);
   // tmp_gen->Clear();
 
   TString histname;

   for(int i=0; i < 7; i++) {

      histname.Clear();
      histname = histoName;
      histname += trigger[i];
      histname += suffix;

      TH1F *tmp_help_data1 = new TH1F();
      TH1F *tmp_help_data2 = new TH1F();
      TH1F *tmp_help_mc = new TH1F();
      //  TH1F *tmp_help_gen = new TH1F();

      if(i < 5){
         jet_file->cd();  
         tmp_help_data1 = 0;
         tmp_help_data1 = (TH1F*) gDirectory->FindObjectAny(histname);
              
         jetmon_file->cd(); 
         tmp_help_data2 = 0;
         tmp_help_data2 = (TH1F*) gDirectory->FindObjectAny(histname);
      }
      
      if(i >=5){
         jet_file->cd();  
         tmp_help_data1 = 0;
         tmp_help_data1 = (TH1F*) gDirectory->FindObjectAny(histname);

         jetht_file->cd(); 
         tmp_help_data2 = 0;
         tmp_help_data2 = (TH1F*) gDirectory->FindObjectAny(histname);
      }

      mc_file->cd(); 
      tmp_help_mc = (TH1F*) gDirectory->FindObjectAny(histname);

      //  mc_file->cd(); 
      //  tmp_help_gen = (TH1F*) gDirectory->FindObjectAny("Gen" + histname);

      // cout << "Data: " << i << "   " << tmp_help_data1->GetEntries()+tmp_help_data2->GetEntries() << endl;
      //   cout << "MC: " << i << "   " << tmp_mc->GetEntries() << endl;

      TH1F *tmp_help_data = new TH1F(*tmp_help_data1);
      tmp_help_data->Clear();
      tmp_help_data->Add(tmp_help_data1);
      tmp_help_data->Add(tmp_help_data2); 
    
      tmp_help_mc->Scale(tmp_help_data->Integral()/tmp_help_mc->Integral());

      //  tmp_help_gen->Scale(tmp_help_data->Integral()/tmp_help_gen->Integral());

      tmp_mc->Add(tmp_help_mc);
      tmp_data->Add(tmp_help_data);
      //   tmp_gen->Add(tmp_help_gen);
   }

   if( histname.Contains("Alpha")){
      //   if(  histname.Contains("Bla") ) {
      TCanvas *c = new TCanvas();
      //c->SetLogx();
      c->SetLogy();
      c->SetBottomMargin(0.25 + 0.75*c->GetBottomMargin()-0.25*c->GetTopMargin());
      c->cd();
      //   tmp_data->Rebin(2);
      tmp_data->GetXaxis()->SetLabelSize(0);
      tmp_data->GetXaxis()->SetTitle("");
      tmp_data->GetYaxis()->SetTitle("Events");
      tmp_data->GetXaxis()->SetRangeUser(0.0, 0.30 );
      tmp_data->GetYaxis()->SetRangeUser(0.5, 1000000* tmp_data->GetMaximum() );
      tmp_data->SetMarkerStyle(20);
      tmp_data->SetMarkerColor(kBlack);
      tmp_data->GetXaxis()->SetNdivisions(505);
      tmp_data->GetYaxis()->SetNdivisions(505);
      tmp_data->Draw();
      //  tmp_mc->Rebin(2);
      tmp_mc->SetLineColor(38);
      tmp_mc->SetFillColor(38);
      tmp_mc->Draw("samehist");
      tmp_data->Draw("same");
      gPad->RedrawAxis();
      
      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2, 0.6, 0.01);
      label->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
      label->AddText("0.0 #leq |#eta| #leq 5.0, #Delta#phi > 2.7");
      label->Draw("same");
      
      TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(2,0.6,0.11);
      leg1->AddEntry(tmp_data,"Data","PL");
      leg1->AddEntry(tmp_mc,"Simulation","LF");
      leg1->Draw();

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
      r->SetYTitle("Data/Sim.");
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
         f->SetLineColor( kRed );
         f->SetLineWidth( 1 );
         f->Draw("same");
      }     

      TF1 *f2 = new TF1("func", "0.5*[0]*(TMath::Erf([1]*x-[2])+1)", 0, 0.3);
      f2->SetParameter(0, 1.15);
      f2->SetParameter(1, 15);
      f2->SetParameter(2, 0.1);

      f2->SetLineColor(kRed);
      //  f2->Draw("same");

      cout << "Par 0: " << fit->GetParameter(0) << endl; 
      cout << "Par 1: " << fit->GetParameter(1) << endl; 
      cout << "Par 2: " << fit->GetParameter(2) << endl; 

      TLine l;
      l.DrawLine(xMin1, 1., 0.3, 1.);
      c->cd();
      
      TString PlotsName = histoName;
      PlotsName += suffix;
      c->Print("ControlPlots/" + PlotsName + ".eps");
      // c->Print("ControlPlots/Herwig" + PlotsName + ".eps");
      //c->Print("ControlPlots/AfterReweight_"+ PlotsName + ".eps");
      //c->Print("ControlPlots/ForwardExtension" + PlotsName + ".eps");

   }

   else {
      TCanvas *c = new TCanvas();
      c->SetLogx();
      c->SetLogy();
      // c->SetBottomMargin(0.25 + 0.75*c->GetBottomMargin()-0.25*c->GetTopMargin());
      //c->cd();
      tmp_data->Rebin(4);
      tmp_data->GetXaxis()->SetRange(tmp_data->GetXaxis()->GetXmin(), tmp_data->GetXaxis()->GetXmax());
      if(histname.Contains("PtAve")) tmp_data->GetXaxis()->SetRangeUser(62., tmp_data->GetXaxis()->GetXmax());
      tmp_data->GetXaxis()->SetTitle(xTitle);
      tmp_data->GetYaxis()->SetTitle("Events");
      tmp_data->GetYaxis()->SetRangeUser(1.0, 1000* tmp_data->GetMaximum() );
      tmp_data->SetMarkerStyle(20);
      tmp_data->SetMarkerColor(kBlack);
      //  tmp_data->GetXaxis()->SetNdivisions(505);
      tmp_data->Draw();
      tmp_mc->Rebin(4);
      tmp_mc->SetLineColor(38);
      tmp_mc->SetFillColor(38);
      tmp_mc->Draw("samehist");
      tmp_data->Draw("same");
      gPad->RedrawAxis();
      
      TPaveText *label = util::LabelFactory::createPaveTextWithOffset(2, 0.6, 0.01);
      label->AddText("Anti-k_{T} (R=0.5) PFCHS Jets");
      label->AddText("0.0 #leq |#eta| #leq 5.0, #Delta#phi > 2.7");
      label->Draw("same");
      
      TLegend* leg1 = util::LabelFactory::createLegendColWithOffset(2,0.6,0.11);
      leg1->AddEntry(tmp_data,"Data","PL");
      leg1->AddEntry(tmp_mc,"Simulation","LF");
      leg1->Draw();
      
      TString PlotsName = histoName;
      PlotsName += suffix;
      //c->Print("ControlPlots/" + PlotsName + ".eps");
      c->Print("ControlPlots/Herwig" + PlotsName + ".eps");
      // c->Print("ControlPlots/ForwardExtension" + PlotsName + ".eps");
   }

   if(histname.Contains("PtAve")) {
      std::vector<int> pt_bins;
      pt_bins.push_back(62);
      pt_bins.push_back(107);
      pt_bins.push_back(175);
      pt_bins.push_back(205);
      pt_bins.push_back(242);
      pt_bins.push_back(270);
      pt_bins.push_back(310);
      pt_bins.push_back(335);
      pt_bins.push_back(379);
      pt_bins.push_back(410);
      pt_bins.push_back(467);
      pt_bins.push_back(600);
      pt_bins.push_back(1000);
      pt_bins.push_back(2000);
      for(int i = 0; i < 13; i++) {
         tmp_data->GetXaxis()->SetRangeUser(pt_bins.at(i), pt_bins.at(i+1));
         std::cout << "data pt: " << pt_bins.at(i) << " - " << pt_bins.at(i+1) << " Mean =  " << tmp_data->GetMean() << std::endl;
         //std::cout << "data events: " << pt_bins.at(i) << " - " << pt_bins.at(i+1) << " # ev. =  " << tmp_data->GetEntries() << std::endl;
         tmp_mc->GetXaxis()->SetRangeUser(pt_bins.at(i), pt_bins.at(i+1));
         std::cout << "mc pt: " << pt_bins.at(i) << " - " << pt_bins.at(i+1) << " Mean =  " << tmp_mc->GetMean() << std::endl;
         //   std::cout << "mc events: " << pt_bins.at(i) << " - " << pt_bins.at(i+1) << " # ev. =  " << tmp_mc->GetEntries() << std::endl;
         //    tmp_gen->GetXaxis()->SetRangeUser(pt_bins.at(i), pt_bins.at(i+1));
         //  std::cout << "gen pt: " << pt_bins.at(i) << " - " << pt_bins.at(i+1) << " Mean =  " << tmp_gen->GetMean() << std::endl;
      }
   }
}
  

void ControlPlots()
{
   setTDRStyle();
   gROOT->ForceStyle();
   
   // input files
   /*  TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_final.root", "READ");
   TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD_final.root", "READ");
   TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD_final.root", "READ");
   TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final.root", "READ");*/


   // TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_AddAngularHistos_final.root", "READ");
   // TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoA_AddAngularHistos_final.root", "READ");
   // TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoA_AddAngularHistos_final.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_AddAngularHistos_v3.root", "READ");
   //  TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_AddAngularHistos_v3.root", "READ");
   

   /*TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_ForwardExtension_final_v2.root", "READ");
   TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoA_ForwardExtension_final_v2.root", "READ");
   TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoA_ForwardExtension_final_v2.root", "READ");
   TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_v2.root", "READ");*/

   TFile* jet_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/Jet_ReRecoA_nominal_v4.root", "READ");
   TFile* jetht_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetHT_ReRecoBToD_nominal_v4.root", "READ");
   TFile* jetmon_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/JetMon_ReRecoBToD_nominal_v4.root", "READ");
   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_v4.root", "READ");
   TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_v4.root", "READ");


   // TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_v2.root", "READ");

   //TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_SmearedWithMeasuredValues.root", "READ");
   //TFile* mc_file = new TFile("/afs/desy.de/user/k/kheine/zunzuncito/zz-out/MC_QCD_Pt-15to3000_TuneZ2_Flat_final_ReweightAlphaSpectrum.root", "READ");
   
   //MakePlot_Trig(jet_file, jetht_file, jetmon_file, mc_file, "N_{Vtx}", "Before PU Reweighting", "NVtx_", "_AfterTriggerSelection");
   //MakePlot_Trig(jet_file, jetht_file, jetmon_file, mc_file, "N_{Vtx}", "After PU Reweighting", "NVtx_", "_AfterPUReweighting");
   //MakePlot_Trig(jet_file, jetht_file, jetmon_file, mc_file, "N_{Vtx}", "After Final Selection", "NVtx_", "_AfterAsymmHistos");

   MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "p_{T}^{ave} [GeV]", "After Final Selection", "PtAve_", "_AfterAsymmHistos");
   //MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "p_{T, 1} [GeV]", "After Final Selection", "Jet1Pt_", "_AfterAsymmHistos");
   //MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "p_{T, 2} [GeV]", "After Final Selection", "Jet2Pt_", "_AfterAsymmHistos");
   //MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "p_{T, 3} [GeV]", "After Final Selection", "Jet3Pt_", "_AfterAsymmHistos");
   //MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "#alpha", "After Final Selection", "Alpha_", "_AfterAsymmHistos");
   // MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "#Delta #phi", "After Final Selection", "DeltaPhi_", "_AfterAsymmHistos");
   /*MakePlot_Trig(jet_file, jetht_file, jetmon_file, mc_file, "#eta_{jet1}", "After Final Selection", "Jet1Eta_", "_AfterAsymmHistos");
   MakePlot_Trig(jet_file, jetht_file, jetmon_file, mc_file, "#eta_{jet2}", "After Final Selection", "Jet2Eta_", "_AfterAsymmHistos");
   MakePlot_Trig(jet_file, jetht_file, jetmon_file, mc_file, "#eta_{jet3}", "After Final Selection", "Jet3Eta_", "_AfterAsymmHistos");
   MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "#eta_{jet1}", "After Final Selection", "Jet1Eta_", "_AfterAsymmHistos");
   MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "#eta_{jet2}", "After Final Selection", "Jet2Eta_", "_AfterAsymmHistos");
   MakePlot(jet_file, jetht_file, jetmon_file, mc_file, "#eta_{jet3}", "After Final Selection", "Jet3Eta_", "_AfterAsymmHistos");*/


 
  //  for(int ipt = 0; ipt < 13; ++ipt){
//       for(int ieta = 0; ieta < 5; ++ieta){
//          TString ss;
//          ss.Form("AlphaSpectrum_Pt%i_eta%i", ipt, ieta);
//          MakePlot_Alpha(jet_file, jetht_file, jetmon_file, mc_file, "#alpha", "After Final Selection", ss, ipt, ieta);
//       }
//    }

   /*  for(int ipt = 0; ipt < 13; ++ipt){
      for(int ieta = 0; ieta < 5; ++ieta){
         TString ss;
         ss.Form("AlphaProjectionSpectrum_Pt%i_eta%i", ipt, ieta);
         MakePlot_Alpha(jet_file, jetht_file, jetmon_file, mc_file, "#alpha_{proj}", "After Final Selection", ss, ipt, ieta);

         TString ss2;
         ss2.Form("DeltaPhiDijetJ3_Pt%i_eta%i", ipt, ieta);
         MakePlot_Alpha(jet_file, jetht_file, jetmon_file, mc_file, "#Delta#phi(jet_3, dijet-axis)", "After Final Selection", ss2, ipt, ieta);

         //  TString ss3;
         //ss3.Form("Alpha_vs_AlphaProjection_Pt%i_eta%i", ipt, ieta);
         //MakePlot_Alpha(jet_file, jetht_file, jetmon_file, mc_file, "#Delta3phi(jet_3, dijet-axis)", "After Final Selection", ss3, ipt, ieta);
       }
       }*/


        
 
}
