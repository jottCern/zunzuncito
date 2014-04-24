#include "modules.hpp"
#include "event.hpp"
#include "utils.hpp"
#include "processor.hpp"

#include "TFile.h"
#include "TChain.h"

#include <iostream>
#include <stdexcept>

using namespace std;

int main(){
    // 1. define sequence of modules to run on the events:
    processor p;
    // p.add_module(boost::shared_ptr<module>(new eventfilter())); // done on ntuple already
    // p.add_module(boost::shared_ptr<module>(new dumper("before JEC", 5)));
    // p.add_moduleboost::(shared_ptr<module>(new controlplots("NoCuts")));
    p.add_module(boost::shared_ptr<module>(new jec("/afs/desy.de/user/k/kheine/zunzuncito/jetcorrs", 5)));
    // p.add_module(boost::shared_ptr<module>(new dumper("after JEC", 5)));
    // p.add_module(shared_ptr<module>(new controlplots("AfterJEC")));
    // p.add_module(boost::shared_ptr<module>(new jec_uncertainty("/afs/desy.de/user/k/kheine/zunzuncito/jetcorrs/FT_53_V21_AN6_Uncertainty_AK5PFchs.txt", syst_dir::down)));
    // p.add_module(boost::shared_ptr<module>(new dumper("after JEC uncertainty UP", 5)));
    // p.add_module(boost::shared_ptr<module>(new controlplots("AfterJEC")));
    p.add_module(boost::shared_ptr<module>(new reasonablemceventchecker()));
    p.add_module(boost::shared_ptr<module>(new response_histos()));
    // p.add_module(boost::shared_ptr<module>(new controlplots("AfterMCEventChecker")));
    // p.add_module(boost::shared_ptr<module>(new smearmc(5)));
    // p.add_module(boost::shared_ptr<module>(new controlplots("AfterSmearMC")));
    p.add_module(boost::shared_ptr<module>(new checktrigger()));
    p.add_module(boost::shared_ptr<module>(new controlplots("AfterTriggerSelection")));
    p.add_module(boost::shared_ptr<module>(new pureweighting("/afs/desy.de/user/k/kheine/zunzuncito/PUDistributions")));
    // p.add_module(boost::shared_ptr<module>(new dumper("after pureweighting", 5)));
    p.add_module(boost::shared_ptr<module>(new controlplots("AfterPUReweighting")));
    p.add_module(boost::shared_ptr<module>(new eventcuts()));
    p.add_module(boost::shared_ptr<module>(new controlplots("AfterEventCuts")));
    // p.add_module(boost::shared_ptr<module>(new alphareweighting()));
    // p.add_module(boost::shared_ptr<module>(new gluonsplittingreweighting()));
    // p.add_module(boost::shared_ptr<module>(new controlplots("AfterAlphaReweighting")));
    p.add_module(boost::shared_ptr<module>(new asymm_histos()));
    p.add_module(boost::shared_ptr<module>(new controlplots("AfterAsymmHistos")));
      
    // 2. define datasets
    vector<dataset> datasets;
    
    string prefix = "/nfs/dust/cms/user/rathjd/Calibration/Summer2013ReReco_v1/ak5PFCHS_";
  
    // ---- Data ---- //
    // -------------------------------- //
    // dataset jet("Jet_ReRecoA_NoMinPtCutForThirdJet_AddNewAlphaBin_final");
    // dataset jet("Jet_ReRecoA_nominal_exclusive_alpha");
    // dataset jet("Jet_ReRecoA_NoMinPtCutForThirdJet_AddNewAlphaBin_nominal_v3");
    dataset jet("Jet_ReRecoA_nominal_v3");
    // dataset jet("Jet_ReRecoA_ForwardExtension_final_v2");
    // dataset jet("Jet_ReRecoA_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v2");
    // dataset jet("Jet_ReRecoA_ForwardExtensionNextToCentral_final_v1");
    // dataset jet("Jet_ReRecoA_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v1");
    jet.mc = false;
    jet.jetdata = true;
    jet.jethtdata = false;
    jet.jetmondata = false;
    jet.files = glob(prefix + "*_sam0.root");
    //jet.nskip = 2900000;
    //jet.nmax = 100000;
    datasets.push_back(jet);
    // -------------------------------- //
    
    // -------------------------------- //
    // dataset jetht("JetHT_ReRecoBToD_NoMinPtCutForThirdJet_AddNewAlphaBin_final");
    // dataset jetht("JetHT_ReRecoBToD_nominal_exclusive_alpha");
    // dataset jetht("JetHT_ReRecoBToD_NoMinPtCutForThirdJet_AddNewAlphaBin_nominal_v3");
    dataset jetht("JetHT_ReRecoBToD_nominal_v3");
    // dataset jetht("JetHT_ReRecoA_ForwardExtension_final_v2");
    // dataset jetht("JetHT_ReRecoBToD_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v2");
    // dataset jetht("JetHT_ReRecoBToD_ForwardExtensionNextToCentral_final_v1");
    // dataset jetht("JetHT_ReRecoBToD_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v1");
    jetht.mc = false;
    jetht.jetdata = false;
    jetht.jethtdata = true;
    jetht.jetmondata = false;
    jetht.files = glob(prefix + "*_sam[1 3 5].root");
    //jetht.nmax = 100000;
    datasets.push_back(jetht);
    // -------------------------------- //

    // -------------------------------- //
    // dataset jetmon("JetMon_ReRecoBToD_NoMinPtCutForThirdJet_AddNewAlphaBin_final");
    // dataset jetmon("JetMon_ReRecoBToD_nominal_exclusive_alpha");
    // dataset jetmon("JetMon_ReRecoBToD_NoMinPtCutForThirdJet_AddNewAlphaBin_nominal_v3");
    dataset jetmon("JetMon_ReRecoBToD_nominal_v3");
    // dataset jetmon("JetMon_ReRecoA_ForwardExtension_final_v2");
    // dataset jetmon("JetMon_ReRecoBToD_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v2");
    // dataset jetmon("JetMon_ReRecoBToD_ForwardExtensionNextToCentral_final_v1");
    // dataset jetmon("JetMon_ReRecoBToD_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_final_v1");
    jetmon.mc = false;
    jetmon.jetdata = false;
    jetmon.jethtdata = false;
    jetmon.jetmondata = true;
    jetmon.files = glob(prefix + "*_sam[2 4 6].root");
    //jetmon.nmax = 100000;
    // jetmon.nmax = 0;
    datasets.push_back(jetmon); 
    // -------------------------------- //
  
    // ---- MC ---- //
    // -------------------------------- //
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_v2");  
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_MBXS73500_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_JECup_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_JECdown_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_ReweightAlphaSpectrum_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_ClosureFirstHalf_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_ClosureSecondHalf_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_SmearedWithMeasuredValues_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_NoMinPtCutForThirdJet_AddNewAlphaBin_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_GluonSplittingReweighting_v2");

    dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_v3");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_NoMinPtCutForThirdJet_AddNewAlphaBin_v3");  

    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtension_NoMinPtCutForThirdJet_AddNewAlphaBin_v2");

    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionNextToCentral_v1");
    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_final_nominal_ForwardExtensionNextToCentral_NoMinPtCutForThirdJet_AddNewAlphaBin_v1");

    // dataset mc("MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_ForwardExtension_v2");
    // dataset mc("MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_ForwardExtensionNextToCentral_v1");
    // dataset mc("MC_QCD_Pt-15to3000_TuneEE3C_Flat_herwigpp_final_nominal_SmearedWithMeasuredValues_v2");

    // dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_nominal_exclusive_alpha");
    mc.mc = true;
    mc.jetdata = false;
    mc.jethtdata = false;
    mc.jetmondata = false;
    mc.files = glob("/nfs/dust/test/cms/user/rathjd/Calibration/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/merged/ak5PFCHS_*.root");
    // mc.files = glob("/nfs/dust/cms/user/kheine/CalibNTupel/MC/EE3C_herwigpp/Summer12_DR53X_QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp_ak5PFCHS_v2.root");
    // mc.nskip = 4909676;
    //mc.nmax = 4909675;
    // mc.nmax = 100;
    datasets.push_back(mc);
    // -------------------------------- //
    
    // 3. run all modules on all datasets:
    p.run(datasets, "/afs/desy.de/user/k/kheine/zunzuncito/zz-out/");
}
