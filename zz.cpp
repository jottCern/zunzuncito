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
    // p.add_module(shared_ptr<module>(new eventfilter())); // done on ntuple already
    //p.add_module(shared_ptr<module>(new dumper("before JEC", 5)));
    //   p.add_module(shared_ptr<module>(new controlplots("NoCuts")));
    p.add_module(shared_ptr<module>(new jec("/afs/desy.de/user/k/kheine/kheine_dust/zunzuncito/jetcorrs", 5)));
    //p.add_module(shared_ptr<module>(new dumper("after JEC", 5)));
    //  p.add_module(shared_ptr<module>(new controlplots("AfterJEC")));
    //p.add_module(shared_ptr<module>(new jec_uncertainty("/afs/desy.de/user/k/kheine/kheine_dust/zunzuncito/jetcorrs/FT_53_V21_AN6_Uncertainty_AK5PFchs.txt", syst_dir::down)));
    //p.add_module(shared_ptr<module>(new dumper("after JEC uncertainty UP", 5)));
    p.add_module(shared_ptr<module>(new controlplots("AfterJEC")));
    p.add_module(shared_ptr<module>(new reasonablemceventchecker()));
    p.add_module(shared_ptr<module>(new controlplots("AfterMCEventChecker")));
    // p.add_module(shared_ptr<module>(new smearmc(5)));
    //p.add_module(shared_ptr<module>(new controlplots("AfterSmearMC")));
    p.add_module(shared_ptr<module>(new checktrigger()));
    p.add_module(shared_ptr<module>(new controlplots("AfterTriggerSelection")));
    p.add_module(shared_ptr<module>(new pureweighting("/afs/desy.de/user/k/kheine/kheine_dust/PUDistributions")));
    //p.add_module(shared_ptr<module>(new dumper("after pureweighting", 5)));
    p.add_module(shared_ptr<module>(new controlplots("AfterPUReweighting")));
    p.add_module(shared_ptr<module>(new eventcuts()));
    p.add_module(shared_ptr<module>(new controlplots("AfterEventCuts")));
    // p.add_module(shared_ptr<module>(new alphareweighting()));
    //p.add_module(shared_ptr<module>(new controlplots("AfterAlphaReweighting")));
    p.add_module(shared_ptr<module>(new asymm_histos()));
    p.add_module(shared_ptr<module>(new controlplots("AfterAsymmHistos")));
    p.add_module(shared_ptr<module>(new response_histos()));
    p.add_module(shared_ptr<module>(new controlplots("AfterResponseHistos")));
    
    // 2. define datasets
    vector<dataset> datasets;
    
    string prefix = "/nfs/dust/test/cms/user/rathjd/Calibration/Summer2013ReReco_v1/Summer2013ReReco_v1/ak5PFCHS_";
    // string prefix = "/nfs/dust/test/cms/user/rathjd/Calibration/Summer2013ReReco_v1/Summer2013ReReco_v1/ak5Calo_";

  //   dataset jet("Jet_ReRecoA_v2");
//     jet.mc = false;
//     jet.jetdata = true;
//     jet.jethtdata = false;
//     jet.jetmondata = false;
//     jet.files = glob(prefix + "*_sam0.root");
//     //jet.nskip = 2900000;
//     //jet.nmax = 100000;
//     datasets.push_back(jet);
    
//     dataset jetht("JetHT_ReRecoBToD_v2");
//     jetht.mc = false;
//     jetht.jetdata = false;
//     jetht.jethtdata = true;
//     jetht.jetmondata = false;
//     jetht.files = glob(prefix + "*_sam[1 3 5].root");
//     //jetht.nmax = 100000;
//     datasets.push_back(jetht);

//     dataset jetmon("JetMon_ReRecoBToD_v2");
//     jetmon.mc = false;
//     jetmon.jetdata = false;
//     jetmon.jethtdata = false;
//     jetmon.jetmondata = true;
//     jetmon.files = glob(prefix + "*_sam[2 4 6].root");
//     //jetmon.nmax = 100000;
//     // jetmon.nmax = 0;
//     datasets.push_back(jetmon);
    
    //dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_SmearedWithMeasuredValues");
    //dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_JECdown");
    dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_ChangedPLIDefinitionToGenPtAve");
    //dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_MBXS73500");
    //dataset mc("MC_QCD_Pt-15to3000_TuneZ2_Flat_ReweightAlphaSpectrum");
    mc.mc = true;
    mc.jetdata = false;
    mc.jethtdata = false;
    mc.jetmondata = false;
    mc.files = glob("/nfs/dust/test/cms/user/rathjd/Calibration/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/merged/ak5PFCHS_*.root");
    // mc.files = glob("/nfs/dust/test/cms/user/rathjd/Calibration/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/merged/ak5Calo_*.root");
    //  mc.nmax = 100000;
    //mc.nmax = 50;
    datasets.push_back(mc);
    
    // 3. run all modules on all datasets:
    p.run(datasets, "/afs/desy.de/user/k/kheine/kheine_dust/zunzuncito/zz-out/");
}
