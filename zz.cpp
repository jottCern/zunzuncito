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
    p.add_module(shared_ptr<module>(new eventfilter()));
    //p.add_module(shared_ptr<module>(new dumper("before JEC", 5)));
    p.add_module(shared_ptr<module>(new jec("/afs/desy.de/user/o/ottjoc/zunzuncito/jetcorrs", 5)));
    //p.add_module(shared_ptr<module>(new dumper("after JEC", 5)));
    //p.add_module(shared_ptr<module>(new jec_uncertainty("/afs/desy.de/user/o/ottjoc/zunzuncito/jetcorrs/FT_53_V21_AN6_Uncertainty_AK5PFchs.txt", syst_dir::up)));
    //p.add_module(shared_ptr<module>(new dumper("after JEC uncertainty UP", 5)));
    p.add_module(shared_ptr<module>(new asymm_histos()));
    
    // 2. define datasets
    vector<dataset> datasets;
    
    string prefix = "/nfs/dust/test/cms/user/rathjd/Calibration/Summer2013ReReco_v1/Summer2013ReReco_v1/ak5PFCHS_";
    
    dataset runa("runa");
    runa.mc = false;
    runa.files = glob(prefix + "*_sam0.root");
    runa.nmax = 0;
    datasets.push_back(runa);
    
    dataset mc("mc");
    mc.mc = true;
    mc.files = glob("/nfs/dust/test/cms/user/rathjd/Calibration/QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1/merged/ak5PFCHS_0.root");
    mc.nmax = 100000;
    datasets.push_back(mc);
    
    // 3. run all modules on all datasets:
    p.run(datasets, "/nfs/dust/test/cms/user/ottjoc/zz-out/");
}
