#include "processor.hpp"
#include "event.hpp"

#include "TChain.h"
#include "TFile.h"

#include <stdexcept>
#include <iostream>
#include <iomanip>

using namespace std;
using boost::shared_ptr;

namespace {
    
template<typename T>
void connect(TTree & tree, const char * name, T & t){
    tree.SetBranchStatus(name, 1);
    int res = tree.SetBranchAddress(name, &t);
    if(res!=0){
        throw runtime_error(string("SetBranchAddress with name '") + name + "' failed");
    }
}    


}

void processor::add_module(const shared_ptr<module> & m){
    modules.push_back(m);
}

void processor::run(const dataset & d, TFile & outfile){
    if(d.nskip<0){
        throw runtime_error("nskip<0 encounteredl; this is not allowed");
    }
    
    if(verbose)cout << "Starting with dataset " << d.name << endl;
    TChain c("DiJetTree");
    for(size_t i=0; i<d.files.size(); ++i){
        c.Add(d.files[i].c_str());
    }
    size_t n = c.GetEntries();
    if(verbose)cout << "Found " << n << " entries in " << d.files.size() << " files" << endl;
    // connect branches:
    event evt;
    evt.is_mc = d.mc;
    
    c.SetBranchStatus("*", 0);
    connect(c, "RunNumber", evt.RunNumber);
    connect(c, "LumiBlockNumber", evt.LumiBlockNumber);
    connect(c, "EventNumber", evt.EventNumber);
    connect(c, "Weight", evt.Weight);
    
    connect(c, "HltDiPFJetAve40", evt.HltDiPFJetAve40);
    connect(c, "PS_HltDiPFJetAve40", evt.PS_HltDiPFJetAve40);
    connect(c, "HltDiPFJetAve80", evt.HltDiPFJetAve80);
    connect(c, "PS_HltDiPFJetAve80", evt.PS_HltDiPFJetAve80);
    connect(c, "HltDiPFJetAve140", evt.HltDiPFJetAve140);
    connect(c, "PS_HltDiPFJetAve140", evt.PS_HltDiPFJetAve140);
    connect(c, "HltDiPFJetAve260", evt.HltDiPFJetAve260);
    connect(c, "PS_HltDiPFJetAve260", evt.PS_HltDiPFJetAve260);
    connect(c, "HltDiPFJetAve320", evt.HltDiPFJetAve320);
    connect(c, "PS_HltDiPFJetAve320", evt.PS_HltDiPFJetAve320);
    connect(c, "HltDiPFJetAve400", evt.HltDiPFJetAve400);
    connect(c, "PS_HltDiPFJetAve400", evt.PS_HltDiPFJetAve400);
    
    connect(c, "VtxN", evt.VtxN);
    connect(c, "VtxPosX", evt.VtxPosX);
    connect(c, "VtxPosY", evt.VtxPosY);
    connect(c, "VtxPosZ", evt.VtxPosZ);
    connect(c, "VtxNDof", evt.VtxNDof);
    connect(c, "VtxIsFake", evt.VtxIsFake);

    connect(c, "PUMCNumVtx", evt.PUMCNumVtx);
    connect(c, "PUMCNumVtxOOT", evt.PUMCNumVtxOOT);
    connect(c, "PUMCNumTruth", evt.PUMCNumTruth);
    connect(c, "Rho", evt.Rho);
    connect(c, "Rho25", evt.Rho25);
    
    connect(c, "NobjJet", evt.NobjJet);
    connect(c, "JetPt", evt.JetPt);
    connect(c, "JetPhi", evt.JetPhi);
    connect(c, "JetEta", evt.JetEta);
    connect(c, "JetE", evt.JetE);
    connect(c, "JetArea", evt.JetArea);
    connect(c, "JetIDLoose", evt.JetIDLoose);
    connect(c, "JetIDTight", evt.JetIDTight);
    
    connect(c, "JetGenJetDeltaR", evt.JetGenJetDeltaR);
    connect(c, "GenJetPt", evt.GenJetPt);
    connect(c, "GenJetPhi", evt.GenJetPhi);
    connect(c, "GenJetEta", evt.GenJetEta);
    
    size_t ifirst = d.nskip;
    size_t ilast = n;
    if(d.nmax >= 0){
        ilast = d.nskip + d.nmax;
        if(ilast > n){
            cout << "sample contains " << n << " events, but configured nmax=" << d.nmax << " (nskip=" << d.nskip << ")";
            throw runtime_error("configured more events than the sample contains!");
        }
    }
    
    for(size_t im=0; im<modules.size(); ++im){
        modules[im]->start_dataset(d, outfile);
    }
    
    if(verbose){
        cout << "Starting event loop now from index=" << ifirst << " to index < " << ilast << " (" << (ilast - ifirst) << " events)" << endl;
    }
    
    vector<int> cutflow(modules.size() + 1);
    
    for(size_t i=ifirst; i<ilast; ++i){
        if(verbose){
            if(i%10000==0){
                cout << i << endl;
            }
        }
        c.GetEntry(i);
        cutflow[0] += 1;
        for(size_t im=0; im<modules.size(); ++im){           
            if(!modules[im]->process(evt)) break;// and continue with next event ...
            else{
                ++cutflow[im+1];
            }
        }            
    }
    
    outfile.cd();
    outfile.Write();
    outfile.Close();
    
    if(verbose){
        cout << "Processing of dataset '" << d.name << "' done. Cutflow: " << endl;
        cout << "Total events processed: " << cutflow[0] << endl;
        for(size_t i=0; i<modules.size(); ++i){
            cout << "After module " << i << " (" <<  modules[i]->name() << "): " << cutflow[i+1] << endl;
        }
    }
}

void processor::run(const vector<dataset> & datasets, const string & outpath){
    for(size_t i=0; i<datasets.size(); ++i){
        string outfile = outpath + "/" + datasets[i].name + ".root";
        TFile f(outfile.c_str(), "recreate");
        if(!f.IsOpen()) throw runtime_error("could not open output file '" + outfile + "'");
        run(datasets[i], f);
    }
}

