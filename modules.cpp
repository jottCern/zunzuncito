#include "modules.hpp"
#include "utils.hpp"
#include "event.hpp"

#include "JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TFile.h"
#include "TH1F.h"

#include <sstream>

using namespace std;

bool dumper::process(event & evt){
    std::cout << "dumper(" << prefix << "): ";
    dump(evt, nmax, evt.is_mc);
    return true;
}

jec::jec(const std::string & path_, int njetsmax_, const std::string & globaltag_mc, const std::string & globaltag_data): njetsmax(njetsmax_), path(path_), gt_data(globaltag_data), gt_mc(globaltag_mc){
    // make sure path ends with a '/':
    if(path.size() > 0){
        if(path[path.size()-1]!='/'){
            path += '/';
        }
    }
}

void jec::start_dataset(const dataset & d, TFile & outfile){
    std::vector<JetCorrectorParameters> pars;
    const std::string & globaltag = d.mc ? gt_mc : gt_data;
    pars.push_back(path + globaltag + "_L1FastJet_AK5PFchs.txt");
    pars.push_back(path + globaltag + "_L2Relative_AK5PFchs.txt");
    pars.push_back(path + globaltag + "_L3Absolute_AK5PFchs.txt");
    if(!d.mc){
        pars.push_back(path + globaltag + "_L2L3Residual_AK5PFchs.txt");
    }
    fjc.reset(new FactorizedJetCorrector(pars));
}

bool jec::process(event & evt){
    // correct all jets:
    int nmax = evt.NobjJet;
    if(njetsmax > 0 && njetsmax < nmax){
        nmax = njetsmax;
    }
    for(int i=0; i<nmax; ++i){
        fjc->setJetPt(evt.JetPt[i]);
        fjc->setJetEta(evt.JetEta[i]);
        fjc->setJetE(evt.JetE[i]);
        fjc->setJetA(evt.JetArea[i]);
        fjc->setRho(evt.Rho);
        float corr = fjc->getCorrection();
        evt.JetPt[i] *= corr;
        evt.JetE[i] *= corr;
    }
    evt.sort_jets_pt(nmax);
    return true;
}

jec::~jec(){}

jec_uncertainty::jec_uncertainty(const string & path_txtfile, syst_dir::e_syst_dir sdir_, int njetsmax_): sdir(sdir_), njetsmax(njetsmax_){
    jec_unc.reset(new JetCorrectionUncertainty(path_txtfile));
}

bool jec_uncertainty::process(event & evt){
    int nmax = evt.NobjJet;
    if(njetsmax > 0 && njetsmax < nmax){
        nmax = njetsmax;
    }
    for(int i=0; i < nmax; ++i){
        jec_unc->setJetEta(evt.JetEta[i]);
        jec_unc->setJetPt(evt.JetPt[i]);
        float correctionfactor;
        if(sdir == syst_dir::up){
            float unc = jec_unc->getUncertainty(1);
            correctionfactor = (1 + fabs(unc));
        } else {
            float unc = jec_unc->getUncertainty(-1);
            correctionfactor = (1 - fabs(unc));
        }
        evt.JetPt[i] *= correctionfactor;
        evt.JetE[i] *= correctionfactor;
    }
    evt.sort_jets_pt(nmax);
    return true;
}

jec_uncertainty::~jec_uncertainty(){}


bool eventfilter::process(event & evt){
    // 1. we must have a vertex with rho < 2cm, abs(z) < 24cm and ndof >= 4, so reject the event otherwise:
    if(evt.VtxNDof < 4.0f or std::abs(evt.VtxPosZ) > 24.0f or hypotf(evt.VtxPosX, evt.VtxPosY) > 2.0f) return false;
    // 2. TODO: beam-halo filter!
    
    // if it passed both: keep it:
    return true;
}

eventfilter::~eventfilter(){}



int standard_ptave_binning::nbins() const {
    // TODO: number of bins!
    return 1;
}


int standard_ptave_binning::ibin(const event & evt) const{
    if(evt.NobjJet < 2) return -1;
    float ptave = 0.5 * (evt.JetPt[0] + evt.JetPt[1]);
    // TODO: depending on ptave, ask also for trigger, but only for data (use evt.is_mc)
    return 0;
    
}

standard_ptave_binning::~standard_ptave_binning(){}


asymm_histos::asymm_histos(const shared_ptr<pt_binning> & binning, const std::string & dir_): ptbinning(binning), dir(dir_){
}


int asymm_histos::nbins_eta() const{
    return 4;
}

int asymm_histos::ibin_eta(const event & evt) const{
    if(evt.NobjJet < 2) return -1;
    float eta0 = abs(evt.JetEta[0]);
    float eta1 = abs(evt.JetEta[1]);
    const float eta_high[] = {1.0f, 1.5f, 2.3f, 4.7f};
    const int neta = sizeof(eta_high) / sizeof(float);
    for(int ieta = 0; ieta < neta; ++ieta){
        // both the leading and the sub-leading jet must be within the same eta bin extending from eta_lo to eta_hi:
        float eta_hi = eta_high[ieta];
        float eta_lo = ieta == 0? 0.0f : eta_high[ieta-1];
        if(eta0 >= eta_lo and eta0 < eta_hi and eta1 >= eta_lo and eta1 < eta_hi) return ieta;
    }
    return -1;
}
   


int asymm_histos::nbins_alpha() const{
    return 6;
}

int asymm_histos::ibin_alpha(const event & evt) const{
    if(evt.NobjJet < 3) return -1;
    // reject events with small third jet pt. TODO: at which jet pt threshold exactly?
    if(evt.JetPt[2] < 6.0f) return -1;
    float alpha = evt.JetPt[2] / (0.5f * (evt.JetPt[0] + evt.JetPt[1]));
    const float alpha_min[] = {0.1f, 0.12f, 0.14f, 0.16f, 0.18f, 0.20f}; // lower bin borders of alpha binning
    const int nalpha = sizeof(alpha_min) / sizeof(float);
    // find the bin index in alpha for this event. The last bin extends to +infinity implicitly
    for(int ialpha=0; ialpha < nalpha; ++ialpha){
        if(alpha >= alpha_min[ialpha] and (ialpha + 1 == nalpha or alpha < alpha_min[ialpha + 1])){
            return ialpha;
        }
    }
    return -1;
}



void asymm_histos::start_dataset(const dataset & d, TFile & outfile){
    // clear the histograms from previous dataset:
    histos.clear();
    // create all histograms in the output file:
    outfile.cd();
    if(!dir.empty()){
        outfile.mkdir(dir.c_str());
        outfile.cd(dir.c_str());
    }
    // create the asymmetry histograms: histograms from 0 to 1 with 1000 bins:
    const int n_pt = ptbinning->nbins();
    const int n_alpha = nbins_alpha();
    const int n_eta = nbins_eta();
    for(int ipt=0; ipt<n_pt; ++ipt){
        for(int ieta=0; ieta < n_eta; ++ieta){
            for(int ialpha=0; ialpha < n_alpha; ++ialpha){
                stringstream ss;
                ss << "pt" << ipt << "_eta" << ieta << "_alpha" << ialpha;
                histos.push_back(new TH1F(ss.str().c_str(), ss.str().c_str(), 1000, 0.0, 1.0));
            }
        }
    }
}



bool asymm_histos::process(event & evt){
    if(evt.NobjJet < 3) return false;
    int ptbin = ptbinning->ibin(evt);
    if(ptbin < 0) return false;
    assert(ptbin < ptbinning->nbins());
    
    int etabin = ibin_eta(evt);
    if(etabin < 0) return false;
    assert(etabin < nbins_eta());
    
    int alphabin = ibin_alpha(evt);
    if(alphabin < 0) return false;
    assert(alphabin < nbins_alpha());
    
    const int n_pt = ptbinning->nbins();
    const int n_eta = nbins_eta();
    const int n_alpha = nbins_alpha();
    int histo_index = ptbin * (n_eta * n_alpha) + etabin * n_alpha + alphabin;
    assert(histo_index < int(histos.size()));
    float asymmetry = (evt.JetPt[0] - evt.JetPt[1]) / (evt.JetPt[0] + evt.JetPt[1]);
    // asymmetry defined that way should always be positive, if the pt-sorting was right:
    assert(asymmetry >= 0.0f);
    histos[histo_index]->Fill(asymmetry, evt.Weight);
    return true;
}

