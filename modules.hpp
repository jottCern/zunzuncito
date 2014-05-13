#ifndef MODULES_HPP
#define MODULES_HPP

#include "processor.hpp"
#include "fwd.hpp"
#include <string>
#include "TF1.h"
#include "TH2.h"


// call "dump" (from utils) for each event, with a configureable prefix. This is a quick-and-dirty method
// to inspect manually some events.
class dumper: public module{
public:
    dumper(const std::string & prefix_, int nmax_): prefix(prefix_), nmax(nmax_){}
    
    virtual bool process(event & evt);
    
    virtual ~dumper(){}
    
    virtual std::string name(){
        return "dumper";
    }
private:
    std::string prefix;
    int nmax;
};


// apply jet energy corrections on all jets in the event, and re-sort the jets descending in pt.
// NOTE: running the jet energy correction on all 50 jets in the event takes
// quite some time. So one could think about correction only the first 5 or so
// which should be enough to get updated values for ptaverage, alpha, even if there is
// some re-ordering after the correction.
class jec: public module{
public:
    
   // path is the path for the txt files. The name of the text files to use is constructed
   // using globaltag_data and globaltag_mc, resp. The jet type is always AK5PFchs. The correction is always done
   // to level 3 (including the residual correction for data).
   // njetsmax is the number of jets to correct. The default of -1 corrects ALL (50) jets
   // in the event, which is the correct thing to do but costs time.
   explicit jec(const std::string & path_, int njetsmax = -1, const std::string & globaltag_mc = "START53_V27", const std::string & globaltag_data = "FT_53_V21_AN6");
   
   virtual void start_dataset(const dataset & d, TFile & outfile);
   virtual bool process(event & evt);
   
   virtual ~jec();
   
   virtual std::string name(){
       return "jec";
   }
   
private:
    int njetsmax;
    std::string path, gt_data, gt_mc;
    std::auto_ptr<FactorizedJetCorrector> fjc;
};


// enumeration for the direction of a systematic uncertainty:
namespace syst_dir {
enum e_syst_dir{
    up, down
};
}

// apply a jet energy shift according to the jet energy uncertainties on all jets, and re-sort the jets descending in pt.
// This module should be run on *corrected* jets, i.e. after running jec.
class jec_uncertainty: public module{
public:
    
    // path_txtfile is the path to the uncertainty txt file to use
    // syst_dir is either up or down, depending on whether to increase or decrease all jet momenta.
    // njetsmax controls on how many jets the uncertainty is applied, -1 means to do it for all (50) jets in the event, which however takes some time;
    //    see jec module.
    jec_uncertainty(const std::string & path_txtfile, syst_dir::e_syst_dir syst_dir, int njetsmax = -1);
    
    virtual bool process(event & evt);
    
    virtual ~jec_uncertainty();
    
    virtual std::string name(){
       return "jec_uncertainty";
    }
    
private:
    std::auto_ptr<JetCorrectionUncertainty> jec_unc;
    syst_dir::e_syst_dir sdir;
    int njetsmax;
};

// select events that pass the trigger selection
class checktrigger: public module {
public:
    virtual bool process(event & evt);
    virtual ~checktrigger();
    
    virtual std::string name(){
       return "select events passing trigger";
    }
};

// reweight pile-up scenario of MC sample to match PU distribution in data
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
class pureweighting: public module{
public:
    
   // path is the path to the pile-up distributions needed as input: true_mc is folder for true pile-up distributions in mc and pu_data is the measured/observed pu distribution in data per trigger 
   explicit pureweighting(const std::string & path_, const std::string & true_mc_ = "TrueDistributions/", const std::string & pu_data_ = "Cert_2012_190456-208686_ReReco/");
   // explicit pureweighting(const std::string & path_, const std::string & true_mc_ = "TrueDistributions/", const std::string & pu_data_ = "Cert_2012_190456-208686_ReReco_MBXS73500/");
   
   virtual void start_dataset(const dataset & d, TFile & outfile);
   virtual bool process(event & evt);
   
   virtual ~pureweighting();
   
   virtual std::string name(){
       return "pu reweighting";
   }
   
private:
   std::string path, true_mc, pu_data;
   TH1F *h_DiPFJetAve40, *h_DiPFJetAve80, *h_DiPFJetAve140, *h_DiPFJetAve200, *h_DiPFJetAve260, *h_DiPFJetAve320, *h_DiPFJetAve400;
};


// abstract base class for the pt binning of events.
// Implementations could use both e.g. leading jet pt or average
// pt -- together with the corresponding trigger -- as criterion.
class pt_binning{
public:
    // number of bins:
    virtual int nbins() const = 0;
    
    // the bin number for this event (0..nbins-1). The special
    // value -1 means that this event does not belong to any bin.
    virtual int ibin(const event & evt) const = 0;

    virtual ~pt_binning(){}
};


class standard_ptave_binning: public pt_binning{
public:
    virtual int nbins() const;
    virtual int ibin(const event & evt) const;
    virtual ~standard_ptave_binning();
};

// create the asymmetry histograms: one for each bin in (ptave, eta, alpha).
// The eta and alpha-bins are hard-coded; the pt binning is provided by the pt_binning class.
class asymm_histos: public module{
public:
    // dir is the directory in which to create the output histograms
   explicit asymm_histos(const boost::shared_ptr<pt_binning> & binning = boost::shared_ptr<pt_binning>(new standard_ptave_binning), const std::string & dir = "");
    
    virtual void start_dataset(const dataset & d, TFile & outfile);
    virtual bool process(event & evt);
    virtual std::string name(){
        return "asymmetry histograms";
    }
    
private:
    boost::shared_ptr<pt_binning> ptbinning;
    std::string dir;
    std::vector<TH1F*> histos_asymm;
    std::vector<TH1F*> histos_genasymm;
    std::vector<TH1F*> histos_asymm_forward;
    std::vector<TH1F*> histos_genasymm_forward;
    std::vector<TH1F*> histos_alphaspectrum;
      
    // those methods define the binning in eta and alpha. They
    // are called from 'start_dataset' and 'process', which are
    // independent of the pt and eta binning.
    int nbins_eta() const;
    int ibin_eta(const event & evt) const;
    
    int nbins_alpha() const;
    int ibin_alpha(const event & evt) const;
};

// create the response histograms: one for each bin in (ptgen, geneta).
// The eta bins are hard-coded; the pt binning is provided by the pt_binning class.
class response_histos: public module{
public:
    // dir is the directory in which to create the output histograms
    explicit response_histos(const std::string & dir = "");
    
    virtual void start_dataset(const dataset & d, TFile & outfile);
    virtual bool process(event & evt);
    virtual std::string name(){
        return "response histograms";
    }
    
private:
    std::string dir;
    std::vector<TH1F*> histos_response, histos_response_uds, histos_response_c, histos_response_b, histos_response_gluon;
      
    // those methods define the binning in geneta. They
    // are called from 'start_dataset' and 'process', which are
    // independent of the pt and eta binning.
    int nbins_etagen() const;
    int ibin_etagen(const float & eta) const;

    int nbins_ptgen() const;
    int ibin_ptgen(const float & ptgen) const;
};


// apply standard filters: good vertex (ndof>=4, rho<2, |z|<24) and beam-halo veto e.g.
class eventfilter: public module {
public:
    virtual bool process(event & evt);
    virtual ~eventfilter();
    
    virtual std::string name(){
       return "eventfilter";
    }
};

// apply diJet event selection cuts
class eventcuts: public module {
public:
    virtual bool process(event & evt);
    virtual ~eventcuts();
    
    virtual std::string name(){
       return "event selection cuts";
    }
};

// plot some control distributions
class controlplots: public module {
public:
    // suffix is used to identify histos in output, dir is the directory in which to create the output histograms
   explicit controlplots(const std::string & suffix_, const std::string & dir = "");

    virtual void start_dataset(const dataset & d, TFile & outfile);
    virtual bool process(event & evt);
    virtual ~controlplots();
    
    virtual std::string name(){
       return "plot control distributions";
    }

private:
    std::string suffix, dir;
   TH1F* NVtx_HltDiPFJetAve40, *Jet1Pt_HltDiPFJetAve40, *Jet2Pt_HltDiPFJetAve40, *Jet3Pt_HltDiPFJetAve40, *Jet1Eta_HltDiPFJetAve40, *Jet2Eta_HltDiPFJetAve40, *Jet3Eta_HltDiPFJetAve40, *GenJet1Pt_HltDiPFJetAve40, *GenJet2Pt_HltDiPFJetAve40, *GenJet3Pt_HltDiPFJetAve40, *PtAve_HltDiPFJetAve40, *GenPtAve_HltDiPFJetAve40, *Alpha_HltDiPFJetAve40, *GenAlpha_HltDiPFJetAve40, *DeltaPhi_HltDiPFJetAve40;

   TH1F* NVtx_HltDiPFJetAve80, *Jet1Pt_HltDiPFJetAve80, *Jet2Pt_HltDiPFJetAve80, *Jet3Pt_HltDiPFJetAve80, *Jet1Eta_HltDiPFJetAve80, *Jet2Eta_HltDiPFJetAve80, *Jet3Eta_HltDiPFJetAve80, *GenJet1Pt_HltDiPFJetAve80, *GenJet2Pt_HltDiPFJetAve80, *GenJet3Pt_HltDiPFJetAve80, *PtAve_HltDiPFJetAve80, *GenPtAve_HltDiPFJetAve80, *Alpha_HltDiPFJetAve80, *GenAlpha_HltDiPFJetAve80, *DeltaPhi_HltDiPFJetAve80;

   TH1F* NVtx_HltDiPFJetAve140, *Jet1Pt_HltDiPFJetAve140, *Jet2Pt_HltDiPFJetAve140, *Jet3Pt_HltDiPFJetAve140, *Jet1Eta_HltDiPFJetAve140, *Jet2Eta_HltDiPFJetAve140, *Jet3Eta_HltDiPFJetAve140, *GenJet1Pt_HltDiPFJetAve140, *GenJet2Pt_HltDiPFJetAve140, *GenJet3Pt_HltDiPFJetAve140, *PtAve_HltDiPFJetAve140, *GenPtAve_HltDiPFJetAve140, *Alpha_HltDiPFJetAve140, *GenAlpha_HltDiPFJetAve140, *DeltaPhi_HltDiPFJetAve140;

    TH1F* NVtx_HltDiPFJetAve200, *Jet1Pt_HltDiPFJetAve200, *Jet2Pt_HltDiPFJetAve200, *Jet3Pt_HltDiPFJetAve200,*Jet1Eta_HltDiPFJetAve200, *Jet2Eta_HltDiPFJetAve200, *Jet3Eta_HltDiPFJetAve200, *GenJet1Pt_HltDiPFJetAve200, *GenJet2Pt_HltDiPFJetAve200, *GenJet3Pt_HltDiPFJetAve200, *PtAve_HltDiPFJetAve200, *GenPtAve_HltDiPFJetAve200, *Alpha_HltDiPFJetAve200, *GenAlpha_HltDiPFJetAve200, *DeltaPhi_HltDiPFJetAve200;

   TH1F* NVtx_HltDiPFJetAve260, *Jet1Pt_HltDiPFJetAve260, *Jet2Pt_HltDiPFJetAve260, *Jet3Pt_HltDiPFJetAve260, *Jet1Eta_HltDiPFJetAve260, *Jet2Eta_HltDiPFJetAve260, *Jet3Eta_HltDiPFJetAve260, *GenJet1Pt_HltDiPFJetAve260, *GenJet2Pt_HltDiPFJetAve260, *GenJet3Pt_HltDiPFJetAve260, *PtAve_HltDiPFJetAve260, *GenPtAve_HltDiPFJetAve260, *Alpha_HltDiPFJetAve260, *GenAlpha_HltDiPFJetAve260, *DeltaPhi_HltDiPFJetAve260;

   TH1F* NVtx_HltDiPFJetAve320, *Jet1Pt_HltDiPFJetAve320, *Jet2Pt_HltDiPFJetAve320, *Jet3Pt_HltDiPFJetAve320, *Jet1Eta_HltDiPFJetAve320, *Jet2Eta_HltDiPFJetAve320, *Jet3Eta_HltDiPFJetAve320, *GenJet1Pt_HltDiPFJetAve320, *GenJet2Pt_HltDiPFJetAve320, *GenJet3Pt_HltDiPFJetAve320, *PtAve_HltDiPFJetAve320, *GenPtAve_HltDiPFJetAve320, *Alpha_HltDiPFJetAve320, *GenAlpha_HltDiPFJetAve320, *DeltaPhi_HltDiPFJetAve320;

   TH1F* NVtx_HltDiPFJetAve400, *Jet1Pt_HltDiPFJetAve400, *Jet2Pt_HltDiPFJetAve400, *Jet3Pt_HltDiPFJetAve400, *Jet1Eta_HltDiPFJetAve400, *Jet2Eta_HltDiPFJetAve400, *Jet3Eta_HltDiPFJetAve400, *GenJet1Pt_HltDiPFJetAve400, *GenJet2Pt_HltDiPFJetAve400, *GenJet3Pt_HltDiPFJetAve400, *PtAve_HltDiPFJetAve400, *GenPtAve_HltDiPFJetAve400, *Alpha_HltDiPFJetAve400, *GenAlpha_HltDiPFJetAve400, *DeltaPhi_HltDiPFJetAve400;
};

// use only reasonable MC events --> remove events with large weight and large response
class reasonablemceventchecker: public module {
public:
    virtual bool process(event & evt);
    virtual ~reasonablemceventchecker();
    
    virtual std::string name(){
       return "remove mc events with high weight and large response";
    }
};

// smear MC module
class smearmc: public module {
public:
    explicit smearmc(int njetsmax = -1);
    virtual bool process(event & evt);
    virtual ~smearmc();
    
    virtual std::string name(){
       return "smear mc events with certain factor/resolution";
    }

private:
    int njetsmax;
    float getsmearfactor(float eta);
};

// reweight MC events to match alpha spectrum in data 
class alphareweighting: public module {
public:
    virtual void start_dataset(const dataset & d, TFile & outfile);
    virtual bool process(event & evt);
    virtual ~alphareweighting();
    
    virtual std::string name(){
       return "reweight mc events to match alpha spectrum in data";
    }

private:
    TF1 *f;
};

// reweight MC events depending on gluon splitting 
class gluonsplittingreweighting: public module {
public:
    virtual void start_dataset(const dataset & d, TFile & outfile);
    virtual bool process(event & evt);
    virtual ~gluonsplittingreweighting();
    
    virtual std::string name(){
       return "reweight mc events depending on gluon splitting";
    }

private:
   TH1F *splitting_fraction_before, *splitting_fraction_after;
};

#endif
 
