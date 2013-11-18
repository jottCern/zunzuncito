#ifndef MODULES_HPP
#define MODULES_HPP

#include "processor.hpp"
#include "fwd.hpp"
#include <string>


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
    explicit asymm_histos(const shared_ptr<pt_binning> & binning = shared_ptr<pt_binning>(new standard_ptave_binning), const std::string & dir = "");
    
    virtual void start_dataset(const dataset & d, TFile & outfile);
    virtual bool process(event & evt);
    virtual std::string name(){
        return "asymmetry histograms";
    }
    
private:
    shared_ptr<pt_binning> ptbinning;
    std::string dir;
    std::vector<TH1F*> histos;
    
    // those methods define the binning in eta and alpha. They
    // are called from 'start_dataset' and 'process', which are
    // independent of the pt and eta binning.
    int nbins_eta() const;
    int ibin_eta(const event & evt) const;
    
    int nbins_alpha() const;
    int ibin_alpha(const event & evt) const;
};


// apply standard filters: good vertex (ndof>=4, rho<2, |z|<24) and beam-halo veto
class eventfilter: public module {
public:
    virtual bool process(event & evt);
    virtual ~eventfilter();
    
    virtual std::string name(){
       return "eventfilter";
    }
};

#endif
