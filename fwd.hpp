#ifndef FWD_HPP
#define FWD_HPP

// This header declares some classes. It is meant to be included (instead of the full header)
// in other .hpp files in cases when such a forward declaration is enough.

#include <boost/shared_ptr.hpp>

using boost::shared_ptr;

struct event;
class module;
class dataset;

class FactorizedJetCorrector;
class JetCorrectionUncertainty;

class TFile;
class TH1F;

#endif
