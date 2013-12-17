// $Id: utils.cc,v 1.3 2011/09/16 11:49:09 mschrode Exp $

#include "ConfigParser.h"
#include "FileOps.h"
#include "HistOps.h"
#include "LabelFactory.h"
#include "StyleSettings.h"
#include "utils.h"


std::vector<TH1*> util::LabelFactory::hDummies_;
unsigned int util::HistOps::COUNT_ = 0;





