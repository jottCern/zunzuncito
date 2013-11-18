#ifndef PROCESSOR_HPP
#define PROCESSOR_HPP

#include "fwd.hpp"
#include <boost/shared_ptr.hpp>
#include <vector>

// abstract base class for all modules, as used by the processor
class module{
public:
    virtual void start_dataset(const dataset & d, TFile & outfile){}
    virtual bool process(event & evt) = 0;
    virtual std::string name() = 0;
    virtual ~module(){}
};

// This is the class which performs the event loop in its "run" methods.
// It keeps an (ordered) list of modules to run on the events. For each dataset
//   module::start_dataset
// is called for all modules in the list. For each event in the dataset
//    module::process
// is called. If the return value of module::process is false, processing of that event is stopped and
// the next event is read from the input; this allows to implement selections / shortcut processing of
// irrelevant events.
class processor{
public:
    // append the given module to the module list
    void add_module(const boost::shared_ptr<module> & m);
    
    // let all modules run on the given dataset:
    void run(const dataset & d, TFile & outfile);

    // call run(dataset, outfile) for each dataset in datasets, where "outfile" refers to a newly
    // created file with the name
    //    outpath + dataset::name + ".root"
    void run(const std::vector<dataset> & datasets, const std::string & outpath);

    processor(): verbose(true){}

    // set verbosity on / off. If on, more progress information will be printed.
    void set_verbose(bool v){
        verbose = v;
    }
    
private:
    std::vector<boost::shared_ptr<module> > modules;
    bool verbose;
};

#endif
