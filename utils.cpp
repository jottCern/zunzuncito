#include "utils.hpp"
#include "event.hpp"
#include <glob.h>
#include <stdexcept>
#include <iostream>

using namespace std;

std::vector<std::string> glob(const std::string& pattern){
    glob_t glob_result;
    int res = glob(pattern.c_str(), GLOB_TILDE | GLOB_BRACE | GLOB_ERR | GLOB_MARK, NULL, &glob_result);
    if(res != 0){
        string errmsg;
        if(res == GLOB_NOSPACE) errmsg = "out of memory";
        else if(res == GLOB_ABORTED) errmsg = "I/O error";
        else if(res == GLOB_NOMATCH) errmsg = "no match";
        else errmsg = "unknwon glob return value";
        throw runtime_error("glob error for pattern '" + pattern + "': " + errmsg);
    }
    vector<string> ret;
    ret.reserve(glob_result.gl_pathc);
    for(unsigned int i=0; i<glob_result.gl_pathc; ++i){
        ret.push_back(glob_result.gl_pathv[i]);
    }
    globfree(&glob_result);
    return ret;
}


void dump(const event & evt, int njetsmax, bool with_mc_info){
    cout << "event " << evt.RunNumber << ":" << evt.EventNumber << " (ls=" << evt.LumiBlockNumber << ")";
    if(with_mc_info){
        cout << ", weight=" << evt.Weight << endl;
    }
    else{
        cout << endl;
    }
    cout << "  rho=" << evt.Rho << ", rho25=" << evt.Rho25 << ", nvtx=" << evt.VtxN << endl;
    cout << "  vertex (x,y,z) = (" << evt.VtxPosX << ", " << evt.VtxPosY << ", " << evt.VtxPosZ << "); ndof = " << evt.VtxNDof << endl;
    int n = min(evt.NobjJet, njetsmax);
    for(int i=0; i<n; ++i){
        cout << "  jet " << i << ": pt=" << evt.JetPt[i] << ", eta=" << evt.JetEta[i] << ", phi=" << evt.JetPhi[i] << ", area=" << evt.JetArea[i];
        cout << " (loose id: " << (evt.JetIDLoose[i] ? "passed" : "FAILED") << ")" << endl;
        if(with_mc_info){
            cout << "     genjet: dr=" << evt.JetGenJetDeltaR[i] << ", pt=" << evt.GenJetPt[i] << endl;
        }
    }
}

