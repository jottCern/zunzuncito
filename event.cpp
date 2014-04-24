#include "event.hpp"
#include <algorithm> // for std::sort

using namespace std;

namespace{

// Helper class to sort jets by pt; to be used to sort a vector of indices, as third argument to std::sort.
// compare to integers i,j; a pair (i,j) is considered sorted if
// and only if pts[i] > pts[j], where pts oif the vector given in the constructor.
struct index_greater_pt{
    float * pts;
    
    explicit index_greater_pt(float * pts_): pts(pts_){}
    
    bool operator()(size_t i, size_t j)const {
        return pts[i] > pts[j];
    }
};

// do
// collection[i] = collection[indices[i]]
// for all i.
template<typename T>
void apply_permutation(T * collection, const vector<size_t> & indices){
    T coll2[indices.size()];
    std::copy(collection, collection + indices.size(), coll2);
    for(size_t i=0; i<indices.size(); ++i){
        collection[i] = coll2[indices[i]];
    }
}


}

void event::sort_jets_pt(int njets){
    int nmax = min(NobjJet, 50);
    if(njets > 0 and njets < NobjJet){
        nmax = njets;
    }
    
    // first, sort an index-vector to get the permutation which sorts the jet collection:
    vector<size_t> indices(nmax);
    for(size_t i=0; i<indices.size(); ++i){
        indices[i] = i;
    }
    sort(indices.begin(), indices.end(), index_greater_pt(JetPt));
    
    // now using the sorted index-vector, re-arrange ALL event content that uses jet indices.
    apply_permutation(JetPt, indices);
    apply_permutation(JetEta, indices);
    apply_permutation(JetE, indices);
    apply_permutation(JetPhi, indices);
    apply_permutation(JetArea, indices);
    apply_permutation(JetIDLoose, indices);
    apply_permutation(JetIDTight, indices);
    
    apply_permutation(JetGenJetDeltaR, indices);
    apply_permutation(GenJetPt, indices);
    apply_permutation(GenJetEta, indices);
    apply_permutation(GenJetPhi, indices);

    apply_permutation(GenPartPt_algo, indices);
    apply_permutation(GenPartEta_algo, indices);
    apply_permutation(GenPartPhi_algo, indices);
    apply_permutation(GenPartE_algo, indices);
    apply_permutation(GenPartId_algo, indices);
    apply_permutation(GenPartPt_phys, indices);
    apply_permutation(GenPartEta_phys, indices);
    apply_permutation(GenPartPhi_phys, indices);
    apply_permutation(GenPartE_phys, indices);
    apply_permutation(GenPartId_phys, indices);
}

