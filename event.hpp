#ifndef EVENT_HPP
#define EVENT_HPP

#include "TDataType.h"

// information about a dataset to process; needed by the processor.
class dataset{
public:
    // whether or not this dataset is mc. This is used for jet energy corrections, reweighting, etc.
    bool mc;
    bool jetdata;
    bool jethtdata;
    bool jetmondata;
        
    // the filenames that make up this dataset:
    std::vector<std::string> files;
    
    // how many events to process in total, and how many to skip form the beginning.
    // nmex = -1 means to process all events.
    int nmax, nskip;
    
    // name of this dataset; should be unique. This is used e.g. to build the output root filename
    std::string name;
    
    explicit dataset(const std::string & name_): nmax(-1), nskip(0), name(name_){}
};


// structure of an event. This is a stripped-down version of what
// DiJetTree->MakeClass produced.
struct event{
    // general event info:
    UInt_t RunNumber;
    UInt_t LumiBlockNumber;
    UInt_t EventNumber;
    Float_t Weight;
    bool is_mc; // not in the tree; filled via dataset
    bool is_jetdata; // not in the tree; filled via dataset
    bool is_jethtdata; // not in the tree; filled via dataset
    bool is_jetmondata; // not in the tree; filled via dataset
    
    // triggers:
    Bool_t HltDiPFJetAve40;
    Int_t  PS_HltDiPFJetAve40;
    Bool_t HltDiPFJetAve80;
    Int_t  PS_HltDiPFJetAve80;
    Bool_t HltDiPFJetAve140;
    Int_t  PS_HltDiPFJetAve140;
    Bool_t HltDiPFJetAve200;
    Int_t  PS_HltDiPFJetAve200;
    Bool_t HltDiPFJetAve260;
    Int_t  PS_HltDiPFJetAve260;
    Bool_t HltDiPFJetAve320;
    Int_t  PS_HltDiPFJetAve320;
    Bool_t HltDiPFJetAve400;
    Int_t  PS_HltDiPFJetAve400;
    
    // vertex info:
    Int_t           VtxN;
    Float_t         VtxPosX;
    Float_t         VtxPosY;
    Float_t         VtxPosZ;
    Float_t         VtxNDof;
    Bool_t          VtxIsFake;
    
    // pileup-related:
    Int_t           PUMCNumVtx;
    Int_t           PUMCNumVtxOOT;
    Float_t         PUMCNumTruth;
    Float_t         Rho;
    Float_t         Rho25;
    
    // jets:
    Int_t           NobjJet;
    Float_t         JetPt[50];   //[NobjJet]
    Float_t         JetPhi[50];   //[NobjJet]
    Float_t         JetEta[50];   //[NobjJet]
    Float_t         JetE[50];   //[NobjJet]
    Float_t         JetArea[50];   //[NobjJet]
    Bool_t          JetIDLoose[50];   //[NobjJet]
    Bool_t          JetIDTight[50];   //[NobjJet]
    
    Float_t         JetGenJetDeltaR[50];   //[NobjJet]
    Float_t         GenJetPt[50];   //[NobjJet]
    Float_t         GenJetPhi[50];   //[NobjJet]
    Float_t         GenJetEta[50];   //[NobjJet]

    Float_t         GenEvtScale;

    Float_t         GenPartPt_algo[50];
    Float_t         GenPartEta_algo[50]; 
    Float_t         GenPartPhi_algo[50];
    Float_t         GenPartE_algo[50];
    Int_t           GenPartId_algo[50];
    Float_t         GenPartPt_phys[50];
    Float_t         GenPartEta_phys[50]; 
    Float_t         GenPartPhi_phys[50];
    Float_t         GenPartE_phys[50];
    Int_t           GenPartId_phys[50];


    // to add events info:
    // 1. add the corresponding data member here
    // 2. make sure to make the corresponding SetBranchAddress call by calling connect in the processor.cpp
    // 3. if it is a jet-related data (i.e. array of length NobjJet), make sure to modify sort_jets_pt to also re-shuffle this new data member.
    

    // sort all data in the jet-related vectors such that the jets are sorted descending in pt.
    // This method should be used after modifying jets by jet energy corrections, resolution smearing, etc. 
    // If njets!=-1, then only the first njets jets are re-sorted, the rest is leaved untouched. (If njets==-1,
    // all 50 jets are re-sorted).
    void sort_jets_pt(int njets = -1);
};

#endif
