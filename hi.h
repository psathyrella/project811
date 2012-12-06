#ifndef HI
#define HI

#include <iostream>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <math.h>
#include <fstream>

#include "TH1F.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"

using namespace std;
class hi {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event;
   Int_t           npart;
   Float_t         pt[892];   //[npart]
   Float_t         eta[892];   //[npart]
   Float_t         phi[892];   //[npart]
   Int_t           pdg[892];   //[npart]
   Int_t           chg[892];   //[npart]
   Int_t           sta[892];   //[npart]
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_npart;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_chg;   //!
   TBranch        *b_sta;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!

   hi(TTree *tree=0);
   virtual ~hi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};
#endif
