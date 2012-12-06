#include "hi.h"
//----------------------------------------------------------------------------------------
hi::hi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output_test.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("output_test.root:/ana");
      dir->GetObject("hi",tree);

   }
   Init(tree);
}

//----------------------------------------------------------------------------------------
hi::~hi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

//----------------------------------------------------------------------------------------
Int_t hi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
//----------------------------------------------------------------------------------------
Long64_t hi::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

//----------------------------------------------------------------------------------------
void hi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("pdg", pdg, &b_pdg);
   fChain->SetBranchAddress("chg", chg, &b_chg);
   fChain->SetBranchAddress("sta", sta, &b_sta);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   Notify();
}

//----------------------------------------------------------------------------------------
Bool_t hi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void hi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t hi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//----------------------------------------------------------------------------------------
void hi::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);

    for(unsigned ipart=0; ipart<npart; ipart++) {
      // thetaVals.push_back(theta(eta[ipart]));
      // phiVals.push_back(phi[ipart]);
    }
    break;
  }
  // sort(zVals.begin(), zVals.end()); // no z vals in the tree yet
  // sort(thetaVals.begin(), thetaVals.end());
  // sort(phiVals.begin(), phiVals.end());
}
