#include "Simulator.h"

//----------------------------------------------------------------------------------------
Track::Track(float pt, float eta, float phi, float mass, float dzVal)
{
  vec.SetPtEtaPhiM(pt,eta,phi,mass);
  dz = dzVal;
}
//----------------------------------------------------------------------------------------
void Track::dump()
{
  cout
    << setw(12) << vec.Pt()
    << setw(12) << vec.Eta()
    << setw(12) << vec.Phi()
    << setw(12) << vec.M()
    << setw(12) << dz
    << endl;
}
//----------------------------------------------------------------------------------------
Simulator::Simulator(float etaMinVal, float etaMaxVal, float phiMinVal, float phiMaxVal)
{
  event = new Event;
  
  etaMin = etaMinVal;
  etaMax = etaMaxVal;
  phiMin = phiMinVal;
  phiMax = phiMaxVal;

  zVar  = new RooRealVar("zVar","zVar",-20,20);
  zSig  = new RooRealVar("zSig","zSig",3);
  zGaus = new RooGaussian("zGaus","zGaus",*zVar,RooConst(0),*zSig);
  thetaVar = new RooRealVar("thetaVar","thetaVar",theta(etaMax),theta(etaMin));
  etaPdf = new EtaPdf("etaPdf","etaPdf",*thetaVar);
  phiVar  = new RooRealVar("phiVar","phiVar",phiMin,phiMax);
  phiPdf = new RooUniform("phiPdf","phiPdf",*phiVar);
  zData = 0;
  phiData = 0;
  zFrame = 0;
  thetaFrame = 0;
  thetaData = 0;
}
//----------------------------------------------------------------------------------------
Simulator::~Simulator()
{
  delete zVar;
  delete zSig;
  delete zGaus;
  delete thetaVar;
  delete etaPdf;
}
//----------------------------------------------------------------------------------------
void Simulator::generateZVals(int n)
{
  zData = zGaus->generate(*zVar,n);
  for(int ientry=0; ientry<zData->numEntries(); ientry++) {
    float val = ((RooRealVar*)(zData->get(ientry)->find("zVar")))->getVal();
    zVals.push_back(val);
  }
  // sort(zVals.begin(), zVals.end());
}
//----------------------------------------------------------------------------------------
void Simulator::generateThetaVals(int n)
{
  thetaData = etaPdf->generate(*thetaVar,n);
  for(int ientry=0; ientry<thetaData->numEntries(); ientry++) {
    float val = ((RooRealVar*)(thetaData->get(ientry)->find("thetaVar")))->getVal();
    thetaVals.push_back(val);
  }
  // sort(thetaVals.begin(), thetaVals.end());
}
//----------------------------------------------------------------------------------------
void Simulator::generatePhiVals(int n)
{
  phiData = phiPdf->generate(*phiVar,n);
  for(int ientry=0; ientry<phiData->numEntries(); ientry++) {
    float val = ((RooRealVar*)(phiData->get(ientry)->find("phiVar")))->getVal();
    phiVals.push_back(val);
  }
}
//----------------------------------------------------------------------------------------
void Simulator::plotZ()
{
  zFrame = zVar->frame(Title("-"));
  assert(zData);
  zData->plotOn(zFrame);
  zGaus->plotOn(zFrame);
  zFrame->Draw();
}
//----------------------------------------------------------------------------------------
void Simulator::plotTheta()
{
  thetaFrame = thetaVar->frame(Title("-"));
  assert(thetaData);
  thetaData->plotOn(thetaFrame);
  etaPdf->plotOn(thetaFrame);
  thetaFrame->Draw();
}
//----------------------------------------------------------------------------------------
int Simulator::getNTracks()
{
  // around 100 charged particles per unit eta
  return 100*(etaMax - etaMin);
}
//----------------------------------------------------------------------------------------
float Simulator::theta(float eta)
{
  return 2*atan(pow(e,-eta));
}
//----------------------------------------------------------------------------------------
float Simulator::eta(float theta)
{
  return -log(tan(theta/2));
}
//----------------------------------------------------------------------------------------
void Simulator::generate(int n)
{
  generateZVals(n);
  generateThetaVals(n);
  generatePhiVals(n);
  
  for(unsigned itrk=0; itrk<n; itrk++) {
    event->tracks.push_back(Track(5, eta(thetaVals[itrk]), phiVals[itrk], 0, zVals[itrk]));
  }
}
//----------------------------------------------------------------------------------------
void Simulator::readHiFile(TString fname, int n)
{
  // gROOT->Macro("hi.C");
  TFile file(fname);
  file.cd("ana");
  TTree *tree = (TTree*)file.Get("hi");
  hi *lyzer = new hi(tree);

  int nMax = min(lyzer->fChain->GetEntriesFast(),Long64_t(n));
  for (Long64_t jentry=0; jentry<nMax;jentry++) {
    Long64_t ientry = lyzer->LoadTree(jentry);
    if (ientry < 0) break;
    lyzer->fChain->GetEntry(jentry);

    for(unsigned ipart=0; ipart<lyzer->npart; ipart++) {
      cout << lyzer->eta[ipart] << setw(12) << lyzer->phi[ipart] << endl;
      thetaVals.push_back(theta(lyzer->eta[ipart]));
      phiVals.push_back(lyzer->phi[ipart]);
      if(thetaVals.size() >= nMax) break;
    }
    if(thetaVals.size() >= nMax) break;
  }
  // sort(zVals.begin(), zVals.end()); // no z vals in the tree yet
  generateZVals(thetaVals.size());
  // sort(thetaVals.begin(), thetaVals.end());
  // sort(phiVals.begin(), phiVals.end());
  for(unsigned itrk=0; itrk<nMax; itrk++) {
    event->tracks.push_back(Track(5, eta(thetaVals[itrk]), phiVals[itrk], 0, zVals[itrk]));
    cout << "pushing track: ";
    event->tracks.back().dump();
  }
}
