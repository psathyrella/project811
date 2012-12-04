#include "Simulator.h"

//----------------------------------------------------------------------------------------
Simulator::Simulator(float etaMinVal, float etaMaxVal)
{
  etaMin = etaMinVal;
  etaMax = etaMaxVal;
  zVar  = new RooRealVar("zVar","zVar",-20,20);
  zSig  = new RooRealVar("zSig","zSig",5);
  zGaus = new RooGaussian("zGaus","zGaus",*zVar,RooConst(0),*zSig);
  thetaVar = new RooRealVar("thetaVar","thetaVar",theta(etaMax),theta(etaMin));
  thetaVar->Print();
  etaPdf = new EtaPdf("etaPdf","etaPdf",*thetaVar);
  zData = 0;
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
}
//----------------------------------------------------------------------------------------
void Simulator::generateZVals(int n)
{
  zData = zGaus->generate(*zVar,n);
  for(int ientry=0; ientry<zData->numEntries(); ientry++) {
    float val = ((RooRealVar*)(zData->get(ientry)->find("zVar")))->getVal();
    zVals.push_back(val);
  }
  sort(zVals.begin(), zVals.end());
}
//----------------------------------------------------------------------------------------
void Simulator::generateThetaVals(int n)
{
  thetaData = etaPdf->generate(*thetaVar,n);
  for(int ientry=0; ientry<thetaData->numEntries(); ientry++) {
    float val = ((RooRealVar*)(thetaData->get(ientry)->find("thetaVar")))->getVal();
    thetaVals.push_back(val);
  }
  sort(thetaVals.begin(), thetaVals.end());
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
