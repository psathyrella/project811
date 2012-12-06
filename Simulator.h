#ifndef SIMULATOR
#define SIMULATOR

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

#include "RooPlot.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooConstVar.h"

#include "EtaPdf.h"
#include "hi.h"

#define e 2.71828182846

using namespace std;
using namespace RooFit;

//========================================================================================
class Track
{
public:
  Track(float pt, float eta, float phi, float mass, float dzVal);
  void dump();
  TLorentzVector vec;
  float dz;
};
//========================================================================================
class Event
{
public:
  vector<Track> tracks;
};
//========================================================================================
class Simulator
{
public:
  Simulator(float etaMinVal=0, float etaMaxVal=5, float phiMinVal=0, float phiMaxVal=0.5);
  ~Simulator();
  void generateZVals(int n);
  void generateThetaVals(int n);
  void generatePhiVals(int n);
  void generate(int n);
  void readHiFile(TString fname, int n=99999999);
  void plotZ();
  void plotTheta();
  void plotPhi();
  int getNTracks();
  float theta(float eta);
  float eta(float theta);

  Event *event;
  float etaMin,etaMax,phiMin,phiMax;
  RooRealVar *zVar,*zSig,*thetaVar;
  RooGaussian *zGaus;
  EtaPdf *etaPdf;
  RooDataSet *zData,*thetaData;
  RooPlot *zFrame,*thetaFrame;
  vector<float> zVals,thetaVals,phiVals;
};
#endif
