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
#include "Util.h"

using namespace std;
using namespace RooFit;
//========================================================================================
class Simulator
{
public:
  Simulator(float etaMinVal, float etaMaxVal);
  ~Simulator();
  void generateZVals(int n);
  void plotZ();
  void generateThetaVals(int n);
  void plotTheta();
  int getNTracks();
  float theta(float eta);

  float etaMin,etaMax;
  RooRealVar *zVar,*zSig,*thetaVar;
  RooGaussian *zGaus;
  EtaPdf *etaPdf;
  RooDataSet *zData,*thetaData;
  RooPlot *zFrame,*thetaFrame;
  vector<float> zVals,thetaVals;
};
#endif
