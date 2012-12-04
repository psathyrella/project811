#ifndef UTIL
#define UTIL

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

using namespace std;
using namespace RooFit;

#define e 2.71828182846

//----------------------------------------------------------------------------------------
//========================================================================================
class Track
{
public:
  TLorentzVector vec;
  float dz;
};
//========================================================================================
class Event
{
public:
  vector<Track> tracks;
};
#endif
