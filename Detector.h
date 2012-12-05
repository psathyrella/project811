#ifndef DETECTOR
#define DETECTOR

#include <iostream>
#include <cassert>
#include <sstream>
#include <utility>
#include <math.h>
#include <fstream>

#include "TH1F.h"
#include "TROOT.h"
#include "TMarker.h"
#include "TLine.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH2.h"
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
#include "Simulator.h"

using namespace std;
using namespace RooFit;

class Detector
{
public:
  Detector(TString txtFile);
  void drawAscii(vector<Track> *trks=0);
  void draw(vector<Track> *trks=0);
  pair<float,float> findRPhi(float zVal, Track trk);   // find (r,phi) at the given z for this track
  pair<unsigned,unsigned> findClosestChar(float rVal, float zVal);
  int findClosestPixel(float r, float phi, float z);
  void propagateTrack(Track trk);      // figure out which pixels the track hits
  
  vector<float> rVals,phiVals,zVals;   // r,phi coordinates of the centers of the pixels
  vector<float> zLayers,phiLayers;     // z/phi coords of each layer
  vector<bool> isHit;                  // did this pixel get a hit in this event?
  map<TString,float> minVals,maxVals;  // physical extent of the detector
  int nZ,nR;                           // points for drawing ascii art
  vector<vector<char> > rast;
};  
#endif
