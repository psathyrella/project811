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
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TH3F.h"
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
#include "LineFit.h"

using namespace std;
using namespace RooFit;

class Detector
//
// note: because of the cylindrical -> cartesian transition, z is the first coordinate when
// plotting or interfacing to LineFit, i.e. it goes (z,y,x)
//
{
public:
  Detector(TString txtFile);
  void drawAscii(vector<Track> *trks=0);
  void draw(vector<Track> *trks=0);
  // void draw3d(vector<Track> *trks=0, float zMinDraw=-1, float zMaxDraw=-1, float rMinDraw=-1, float rMaxDraw=-1, TString option="");
  void draw3d(vector<Track> *trks=0, float zMinDraw=-1, float zMaxDraw=-1, float xMinDraw=-1, float xMaxDraw=-1, float yMinDraw=-1, float yMaxDraw=-1, TString option="");
  pair<float,float> findRPhi(float zVal, Track trk);   // find (r,phi) at the given z for this track
  pair<unsigned,unsigned> findClosestChar(float rVal, float zVal);
  int findClosestPixel(float r, float phi, float z);
  void propagateTrack(Track trk);      // figure out which pixels the track hits
  vector<float> chooseHits();          // pick a set of hits to pass to the fitter
  void fitTrack(vector<float> hits);
  
  vector<float> rVals,phiVals,zVals;   // r,phi coordinates of the centers of the pixels
  vector<short> markerStyles,markerSizes;
  vector<float> zLayers,phiLayers;     // z/phi coords of each layer
  vector<bool> isHit,isEdge;                  // did this pixel get a hit in this event?
  map<TString,float> minVals,maxVals;  // physical extent of the detector
  float dR,dS,dZ;                      // pixel-spacings. dR and dS = r*dPhi are the resolutions
  int nZ,nR;                           // points for drawing ascii art
  vector<vector<char> > rast;
  int bkgColor,hitColor;  // for plotting
  vector<LineFit*> lines; // collection of fitted lines
};  
#endif
