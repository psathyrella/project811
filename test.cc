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

#include "MitStyleRemix.h"
#include "EtaPdf.h"
#include "Simulator.h"
#include "Detector.h"
#include "LineFit.h"

float unit;
using namespace std;
using namespace RooFit;
//----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  RooRandom::randomGenerator()->SetSeed(getpid());

  TCanvas can("can","",700,900);
  // can.SetPhi(can.GetPhi()+180);
  // can.SetTheta(0);
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.15);
  SetStyle();
  cout << can.GetPhi() << " " << can.GetTheta() << endl;

  Detector trk("tracker.txt");
  // Detector tof("tof.txt");

  // Simulator sim(1.5,5,0,0.5);
  Simulator sim(2.5,5,trk.minVals["phi"],trk.maxVals["phi"]);
  // sim.generate(sim.getNTracks());
  sim.generate(1);
  // sim.readHiFile("output_test.root",5);

  for(unsigned itrk=0; itrk<sim.event->tracks.size(); itrk++) {
    sim.event->tracks[itrk].dump();
    trk.propagateTrack(sim.event->tracks[itrk]);
    // tof.propagateTrack(sim.event->tracks[itrk]);
  }
  trk.fitTrack(trk.chooseHits());


  cout << "drawing" << endl;
  float zMin(0),zMax(320),rMin(0),rMax(20),xMin(0),xMax(20),yMin(0),yMax(20);
  // tof.draw3d(&sim.event->tracks,zMin,zMax,rMin,rMax);
  trk.draw3d(&sim.event->tracks,xMin,xMax,yMin,yMax,zMin,zMax);
  can.SaveAs("/afs/cern.ch/user/d/dkralph/www/foo.png");

}
