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

float unit;
using namespace std;
using namespace RooFit;
//----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  RooRandom::randomGenerator()->SetSeed(getpid());

  TCanvas can("can","",900,700);
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.15);
  SetStyle();

  Simulator sim(2.5,5,0,.5);
  // sim.generate(sim.getNTracks());
  sim.generate(20);
  // sim.readHiFile("output_test.root",5);

  Detector trk("tracker.txt");
  Detector tof("tof.txt");
  for(unsigned itrk=0; itrk<sim.event->tracks.size(); itrk++) {
    trk.propagateTrack(sim.event->tracks[itrk]);
    tof.propagateTrack(sim.event->tracks[itrk]);
  }

  cout << "drawing" << endl;
  // TH3F hist("hist",";z [cm];x [cm];y [cm]",100,0,200,100,0*cos(0),20*cos(0.5),100,0*sin(0),20*sin(0.5));
  // hist.Draw();
  trk.draw3d(&sim.event->tracks,0,200,0,20);
  tof.draw3d(&sim.event->tracks,0,200,0,20,"same");
  can.SaveAs("/afs/cern.ch/user/d/dkralph/www/foo.png");

}
