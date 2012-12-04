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

float unit;
using namespace std;
using namespace RooFit;
//----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  TCanvas can("can","",900,700);
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.15);
  SetStyle();

  Simulator sim(0,5,0,1);
  sim.generate(sim.getNTracks());

  TString plotDir(".");
  sim.plotZ();
  can.SaveAs(plotDir+"/z.png");

  sim.plotTheta();
  can.SaveAs(plotDir+"/theta.png");
}
