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

#include "/afs/cern.ch/user/d/dkralph/cms/cmssw/029/CMSSW_5_3_2_patch4/src/MitHzz4l/Util/interface/MitStyleRemix.h"

float unit;
using namespace std;
using namespace RooFit;
//----------------------------------------------------------------------------------------
vector<float> getZvals()
// generate a set of gaussian-distributed z values, put them in a vector, sort, return
{
  RooRealVar zvar("zvar","zvar",-20,20);
  RooRealVar sig("sig","sig",5);

  RooGaussian gaus("gaus","gaus",zvar,RooConst(0),sig);
  RooDataSet *data = gaus.generate(zvar,200);

  vector<float> vals;
  for(int ientry=0; ientry<data->numEntries(); ientry++) {
    float val = ((RooRealVar*)(data->get(ientry)->find("zvar")))->getVal();
    vals.push_back(val);
  }
  sort(vals.begin(), vals.end());
  return vals;
}
//----------------------------------------------------------------------------------------
void fillCloseHist(TH1F &hist, vector<float> vals)
// fill histogram with the distance to the nearest neighbor for each of the z values
{
  for(unsigned ival=0; ival<vals.size(); ival++) {
    float val(vals[ival]),closest;
    if(ival==0)
      closest = vals[ival+1];
    else if(ival==vals.size()-1)
      closest = vals[ival-1];
    else
      closest = min(fabs(vals[ival-1] - val), fabs(vals[ival+1] - val));
    hist.Fill(closest);
  }
}
//----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  TCanvas can("can","",900,700);
  cout << "margin: " << gPad->GetLeftMargin() << endl;
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.15);
  cout << "margin: " << gPad->GetLeftMargin() << endl;
  SetStyle();

  gStyle->SetOptStat("m");
  TH1F hist("hist",";d_{near};entries",100,0,.2);

  for(int itoy=0; itoy<100; itoy++)
    fillCloseHist(hist, getZvals());

  hist.Draw("hist");
  // RooPlot *frame = zvar.frame(Title("X"));
  // data->plotOn(frame);
  // gaus.plotOn(frame);
  // frame->Draw();
  can.SaveAs("$HOME/www/foo.png");
}
