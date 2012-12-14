#include <iostream>
#include <cassert>
#include <sstream>
#include <math.h>
#include <fstream>

#include "TH1F.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
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

float unit;
using namespace std;
using namespace RooFit;
//----------------------------------------------------------------------------------------
class SignalPdf
{
public:
  SignalPdf(TString name, RooRealVar &timePar, float mean, float mean_lo, float mean_hi, float sigma, float sigma_lo, float sigma_hi, bool excludeTop=false);
  // box
  TH1F *hbox;
  RooDataHist *dbox;
  RooHistPdf  *box;
  RooRealVar *mean,*sigma,*time;
  RooGaussian *gaus;
  // box (X) gaussian
  RooFFTConvPdf *sig;
};
//----------------------------------------------------------------------------------------
SignalPdf::SignalPdf(TString name, RooRealVar &timePar, float meanVal, float mean_lo, float mean_hi, float sigmaVal, float sigma_lo, float sigma_hi, bool excludeTop)
{
  time = &timePar;
  // box
  if(excludeTop)
    hbox = new TH1F("hbox_"+name,"",2,-1*unit,1*unit);
  else
    hbox = new TH1F("hbox_"+name,"",3,-1*unit,2*unit);
  hbox->Fill(0.5);
  dbox = new RooDataHist("dbox_"+name, "dbox_"+name, *time, hbox);
  const RooArgSet *args = dbox->get();
  args->Print();
  box  = new RooHistPdf("box_"  +name, "box_"+ name, *time, *dbox);
  // gaussian  
  mean  = new RooRealVar("mean_"+name, "mean_"+name, meanVal,  mean_lo, mean_hi);
  sigma = new RooRealVar("sigma_"+name,"sigma_"+name,sigmaVal, sigma_lo,sigma_hi);
  gaus  = new RooGaussian("gaus_"+name,"gaus_"+name, *time,*mean,*sigma);
  // box (X) gaussian
  sig = new RooFFTConvPdf("signal_"+name,"signal_"+name, *time, *box, *gaus, 2);
}
//----------------------------------------------------------------------------------------
TH1F *getToyHistFromTxt(TString fname)
{
  TH1F *htoys = new TH1F("htoys","",30,-100,100);
  
  ifstream ifs(fname);
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    stringstream ss(line);
    float val;
    if(val==3150) continue; // initial value
    ss >> val;
    htoys->Fill(val);
  }
  return htoys;
}
//----------------------------------------------------------------------------------------
double calcLikelihood(RooDataSet *data, SignalPdf *pdf, double deltaT)
{
  pdf->mean->setVal(deltaT);
  double loVal(1000),hiVal(8000);
  
  double totalLogLik(0);
  const TTree *tree = data->tree();
  for(unsigned ientry=0; ientry<data->numEntries(); ientry++) {
    double val = ((RooRealVar*)(data->get(ientry)->find("time")))->getVal();
    if(val>loVal && val<hiVal) continue;
    pdf->time->setVal(val);
    double logLik_j = log(pdf->sig->getVal());
    totalLogLik += logLik_j;
  }
  ofstream txtOutFile("lik-range.txt",ios_base::app);
  txtOutFile << setw(12) << deltaT << setw(12) << totalLogLik << endl;
  txtOutFile.close();
  return 1;
}
//----------------------------------------------------------------------------------------
void makeGraph(TString fname)
{
  vector<double> mHs,limits;
  ifstream ifs(fname);
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    stringstream ss(line);
    double mH;
    double limit;
    ss >> mH >> limit;
    mHs.push_back(mH);
    limits.push_back(limit);
  }

  Double_t xArr[mHs.size()];
  Double_t yArr[limits.size()];
  for(unsigned i=0; i<mHs.size(); i++) {
    xArr[i] = mHs[i];
    yArr[i] = limits[i];
  }
  
  TGraphAsymmErrors txtGr(mHs.size(),xArr,yArr);
  TCanvas can;
  txtGr.SetMarkerStyle(22);
  txtGr.SetTitle(";#delta t;log likelihood");
  TF1 poly("poly","pol2",-.9,.9);
  txtGr.Fit("poly");
  txtGr.Draw("ap");
  can.SaveAs("$HOME/www/lik-range.png");
}
//----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  assert(argc>1);
  bool generateToys = atoi(argv[1]);
  
  TCanvas can;
  unit = 10500;

  if(generateToys) {
    RooRealVar time("time","time",-5*unit,5*unit);
    time.setRange("fitrange",-8000,8000);
    // time.setBins(1000,"cache");   // nbins for fft (wtf? does something weird)
    time.setBins(1000);

    float mean_lo(time.getMin()),mean_hi(time.getMax());
    SignalPdf signal("signal", time, 0*unit, mean_lo, mean_hi, .06*unit, 0*unit, 1*unit);
    
    RooRandom::randomGenerator()->SetSeed(getpid());
    RooDataSet *dsignal;
    for(int i=0; i<1; i++) {
      dsignal = signal.sig->generate(time,15200);
      system("rm -fv lik-range.txt");
      for(double deltaT=-75; deltaT<76; deltaT+=1)
	calcLikelihood(dsignal, &signal, deltaT);
      assert(0);
      SignalPdf signalForFit("signalForFit", time, 0*unit, mean_lo, mean_hi, 0.08*unit, 0*unit, 1*unit);
      RooFitResult *fitres = signalForFit.sig->fitTo(*dsignal,Range("fitrange"));
      cout << "final: " << signalForFit.mean->getVal() << endl;

      RooPlot *frame = time.frame(Title("."),Range("fitrange"));
      dsignal->plotOn(frame,Binning(100));
      signalForFit.sig->plotOn(frame);
      frame->Draw();
      can.SaveAs("$HOME/www/foo.png");
    }
  
  }

  if(!generateToys) {
    makeGraph("lik-range.txt");
    
    // TH1F *htoys = getToyHistFromTxt("toys.txt");
    // htoys->SetTitle(";time [ns];N Toys");
    // htoys->Draw();
    // can.SaveAs("$HOME/www/toys.png");
  }
}
