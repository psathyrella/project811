#ifndef LINEFIT
#define LINEFIT
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
// #include "Detector.h"

using namespace ROOT::Math;
using namespace std;

void lineFcn(double t, double *p, double &x, double &y, double &z);
double distance2(double x,double y,double z, double *p);
void SumDistance2(int &, double *, double & sum, double * par, int );

//----------------------------------------------------------------------------------------
class LineFit
{
public:
  LineFit(vector<float> vals);              // list of x,y,z coords
  void fit();
  float getChiSquare(float dR, float dS);
  float closestApproachToZAxis(float zMin=-50, float zMax=400, float dz=10);

  TGraph2D gr;
  TVirtualFitter *min;
  TPolyLine3D *line;
  double parFit[4];     // fitted parameter values
};
#endif
