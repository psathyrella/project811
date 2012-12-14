#include "LineFit.h"

//----------------------------------------------------------------------------------------
void lineFcn(double t, double *p, double &x, double &y, double &z)
// set (x,y,z) given a z (=t) value and the parameters p[0-3]
{ 
  // a parameteric line is define from 6 parameters but 4 are independent
  // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
  // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
  x = p[0] + p[1]*t; 
  y = p[2] + p[3]*t;
  z = t; 
} 
//----------------------------------------------------------------------------------------
double distance2(double x,double y,double z, double *p)
{ 
  // calculate distance line-point 
  // distance line point is D= | (xp-x0) cross  ux | 
  // where ux is direction of line and x0 is a point in the line (like t = 0) 
  XYZVector xp(x,y,z); 
  XYZVector x0(p[0], p[2], 0. ); 
  XYZVector x1(p[0] + p[1], p[2] + p[3], 1. ); 
  XYZVector u = (x1-x0).Unit(); 
  double d2 = ((xp-x0).Cross(u)) .Mag2(); 
  return d2; 
}
//----------------------------------------------------------------------------------------
void SumDistance2(int &, double *, double & sum, double * par, int )
{ 
  // function to be minimized 
  TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
  assert(gr != 0);
  double * x = gr->GetX();
  double * y = gr->GetY();
  double * z = gr->GetZ();
  int npoints = gr->GetN();
  sum = 0;
  for(int i =0; i<npoints; i++) { 
    double d = distance2(x[i],y[i],z[i],par); 
    sum += d;
  }
}
//----------------------------------------------------------------------------------------
LineFit::LineFit(vector<float> vals) // x1,y1,z1,x2,y2,z2,...
{
  parFit[0] = parFit[1] = parFit[2] = parFit[3] = 0;
  for(unsigned iVal=0; iVal<vals.size()-2; iVal+=3) {
    gr.SetPoint(iVal/3, vals[iVal], vals[iVal+1], vals[iVal+2]);
  }
}
//----------------------------------------------------------------------------------------
void LineFit::fit()
{
  min = TVirtualFitter::Fitter(0,4);
  min->SetObjectFit(&gr);
  void (*fcn)(int&, double*, double&, double*, int);
  fcn = &SumDistance2;
  min->SetFCN(fcn);

  Double_t arglist[10];
  arglist[0] = 3;
  min->ExecuteCommand("SET PRINT",arglist,-1);
  
  double pStart[4] = {1,1,1,1};
  min->SetParameter(0,"x0",pStart[0],0.01,0,0);
  min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
  min->SetParameter(2,"y0",pStart[2],0.01,0,0);
  min->SetParameter(3,"Ay",pStart[3],0.01,0,0);
    
  arglist[0] = 1000; // number of function calls 
  arglist[1] = 0.001; // tolerance 
  min->ExecuteCommand("MIGRAD",arglist,2);
  // if (minos) min->ExecuteCommand("MINOS",arglist,0);

  int nvpar,nparx; 
  double amin,edm, errdef;
  min->GetStats(amin,edm,errdef,nvpar,nparx);
  // min->PrintResults(1,amin);

  // get fit parameters
  for(int i=0; i<4; i++) 
    parFit[i] = min->GetParameter(i);
   
  // cout << "minimized fcn: " << amin << endl;

  // draw the fitted line
  int nPts = 10;
  double t0 = 0;
  double dt = 1.9*(gr.GetZ())[gr.GetN()-1];
  line = new TPolyLine3D(nPts);
  for(int iPt=0; iPt<nPts; iPt++) {
    double t = t0 + dt*iPt/nPts;
    double x,y,z;
    lineFcn(t,parFit,x,y,z);
    line->SetPoint(iPt,x,y,z);
  }
}
//----------------------------------------------------------------------------------------
float LineFit::getChiSquare(float dR, float dS)
{
  //
  // NOTE: I don't think this is actually right... but it at least scales with the quality of the fit
  // in the proper way, and I just need to know which of the possible lines has the best fit, so... fuck it, for the moment.
  //
  float sum(0);
  for(unsigned ih=0; ih<gr.GetN(); ih++) {
    // float diff2 = distance2((gr.GetX())[ih], (gr.GetY())[ih], (gr.GetZ())[ih], parFit);
    // cout << ih << setw(12) << sqrt(diff2) << endl;
    // sum += diff2 / (dR*dR + dS*dS);
    double xHit((gr.GetX())[ih]),yHit((gr.GetY())[ih]),zHit((gr.GetZ())[ih]); // coords of the hit
    double xLine,yLine,zLine;
    lineFcn(zHit, parFit, xLine, yLine, zLine); // coords of the line at the z value of the hit
    double dX2 = (xHit - xLine)*(xHit - xLine);
    double dY2 = (yHit - yLine)*(yHit - yLine);
    double delta = sqrt(dR*dR + dS*dS);
    double dZ2 = delta;  // not really... but ok for now
    // cout << "hits: "
    // 	 << setw(12) << xHit
    // 	 << setw(12) << xLine
    // 	 << setw(12) << yHit
    // 	 << setw(12) << yLine
    // 	 << setw(12) << zHit
    // 	 << setw(12) << zLine
    // 	 << "   " << ih << setw(12) << sqrt(dX2) << setw(12) << sqrt(dY2) << setw(12) << sqrt(dZ2) << " delta: " << delta << endl;

    sum += dX2 / (delta*delta);
    sum += dY2 / (delta*delta);
    sum += dZ2 / (delta*delta);
  }
  int ndof(3*gr.GetN() - 4 - 1);
  // cout << "ndof: " << ndof << endl;
  return sqrt(sum) / ndof;
}
//----------------------------------------------------------------------------------------
float LineFit::closestApproachToZAxis(float zMin, float zMax, float dz)
// return the z at which the line passes closest to the z axis
{
  float nDzStart(10);         // divide the interval into this many pieces
  float lastDz(zMax-zMin);    // the last step size that was tried
  int nTries(0),nTriesMax(3);
  float minDistance(99999999999999);
  float bestZ(zMin);
  float zMinTry(zMin),zMaxTry(zMax);
  while(nTries<nTriesMax && lastDz>dz) { // give up when the step size gets smaller than dz or we've already tried nTriesMax times
    float tmpMinDistance(minDistance);
    float tmpBestZ(bestZ);
    for(float zTry=zMinTry; zTry<zMaxTry+(zMaxTry-zMinTry)/nDzStart; zTry += (zMaxTry-zMinTry)/nDzStart) {
      float distance(sqrt(distance2(0, 0, zTry, parFit))); // minimum distance from the z-axis at zTry to the line
      // cout << "closestApproachToZAxis   trying: " << zTry << setw(12) << distance << endl;
      if(distance < tmpMinDistance) {
	tmpMinDistance = distance;
	tmpBestZ = zTry;
      }
    }
    // cout << "closestApproachToZAxis   best for this loop: " << tmpMinDistance << setw(12) << tmpBestZ << endl;
    if(tmpMinDistance <= minDistance) {
      minDistance = tmpMinDistance;
      bestZ = tmpBestZ;
      zMinTry = bestZ - 1*(zMaxTry-zMinTry)/nDzStart; // set new limits to +/- two intervals
      zMaxTry = bestZ + 1*(zMaxTry-zMinTry)/nDzStart;
      // cout << "closestApproachToZAxis   setting best" << setw(12) << zMinTry << setw(12) << zMaxTry << endl;
    }

    nTries++;
    lastDz = (zMax-zMin)/nDzStart;
  }
  // cout << "closestApproachToZAxis returning: " << setw(12) << bestZ << setw(12) << minDistance << endl;
  return bestZ;
}
    
