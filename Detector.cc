#include "Detector.h"

Detector::Detector(TString txtFile)
{
  ifstream ifs(txtFile);
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;

    stringstream ss(line);
    if(line[0]=='$') {
      TString dollarSign,var;
      float min,max;
      ss >> dollarSign >> var >> min >> max;
      minVals[var] = min;
      maxVals[var] = max;
    } else if(line[0]=='d') {
      // make evenly-spaced pixels
      TString tmp;
      // note: these deltas are just approximate
      float dR,dS,dZ; // dS = r*dPhi
      ss >> tmp >> dR >> dS >> dZ;

      int nR = int((maxVals["r"] - minVals["r"]) / dR);
      float realDR((maxVals["r"] - minVals["r"]) / nR);
      vector<float> tmpRVals;
      // cout << "pushing rs: " << endl;
      for(unsigned iR=0; iR<nR; iR++) {
	// cout << "   " << minVals["r"] + 0.5*realDR + iR*realDR << endl;
	tmpRVals.push_back(minVals["r"] + 0.5*realDR + iR*realDR);
      }
      int nZ = int((maxVals["z"] - minVals["z"]) / dZ);
      float realDZ((maxVals["z"] - minVals["z"]) / nZ);
      vector<float> tmpZVals;
      // cout << "pushing zs: " << endl;
      for(unsigned iZ=0; iZ<nZ; iZ++) {
	// cout << "   " << minVals["z"] + 0.5*realDZ + iZ*realDZ << endl;
	tmpZVals.push_back(minVals["z"] + 0.5*realDZ + iZ*realDZ);
      }

      for(unsigned iZ=0; iZ<tmpZVals.size(); iZ++) {
	for(unsigned iR=0; iR<tmpRVals.size(); iR++) {
	  int nPhi = int(tmpRVals[iR]*(maxVals["phi"] - minVals["phi"]) / dS);
	  nPhi = max(1,nPhi);
	  float realDPhi((maxVals["phi"] - minVals["phi"]) / nPhi);
	  // cout << "nPhi: " << nPhi << "realDPhi: " << setw(12) << realDPhi << "realDS: " << setw(12) << tmpRVals[iR]*realDPhi << endl;
	  // cout << "pushing phis: " << endl;
	  for(unsigned iPhi=0; iPhi<nPhi; iPhi++) {
	    float phiVal(minVals["phi"] + 0.5*realDPhi + iPhi*realDPhi);
	    // cout << "   " << phiVal << endl;
	    rVals.push_back(tmpRVals[iR]);
	    zVals.push_back(tmpZVals[iZ]);
	    phiVals.push_back(phiVal);
	    isHit.push_back(false);
	    if(zLayers.size()==0 || tmpZVals[iZ] != zLayers.back())
	      zLayers.push_back(tmpZVals[iZ]);
	    if(phiLayers.size()==0 || phiVal != phiLayers.back())
	      phiLayers.push_back(phiVal);
	  }
	}
      }
	
      
    } else {
      // read each pixel from text file
      float r,phi,z;
      ss >> r >> phi >> z;
      rVals.push_back(r);
      phiVals.push_back(phi);
      zVals.push_back(z);
      if(zLayers.size()==0 || z != zLayers.back())
	zLayers.push_back(z);
      if(phiLayers.size()==0 || phi != phiLayers.back())
	phiLayers.push_back(phi);
      isHit.push_back(false);
    }
  }

  assert(minVals.find("r") != minVals.end());
  assert(minVals.find("phi") != minVals.end());
  assert(minVals.find("z") != minVals.end());

  nR=55;
  nZ=160;
}
//----------------------------------------------------------------------------------------
int Detector::findClosestPixel(float r, float phi, float z)
{
  if(r < minVals["r"] || r > maxVals["r"])		return -1;
  if(phi < minVals["phi"] || phi > maxVals["phi"])	return -1;
  if(z < minVals["z"] || z > maxVals["z"])		return -1;
  float close(999999);
  int iClose(-1);
  for(unsigned iPix=0; iPix<rVals.size(); iPix++) {
    float dR   = r	- rVals[iPix];
    float dPhi = phi	- phiVals[iPix];
    float dZ   = z	- zVals[iPix];
    float delta(sqrt(dR*dR + dPhi*dPhi + dZ*dZ));
    if(delta < close) {
      iClose = iPix;
      close = delta;
    }
  }
  return iClose;
}
//----------------------------------------------------------------------------------------
pair<unsigned,unsigned> Detector::findClosestChar(float rVal, float zVal)
{
  unsigned iRclose(0),iZclose(0);
  float close(9999999);
  // cout << "pix: " << rVal << setw(12) << zVal << endl;
  for(unsigned iR=0; iR<nR; iR++) {
    for(unsigned iZ=0; iZ<nZ; iZ++) {
      float rPos = minVals["r"] + (maxVals["r"] - minVals["r"]) * iR / nR;
      float zPos = minVals["z"] + (maxVals["z"] - minVals["z"]) * iZ / nZ;
      float dR = rPos - rVal;
      float dZ = zPos - zVal;

      if(sqrt(dR*dR + dZ*dZ) < close) {
	iRclose = iR;
	iZclose = iZ;
	close = sqrt(dR*dR + dZ*dZ);
	// cout << "  setting with: " << setw(12) << iRclose << setw(12) << iZclose << setw(12) << rPos << setw(12) << zPos << "( " << close << ")" << endl;
      }
    }
  }

  return pair<unsigned,unsigned> (iRclose,iZclose);
}
//----------------------------------------------------------------------------------------
void Detector::draw(vector<Track> *trks)
{
  TCanvas can("can","can",900,600);
  TH2F hist("hist",";z [cm];r [cm]",100,minVals["z"],maxVals["z"],100,minVals["r"],maxVals["r"]);
  hist.Draw("colz");

  vector<int> markerStyles;
  vector<int> markerSizes;
  markerStyles.push_back(20);   markerSizes.push_back(1);
  markerStyles.push_back(24);   markerSizes.push_back(2.5);
  markerStyles.push_back(23);   markerSizes.push_back(1);
  markerStyles.push_back(23);   markerSizes.push_back(1);
  markerStyles.push_back(23);   markerSizes.push_back(1);
  markerStyles.push_back(23);   markerSizes.push_back(1);
  markerStyles.push_back(23);   markerSizes.push_back(1);
  int iMarkStyle(0);
  for(unsigned iPix=0; iPix<rVals.size(); iPix++) {
    cout << "mark with: " << setw(12) << rVals[iPix] << setw(12) << zVals[iPix] << endl;
    if(iPix>0 && phiVals[iPix] != phiVals[iPix-1]) iMarkStyle++;
    TMarker *mark = new TMarker(zVals[iPix],rVals[iPix],markerStyles[iMarkStyle]);
    mark->SetMarkerSize(markerSizes[iMarkStyle]);
    if(isHit[iPix])
      mark->SetMarkerColor(kRed);
    mark->Draw();
  }

  if(trks) {
    for(unsigned itrk=0; itrk<trks->size(); itrk++) {
      pair<float,float> loCoords(findRPhi(minVals["z"], (*trks)[itrk]));
      pair<float,float> hiCoords(findRPhi(maxVals["z"], (*trks)[itrk]));
      TLine *line = new TLine(minVals["z"], loCoords.first, maxVals["z"], hiCoords.first);
      line->Draw();
    }
  }

  // can.SaveAs("/afs/cern.ch/user/d/dkralph/www/foo.png");
}  
//----------------------------------------------------------------------------------------
void Detector::drawAscii(vector<Track> *trks)
{
  for(unsigned iR=0; iR<nR; iR++) {
    vector<char> tmpVec;
    rast.push_back(tmpVec);
    for(unsigned iZ=0; iZ<nZ; iZ++) {
      rast.back().push_back(' ');
    }
  }

  for(unsigned iPix=0; iPix<rVals.size(); iPix++) {
    pair<unsigned,unsigned> coords(findClosestChar(rVals[iPix], zVals[iPix]));

    assert(rast.size() >= coords.first);
    assert(rast[coords.first].size() >= coords.second);
    if(isHit[iPix])
      rast[coords.first][coords.second] = 'o';
    else
      rast[coords.first][coords.second] = 'x';
  }

  if(trks) {
    for(unsigned itrk=0; itrk<trks->size(); itrk++) {
      for(unsigned iZ=0; iZ<nZ; iZ++) {
	float zPos = minVals["z"] + (maxVals["z"] - minVals["z"]) * iZ / nZ;
	pair<float,float> trkCoords(findRPhi(zPos, (*trks)[itrk]));
	pair<unsigned,unsigned> charCoords(findClosestChar(trkCoords.first, zPos));
	assert(rast.size() >= charCoords.first);
	assert(rast[charCoords.first].size() >= charCoords.second);
	if(rast[charCoords.first][charCoords.second] == ' ')
	  rast[charCoords.first][charCoords.second] = '.';
      }
    }
  }

  cout << "z = ";
  for(unsigned iZ=0; iZ<nZ; iZ++) {
    float zPos = minVals["z"] + (maxVals["z"] - minVals["z"]) * iZ / nZ;
    if((iZ%10) == 0)
      cout << setw(10) << setprecision(3) << zPos;
  }
  cout << endl << endl;

  for(unsigned iR=0; iR<nR; iR++) {
    for(unsigned iZ=0; iZ<nZ; iZ++) {
      float rPos = minVals["r"] + (maxVals["r"] - minVals["r"]) * iR / nR;
      if(iZ==0)
	cout << "r = " << setw(5) << setprecision(3) << rPos << " ";
      cout << rast[iR][iZ];
    }
    cout << endl;
  }
  
}
//----------------------------------------------------------------------------------------
pair<float,float> Detector::findRPhi(float zVal, Track trk)
{
  pair<float,float> coords;
  assert(zVal > trk.dz);
  coords.first = (zVal - trk.dz)*tan(trk.vec.Theta());
  coords.second = trk.vec.Phi();
  return coords;
}
//----------------------------------------------------------------------------------------
void Detector::propagateTrack(Track trk)
{
  for(unsigned ilayer=0; ilayer<zLayers.size(); ilayer++) {
    // cout << "layer: " << ilayer << endl;
    pair<float,float> trkCoords(findRPhi(zLayers[ilayer],trk));
    // cout << "  looking for trk with " << setw(12) << trkCoords.first << setw(12) << trkCoords.second << endl;
    int iHit(findClosestPixel(trkCoords.first, trkCoords.second, zLayers[ilayer]));
    if(iHit<0)
      ;// cout << "     not in detector" << endl;
    else {
      // cout << "  using pixel: " << setw(12) << iHit << setw(12) << rVals[iHit] << setw(12) << phiVals[iHit] << setw(12) << zVals[iHit] << endl;
      isHit[iHit] = true;
    }
  }
}
		    
      
