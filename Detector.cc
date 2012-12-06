#include "Detector.h"

Detector::Detector(TString txtFile)
{
  ifstream ifs(txtFile);
  assert(ifs.is_open());
  string line;
  // note: these deltas are just approximate
  float dR,dS,dZ; // dS = r*dPhi
  vector<float> tmpRVals,tmpZVals;
  while(getline(ifs,line)) {
    if(line[0]=='#' || line[0]==' ' || line.size()==0) continue;

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
      ss >> tmp >> dR >> dS >> dZ;
    } else if(line[0]=='z') {
      // read in explicit layers
      TString var;
      ss >> var;
      float zVal;
      while(!ss.eof()) {
	ss >> zVal;
	tmpZVals.push_back(zVal);
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

  int nR = int((maxVals["r"] - minVals["r"]) / dR);
  float realDR((maxVals["r"] - minVals["r"]) / nR);
  // cout << "pushing rs: " << endl;
  for(unsigned iR=0; iR<nR; iR++) {
    // cout << "   " << minVals["r"] + 0.5*realDR + iR*realDR << endl;
    tmpRVals.push_back(minVals["r"] + 0.5*realDR + iR*realDR);
  }
  if(dZ != -1) {
    int nZ = int((maxVals["z"] - minVals["z"]) / dZ);
    float realDZ((maxVals["z"] - minVals["z"]) / nZ);
    // cout << "pushing zs: " << endl;
    for(unsigned iZ=0; iZ<nZ; iZ++) {
      // cout << "   " << minVals["z"] + 0.5*realDZ + iZ*realDZ << endl;
      tmpZVals.push_back(minVals["z"] + 0.5*realDZ + iZ*realDZ);
    }
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
	if(// iZ==0 || iZ==tmpZVals.size()-1 || 
	   iR==0 || iR==tmpRVals.size()-1 || iPhi==0 || iPhi==nPhi-1)
	  // if(iPhi==0 || iPhi==nPhi-1)
	  isEdge.push_back(true);
	else
	  isEdge.push_back(false);
	markerStyles.push_back( (iPhi==0) ? 20 : 24); // closed : open circles
	markerSizes.push_back(1 + iPhi);
	if(zLayers.size()==0 || tmpZVals[iZ] != zLayers.back())
	  zLayers.push_back(tmpZVals[iZ]);
	if(phiLayers.size()==0 || phiVal != phiLayers.back())
	  phiLayers.push_back(phiVal);
      }
    }
  }
	
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
void Detector::draw3d(vector<Track> *trks, float zMinDraw, float zMaxDraw, float rMinDraw, float rMaxDraw, TString option)
{
  if(zMinDraw==zMaxDraw) {
    zMinDraw = minVals["z"];
    zMaxDraw = maxVals["z"];
  }
  if(rMinDraw==rMaxDraw) {
    rMinDraw = minVals["r"];
    rMaxDraw = maxVals["r"];
  }
  
  // TCanvas can("can","can",2000,1000);
  if(!option.Contains("same",TString::kIgnoreCase)) {
    TH3F *hist = new TH3F("hist",";z [cm];x [cm];y [cm]",100,zMinDraw,zMaxDraw,100,rMinDraw*cos(minVals["phi"]),rMaxDraw*cos(maxVals["phi"]),100,rMinDraw*sin(minVals["phi"]),rMaxDraw*sin(maxVals["phi"]));
    hist->SetDirectory(0);
    hist->Draw();
  }
  for(unsigned iPix=0; iPix<rVals.size(); iPix++) {
    if(!isEdge[iPix] && !isHit[iPix]) continue;
    float *coords = new float[3];
    coords[0] = zVals[iPix];
    // coords[1] = phiVals[iPix];
    // coords[2] = rVals[iPix];
    coords[1] = rVals[iPix]*cos(phiVals[iPix]);
    coords[2] = rVals[iPix]*sin(phiVals[iPix]);
    TPolyMarker3D *mark = new TPolyMarker3D(3,coords,20);
    if(isHit[iPix])
      mark->SetMarkerColor(kRed);
    mark->Draw();
  }
  if(trks) {
    for(unsigned itrk=0; itrk<trks->size(); itrk++) {
      // pair<float,float> loPair(findRPhi(zMinDraw, (*trks)[itrk]));
      pair<float,float> loPair(findRPhi((*trks)[itrk].dz, (*trks)[itrk]));
      pair<float,float> hiPair(findRPhi(zMaxDraw, (*trks)[itrk]));
      float *coords = new float[3];
      // coords[0] = zMinDraw;
      coords[0] = (*trks)[itrk].dz;
      coords[1] = loPair.first*cos(loPair.second);
      coords[2] = loPair.first*sin(loPair.second);
      coords[3] = zMaxDraw;
      coords[4] = hiPair.first*cos(hiPair.second);
      coords[5] = hiPair.first*sin(hiPair.second);
      // cout
      // 	<< "drawing line with: "
      // 	<< setw(12) << coords[0]
      // 	<< setw(12) << coords[1]
      // 	<< setw(12) << coords[2]
      // 	<< setw(12) << coords[3]
      // 	<< setw(12) << coords[4]
      // 	<< setw(12) << coords[5]
      // 	<< endl;
      TPolyLine3D *line = new TPolyLine3D(2,coords);
      line->SetLineColor(kBlue);
      line->Draw();
    }
  }
  // cout << "saving..." << endl;
  // can.SaveAs("/afs/cern.ch/user/d/dkralph/www/foo.png");
  // cout << "saved..." << endl;
}
  
//----------------------------------------------------------------------------------------
void Detector::draw(vector<Track> *trks)
{
  TCanvas can("can","can",900,600);
  TH2F hist("hist",";z [cm];r [cm]",100,minVals["z"],maxVals["z"],100,minVals["r"],maxVals["r"]);
  hist.Draw("colz");

  // cout << "mark with: " << endl;
  for(unsigned iPix=0; iPix<rVals.size(); iPix++) {
    // cout 
    //   << setw(12) << rVals[iPix]
    //   << setw(12) << phiVals[iPix]
    //   << setw(12) << zVals[iPix]
    //   << setw(12) << markerStyles[iPix]
    //   << setw(12) << markerSizes[iPix]
    //   << endl;
    TMarker *mark = new TMarker(zVals[iPix],rVals[iPix],markerStyles[iPix]);
    mark->SetMarkerSize(markerSizes[iPix]);
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

  can.SaveAs("/afs/cern.ch/user/d/dkralph/www/foo.png");
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
  // cout << "finding r,phi for z: " << zVal << " and ";
  // trk.dump();
  // assert(zVal > trk.dz);
  coords.first = (zVal - trk.dz)*tan(trk.vec.Theta());
  coords.second = trk.vec.Phi();
  // cout << "coords: " << coords.first << " " << coords.second << endl;
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
		    
      
