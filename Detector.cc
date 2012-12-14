#include "Detector.h"

Detector::Detector(TString txtFile) : bkgColor(1),hitColor(632),dR(-1),dS(-1),dZ(-1)
{
  ifstream ifs(txtFile);
  assert(ifs.is_open());
  string line;
  // note: these deltas are just approximate
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
    } else if(line[0]=='s') {
      // set detector's plotting style
      TString tmp;
      ss >> tmp >> bkgColor >> hitColor;
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
  if(dR==-1 && rVals.size()>1)       dR = fabs(rVals[1]   - rVals[0]);
  if(dS==-1 && phiVals.size()>1)     dS = fabs(rVals[1]*phiVals[1] - rVals[0]*phiVals[0]);
  if(dZ==-1 && zVals.size()>1)       dZ = fabs(zVals[1]   - zVals[0]);

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
	// cout << "pushing: " << setw(12) << rVals.back() << setw(12) << zVals.back() << setw(12) << phiVals.back() << endl;
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
void Detector::draw3d(vector<Track> *trks, float xMinDraw, float xMaxDraw, float yMinDraw, float yMaxDraw, float zMinDraw, float zMaxDraw, TString option)
{
  if(zMinDraw==zMaxDraw) {
    zMinDraw = minVals["z"];
    zMaxDraw = maxVals["z"];
  }
  // if(rMinDraw==rMaxDraw) {
  //   rMinDraw = minVals["r"];
  //   rMaxDraw = maxVals["r"];
  // }
  
  if(!option.Contains("same",TString::kIgnoreCase)) {
    // float xMin = min(rMinDraw*cos(minVals["phi"]),rMinDraw*cos(maxVals["phi"]));
    // float xMax = max(rMaxDraw*cos(minVals["phi"]),rMaxDraw*cos(maxVals["phi"]));
    // float yMin = min(rMinDraw*sin(minVals["phi"]),rMinDraw*sin(maxVals["phi"]));
    // float yMax = max(rMaxDraw*sin(minVals["phi"]),rMaxDraw*sin(maxVals["phi"]));
    TH3F *hist = new TH3F("hist",";x [cm];y [cm];z [cm]",100,xMinDraw,xMaxDraw,100,yMinDraw,yMaxDraw,100,zMinDraw,zMaxDraw);
    hist->SetDirectory(0);
    hist->Draw();
  }
  for(unsigned iPix=0; iPix<rVals.size(); iPix++) {
    if(!isEdge[iPix] && !isHit[iPix]) continue;
    float *coords = new float[3];
    coords[0] = rVals[iPix]*cos(phiVals[iPix]);
    // cout << setw(12) << rVals[iPix] << setw(12) << phiVals[iPix] << setw(12) << rVals[iPix]*cos(phiVals[iPix])
    // 	 << setw(12) << rVals[iPix]*sin(phiVals[iPix]) << endl;
    coords[1] = rVals[iPix]*sin(phiVals[iPix]);
    coords[2] = zVals[iPix];
    TPolyMarker3D *mark = new TPolyMarker3D(3,coords,20);
    if(isHit[iPix])
      mark->SetMarkerColor(hitColor);
    else
      mark->SetMarkerColor(bkgColor);
    mark->Draw();
  }
  if(trks) {
    for(unsigned itrk=0; itrk<trks->size(); itrk++) {
      pair<float,float> loPair(findRPhi((*trks)[itrk].dz, (*trks)[itrk]));
      pair<float,float> hiPair(findRPhi(zMaxDraw, (*trks)[itrk]));
      float *coords = new float[3];
      coords[0] = loPair.first*cos(loPair.second);
      coords[1] = loPair.first*sin(loPair.second);
      coords[2] = (*trks)[itrk].dz;
      coords[3] = hiPair.first*cos(hiPair.second);
      coords[4] = hiPair.first*sin(hiPair.second);
      coords[5] = zMaxDraw;
      TPolyLine3D *line = new TPolyLine3D(2,coords);
      line->SetLineColor(kBlue);
      line->SetLineWidth(2);
      line->Draw();
    }
  }
  for(unsigned iline=0; iline<lines.size(); iline++) {
    lines[iline]->line->SetLineColor(kRed+1);
    lines[iline]->line->SetLineWidth(2);
    lines[iline]->line->Draw("same");
  }
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
// find r,phi coordinates of a track for a given z (note: may not be correct for negative z)
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
// figure out which pixels in each z layer the track hits
{
  for(unsigned ilayer=0; ilayer<zLayers.size(); ilayer++) {
    pair<float,float> trkCoords(findRPhi(zLayers[ilayer],trk));
    int iHit(findClosestPixel(trkCoords.first, trkCoords.second, zLayers[ilayer]));
    if(iHit<0)
      ;
    else {
      isHit[iHit] = true;
      hitIndices.push_back(iHit);
    }
  }
}
//----------------------------------------------------------------------------------------
vector<int> Detector::chooseHits(int iStartLayer, bool useBeamSpotConstraint)
{
  vector<LineFit> lfs;
  vector<int> hit1s,hit2s,hit3s; // hits in the first, second, and third layers for each LineFit
  // loop over all the hits in iStartLayer
  for(unsigned ihit=0; ihit<hitIndices.size(); ihit++) {
    int hit1(hitIndices[ihit]);
    if(zVals[hit1] != zLayers[iStartLayer]) continue;
    // loop over all the hits in the next layer
    for(unsigned jhit=0; jhit<hitIndices.size(); jhit++) {
      int hit2(hitIndices[jhit]);
      if(zVals[hitIndices[jhit]] != zLayers[iStartLayer+1]) continue;

      vector<float> hits;
      if(useBeamSpotConstraint) {
	hits.push_back(0);  // hacky beam-spot constraint
	hits.push_back(0);
	hits.push_back(0);
      }
      hits.push_back(rVals[hit1]*cos(phiVals[hit1]));  // x
      hits.push_back(rVals[hit1]*sin(phiVals[hit1]));  // y
      hits.push_back(zVals[hit1]);                     // z
      hits.push_back(rVals[hit2]*cos(phiVals[hit2]));  // x
      hits.push_back(rVals[hit2]*sin(phiVals[hit2]));  // y
      hits.push_back(zVals[hit2]);                     // z

      int nToAssert(useBeamSpotConstraint ? 9 : 6);
      if(iStartLayer > zLayers.size()-3) { // if there's only two layers left
	assert(hits.size()==nToAssert);
	LineFit lf(hits);
	lf.fit();
	lfs.push_back(lf);
	hit1s.push_back(hit1);
	hit2s.push_back(hit2);
	hit3s.push_back(-1);
      } else { // if there's three or more layers left
	// loop over all the hits in the next next layer
	for(unsigned khit=0; khit<hitIndices.size(); khit++) {
	  int hit3(hitIndices[khit]);
	  if(zVals[hitIndices[khit]] != zLayers[iStartLayer+2]) continue;
	  
	  if(hits.size() > nToAssert) hits.pop_back();
	  if(hits.size() > nToAssert) hits.pop_back();
	  if(hits.size() > nToAssert) hits.pop_back();
	  assert(hits.size() == nToAssert);
	  hits.push_back(rVals[hit3]*cos(phiVals[hit3]));  // x
	  hits.push_back(rVals[hit3]*sin(phiVals[hit3]));  // y
	  hits.push_back(zVals[hit3]);                     // z

	  LineFit lf(hits);
	  lf.fit();
	  lfs.push_back(lf);
	  hit1s.push_back(hit1);
	  hit2s.push_back(hit2);
	  hit3s.push_back(hit3);
	}
      }
    }
  }

  // find the best line
  double minChi2(99999);
  int iBestLine(-1);
  for(unsigned ilf=0; ilf<lfs.size(); ilf++) {
    
    // require closest distance to z axis to be small
    float zAtClosest = lfs[ilf].closestApproachToZAxis();
    float dzMax = (iStartLayer > zLayers.size()-3) ? 50 : 15; // use looser constraint if we've only got two hits to work with
    if(fabs(zAtClosest) > dzMax) {
      cout << "  dz too large: " << zAtClosest << endl;
      continue;
    }

    // require dxy at z=0 to be small
    double x,y,z;
    float dxyMax = (iStartLayer > zLayers.size()-3) ? 8 : 1.5; // use looser constraint if we've only got two hits to work with
    lineFcn(0, lfs[ilf].parFit, x, y, z);
    if(sqrt(x*x + y*y) > 2.5) {
      cout << " dxy too large: " << sqrt(x*x + y*y) << endl;
      continue;
    }

    double chi2(lfs[ilf].getChiSquare(dR,dS));
    if(chi2 < minChi2) {
      iBestLine = ilf;
      minChi2 = chi2;
    }
  }

  vector<int> hitsOfBestLine;
  if(iBestLine>=0) {
    LineFit *lfBest = new LineFit(lfs[iBestLine]);
    cout << "smallest dz: " << lfBest->closestApproachToZAxis() << endl;
    lines.push_back(lfBest);
    hitsOfBestLine.push_back(hit1s[iBestLine]);
    hitsOfBestLine.push_back(hit2s[iBestLine]);
    if(hit3s[iBestLine] >= 0)
      hitsOfBestLine.push_back(hit3s[iBestLine]);
  }

  return hitsOfBestLine;
}
//----------------------------------------------------------------------------------------
void Detector::fitTrack(vector<float> hits)
{
  LineFit *lf = new LineFit(hits);
  lf->fit();
  lines.push_back(lf);
}
//----------------------------------------------------------------------------------------
void Detector::findAllTracks()
{
  int iStartLayer(0);
  while(hitIndices.size() > 0) {
    cout << "layer: " << iStartLayer << " hitIndices before: " << hitIndices.size() << endl;
    vector<int> hitsOfBestLine(chooseHits(iStartLayer)); // find the three hits that fit to the best line
    
    if(hitsOfBestLine.size() != 3) {
      if(iStartLayer < zLayers.size() - 1) {
	iStartLayer++;
	continue;
      } else {
	cout << "no more lines!" << endl;
	break;
      }
    }
	
    // remove them from hitIndices
    assert(hitIndices.size()>2);
    for(unsigned ih=0; ih<hitIndices.size(); ih++)
      if(hitIndices[ih] == hitsOfBestLine[0])
	hitIndices.erase(hitIndices.begin()+ih);

    for(unsigned ih=0; ih<hitIndices.size(); ih++)
      if(hitIndices[ih] == hitsOfBestLine[1])
	hitIndices.erase(hitIndices.begin()+ih);

    for(unsigned ih=0; ih<hitIndices.size(); ih++)
      if(hitIndices[ih] == hitsOfBestLine[2])
	hitIndices.erase(hitIndices.begin()+ih);

  }
}
//----------------------------------------------------------------------------------------
void Detector::calcResolution(vector<Track> tracks)
{
  vector<int> foundTracks,usedLines; // find the best match, then remove this line and track, the next best...
  vector<float> distances;
  while(foundTracks.size() < tracks.size()) {
    float minDistance(999999);
    int iMatchTrk(-1),iMatchLine(-1);
    for(unsigned itrk=0; itrk<tracks.size(); itrk++) { // find the fitted line that comes closest to the track's origin
      bool alreadyFound(false);
      for(unsigned iFoundTrk=0; iFoundTrk<foundTracks.size(); iFoundTrk++)
	if(itrk == foundTracks[iFoundTrk])
	  alreadyFound = true;
      if(alreadyFound)
	continue;
	    
      for(unsigned iline=0; iline<lines.size(); iline++) {
	bool alreadyUsed(false);
	for(unsigned iUsedLine=0; iUsedLine<usedLines.size(); iUsedLine++)
	  if(iline == usedLines[iUsedLine])
	    alreadyUsed = true;
	if(alreadyUsed)
	  continue;

	float distance(sqrt(distance2(0, 0, tracks[itrk].dz, lines[iline]->parFit)));
	if(distance < minDistance) {
	  minDistance = distance;
	  iMatchTrk = itrk;
	  iMatchLine = iline;
	  // cout << "  set min to: " << setw(12) << minDistance << setw(12) << iMatchTrk << setw(12) << iMatchLine << endl;
	}
      }
    }
    // assert(iMatchTrk >= 0);
    // assert(iMatchLine >= 0);
    if(iMatchTrk < 0 || iMatchLine < 0) {
      cout << "out of tracks or lines" << endl;
      break;
    }
    foundTracks.push_back(iMatchTrk);
    usedLines.push_back(iMatchLine);
    distances.push_back(minDistance);
    // cout << "pushed back: " << distances.back() << setw(12) << iMatchTrk << setw(12) << iMatchLine << endl;
  }

  cout << "distances of " << distances.size() << " tracks" << endl;
  double sum(0);
  for(unsigned idist=0; idist<distances.size(); idist++) {
    cout << setw(12) << distances[idist] << endl;
    sum += distances[idist];
  }
  cout << "  mean: " << sum / distances.size() << endl;
}
