#include "LLPFinder.h"
#include "marlin/Global.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "EVENT/TrackState.h"
#include "IMPL/TrackStateImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/VertexImpl.h"
#include <IMPL/LCRelationImpl.h>
#include "UTIL/Operators.h"
#include <math.h>

#include <DD4hep/Detector.h>
#include <DD4hep/DD4hepUnits.h>
#include <DDRec/Vector3D.h>


#include "HelixClassT.h"

using namespace lcio ;
using namespace marlin ;


LLPFinder aLLPFinder ;


LLPFinder::LLPFinder() : Processor("LLPFinder") {

  _description = "LLP Finder Processor " ;

  registerInputCollection(LCIO::TRACK,
			  "TrackCollection",
			  "Name of input collection of tracks",
			  _trackColName,
			  std::string("LDCTracks"));
/*
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			   "RecoParticleCollection",
			   "Name of output collection of reconstructed particles",
			   _recoPartColName,
			   std::string("LLPRecoParticles"));
*/
  registerOutputCollection(LCIO::VERTEX,
			   "VertexCollection",
			   "Name of output collection of neutral vertices",
			   _vertexColName,
			   std::string("LLPVertices"));

  registerOutputCollection(LCIO::LCRELATION,
			    "VertexTracksLinkName" ,
			    "Name of the VertexTracksLink output collection"  ,
			    _vertexTracksLinkName ,
			    std::string("VertexTracksLink") );

  registerProcessorParameter("CutOnHelixDistance",
			     "Cut on maximum closest distance between first/last hit of two tracks",
			     _rHelixCut,
			     float(100.0));

  registerProcessorParameter("UseHelixDistance",
			     "Use helix distance method to find vertex",
			     _useHelixDist,
			     bool(true));

  registerProcessorParameter("UseCutOnHelixDistance",
			     "Use cut on hit distance as cut on helix distance",
			     _useHelixDistCut,
			     bool(true));

  registerProcessorParameter("CutOnTrackCurvatureRatio",
			     "Cut on maximum ratio of the curvatures of the considered tracks",
			     _curvatureRatioCut,
			     float(0.95));

  registerProcessorParameter("CutOnTrackOpeningAngle",
			     "Cut on maximum cosine of opening angle between considered tracks",
			     _cosOpeningAngleCut,
			     float(0.995));

  registerProcessorParameter("CutOnRefPointDistance",
			     "Maximum distance between the reference points of considered TrackStates",
			     _refPointDistCut,
			     float(600));

  registerProcessorParameter("CutOnPtAtVtx",
			     "Minimum pT of sum of the track momenta",
			     _vtxPtCut,
			     float(0.3));

  registerProcessorParameter("CutOnTrackNdf",
			     "Reject vtx candidate if a track has smaller Ndf",
			     _trkNdfCut,
			     int(10));

  registerProcessorParameter("CutOnLeadingTrackPt",
			     "Cut on the minimum pT of the leading track",
			     _leadingTrackPtCut,
			     float(0.4));

  registerProcessorParameter("CutOnSecondTrackPt",
			     "Cut on the minimum pT of the second track",
			     _secondTrackPtCut,
			     float(0.4));

  registerProcessorParameter("CutOnShortTrackPt",
			     "Minimum pT of the artificial short high-pT tracks",
			     _shortTrackPtCut,
			     float(4));

  registerProcessorParameter("CutOnShortTrackNdf",
			     "Maximum Ndf of the artificial short high-pT tracks",
			     _shortTrackNdfCut,
			     int(60));

  registerProcessorParameter("MinDistToUseZCut",
			     "Minimum track length in Z to use MaxRefToArcZRatio cut",
			     _minLenToCutOnZ,
			     float(100));

  registerProcessorParameter("MaxRefToArcPhiRatio",
			     "Maximum ratio of ref. point phi to arc phi w.r.t. vtx phi",
			     _refToArcPhiCut,
			     float(0.05));

  registerProcessorParameter("MaxRefToArcZRatio",
			     "Minimum ratio of ref. point Z to track length in Z w.r.t. vtx Z",
			     _refToLenZCut,
			     float(-0.05));

  registerProcessorParameter("CutOnRefPointsAndCirclesDist",
			     "Cut on combined distances between ref. points and helix circle centres",
			     _refPointsCirclesDistCut,
			     float(-2000));



}

void LLPFinder::init() {

  dd4hep::Detector& theDet = dd4hep::Detector::getInstance();
  double bfieldV[3] ;
  theDet.field().magneticField( { 0., 0., 0. }  , bfieldV  ) ;
  _bField = bfieldV[2]/dd4hep::tesla ;

  _nRun = -1;
  _nEvt = 0;

}


void LLPFinder::processRunHeader( LCRunHeader* ) {

  _nRun++ ;
  _nEvt = 0;

}

void LLPFinder::processEvent( LCEvent * evt ) {


  try {

    LCCollection * col = evt->getCollection( _trackColName.c_str() );

    int nelem = col->getNumberOfElements();

    TrackPairVec  trkPairs;
    trkPairs.clear();

    std::map<Track*,int> trackUsed;

    for (int i=0;i<nelem;++i) {
      Track * trk = dynamic_cast<Track*>(col->getElementAt(i));
      trackUsed[trk] = 0;
    }

    for (int i=0;i<nelem-1;++i) {
      Track * firstTrack = dynamic_cast<Track*>(col->getElementAt(i));

      TrackState *firstTrackFHTS = (TrackState*) firstTrack->getTrackState(lcio::TrackState::AtFirstHit);
      TrackState *firstTrackLHTS = (TrackState*) firstTrack->getTrackState(lcio::TrackState::AtLastHit);

      float ftFirstHitLoc[3];
      float ftLastHitLoc[3];
      for (int idx=0;idx<3;idx++){
	      ftFirstHitLoc[idx] = *(firstTrackFHTS->getReferencePoint()+idx);
	      ftLastHitLoc[idx] = *(firstTrackLHTS->getReferencePoint()+idx);
      }

      for (int j=i+1;j<nelem;++j) {
	      Track * secondTrack = dynamic_cast<Track*>(col->getElementAt(j));

	      streamlog_out( DEBUG0 ) << " ***************** candidate tracks : "
				    << " t" << i << " " << lcshort( firstTrack ) << "\n"
				    << " t" << j << " " << lcshort( secondTrack )
				    << std::endl ;

      	TrackState *secondTrackFHTS = (TrackState*) secondTrack->getTrackState(lcio::TrackState::AtFirstHit);
      	TrackState *secondTrackLHTS = (TrackState*) secondTrack->getTrackState(lcio::TrackState::AtLastHit);

        std::vector<TrackState*> trackStates;
        trackStates.push_back(firstTrackFHTS);
        trackStates.push_back(firstTrackLHTS);
        trackStates.push_back(secondTrackFHTS);
        trackStates.push_back(secondTrackLHTS);

      	float stFirstHitLoc[3];
        float stLastHitLoc[3];
        for (int idx=0;idx<3;idx++){
      		stFirstHitLoc[idx] = *(secondTrackFHTS->getReferencePoint()+idx);
      		stLastHitLoc[idx] = *(secondTrackLHTS->getReferencePoint()+idx);
        }

        // find the smallest distance between tracks first/last hits
        // indices: 0 - ftFirstHit-stFirstHit, 1 - ftFirstHit-stLastHit, 2 - ftLastHit-stFirstHit, 3 - ftLastHit-stLastHit
    	  std::vector<float> hitDists = getHitDists(firstTrack,secondTrack);
    	  int minDistIdx = std::min_element(hitDists.begin(),hitDists.end()) - hitDists.begin();
    	  float minDist = *std::min_element(hitDists.begin(), hitDists.end());

        // get two closest TrackStates
        TrackStateImpl firstTrackState;
        TrackStateImpl secondTrackState;
        if(minDistIdx==0) {
          // firstTrackState = fixTrackStateDirection(firstTrackFHTS,firstTrackLHTS);
          // secondTrackState = fixTrackStateDirection(secondTrackFHTS,secondTrackLHTS);
          firstTrackState = fixTrackStateDirection(trackStates.at(0),trackStates.at(1));
          secondTrackState = fixTrackStateDirection(trackStates.at(2),trackStates.at(3));
        }
        if(minDistIdx==1) {
          // firstTrackState=fixTrackStateDirection(firstTrackFHTS,firstTrackLHTS);
          // secondTrackState=fixTrackStateDirection(secondTrackLHTS,secondTrackFHTS);
          firstTrackState=fixTrackStateDirection(trackStates.at(0),trackStates.at(1));
          secondTrackState=fixTrackStateDirection(trackStates.at(3),trackStates.at(2));
        }
        if(minDistIdx==2) {
          // firstTrackState=fixTrackStateDirection(firstTrackLHTS,firstTrackFHTS);
          // secondTrackState=fixTrackStateDirection(secondTrackFHTS,secondTrackLHTS);
          firstTrackState=fixTrackStateDirection(trackStates.at(1),trackStates.at(0));
          secondTrackState=fixTrackStateDirection(trackStates.at(2),trackStates.at(3));
        }
        if(minDistIdx==3) {
          // firstTrackState=fixTrackStateDirection(firstTrackLHTS,firstTrackFHTS);
          // secondTrackState=fixTrackStateDirection(secondTrackLHTS,secondTrackFHTS);
          firstTrackState=fixTrackStateDirection(trackStates.at(1),trackStates.at(0));
          secondTrackState=fixTrackStateDirection(trackStates.at(3),trackStates.at(2));
        }
        // get track states at last hits
        TrackStateImpl firstTrackLastHitTS = getOtherTrackState(firstTrack, firstTrackState);
        TrackStateImpl secondTrackLastHitTS = getOtherTrackState(secondTrack, secondTrackState);

        // get helices at fixed TrackStates
        std::vector<HelixClass> fixedHelices; // first helix on 0, 1 indices
        for (int iTS=0; iTS<4; iTS+=2){
          TrackStateImpl fixedTsFh = fixTrackStateDirection(trackStates[iTS],trackStates[iTS+1]);
          TrackStateImpl fixedTsLh = fixTrackStateDirection(trackStates[iTS+1],trackStates[iTS]);
          fixedHelices.push_back( getHelix(fixedTsFh) );
          fixedHelices.push_back( getHelix(fixedTsLh) );
        }
        // check charge products of helices
        std::vector<int> chargeProds;
        std::map<float,int> distCharge;
        for (int iHel=0; iHel<2; iHel++) {
          for (int jHel=2; jHel<4; jHel++) {
            int chargeProd = fixedHelices[iHel].getCharge() * fixedHelices[jHel].getCharge();
            chargeProds.push_back(chargeProd);
          }
        }
        float minLegalDist = 50000;
        for (int iCh=0; iCh<4; iCh++) {
          streamlog_out(DEBUG1) << " Charge product " << chargeProds[iCh] << " at dist " << iCh << std::endl;
          if(chargeProds[iCh] < 0) {
            distCharge[hitDists[iCh]] = 1;
            streamlog_out(DEBUG1) << " Hit dist " << hitDists[iCh] << " at dist " << iCh << std::endl;
            if (hitDists[iCh] < minLegalDist) minLegalDist = hitDists[iCh];
          }
          else distCharge[hitDists[iCh]] = 0;
        }
        // find index of the smallest "legal" distance
        std::vector<float>::iterator itr = std::find(hitDists.begin(), hitDists.end(), minLegalDist);
        int minLegalDistIdx = std::distance(hitDists.begin(), itr);

        minDist = minLegalDist;
        minDistIdx = minLegalDistIdx;

/*
        TrackStateImpl fixedFtFh = fixTrackStateDirection(firstTrackFHTS,firstTrackLHTS);
        TrackStateImpl fixedFtLh = fixTrackStateDirection(firstTrackLHTS,firstTrackFHTS);
        TrackStateImpl fixedStFh = fixTrackStateDirection(secondTrackFHTS,secondTrackLHTS);
        TrackStateImpl fixedStLh = fixTrackStateDirection(secondTrackLHTS,secondTrackFHTS);

        HelixClass firstHelixFH = getHelix(fixedFtFh);
        HelixClass firstHelixLH = getHelix(fixedFtLh);
        HelixClass secondHelixFH = getHelix(fixedStFh);
        HelixClass secondHelixLH = getHelix(fixedStLh);
*/
        float r1 = firstTrack->getRadiusOfInnermostHit();
      	float r2 = secondTrack->getRadiusOfInnermostHit();

        HelixClass firstHelix, secondHelix, firstHelix_LH, secondHelix_LH;
        if(minDistIdx<2) {
          firstHelix = fixedHelices[0]; // for minDistIdx<2 we always have first hit in first helix
          secondHelix = fixedHelices[minDistIdx+2];
          firstHelix_LH = fixedHelices[1];
          secondHelix_LH = fixedHelices[((minDistIdx==0)?3:2)];
        }
        else {
          firstHelix = fixedHelices[1];
          secondHelix = fixedHelices[minDistIdx];
          firstHelix_LH = fixedHelices[0];
          secondHelix_LH = fixedHelices[((minDistIdx==2)?3:2)];
        }

        float distLLP;
        float distLLP1;
        float distLLP2;
      	float distHits;
        float momentum[3];
        float momentum1[3];
        float momentum2[3];
        float vertex[3];
        float vertex1[3];
        float vertex2[3];
      	float llp_vtx[3];

        // get distance between helices
        distLLP1 = firstHelix.getDistanceToHelix(&secondHelix, vertex1, momentum1);
        distLLP2 = secondHelix.getDistanceToHelix(&firstHelix, vertex2, momentum2);

        float px1 = firstHelix.getMomentum()[0];
      	float py1 = firstHelix.getMomentum()[1];
      	float pz1 = firstHelix.getMomentum()[2];
      	float pp1 = sqrt(px1*px1+py1*py1+pz1*pz1);
        float px2 = secondHelix.getMomentum()[0];
      	float py2 = secondHelix.getMomentum()[1];
      	float pz2 = secondHelix.getMomentum()[2];
      	float pp2 = sqrt(px2*px2+py2*py2+pz2*pz2);

        // getDistanceToHelix() is "asymmetric" and works better called on a track with higher momentum
        if (pp1 > pp2) {
          distLLP = distLLP1;
          vertex[0] = vertex1[0];
          vertex[1] = vertex1[1];
          vertex[2] = vertex1[2];
          momentum[0] = momentum1[0];
          momentum[1] = momentum1[1];
          momentum[2] = momentum1[2];
        }
        else {
          distLLP = distLLP2;
          vertex[0] = vertex2[0];
          vertex[1] = vertex2[1];
          vertex[2] = vertex2[2];
          momentum[0] = momentum2[0];
          momentum[1] = momentum2[1];
          momentum[2] = momentum2[2];
        }

        // check if TS first hits are indeed closer to vertex than last hits
        float fTsToVtx = getDistance(firstTrackState.getReferencePoint(), vertex);
        float sTsToVtx = getDistance(secondTrackState.getReferencePoint(), vertex);
        float ftLhToVtx = getDistance(firstTrackLastHitTS.getReferencePoint(), vertex);
        float stLhToVtx = getDistance(secondTrackLastHitTS.getReferencePoint(), vertex);
        if ( (ftLhToVtx+stLhToVtx) < (fTsToVtx+sTsToVtx) ) { // if not, use the helices on the opposite ends
          streamlog_out( DEBUG7 ) << " *** Reconstructed vertex at (" 
          << vertex[0] << ", " << vertex[1] << ", " << vertex[2] 
          << ") closer to the last hits! Inverting helices... " << std::endl;

          firstHelix = firstHelix_LH;
          secondHelix = secondHelix_LH;
          distLLP1 = firstHelix.getDistanceToHelix(&secondHelix, vertex1, momentum1);
          distLLP2 = secondHelix.getDistanceToHelix(&firstHelix, vertex2, momentum2);

          // swap also trackStates
          TrackStateImpl tempTS1 = firstTrackState;
          TrackStateImpl tempTS2 = secondTrackState;
          firstTrackState = firstTrackLastHitTS;
          secondTrackState = secondTrackLastHitTS;
          firstTrackLastHitTS = tempTS1;
          secondTrackLastHitTS = tempTS2;

          // reassign momenta
          px1 = firstHelix.getMomentum()[0];
          py1 = firstHelix.getMomentum()[1];
          pz1 = firstHelix.getMomentum()[2];
          pp1 = sqrt(px1*px1+py1*py1+pz1*pz1);
          px2 = secondHelix.getMomentum()[0];
          py2 = secondHelix.getMomentum()[1];
          pz2 = secondHelix.getMomentum()[2];
          pp2 = sqrt(px2*px2+py2*py2+pz2*pz2);
        }

      	
      	float pt1 = sqrt(px1*px1+py1*py1);
        float xC1 = firstHelix.getXC();
        float yC1 = firstHelix.getYC();
      	float pt2 = sqrt(px2*px2+py2*py2);
        float xC2 = secondHelix.getXC();
        float yC2 = secondHelix.getYC();

        float charge1 = firstHelix.getCharge();
      	float charge2 = secondHelix.getCharge();
      	float prodCharge = charge1*charge2;
        float omega1 = firstHelix.getOmega();
      	float omega2 = secondHelix.getOmega();

        int Ndf1 = firstTrack->getNdf();
        int Ndf2 = secondTrack->getNdf();

        // getDistanceToHelix() is "asymmetric" and works better called on a track with higher momentum
        if (pp1 > pp2) {
          distLLP = distLLP1;
          vertex[0] = vertex1[0];
          vertex[1] = vertex1[1];
          vertex[2] = vertex1[2];
          momentum[0] = momentum1[0];
          momentum[1] = momentum1[1];
          momentum[2] = momentum1[2];
        }
        else {
          distLLP = distLLP2;
          vertex[0] = vertex2[0];
          vertex[1] = vertex2[1];
          vertex[2] = vertex2[2];
          momentum[0] = momentum2[0];
          momentum[1] = momentum2[1];
          momentum[2] = momentum2[2];
        }

        float leadingPt, secondPt;
        int leadingNdf, secondNdf;
         if (pt1 > pt2) {
           leadingPt = pt1;
           secondPt = pt2;
           leadingNdf = Ndf1;
           secondNdf = Ndf2;
         }
         else {
           leadingPt = pt2;
           secondPt = pt1;
           leadingNdf = Ndf2;
           secondNdf = Ndf1;
         }

        // calculate variables for cuts
        float ptAtVtx = std::hypot( (px1+px2),(py1+py2) );
        float thetaAtVtx = atan2(ptAtVtx,(pz1+pz2));
        float cosOpenAngle = (px1*px2 + py1*py2 + pz1*pz2) / (pp1*pp2) ;
        float circlesDist = std::hypot( (xC2-xC1),(yC2-yC1) );
        float refPointsCirclesDist = 2.2 * minDist - circlesDist ; // hardcoded dependence tested for ILD-TPC!
        float curvRatio = fabs(omega1/omega2);

        // angles and arcs for the cut on random intersections
        float phiVtx1 = getPhiOnHelix(vertex,xC1,yC1);
        float phiVtx2 = getPhiOnHelix(vertex,xC2,yC2);
        float phiRef1 = getPhiOnHelix(firstTrackState.getReferencePoint(),xC1,yC1);
        float phiRef2 = getPhiOnHelix(secondTrackState.getReferencePoint(),xC2,yC2);
        float phiLH1 = getPhiOnHelix(firstTrackLastHitTS.getReferencePoint(),xC1,yC1);
        float phiLH2 = getPhiOnHelix(secondTrackLastHitTS.getReferencePoint(),xC2,yC2);
        float arc1 = getArcLength(phiRef1, phiLH1, charge1);
        float arc2 = getArcLength(phiRef2, phiLH2, charge2);

        // get angles w.r.t. to vertex phi
        phiRef1 = rotateAngle(phiRef1, -phiVtx1);
        phiRef2 = rotateAngle(phiRef2, -phiVtx2);
        phiLH1 = rotateAngle(phiLH1, -phiVtx1);
        phiLH2 = rotateAngle(phiLH2, -phiVtx2);

        bool cutPhi1, cutPhi2;
        if (charge1 < 0) {
          cutPhi1 = charge1 * phiRef1 / arc1 < -_refToArcPhiCut;
          cutPhi2 = charge2 * ( 2 * M_PI - phiRef2 ) / arc2 > _refToArcPhiCut;
        }
        else {
          cutPhi1 = charge1 * ( 2 * M_PI - phiRef1 ) / arc1 > _refToArcPhiCut;
          cutPhi2 = charge2 * phiRef2 / arc2 < -_refToArcPhiCut;
        }

        float vtx_refpointZ1 = firstTrackState.getReferencePoint()[2] - vertex[2];
        float vtx_lhZ1 = firstTrackLastHitTS.getReferencePoint()[2] - vertex[2];
        float vtx_refpointZ2 = secondTrackState.getReferencePoint()[2] - vertex[2];
        float vtx_lhZ2 = secondTrackLastHitTS.getReferencePoint()[2] - vertex[2];
        float refpoint_lhZ1 = fabs( firstTrackState.getReferencePoint()[2] - firstTrackLastHitTS.getReferencePoint()[2] );
        float refpoint_lhZ2 = fabs( secondTrackState.getReferencePoint()[2] - secondTrackLastHitTS.getReferencePoint()[2] );
        int sign_pz1 = (pz1 > 0) ? 1 : ((pz1 < 0) ? -1 : 0);
        int sign_pz2 = (pz2 > 0) ? 1 : ((pz2 < 0) ? -1 : 0);

        bool cutZ1 = sign_pz1 * vtx_refpointZ1 / refpoint_lhZ1 < _refToLenZCut; 
        bool cutZ2 = sign_pz2 * vtx_refpointZ2 / refpoint_lhZ2 < _refToLenZCut;

        if (pt1>pt2) {
          if(minDistIdx==0 || minDistIdx==1) {for(int idx=0;idx<3;idx++) llp_vtx[idx] = ftFirstHitLoc[idx];}
          else {for(int idx=0;idx<3;idx++) llp_vtx[idx] = ftLastHitLoc[idx];}
        }
        else {
          if(minDistIdx==0 || minDistIdx==2) {for(int idx=0;idx<3;idx++) llp_vtx[idx] = stFirstHitLoc[idx];}
          else {for(int idx=0;idx<3;idx++) llp_vtx[idx] = stLastHitLoc[idx];}
        }

        if (!_useHelixDist) {
          vertex[0] = llp_vtx[0]; vertex[1] = llp_vtx[1]; vertex[2] = llp_vtx[2];
        }

        // float radLLP = sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]+vertex[2]*vertex[2]);
        float radLLP = sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]);

        streamlog_out( DEBUG0 ) << " ***************** the vertex for tracks : " << dd4hep::rec::Vector3D( (const float*) vertex ) << "\n"
                  << " with helix1 ref. point : " << dd4hep::rec::Vector3D( (const float*) firstHelix.getReferencePoint() )
                  << " and helix2 ref. point : " << dd4hep::rec::Vector3D( (const float*) secondHelix.getReferencePoint() ) << "\n"
                  << " with helix1 momentum : " << dd4hep::rec::Vector3D( (const float*) firstHelix.getMomentum() )
                  << " and helix2 momentum : " << dd4hep::rec::Vector3D( (const float*) secondHelix.getMomentum() ) << "\n"
                << " t1 phi " << firstTrackState.getPhi() << "\n"
                << " t2 phi " << secondTrackState.getPhi() << "\n"
                << " distLLP " << distLLP << "\n"
                << " curvatureRatio " << curvRatio << "\n"
                << " cosOpenAngle " << cosOpenAngle
                << std::endl ;
        streamlog_out( DEBUG0 ) << " ***************** first/last hits distances: ";
        for (int idx=0; idx<4; idx++) streamlog_out( DEBUG0 ) << hitDists[idx] << " ";
        streamlog_out( DEBUG0 ) << " min distance: " << minDist << " at index " << minDistIdx << std::endl;

        // cuts

        bool distCut = minDist < _rHelixCut;
        if (_useHelixDistCut) distCut = distLLP < _rHelixCut;

        if( distCut ) { // cut on distance between helices
          if (!_useHelixDist) distLLP = minDist;

          streamlog_out( DEBUG7 ) << " ### Passed cut on distance between helices: " << distLLP << " < " << _rHelixCut << std::endl;

          // reject tracks that are just randomly intersecting
          if ((refpoint_lhZ1 > _minLenToCutOnZ && cutZ1) || (refpoint_lhZ2 > _minLenToCutOnZ && cutZ2))
            continue;
          if ((refpoint_lhZ1 < _minLenToCutOnZ && cutPhi1) || (refpoint_lhZ2 < _minLenToCutOnZ && cutPhi2))
            continue;

          streamlog_out( DEBUG7 ) << " ### Passed cut on random intersections... " << std::endl;

          // cut on ratios of track curvatures (pTs)
          if (fabs(1.-curvRatio) > fabs(1.-_curvatureRatioCut) && cosOpenAngle > _cosOpeningAngleCut) {

            streamlog_out( DEBUG7 ) << " ### Passed cut on track curvatures (" << fabs(1.-curvRatio) << " > " << fabs(1.-_curvatureRatioCut)
            << ") and opening angle (" << cosOpenAngle << " > " << _cosOpeningAngleCut << ")" << std::endl;

            if (refPointsCirclesDist < _refPointsCirclesDistCut) { // cut on distance between two reference points and helix circles

            streamlog_out( DEBUG7 ) << " ### Passed cut on distance between two reference points and helix circles: " 
                                    << refPointsCirclesDist << " < " << _refPointsCirclesDistCut << std::endl;

              if (ptAtVtx > _vtxPtCut) {

                streamlog_out( DEBUG7 ) << " ### Passed cut on pT at vertex: " << ptAtVtx << " > " << _vtxPtCut << std::endl;

                if ( Ndf2 > _trkNdfCut ) {

                  streamlog_out( DEBUG7 ) << " ### Passed cut on Ndf: " << Ndf2 << " > " << _trkNdfCut << std::endl;

                  streamlog_out( DEBUG7 ) << "### Accepted theta : "
                  << thetaAtVtx << " at event " << _nEvt << std::endl;
                  streamlog_out( DEBUG7 ) << "### ptAtVtx : "
                  << ptAtVtx << " pz1: " << pz1 << " pz2 " << pz2 << std::endl;

                  if (!(leadingPt < _leadingTrackPtCut && secondPt < _secondTrackPtCut)) {
                  // if (secondPt > _leadingTrackPtCut*leadingPt + _secondTrackPtCut ) {

                    streamlog_out( DEBUG7 ) << " ### Passed cut on leading and second tracks pTs " << std::endl;
                    
                    if (!(leadingPt > _shortTrackPtCut && leadingNdf < _shortTrackNdfCut)) {

                      streamlog_out( DEBUG7 ) << " ### Passed cut on artificial, short tracks " << std::endl;

                      streamlog_out( DEBUG2 ) << "Accepted min. distance between refPoints: "
                      << minDist << " at index " << minDistIdx << std::endl;

                      TrackPair * trkPair = new TrackPair();
                      trkPair->setFirstTrack( firstTrack );
                      trkPair->setSecondTrack( secondTrack );
                      trkPair->setDistance( distLLP );
                      trkPair->setVertex( vertex );
                      trkPair->setMomentum( momentum );
                      //trkPair->setCode( code );
                      trkPairs.push_back( trkPair );
                      streamlog_out( DEBUG0 ) << firstTrack << " " << secondTrack << std::endl;
                      streamlog_out( DEBUG0 ) << "!!!!   Track Pair saved !!!!" << std::endl;


                      //streamlog_out( DEBUG0 ) << "Code = " << code << std::endl;
                      streamlog_out( DEBUG0 ) << "Vertex = " << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
                      streamlog_out( DEBUG0 ) << "Momentum = " << momentum[0] << " " << momentum[1] << " " << momentum[2] << std::endl;
                  
                    }
                  }
                }
              }

            }

	        }
        }

      }
    } //MESSAGE ------

    streamlog_out( DEBUG0 ) << std::endl;

    // Sorting of all vertices in ascending order of the track misdistance

    int nTrkPairs = int(trkPairs.size());
    streamlog_out( DEBUG0 ) << "Number of track pairs = " << nTrkPairs << std::endl;

    if (nTrkPairs>0) { // LLPs are present in event

      Sorting( trkPairs );

      // Declaration of the output collections
      //LCCollectionVec * colRecoPart = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
      LCCollectionVec * colVertex   = new LCCollectionVec(LCIO::VERTEX);
      LCCollectionVec * colVtxToTrkRel   = new LCCollectionVec(LCIO::LCRELATION);

      for (int iTrkP=0;iTrkP<nTrkPairs;++iTrkP) {
      	TrackPair * pair = trkPairs[iTrkP];
      	Track * firstTrack = pair->getFirstTrack();
      	Track * secondTrack = pair->getSecondTrack();
      	if (trackUsed[firstTrack]==0&&trackUsed[secondTrack]==0) {

      	  //ReconstructedParticleImpl * part = new ReconstructedParticleImpl();
      	  VertexImpl * vtx = new VertexImpl();

      	  float vertex[3];
      	  //float momentum[3];
      	  //int code = pair->getCode();
      	  for (int iC=0;iC<3;++iC) {
      	    vertex[iC] = pair->getVertex()[iC];
      	    //momentum[iC] = pair->getMomentum()[iC];
      	  }

      	  float distance = pair->getDistance();
      	  vtx->setPosition( vertex );
          vtx->addParameter( distance );

      	  //part->setMomentum( momentum );
      	  //part->setType( code );

       	  //streamlog_out( DEBUG0 ) << "Code = " << code << "  Distance = " << distance << std::endl;
       	  streamlog_out( DEBUG0 ) << "Vertex = ("
       		    << vertex[0] << ","
       		    << vertex[1] << ","
       		    << vertex[2] << "), "
              << "Distance = " << distance << std::endl;
          /*
       	  streamlog_out( DEBUG0 ) << "Momentum = ("
       		    << momentum[0] << ","
       		    << momentum[1] << ","
       		    << momentum[2] << ")" << std::endl;
          */
       	  streamlog_out( DEBUG0 ) << firstTrack << " " << secondTrack << std::endl;

          /*
      	  part->setMass( mass );
      	  vtx->setAssociatedParticle( part );
      	  part->setStartVertex( vtx );
      	  part->addTrack( firstTrack );
      	  part->addTrack( secondTrack );

      	  colRecoPart->addElement( part );
          */

      	  colVertex->addElement( vtx );

      	  trackUsed[firstTrack] = 1;
      	  trackUsed[secondTrack] = 1;

          // set relations from vertex to corresponding tracks
          LCRelationImpl *vtxTracksRelNav1 = new LCRelationImpl( vtx , firstTrack, 1.0 ) ;
          LCRelationImpl *vtxTracksRelNav2 = new LCRelationImpl( vtx , secondTrack , 1.0 ) ;

          // vtxTracksRelNav1.addRelation(   vtx , firstTrack, 1.0 ) ;
          // vtxTracksRelNav2.addRelation(   vtx , secondTrack , 1.0 ) ;

          colVtxToTrkRel->addElement( vtxTracksRelNav1 );
          colVtxToTrkRel->addElement( vtxTracksRelNav2 );

          streamlog_out( DEBUG3 ) << " =============== tracks : "
  				    << " t" << iTrkP << " " << lcshort( firstTrack ) << "\n"
  				    << " t" << iTrkP << " " << lcshort( secondTrack ) << "\n"
  				    << " =============== saved !!! " << std::endl ;
      	}
      }

      //evt->addCollection( colRecoPart,_recoPartColName.c_str() );
      evt->addCollection( colVertex, _vertexColName.c_str() );
      evt->addCollection( colVtxToTrkRel, _vertexTracksLinkName.c_str() );

      streamlog_out( DEBUG0 ) << "!!!!   Collections saved !!!!" << std::endl;

    }

    // Clean up memory
    for (int iTrkP=0;iTrkP<nTrkPairs;++iTrkP) {
      TrackPair * trkPair = trkPairs[iTrkP];
      delete trkPair;
    }
    trkPairs.clear();

    //    getchar();

  }
  catch(DataNotAvailableException &e) {}



  _nEvt++;

}


void LLPFinder::check( LCEvent * ) { }

void LLPFinder::end(){ }

void LLPFinder::Sorting( TrackPairVec & trkPairVec ) {

  int sizeOfVector = int(trkPairVec.size());
  TrackPair *one,*two,*Temp;

  for (int i = 0 ; i < sizeOfVector-1; i++)
    for (int j = 0; j < sizeOfVector-i-1; j++) {
    	one = trkPairVec[j];
    	two = trkPairVec[j+1];
    	if( one->getDistance() > two->getDistance() ) {
  	    Temp = trkPairVec[j];
  	    trkPairVec[j] = trkPairVec[j+1];
  	    trkPairVec[j+1] = Temp;
  	  }
    }

}

float LLPFinder::Rmin( Track* track ) {

   // find track extrema

  float rmin = 1000000.;
  TrackerHitVec hitvec = track->getTrackerHits();
  int nhits = (int)hitvec.size();
  float zmax =-99999.;
  float zmin =99999.;
  for(int ih =0;ih<nhits;++ih){
    float z = (float)hitvec[ih]->getPosition()[2];
    if(z<zmin)zmin=z;
    if(z>zmax)zmax=z;
  }
  float tanLambda = track->getTanLambda();
  //  streamlog_out( DEBUG0 ) << " V0 Check : " << tanLambda << " z = " << zmin << " - " << zmax << std::endl;
  float zzz = zmin;
  if(tanLambda<0)zzz=zmax;

  float zstart = 0;
  if(fabs(zmin)<fabs(zmax))zstart = zmin;
  if(fabs(zmax)<fabs(zmin))zstart = zmax;
  //streamlog_out( DEBUG0 ) << " V0 Check " << zstart << " - " << zzz << std::endl;
  for(int ih =0;ih<nhits;++ih){
    float z = (float)hitvec[ih]->getPosition()[2];
    if(fabs(z-zstart)<250){
      float x = (float)hitvec[ih]->getPosition()[0];
      float y = (float)hitvec[ih]->getPosition()[1];
      float r2 = x*x+y*y;
      float  r = sqrt(r2);
      if(r<rmin)rmin = r;
    }

  }

  return rmin;

}

float LLPFinder::getDistance(const float *p1, const float *p2) {
  return sqrt( (p1[0]-p2[0])*(p1[0]-p2[0])
              +(p1[1]-p2[1])*(p1[1]-p2[1])
              +(p1[2]-p2[2])*(p1[2]-p2[2]) );
}

std::vector<float> LLPFinder::getHitDists( Track* firstTrack, Track* secondTrack ) {

      std::vector<float> dists;

      TrackStateImpl firstTrackFHTS = *(firstTrack->getTrackState(lcio::TrackState::AtFirstHit));
      TrackStateImpl firstTrackLHTS = *(firstTrack->getTrackState(lcio::TrackState::AtLastHit));

      float ftFirstHitLoc[3];
      float ftLastHitLoc[3];

      for (int idx=0;idx<3;idx++){
	      ftFirstHitLoc[idx] = *(firstTrackFHTS.getReferencePoint()+idx);
	      ftLastHitLoc[idx] = *(firstTrackLHTS.getReferencePoint()+idx);
      }

      streamlog_out( DEBUG0 ) << " ------- First Track: \n";
      streamlog_out( DEBUG0 ) << " - first hit: ";
      for (int idx=0;idx<3;idx++){
	      streamlog_out( DEBUG0 ) << ftFirstHitLoc[idx] << " ";
      }
      streamlog_out( DEBUG0 ) << "\n";
      streamlog_out( DEBUG0 ) << " - last hit: ";
      for (int idx=0;idx<3;idx++){
	      streamlog_out( DEBUG0 ) << ftLastHitLoc[idx] << " ";
      }
      streamlog_out( DEBUG0 ) << "\n";


      TrackStateImpl secondTrackFHTS = *(secondTrack->getTrackState(lcio::TrackState::AtFirstHit));
      TrackStateImpl secondTrackLHTS = *(secondTrack->getTrackState(lcio::TrackState::AtLastHit));

      float stFirstHitLoc[3];
      float stLastHitLoc[3];
      for (int idx=0;idx<3;idx++){
	stFirstHitLoc[idx] = *(secondTrackFHTS.getReferencePoint()+idx);
	stLastHitLoc[idx] = *(secondTrackLHTS.getReferencePoint()+idx);
      }

      streamlog_out( DEBUG0 ) << " ------- Second Track: \n";
      streamlog_out( DEBUG0 ) << " - first hit: ";
      for (int idx=0;idx<3;idx++){
	      streamlog_out( DEBUG0 ) << stFirstHitLoc[idx] << " ";
      }
      streamlog_out( DEBUG0 ) << "\n";
      streamlog_out( DEBUG0 ) << " - last hit: ";
      for (int idx=0;idx<3;idx++){
	      streamlog_out( DEBUG0 ) << stLastHitLoc[idx] << " ";
      }
      streamlog_out( DEBUG0 ) << "\n";

      dists.push_back( sqrt( pow((ftFirstHitLoc[0]-stFirstHitLoc[0]),2) + pow((ftFirstHitLoc[1]-stFirstHitLoc[1]),2) + pow((ftFirstHitLoc[2]-stFirstHitLoc[2]),2) ) );
      dists.push_back( sqrt( pow((ftFirstHitLoc[0]-stLastHitLoc[0]),2) + pow((ftFirstHitLoc[1]-stLastHitLoc[1]),2) + pow((ftFirstHitLoc[2]-stLastHitLoc[2]),2) ) );
      dists.push_back( sqrt( pow((ftLastHitLoc[0]-stFirstHitLoc[0]),2) + pow((ftLastHitLoc[1]-stFirstHitLoc[1]),2) + pow((ftLastHitLoc[2]-stFirstHitLoc[2]),2) ) );
      dists.push_back( sqrt( pow((ftLastHitLoc[0]-stLastHitLoc[0]),2) + pow((ftLastHitLoc[1]-stLastHitLoc[1]),2) + pow((ftLastHitLoc[2]-stLastHitLoc[2]),2) ) );

      return dists;

}

TrackStateImpl LLPFinder::flipTrackState( TrackState* ts ) {

  int loc = ts->getLocation();
	float d0 = ts->getD0();
	float z0 = ts->getZ0();
	float phi0 = ts->getPhi();
	float tanLambda = ts->getTanLambda();
	float omega = ts->getOmega();

  float refPoint[3];
	for (int iref=0; iref<3; iref++) refPoint[iref] = ts->getReferencePoint()[iref];

  FloatVec covMatrix;
  int covSize = ts->getCovMatrix().size();
  covMatrix.resize(covSize);
  for (int icov=0; icov<covSize; icov++) covMatrix[icov]=ts->getCovMatrix()[icov];

  float phi = phi0;
	if (phi0 > 0) phi = -1*(M_PI - phi);
	else phi = M_PI + phi;

  TrackStateImpl flipped_ts;

  flipped_ts.setLocation(loc);
	flipped_ts.setD0(-d0);
	flipped_ts.setPhi(phi);
	flipped_ts.setOmega(-omega);
	flipped_ts.setZ0(z0);
	flipped_ts.setTanLambda(-tanLambda);
	flipped_ts.setReferencePoint(refPoint);
	flipped_ts.setCovMatrix(covMatrix);

  return flipped_ts;

}

TrackStateImpl LLPFinder::fixTrackStateDirection( TrackState* ts, TrackState* other_ts ) {

  float d0 = ts->getD0();
  float z0 = ts->getZ0();
  float phi = ts->getPhi();
  float tanLambda = ts->getTanLambda();
  float omega = ts->getOmega();
  float ref[3];
  for (int i=0; i<3; i++) {ref[i] = ts->getReferencePoint()[i];}

  HelixClass helix;
  helix.Initialize_Canonical(phi,d0,z0,omega,tanLambda,_bField,ref);

  int direction = 1;
	if (ts->getReferencePoint()[2] > other_ts->getReferencePoint()[2])
		direction = -1;

  TrackStateImpl fixed_ts;

  if (direction * helix.getMomentum()[2] > 0)
		fixed_ts = *ts;
	else if (direction * helix.getMomentum()[2] < 0)
		fixed_ts = flipTrackState(ts);
	else
		streamlog_out(ERROR0) << " dir * p_z = " << direction * helix.getMomentum()[2] << std::endl;

  return fixed_ts;
}

TrackStateImpl LLPFinder::getOtherTrackState(Track *trk, TrackStateImpl ts) {
  int loc = ts.getLocation();
  TrackStateImpl other_ts = (loc == 2) ? *trk->getTrackState(3) : *trk->getTrackState(2);
  return other_ts;
}

HelixClass LLPFinder::getHelix(Track* tr) {

  float d0 = tr->getD0();
  float z0 = tr->getZ0();
  float phi = tr->getPhi();
  float tanLambda = tr->getTanLambda();
  float omega = tr->getOmega();
  float ref[3];
  for (int idx=0; idx<3; idx++) {ref[idx] = tr->getReferencePoint()[idx];}

  HelixClass helix;
  helix.Initialize_Canonical(phi,d0,z0,omega,tanLambda,_bField,ref);

  return helix;

}

HelixClass LLPFinder::getHelix(TrackStateImpl ts) {

  float d0 = ts.getD0();
  float z0 = ts.getZ0();
  float phi = ts.getPhi();
  float tanLambda = ts.getTanLambda();
  float omega = ts.getOmega();
  float ref[3];
  for (int idx=0; idx<3; idx++) {ref[idx] = ts.getReferencePoint()[idx];}

  HelixClass helix;
  helix.Initialize_Canonical(phi,d0,z0,omega,tanLambda,_bField,ref);

  return helix;

}

float LLPFinder::getPhiOnHelix(const float *vertex, float xC, float yC) {
  float phi = atan2( vertex[1]-yC, vertex[0]-xC );
  if (phi < 0) phi += 2*M_PI;
  return phi;
}

float LLPFinder::rotateAngle(float angle, float rotation) {
	float rotatedAngle = angle + rotation;
	if (rotatedAngle < 0) rotatedAngle += 2*M_PI;
	return rotatedAngle;
}

float LLPFinder::getArcLength(float phi1, float phi2, int q) {
  float dPhi;
	if (q < 0) {
    dPhi = phi2 - phi1;
    if (dPhi < 0) dPhi += 2*M_PI;
  }
  else {
    dPhi = phi1 - phi2;
    if (dPhi < 0) dPhi += 2*M_PI;
  }
	return dPhi;
}