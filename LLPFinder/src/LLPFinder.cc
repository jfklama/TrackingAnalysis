#include "LLPFinder.h"
#include "marlin/Global.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "EVENT/TrackState.h"
#include "IMPL/TrackStateImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/VertexImpl.h"
#include "UTIL/Operators.h"
#include <math.h>

#include <DD4hep/Detector.h>
#include <DD4hep/DD4hepUnits.h>
#include <DDRec/Vector3D.h>


#include "HelixClass.h"

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

  registerProcessorParameter("CutOnHitDistance",
			     "Cut on closest distance between first/last hit of two tracks",
			     _rHitCut,
			     float(100.0));

  registerProcessorParameter("UseHelixDistance",
			     "Use helix distance method to find vertex",
			     _useHelixDist,
			     bool(false));

  registerProcessorParameter("UseCutOnHelixDistance",
			     "Use cut on hit distance as cut on helix distance",
			     _useHelixDistCut,
			     bool(false));



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

      	float stFirstHitLoc[3];
        float stLastHitLoc[3];
        for (int idx=0;idx<3;idx++){
      		stFirstHitLoc[idx] = *(secondTrackFHTS->getReferencePoint()+idx);
      		stLastHitLoc[idx] = *(secondTrackLHTS->getReferencePoint()+idx);
        }

        // find the smallest distance between tracks first/last hits
    	  std::vector<float> hitDists = getHitDists(firstTrack,secondTrack);
    	  int minDistIndex = std::min_element(hitDists.begin(),hitDists.end()) - hitDists.begin();
    	  float minDist = *std::min_element(hitDists.begin(), hitDists.end());

        TrackStateImpl firstTrackState;
        TrackStateImpl secondTrackState;
        if(minDistIndex==0) {
          firstTrackState = fixTrackStateDirection(firstTrackFHTS,firstTrackLHTS);
          secondTrackState = fixTrackStateDirection(secondTrackFHTS,secondTrackLHTS);
        }
        if(minDistIndex==1) {
          firstTrackState=fixTrackStateDirection(firstTrackFHTS,firstTrackLHTS);
          secondTrackState=fixTrackStateDirection(secondTrackLHTS,secondTrackFHTS);
        }
        if(minDistIndex==2) {
          firstTrackState=fixTrackStateDirection(firstTrackLHTS,firstTrackFHTS);
          secondTrackState=fixTrackStateDirection(secondTrackFHTS,secondTrackLHTS);
        }
        if(minDistIndex==3) {
          firstTrackState=fixTrackStateDirection(firstTrackLHTS,firstTrackFHTS);
          secondTrackState=fixTrackStateDirection(secondTrackLHTS,secondTrackFHTS);
        }

        float r1 = firstTrack->getRadiusOfInnermostHit();
      	float r2 = secondTrack->getRadiusOfInnermostHit();

        float d01 = firstTrack->getD0();
        float z01 = firstTrack->getZ0();
        float phi1 = firstTrack->getPhi();
        float tanLambda1 = firstTrack->getTanLambda();
        float omega1 = firstTrack->getOmega();

        HelixClass firstHelix;
        firstHelix.Initialize_Canonical(phi1,d01,z01,omega1,tanLambda1,_bField);
        float charge1 = firstHelix.getCharge();

      	float d02 = secondTrack->getD0();
      	float z02 = secondTrack->getZ0();
      	float phi2 = secondTrack->getPhi();
      	float tanLambda2 = secondTrack->getTanLambda();
      	float omega2 = secondTrack->getOmega();

      	HelixClass secondHelix;
      	secondHelix.Initialize_Canonical(phi2,d02,z02,omega2,tanLambda2,_bField);
      	float charge2 = secondHelix.getCharge();

      	float prodCharge = charge1*charge2;

      	//if (prodCharge<0) { // two tracks with opposite charges
      	if (1) { // FIXME requirement for opposite charges nontrivial for LLPs

      	  float px1 = firstHelix.getMomentum()[0];
      	  float py1 = firstHelix.getMomentum()[1];
      	  float pz1 = firstHelix.getMomentum()[2];
      	  float pp1 = sqrt(px1*px1+py1*py1+pz1*pz1);
      	  float pt1 = sqrt(px1*px1+py1*py1);

      	  float px2 = secondHelix.getMomentum()[0];
      	  float py2 = secondHelix.getMomentum()[1];
      	  float pz2 = secondHelix.getMomentum()[2];
      	  float pp2 = sqrt(px2*px2+py2*py2+pz2*pz2);
      	  float pt2 = sqrt(px2*px2+py2*py2);

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

          /*
          if (pp1>pp2) {
          //if (pt1>pt2) {
      	    distLLP = firstHelix.getDistanceToHelix(&secondHelix, vertex, momentum);

      	  }
      	  else {
      	    distLLP = secondHelix.getDistanceToHelix(&firstHelix, vertex, momentum);
      	  }
          */

      	  distLLP1 = firstHelix.getDistanceToHelix(&secondHelix, vertex1, momentum1);
      	  distLLP2 = secondHelix.getDistanceToHelix(&firstHelix, vertex2, momentum2);

          if (distLLP1 < distLLP2) {
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

      	  if (pt1>pt2) {
      	    if(minDistIndex==0 || minDistIndex==1) {for(int idx=0;idx<3;idx++) llp_vtx[idx] = ftFirstHitLoc[idx];}
      	    else {for(int idx=0;idx<3;idx++) llp_vtx[idx] = ftLastHitLoc[idx];}
      	  }
      	  else {
      	    if(minDistIndex==0 || minDistIndex==2) {for(int idx=0;idx<3;idx++) llp_vtx[idx] = stFirstHitLoc[idx];}
      	    else {for(int idx=0;idx<3;idx++) llp_vtx[idx] = stLastHitLoc[idx];}
      	  }

          if (!_useHelixDist) {
            vertex[0] = llp_vtx[0];
            vertex[1] = llp_vtx[1];
            vertex[2] = llp_vtx[2];
            //distLLP = minDist;
          }

      	  float radLLP = sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]);

          streamlog_out( DEBUG0 ) << " ***************** the vertex for tracks : " << dd4hep::rec::Vector3D( (const float*) vertex ) << "\n"
                    << " with Ref. Point t1 : " << dd4hep::rec::Vector3D( (const float*) firstTrackState.getReferencePoint() )
                    << " and Ref. Point t2 : " << dd4hep::rec::Vector3D( (const float*) secondTrackState.getReferencePoint() ) << "\n"
      				    << " t1 " << lcshort( firstTrack ) << "\n"
      				    << " t2 " << lcshort( secondTrack )
      				    << " distLLP " << radLLP
      				    << std::endl ;
          streamlog_out( DEBUG0 ) << " ***************** first/last hits distances: ";
          for (int idx=0; idx<4; idx++) streamlog_out( DEBUG0 ) << hitDists[idx] << " ";
          streamlog_out( DEBUG0 ) << " min distance: " << minDist << " at index " << minDistIndex << std::endl;


          bool distCut = minDist < _rHitCut;
          if (_useHelixDistCut) distCut = distLLP < _rHitCut;

          if( distCut ) { // cut on distance between first/last hits (or between helices)

            if (!_useHelixDist) distLLP = minDist;

      	    TrackPair * trkPair = new TrackPair();
      	    trkPair->setFirstTrack( firstTrack );
      	    trkPair->setSecondTrack( secondTrack );
            trkPair->setDistance( distLLP );
            trkPair->setVertex( vertex );
      	    trkPair->setMomentum( momentum );
      	    //trkPair->setCode( code );
      	    trkPairs.push_back( trkPair );
      	    streamlog_out( DEBUG0 ) << "!!!!   Track Pair saved !!!!" << std::endl;


       	    //streamlog_out( DEBUG0 ) << "Code = " << code << std::endl;
            streamlog_out( DEBUG0 ) << "Vertex = " << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
       	    streamlog_out( DEBUG0 ) << "Momentum = " << momentum[0] << " " << momentum[1] << " " << momentum[2] << std::endl;

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
       		    << vertex[2] << ")" << std::endl;
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

      	  //trackUsed[firstTrack] = 1;
      	  //trackUsed[secondTrack] = 1;
      	}
      }

      //evt->addCollection( colRecoPart,_recoPartColName.c_str() );
      evt->addCollection( colVertex, _vertexColName.c_str() );

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
  HelixClass helix;
  helix.Initialize_Canonical(phi,d0,z0,omega,tanLambda,_bField);

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
