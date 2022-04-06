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

  _description = "V0 Finder Processor " ;

  registerInputCollection(LCIO::TRACK,
			  "TrackCollection",
			  "Name of input collection of reconstructed particles",
			  _trackColName,
			  std::string("LDCTracks"));

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			   "RecoParticleCollection",
			   "Name of output collection of reconstructed particles",
			   _recoPartColName,
			   std::string("V0RecoParticles"));

  registerOutputCollection(LCIO::VERTEX,
			   "VertexCollection",
			   "Name of output collection of neutral vertices",
			   _vertexColName,
			   std::string("LLPVertices"));
/*
  registerOutputCollection(LCIO::VERTEX,
			   "llpVertexCollection",
			   "Name of output collection of neutral vertices",
			   _llpVertexColName,
			   std::string("LLPVertices"));
*/
//   std::vector<float> rVertCut;
//   rVertCut.push_back(14.);
//   rVertCut.push_back(60.);
//   rVertCut.push_back(320.);
//   rVertCut.push_back(1600.);

  registerProcessorParameter("CutOnRadius",
			     "Cuts on V0 radius",
			     _rVertCut,
			     float(5.0));

  registerProcessorParameter("CutOnHitDistance",
			     "Cut on closest distance between first/last hit of two tracks",
			     _rHitCut,
			     float(100.0));

//   std::vector<float> dVertCut;
//   dVertCut.push_back(0.2);
//   dVertCut.push_back(1.0);
//   dVertCut.push_back(1.5);


  registerProcessorParameter("CutOnTrkDistance",
			     "Cut on two track distance",
			     _dVertCut,
			     float(1.5));

  registerProcessorParameter("MinimumTrackHitRatio",
			     "Minimum ratio of inner track hit radius to reconstructed vertex radius",
			     _minTrackHitRatio,
			     float(0.7));

  registerProcessorParameter("MassRangeGamma",
			     "Maximal deviation in mass for photon candidate",
			     _deltaMassGamma,
			     float(0.01));

  registerProcessorParameter("MassRangeK0S",
			     "Maximal deviation in mass for K0S candidate",
			     _deltaMassK0S,
			     float(0.01));

  registerProcessorParameter("MassRangeL0",
			     "Maximal deviation in mass for Lamda0 candidate",
			     _deltaMassL0,
			     float(0.008));

  registerProcessorParameter("RxyCutGamma",
			     "Minimum radius in xy plane for photon candidate",
			     _rxyCutGamma,
			     float(10.0));

  registerProcessorParameter("RxyCutK0S",
			     "Minimum radius in xy plane for K0S candidate",
			     _rxyCutK0S,
			     float(30.0));

  registerProcessorParameter("RxyCutLambda",
			     "Minimum radius in xy plane for Lambda0 candidate",
			     _rxyCutLambda,
			     float(50.0));



}

void LLPFinder::init() {

  MASSProton  = 0.93827203;
  MASSPion    = 0.13957018;
  MASSLambda0 = 1.115683;
  MASSK0S     = 0.497648;
  MASSGamma   = 0;

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

        // find smallest distance between tracks first/last hits
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

  float d01 = firstTrackState.getD0();
  float z01 = firstTrackState.getZ0();
  float phi1 = firstTrackState.getPhi();
  float tanLambda1 = firstTrackState.getTanLambda();
  float omega1 = firstTrackState.getOmega();
  HelixClass firstHelix;
  firstHelix.Initialize_Canonical(phi1,d01,z01,omega1,tanLambda1,_bField);
  float charge1 = firstHelix.getCharge();

  float r1 = firstTrack->getRadiusOfInnermostHit();

	float r2 = secondTrack->getRadiusOfInnermostHit();

	float d02 = secondTrackState.getD0();
	float z02 = secondTrackState.getZ0();
	float phi2 = secondTrackState.getPhi();
	float tanLambda2 = secondTrackState.getTanLambda();
	float omega2 = secondTrackState.getOmega();
	HelixClass secondHelix;
	secondHelix.Initialize_Canonical(phi2,d02,z02,omega2,tanLambda2,_bField);
	float charge2 = secondHelix.getCharge();
	float prodCharge = charge1*charge2;
	if (prodCharge<0) { // two tracks with opposite charges
	//if (1) { // FIXME requirement for opposite charges nontrivial for LLPs

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

	  float distV0;
	  float distHits;
	  float momentum[3];
	  float vertex[3];
	  float llp_vtx[3];

    if (pp1>pp2) {
    //if (pt1>pt2) {
	    distV0 = firstHelix.getDistanceToHelix(&secondHelix, vertex, momentum);

	  }
	  else {
	    distV0 = secondHelix.getDistanceToHelix(&firstHelix, vertex, momentum);
	  }

	  if (pt1>pt2) {
	    if(minDistIndex==0 || minDistIndex==1) {for(int idx=0;idx<3;idx++) llp_vtx[idx] = ftFirstHitLoc[idx];}
	    else {for(int idx=0;idx<3;idx++) llp_vtx[idx] = ftLastHitLoc[idx];}
	  }
	  else {
	    if(minDistIndex==0 || minDistIndex==2) {for(int idx=0;idx<3;idx++) llp_vtx[idx] = stFirstHitLoc[idx];}
	    else {for(int idx=0;idx<3;idx++) llp_vtx[idx] = stLastHitLoc[idx];}
	  }

	  float radV0 = sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]);
	  float radLLP = sqrt(llp_vtx[0]*llp_vtx[0]+llp_vtx[1]*llp_vtx[1]);

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


	  // check to ensure there are no hits on tracks at radii significantly smaller than reconstructed vertex
	  // TO DO: should be done more precisely using helices
	  //if(r1/radV0<_minTrackHitRatio)continue;
	  //if(r2/radV0<_minTrackHitRatio)continue;

	 //streamlog_out( DEBUG0 ) << "                       V0 : " << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;

	  //	  if (distV0 < _dVertCut && radV0 > _rVertCut ) { // cut on vertex radius and track misdistance
	  if (radV0 > _rVertCut  ) {

	    //streamlog_out( DEBUG0 ) << " ***************** found vertex for tracks : " << dd4hep::rec::Vector3D( (const float*) vertex ) << "\n"
		//		    << " t1 " << lcshort( firstTrack ) << "\n"
		//		    << " t2 " << lcshort( secondTrack )
		///	    << " distV0 " << distV0
		//		    << std::endl ;

    //if( distV0 < _dVertCut ) { // cut on vertex radius and track misdistance
    if( minDist < _rHitCut ) { // cut distance between first/last hits


	    streamlog_out( DEBUG0 ) << "  ***** testing various hypotheses " << std::endl ;

	    // Testing K0 hypothesis
	    float energy1 = sqrt(pp1*pp1+MASSPion*MASSPion);
	    float energy2 = sqrt(pp2*pp2+MASSPion*MASSPion);
	    float energyV0 = energy1 + energy2;
	    float massK0 = sqrt(energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2]);

	    // Testing L0 hypothesis
	    if (charge1<0) {
	      energy1 = sqrt(pp1*pp1+MASSPion*MASSPion);
	      energy2 = sqrt(pp2*pp2+MASSProton*MASSProton);
	    }
	    else {
	      energy1 = sqrt(pp1*pp1+MASSProton*MASSProton);
	      energy2 = sqrt(pp2*pp2+MASSPion*MASSPion);
	    }
	    energyV0 = energy1 + energy2;
	    float massL0 = sqrt(energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2]);

           // Testing L0bar hypothesis
            if (charge1>0) {
              energy1 = sqrt(pp1*pp1+MASSPion*MASSPion);
              energy2 = sqrt(pp2*pp2+MASSProton*MASSProton);
            }
            else {
              energy1 = sqrt(pp1*pp1+MASSProton*MASSProton);
              energy2 = sqrt(pp2*pp2+MASSPion*MASSPion);
            }
            energyV0 = energy1 + energy2;
            float massL0bar = sqrt(energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2]);



	    // Testing photon hypothesis
	    energyV0 = pp1 + pp2;
	    float massGamma = sqrt(energyV0*energyV0-momentum[0]*momentum[0]-momentum[1]*momentum[1]-momentum[2]*momentum[2]);

	    float deltaK0 = fabs(massK0 - MASSK0S);
	    float deltaL0 = fabs(massL0 - MASSLambda0);
	    float deltaGm = fabs(massGamma - MASSGamma);
	    float deltaL0bar = fabs(massL0bar - MASSLambda0);
	    if(radV0<_rxyCutGamma )deltaGm    = 100000.;
	    if(radV0<_rxyCutK0S   )deltaK0    = 100000.;
	    if(radV0<_rxyCutLambda)deltaL0    = 100000.;
	    if(radV0<_rxyCutLambda)deltaL0bar = 100000.;

	    int code = 22;
	    bool massCondition = false;

           if (deltaGm<deltaL0&&deltaGm<deltaK0&&deltaGm<deltaL0bar) {
              code = 22;
              massCondition = deltaGm < _deltaMassGamma;
            }
            else if (deltaK0<deltaL0 && deltaK0<deltaL0bar) {
              code = 310;
              massCondition = deltaK0 < _deltaMassK0S;
            }
            else{
              if (deltaL0<deltaL0bar ) {
                code = 3122;
                massCondition = deltaL0 < _deltaMassL0;
              }else{
                code = -3122;
                massCondition = deltaL0bar < _deltaMassL0;
              }
            }

	   streamlog_out( DEBUG0 ) << "  ***** mass condition :  " <<  massCondition
				  << "  code : " << code  << std::endl ;

	    // if (massCondition) {
	    //if (true) { // always save for now
	      bool ok = true;
	      if(r1/radV0<_minTrackHitRatio|| r2/radV0<_minTrackHitRatio){
		r1 = this->Rmin(firstTrack);
		r2 = this->Rmin(secondTrack);
		if(r1/radV0<_minTrackHitRatio || r2/radV0<_minTrackHitRatio)ok = false;
		//streamlog_out( DEBUG0 ) << " V0X: " << ok << " r = " << radV0 << " r1 = " << r1 << " r2 = " << r2 << std::endl;
	      }
	      //if(!ok)continue;
	      TrackPair * trkPair = new TrackPair();
	      trkPair->setFirstTrack( firstTrack );
	      trkPair->setSecondTrack( secondTrack );
        //trkPair->setDistance( distV0 );
        trkPair->setDistance( minDist );
        //trkPair->setVertex( vertex );
        trkPair->setVertex( llp_vtx );
	      trkPair->setMomentum( momentum );
	      trkPair->setCode( code );
	      trkPairs.push_back( trkPair );
	      streamlog_out( DEBUG0 ) << "!!!!   Track Pair saved !!!!" << std::endl;

	    //}
	    /* else {
 	      streamlog_out( DEBUG0 ) << "Rejected vertex : V = ("
 			<< vertex[0] << ","
 			<< vertex[1] << ","
 			<< vertex[2] << ")" << std::endl;
	    } */

 	    streamlog_out( DEBUG0 ) << "Code = " << code << std::endl;
      //streamlog_out( DEBUG0 ) << "Vertex = " << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
      streamlog_out( DEBUG0 ) << "Vertex = " << llp_vtx[0] << " " << llp_vtx[1] << " " << llp_vtx[2] << std::endl;
 	    streamlog_out( DEBUG0 ) << "Momentum = " << momentum[0] << " " << momentum[1] << " " << momentum[2] << std::endl;


	  }
	}
      }
    }

    }//MESSAGE ------

     streamlog_out( DEBUG0 ) << std::endl;

    // Sorting of all vertices in ascending order of the track misdistance

    int nTrkPairs = int(trkPairs.size());
    streamlog_out( DEBUG0 ) << "Number of track pairs = " << nTrkPairs << std::endl;

    if (nTrkPairs>0) { // V0s are present in event

            //streamlog_out( DEBUG0 ) << "Number of track pairs = " << nTrkPairs << std::endl;

      Sorting( trkPairs );

      // Declaration of the output collections
      LCCollectionVec * colRecoPart = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
      LCCollectionVec * colVertex   = new LCCollectionVec(LCIO::VERTEX);

      for (int iTrkP=0;iTrkP<nTrkPairs;++iTrkP) {
	TrackPair * pair = trkPairs[iTrkP];
	Track * firstTrack = pair->getFirstTrack();
	Track * secondTrack = pair->getSecondTrack();
	if (trackUsed[firstTrack]==0&&trackUsed[secondTrack]==0) {

	  ReconstructedParticleImpl * part = new ReconstructedParticleImpl();
	  VertexImpl * vtx = new VertexImpl();

	  float vertex[3];
	  float momentum[3];
	  int code = pair->getCode();
	  for (int iC=0;iC<3;++iC) {
	    vertex[iC] = pair->getVertex()[iC];
	    momentum[iC] = pair->getMomentum()[iC];
	  }

	  float distance = pair->getDistance();
	  vtx->setPosition( vertex );
	  vtx->addParameter( distance );

	  part->setMomentum( momentum );
	  part->setType( code );

 	  streamlog_out( DEBUG0 ) << "Code = " << code << "  Distance = " << distance << std::endl;
 	  streamlog_out( DEBUG0 ) << "Vertex = ("
 		    << vertex[0] << ","
 		    << vertex[1] << ","
 		    << vertex[2] << ")" << std::endl;

 	  streamlog_out( DEBUG0 ) << "Momentum = ("
 		    << momentum[0] << ","
 		    << momentum[1] << ","
 		    << momentum[2] << ")" << std::endl;
 	  streamlog_out( DEBUG0 ) << firstTrack << " " << secondTrack << std::endl;


	  float mass = 0;
	  if ( code == 22)
	    mass = 0;
	  else if ( code == 310 )
	    mass = MASSK0S;
	  else
	    mass = MASSLambda0;

	  part->setMass( mass );
	  vtx->setAssociatedParticle( part );
	  part->setStartVertex( vtx );
	  part->addTrack( firstTrack );
	  part->addTrack( secondTrack );

	  colRecoPart->addElement( part );
	  colVertex->addElement( vtx );

	  //trackUsed[firstTrack] = 1;
	  //trackUsed[secondTrack] = 1;
	}
      }

      evt->addCollection( colRecoPart,_recoPartColName.c_str() );
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
    for (int j = 0; j < sizeOfVector-i-1; j++)
      {
	one = trkPairVec[j];
	two = trkPairVec[j+1];
	if( one->getDistance() > two->getDistance() )
	  {
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
