import math
from pyLCIO import IOIMPL, EVENT, IMPL, ROOT, IOException 
import utils

def vertex_analysis(vertexCollection, trackCollection, \
                    vertexTracksRelations, refitTrackToMCLinkCollection, \
                    llp_endpoint, counter, histograms):
    dist0 = 50000
    helixDist0 = 50000
    bestVtx = IMPL.VertexImpl()
    closest_dists = get_closest_dists(trackCollection)
    trueR = utils.R(llp_endpoint[0], llp_endpoint[1])
    trueZ = llp_endpoint[2]
    trueL = utils.spacialDistance(llp_endpoint, [0,0,0])
    in_acceptance = (trueR <= 1808) and abs(trueZ) <= 2350
    in_TPC = in_acceptance and trueR >= 329

    for vtx in vertexCollection:

        helixDist = vtx.getParameters()[0]
        histograms.helixDistances.Fill(helixDist )
        if len(closest_dists) > 0:
            closest_dist = min(closest_dists)
            histograms.hitsDistances.Fill( closest_dist )
        pos = vtx.getPosition()
        histograms.vtxVsR_reco.Fill( utils.R(pos[0],pos[1]) )
        histograms.vtxVsL_reco.Fill( utils.R3(pos[0],pos[1],pos[2]) )
        dist = utils.spacialDistance(llp_endpoint, pos )
        #vtxMisdistances.Fill(dist)
        # if dist < 30:
        # 	print iEvt, dist0, dist, helixDist

        if dist < dist0:
            dist0 = dist
            pos0 = pos
            helixDist0 = helixDist
            bestVtx = vtx
        if dist > 30:
            histograms.fakeVertices.Fill( pos[2], utils.R(pos[0],pos[1]) )

        '''
        if helixDist < helixDist0:
            dist0 = dist
            pos0 = pos
            helixDist0 = helixDist
        '''
    histograms.vtxMisdistances.Fill(dist0)
    if dist0 < 30:
        histograms.recoVertices.Fill( pos0[2], utils.R(pos0[0],pos0[1]) )
        histograms.vtxEffRVsZ_reco.Fill(trueZ,trueR)
        histograms.vtxEffVsL_reco.Fill(trueL)
        # 	_matchingVerticesTPC+=1
        # _nFakeVertices += vertexCollection.getNumberOfElements() - 1
    # else:
        #print 'no good vtx in ev. ' + str(iEvt) + ' with ' + str(vertexCollection.getNumberOfElements()) + ' vertices'
        # _nFakeVertices += vertexCollection.getNumberOfElements()
    counter.increment("_nMatchingVertices", dist0 < 30, 1)
    counter.increment("_matchingVerticesTPC", dist0 < 30 and in_TPC, 1)
    counter.increment("_nFakeVertices", dist0 < 30, vertexCollection.getNumberOfElements() - 1)
    counter.increment("_nFakeVertices", dist0 >= 30, vertexCollection.getNumberOfElements())

    trksVec = []
    verticesVec = []
    for rel in vertexTracksRelations:
        if rel.getFrom() == bestVtx and dist0 < 30:
            trksVec.append(rel.getTo())
            verticesVec.append(rel.getFrom())

    # find MC particles associated to tracks from vtx
    mcP4 = ROOT.TLorentzVector()
    for trk in trksVec:
        for rel in refitTrackToMCLinkCollection:
            if rel.getFrom()==trk and rel.getWeight() > 0.8:
                particle = rel.getTo()
                e = particle.getEnergy()
                px = particle.getMomentum()[0]
                py = particle.getMomentum()[1]
                pz = particle.getMomentum()[2]
                mcP4.SetPxPyPzE(px,py,pz,e)
                histograms.vtxTrksTruePDG.Fill(abs(particle.getPDG()))

    iTrk=0
    while iTrk < len(trksVec):

        trkSt1, trkSt2 = utils.findClosestLegalTrackStates(trksVec[iTrk],trksVec[iTrk+1])
        vtx_cand = verticesVec[iTrk]

        omega1 = trkSt1.getOmega()
        omega2 = trkSt2.getOmega()
        q1 = omega1/abs(omega1)
        q2 = omega2/abs(omega2)
        ratio = abs(omega1/omega2)
        if ratio > 1.:
            ratio = abs(omega2/omega1)
        curlerOmega = 0.002
        histograms.vtxTrksAreBothCurlers.Fill( abs(omega1) > curlerOmega and abs(omega2) > curlerOmega )

        trkStates = [trkSt1, trkSt2]
        momenta = []

        if  utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 330 and utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) < 1800:
        # if True:

            # print iEvt, llp_endpoint[0],llp_endpoint[1], llp_endpoint[2]
            # print iEvt, vtx_cand.getPosition()[0],vtx_cand.getPosition()[1], vtx_cand.getPosition()[2], vtx_cand.getParameters()[0]
            # print trkSt1.getD0(), trkSt1.getZ0(), omega1, trkSt1.getTanLambda(), trkSt1.getPhi()
            # print trkSt2.getD0(), trkSt2.getZ0(), omega2, trkSt2.getTanLambda(), trkSt2.getPhi()
            # print trkSt1.getReferencePoint()[0],trkSt1.getReferencePoint()[1],trkSt1.getReferencePoint()[2]
            # print trkSt2.getReferencePoint()[0],trkSt2.getReferencePoint()[1],trkSt2.getReferencePoint()[2]

            p_vtx = ROOT.TVector3()
            p_vtx.SetXYZ(0,0,0)

            for trkSt in trkStates:
                mom = utils.getTrackMomentum(trkSt)
                p = ROOT.TVector3()
                p.SetXYZ(mom[0],mom[1],mom[2])
                momenta.append(p)
                p_vtx += p

                # print iEvt, mom[0],mom[1],mom[2]

                histograms.vtxTracksP.Fill(p.Mag())
                histograms.vtxTracksPt.Fill(p.Pt())
                histograms.vtxTracksTheta.Fill(p.Theta())
                histograms.vtxTracksPhi.Fill(p.Phi())
                vtx_refpoint = utils.spacialDistance( trkSt.getReferencePoint(),vtx_cand.getPosition() )
                histograms.vtxRefPointDist.Fill(vtx_refpoint)

                histograms.tracksVtxPtVsTheta.Fill(p.Theta(),p.Pt())
                histograms.tracksVtxPvsTheta.Fill(p.Theta(),p.Mag())
            # print ''
            if dist0 < 30:
                histograms.vtxP.Fill( p_vtx.Mag() )
                histograms.vtxPt.Fill( p_vtx.Pt() )
                histograms.vtxTheta.Fill( p_vtx.Theta() )
                histograms.vtxPhi.Fill( p_vtx.Phi() )

            if momenta[0].Pt() > momenta[1].Pt():
                leadingP = momenta[0]
                secondP = momenta[1]
                leadingTrack = trksVec[iTrk]
                secondTrack = trksVec[iTrk+1]
            else:
                secondP = momenta[0]
                leadingP = momenta[1]
                secondTrack = trksVec[iTrk]
                leadingTrack = trksVec[iTrk+1]

            histograms.vtxTracksPtSum.Fill(leadingP.Pt()+secondP.Pt())

            vtx_refpoint1 = utils.spacialDistance( trkSt1.getReferencePoint(),vtx_cand.getPosition() )
            vtx_refpoint2 = utils.spacialDistance( trkSt2.getReferencePoint(),vtx_cand.getPosition() )
            histograms.vtxRefPointDist1VsDist2.Fill(vtx_refpoint1,vtx_refpoint2)

            # if vtx_refpoint1 > 200 and vtx_refpoint2 > 200:
            # 	print iEvt, vtx_cand.getPosition()[0],vtx_cand.getPosition()[1], vtx_cand.getPosition()[2]

            histograms.vtxTrksRefPointZ.Fill( abs(trkSt1.getReferencePoint()[2]-trkSt2.getReferencePoint()[2]) )
            distRef12 = utils.spacialDistance( trkSt1.getReferencePoint(), trkSt2.getReferencePoint() )
            histograms.vtxTrksRefPointsDist.Fill( distRef12 )
            # print iEvt, vtx_cand.getPosition()[0], vtx_cand.getPosition()[1], vtx_cand.getPosition()[2], vtx_cand.getParameters()[0], distRef12
            # print iEvt, trkSt1.getReferencePoint()[0], trkSt1.getReferencePoint()[1], trkSt1.getReferencePoint()[2],trksVec[iTrk].getNdf()
            # print iEvt, trkSt2.getReferencePoint()[0], trkSt2.getReferencePoint()[1], trkSt2.getReferencePoint()[2],trksVec[iTrk+1].getNdf()
            # print ''

            vtx_refpoint1 = utils.spacialDistance( trkSt1.getReferencePoint(),vtx_cand.getPosition() )
            lh_loc1 = 3 if trkSt1.getLocation() == 2 else 2
            lastHit1 = trksVec[iTrk].getTrackState(lh_loc1)
            vtx_lh1 = utils.spacialDistance( lastHit1.getReferencePoint(), vtx_cand.getPosition() )
            refpoint1_lh1 = utils.spacialDistance( trkSt1.getReferencePoint(), lastHit1.getReferencePoint() )

            vtx_refpoint2 = utils.spacialDistance( trkSt2.getReferencePoint(),vtx_cand.getPosition() )
            lh_loc2 = 3 if trkSt2.getLocation() == 2 else 2
            lastHit2 = trksVec[iTrk+1].getTrackState(lh_loc2)
            vtx_lh2 = utils.spacialDistance( lastHit2.getReferencePoint(), vtx_cand.getPosition() )
            refpoint2_lh2 = utils.spacialDistance( trkSt2.getReferencePoint(), lastHit2.getReferencePoint() )

            vtx_refpointZ1 = trkSt1.getReferencePoint()[2] - vtx_cand.getPosition()[2]
            vtx_lhZ1 = lastHit1.getReferencePoint()[2] - vtx_cand.getPosition()[2]
            vtx_refpointZ2 = trkSt2.getReferencePoint()[2] - vtx_cand.getPosition()[2]
            vtx_lhZ2 = lastHit2.getReferencePoint()[2] - vtx_cand.getPosition()[2]
            refpoint_lhZ1 = abs( trkSt1.getReferencePoint()[2] - lastHit1.getReferencePoint()[2] )
            refpoint_lhZ2 = abs( trkSt2.getReferencePoint()[2] - lastHit2.getReferencePoint()[2] )
            sign_pz1 = momenta[0].Pz() / abs(momenta[0].Pz())
            sign_pz2 = momenta[1].Pz() / abs(momenta[1].Pz())

            xC1, yC1 = utils.getHelixXYC(trkSt1)
            xC2, yC2 = utils.getHelixXYC(trkSt2)
            centresDist = math.sqrt( (xC1-xC2)**2 + (yC1-yC2)**2 )
            centresDistScaled = centresDist / math.sqrt( leadingP.Pt() * secondP.Pt() )
            centresRefsDists = 2.2 * distRef12 - centresDist
            logCentresRefsDists = 300 * math.log10(distRef12) - centresDist
            cosOpenAngle = math.cos(momenta[0].Angle(momenta[1]))

            phiVtx1 = math.atan2( vtx_cand.getPosition()[1]-yC1, vtx_cand.getPosition()[0]-xC1 )
            phiVtx2 = math.atan2( vtx_cand.getPosition()[1]-yC2, vtx_cand.getPosition()[0]-xC2 )
            if phiVtx1 < 0:
                phiVtx1 += 2*math.pi
            if phiVtx2 < 0:
                phiVtx2 += 2*math.pi
            phiRef1 = math.atan2( trkSt1.getReferencePoint()[1]-yC1, trkSt1.getReferencePoint()[0]-xC1 )
            if phiRef1 < 0:
                phiRef1 += 2*math.pi
            # phiRef1 = phiRef1 - phiVtx1
            phiRef2 = math.atan2( trkSt2.getReferencePoint()[1]-yC2, trkSt2.getReferencePoint()[0]-xC2 )
            if phiRef2 < 0:
                phiRef2 += 2*math.pi
            # phiRef2 = phiRef2 - phiVtx2
            phiLH1 = math.atan2( lastHit1.getReferencePoint()[1]-yC1, lastHit1.getReferencePoint()[0]-xC1 )
            if phiLH1 < 0:
                phiLH1 += 2*math.pi
            # phiLH1 = phiLH1 - phiVtx1
            phiLH2 = math.atan2( lastHit2.getReferencePoint()[1]-yC2, lastHit2.getReferencePoint()[0]-xC2 )
            if phiLH2 < 0:
                phiLH2 += 2*math.pi
            arc1 = utils.getArcLength(phiRef1, phiLH1, q1)
            arc2 = utils.getArcLength(phiRef2, phiLH2, q2)
            phiRef1 = utils.rotateAngle(phiRef1, -phiVtx1)
            phiRef2 = utils.rotateAngle(phiRef2, -phiVtx2)
            phiLH1 = utils.rotateAngle(phiLH1, -phiVtx1)
            phiLH2 = utils.rotateAngle(phiLH2, -phiVtx2)
            if arc1 < 0 or arc2 < 0:
                print(arc1, arc2)

            cutPhi1 = True
            cutPhi2 = True
            if q1 < 0:
                cutPhi1 = q1*phiRef1/arc1 < -0.05
                cutPhi2 = q2*(2*math.pi-phiRef2)/arc2 > 0.05
            else:
                cutPhi1 = q1*(2*math.pi-phiRef1)/arc1 > 0.05
                cutPhi2 = q2*phiRef2/arc2 < -0.05

            # if ( (refpoint_lhZ1 > 100 and sign_pz1*vtx_refpointZ1/refpoint_lhZ1 < -0.05) \
            # or (refpoint_lhZ2 > 100 and sign_pz2*vtx_refpointZ2/refpoint_lhZ2 < -0.05) ):
            # 	iTrk += 2
            # 	return
            # if ( (refpoint_lhZ1 < 100 and cutPhi1) \
            # or (refpoint_lhZ2 < 100 and cutPhi2) ):
            # 	iTrk += 2
            # 	return

            histograms.vtxRefZVsLhZ.Fill(sign_pz1*vtx_refpointZ1/refpoint_lhZ1, sign_pz2*vtx_refpointZ2/refpoint_lhZ2)
            histograms.vtxRefZVsLenZ.Fill(sign_pz1*vtx_refpointZ1/refpoint_lhZ1, refpoint_lhZ1)
            histograms.vtxRefZVsLenZ.Fill(sign_pz2*vtx_refpointZ2/refpoint_lhZ2, refpoint_lhZ2)

            if q1 < 0:
                if refpoint_lhZ1 < 100000 and refpoint_lhZ2 < 100000:
                    histograms.vtxRefPhiVsLhPhi.Fill(q1*phiRef1/arc1, q2*(2*math.pi-phiRef2)/arc2)
                    histograms.vtxRefPhiVsArcPhi.Fill(q1*phiRef1/arc1, arc1)
                if refpoint_lhZ2 < 100000:
                    # histograms.vtxRefPhiVsLhPhi.Fill(q2*(2*math.pi-phiRef2)/arc2, q2*(2*math.pi-phiLH2)/arc2)
                    histograms.vtxRefPhiVsArcPhi.Fill(q2*(2*math.pi-phiRef2)/arc2,arc2)
            else:
                if refpoint_lhZ1 < 100000 and refpoint_lhZ2 < 100000:
                    histograms.vtxRefPhiVsLhPhi.Fill(q1*(2*math.pi-phiRef1)/arc1, q2*phiRef2/arc2)
                    histograms.vtxRefPhiVsArcPhi.Fill(q1*(2*math.pi-phiRef1)/arc1, arc1)
                if refpoint_lhZ2 < 100000:
                    # histograms.vtxRefPhiVsLhPhi.Fill(q2*phiRef2/arc2, q2*phiLH2/arc2)
                    histograms.vtxRefPhiVsArcPhi.Fill(q2*phiRef2/arc2, arc2)
            # if not (leadingP.Pt() > 1.5 and leadingTrack.getNdf() < 70):
            if True:
                histograms.vtxTrksCentresDist.Fill(centresDist)
                histograms.cosOpenAngleVtx_reco.Fill(cosOpenAngle)
                histograms.vtxPt1vsPt2.Fill(leadingP.Pt(),secondP.Pt())
                histograms.vtxP1vsP2.Fill(leadingP.Mag(),secondP.Mag())
                histograms.vtxPz1vsPz2.Fill(leadingP.Pz(),secondP.Pz())
                histograms.vtxOmega1vsOmega2.Fill(abs(omega1),abs(omega2))
                histograms.vtxTrksNdf1VsNdf2.Fill(leadingTrack.getNdf(),secondTrack.getNdf())
                histograms.vtxTrksPtVsNdf.Fill(leadingP.Pt(),leadingTrack.getNdf())
                histograms.vtxTrksPtVsNdf.Fill(secondP.Pt(),secondTrack.getNdf())
                histograms.vtxTrksD01VsD02.Fill(leadingTrack.getD0(),secondTrack.getD0())
                histograms.vtxTrksZ01VsZ02.Fill(leadingTrack.getZ0(),secondTrack.getZ0())
                # if 550 * math.log10(distRef12) + 550 < centresDist:
                # if 2.2 * distRef12 + 900 < centresDist:
                histograms.vtxTrksCentresRefPointsDist.Fill(distRef12, centresDist)
                histograms.vtxTrksCentresRefPointsDist_scaled.Fill(distRef12, centresDistScaled)
                histograms.vtxRefPointDistsVsCentres.Fill(vtx_refpoint1,centresDist)
                histograms.vtxRefPointDistsVsCentres.Fill(vtx_refpoint2,centresDist)
                # if utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 330 and utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) < 1800:
                # if p_vtx.Pt() > 1.9:
                    # if centresRefsDists < -2000:
                histograms.corrCentresRefDistsVsPt.Fill(centresRefsDists, p_vtx.Pt())
                # if logCentresRefsDists < -2000:
                # 	if p_vtx.Pt() > 1.9:
                histograms.corrLogCentresRefDistsVsPt.Fill(logCentresRefsDists, p_vtx.Pt())
                if utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 350:
                    histograms.corrCentresRefDistsVsR.Fill(centresRefsDists, utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]))
                histograms.vtxTrksOmega.Fill(abs(omega1))
                histograms.vtxTrksOmega.Fill(abs(omega2))
                histograms.vtxTrksCurvRatio.Fill(ratio)
                histograms.cosOpenAngleVsCurvRatio.Fill(cosOpenAngle, ratio)
                histograms.vtxTrksChi2.Fill( trksVec[iTrk].getChi2() )
                histograms.vtxTrksChi2.Fill( trksVec[iTrk+1].getChi2() )
                histograms.vtxTrksNdf.Fill( trksVec[iTrk].getNdf() )
                histograms.vtxTrksNdf.Fill( trksVec[iTrk+1].getNdf() )
                if trksVec[iTrk].getNdf() > 0:
                    histograms.vtxTrksChi2Ndf.Fill( trksVec[iTrk].getChi2() / trksVec[iTrk].getNdf() )
                if trksVec[iTrk+1].getNdf() > 0:
                    histograms.vtxTrksChi2Ndf.Fill( trksVec[iTrk+1].getChi2() / trksVec[iTrk+1].getNdf() )
            histograms.vtxTrksD0AtIP.Fill( trksVec[iTrk].getD0() )
            histograms.vtxTrksD0AtIP.Fill( trksVec[iTrk+1].getD0() )
            histograms.vtxTrksZ0AtIP.Fill( trksVec[iTrk].getZ0() )
            histograms.vtxTrksZ0AtIP.Fill( trksVec[iTrk+1].getZ0() )
            histograms.vtxTrksD0.Fill( trkSt1.getD0() )
            histograms.vtxTrksD0.Fill( trkSt2.getD0() )
            histograms.vtxTrksZ0.Fill( trkSt1.getZ0() )
            histograms.vtxTrksZ0.Fill( trkSt2.getZ0() )

        iTrk+=2

def get_closest_dists(trackCollection):
	closest_dists = []
	for idx in range( trackCollection.getNumberOfElements() - 1 ):
		firstTrack = trackCollection.getElementAt(idx)
		for jdx in range(idx+1, trackCollection.getNumberOfElements() ):
			secondTrack = trackCollection.getElementAt(jdx)

			ftFH = firstTrack.getTrackState(2).getReferencePoint()
			ftLH = firstTrack.getTrackState(3).getReferencePoint()
			stFH = secondTrack.getTrackState(2).getReferencePoint()
			stLH = secondTrack.getTrackState(3).getReferencePoint()

			dists = []
			dists.append( utils.spacialDistance(ftFH,stFH) )
			dists.append( utils.spacialDistance(ftFH,stLH) )
			dists.append( utils.spacialDistance(ftLH,stFH) )
			dists.append( utils.spacialDistance(ftLH,stLH) )
			closest_dists.append( min(dists) )
			
	return closest_dists