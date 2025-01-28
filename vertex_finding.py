import math
from pyLCIO import IOIMPL, EVENT, UTIL, IMPL, ROOT, IOException 
import utils, constants

def vertex_analysis(event, cat, vertexCollection, trackCollection, mcCollection, \
                    vertexTracksRelations, refitTrackToMCLinkCollection, \
                    llp_endpoint, mclep1, mclep2, weights, \
                    counter, histograms, iEvt):
    dist0 = 50000
    helixDist0 = 50000
    bestVtx = IMPL.VertexImpl()
    closest_dists = get_closest_dists(trackCollection)
    trueR = utils.R(llp_endpoint[0], llp_endpoint[1])
    trueZ = llp_endpoint[2]
    trueL = utils.spacialDistance(llp_endpoint, [0,0,0])
    in_acceptance = (trueR <= 1808) and abs(trueZ) <= 2350
    in_TPC = in_acceptance and trueR >= 329
    cosOpenAngleTruth = math.cos(mclep1.Angle(mclep2.Vect()))

    vtx_conversions = []
    n_vtx_conversions = match_converted_photons(mcCollection, vertexCollection, \
                                                vtx_conversions)
    
    if vertexCollection.getNumberOfElements() == 0:
        if (trueR > 330 and trueR < 1300) and (trueZ > -1000 and trueZ < 1000):
            print('no good vtx in ev. ' + str(iEvt) + ' with ' + str(vertexCollection.getNumberOfElements()) + ' vertices')
            print('true vertex at: ', llp_endpoint[0],llp_endpoint[1],llp_endpoint[2])
            print('')

    for vtx in vertexCollection:

        pos = vtx.getPosition()
        counter.increment("_nSecVerticesTPC", \
                                (utils.R(pos[0], pos[1]) >= 329) \
                            and (utils.R(pos[0], pos[1]) <= 1808) \
                            and abs(pos[2]) <= 2350, \
                                1)

        # if vertex comes from photon conversion we can continue
        if vtx in vtx_conversions:
            counter.increment("_nPhotonConversionsTPC", \
                                (utils.R(pos[0], pos[1]) >= 329) \
                            and (utils.R(pos[0], pos[1]) <= 1808) \
                            and abs(pos[2]) <= 2350, \
                                1)
            # print('found vertex from conversion in evt ' + str(iEvt))
            # continue

        helixDist = vtx.getParameters()[0]
        histograms.helixDistances.Fill(helixDist )
        if len(closest_dists) > 0:
            closest_dist = min(closest_dists)
            histograms.hitsDistances.Fill( closest_dist )
        
        histograms.vtxVsR_reco.Fill( utils.R(pos[0],pos[1]) )
        histograms.vtxVsL_reco.Fill( utils.R3(pos[0],pos[1],pos[2]) )
        dist = utils.spacialDistance(llp_endpoint, pos )
        #vtxMisdistances.Fill(dist)
        # if dist < 30:
        #     print (iEvt, dist0, dist, helixDist)
        # if trueR > 329:
        #     print (iEvt, 'vertex:', pos[0],pos[1],pos[2], "dist:", dist)

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
        # histograms.vtxEffRVsZ_reco.Fill(trueZ,trueR)
        histograms.vtxEffVsL_reco.Fill(trueL)
        histograms.cosOpenAngleVtx_reco.Fill(cosOpenAngleTruth)
        # 	_matchingVerticesTPC+=1
        # _nFakeVertices += vertexCollection.getNumberOfElements() - 1

    # elif (trueR > 330 and trueR < 1300) and (trueZ > -1000 and trueZ < 1000):
    #     print('no good vtx in ev. ' + str(iEvt) + ' with ' + str(vertexCollection.getNumberOfElements()) + ' vertices')
    #     print('true vertex at: ', llp_endpoint[0],llp_endpoint[1],llp_endpoint[2])
    #     if vertexCollection.getNumberOfElements() > 0:
    #         for fakeVtx in vertexCollection:
    #             print('found vertex at: ', fakeVtx.getPosition()[0],fakeVtx.getPosition()[1],fakeVtx.getPosition()[2])
    #     print('')
        # _nFakeVertices += vertexCollection.getNumberOfElements()

    # if dist0 < 30 and vertexCollection.getNumberOfElements() - n_vtx_conversions - 1 < 0:


    trksVec = []
    verticesVec = []
    for rel in vertexTracksRelations:
        # if trueR < 329:
        #     print (iEvt, 'vertex:', rel.getFrom().getPosition()[0],rel.getFrom().getPosition()[1],rel.getFrom().getPosition()[2])
        # if rel.getFrom() == bestVtx and dist0 < 30:
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

    vtx_conversions = []
    vtx_v0s = []
    n_vtx_conversions = match_converted_photons(mcCollection, vertexCollection, \
                                        vtx_conversions)
    counter.increment("_nPhotonConversions", len(vtx_conversions) > 0, n_vtx_conversions)

    # if 'V0Vertices' in event.getCollectionNames():
    true_v0s = match_v0s(mcCollection, vertexCollection, \
                                            vtx_v0s)
    n_v0s = len(true_v0s)

    iTrk=0
    nVtxInEvent = 0
    nMatchedVtxInEvent = 0
    nRecoVtxInTPC = 0
    while iTrk < len(trksVec):

        nVtxInEvent += 1

        trkSt1, trkSt2 = utils.findClosestLegalTrackStates(trksVec[iTrk],trksVec[iTrk+1],verticesVec[iTrk])
        vtx_cand = verticesVec[iTrk]
        pos = vtx_cand.getPosition()
        dist = utils.spacialDistance(llp_endpoint, pos)


            # if  utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 400 and utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) < 1770:
        if utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 330 and utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) < 1770:
        # if True:
        #     and abs(vtx_cand.getPosition()[2]) < 400:
            
            if utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) < 400 or utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 1500: 
                iTrk+=2
                continue


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

            # print (iEvt, llp_endpoint[0],llp_endpoint[1], llp_endpoint[2])
            # print (iEvt, vtx_cand.getPosition()[0],vtx_cand.getPosition()[1], vtx_cand.getPosition()[2], vtx_cand.getParameters()[0])
            # print (trkSt1.getD0(), trkSt1.getZ0(), omega1, trkSt1.getTanLambda(), trkSt1.getPhi())
            # print (trkSt2.getD0(), trkSt2.getZ0(), omega2, trkSt2.getTanLambda(), trkSt2.getPhi())
            # print (trkSt1.getReferencePoint()[0],trkSt1.getReferencePoint()[1],trkSt1.getReferencePoint()[2])
            # print (trkSt2.getReferencePoint()[0],trkSt2.getReferencePoint()[1],trkSt2.getReferencePoint()[2])


            isV0 = 0
            if 'V0Vertices' in event.getCollectionNames():
                v0Collection = event.getCollection('V0Vertices')
                for v0 in v0Collection:
                    v0ToVtx = utils.spacialDistance(v0.getPosition(),vtx_cand.getPosition())
                    if v0ToVtx < 30:
                        # print i, 'v0: ', v0.getPosition()[0],v0.getPosition()[1],v0.getPosition()[2]
                        # print i, 'llp: ', vtx_cand.getPosition()[0],vtx_cand.getPosition()[1],vtx_cand.getPosition()[2]
                        isV0=1
            if isV0 == 1:
                iTrk+=2
                continue


            # minDistToV0 = 1.e9
            # for mc_v0 in true_v0s:
            #     distToV0 = utils.spacialDistance(mc_v0.getEndpoint(),pos)
            #     if distToV0 < minDistToV0:
            #         minDistToV0 = distToV0
            # if minDistToV0 < 300:
            #     print('v0 at ', [mc_v0.getEndpoint()[xyz] for xyz in range(3)], \
            #             ' close to vtx at ', [pos[xyz] for xyz in range(3)])
                # vtxDistToV0.Fill(minDistToV0)

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

            # for pfo in pfoCollection:
            #     # tracks = pfo.getTracks()
            #     pfo_mom = ROOT.TVector3()
            #     pfo_mom.SetXYZ(pfo.getMomentum()[0],pfo.getMomentum()[1],pfo.getMomentum()[2])
            #     print(pfo_mom[0],pfo_mom[1],pfo_mom[2], 'PID: ',[id.getPDG() for id in pfo.getParticleIDs()])
            #     _PIDMethod = "TMVA_BDT_MC_12bins.S"
            #     PIDHan = UTIL.PIDHandler(pfoCollection)
            #     PDG = PIDHan.getParticleID(pfo,PIDHan.getAlgorithmID(_PIDMethod)).getPDG()
            #     print(PDG)
            #     if utils.spacialDistance(pfo.getReferencePoint(), trkSt1.getReferencePoint() ) < 1.e-3: 
            #         # print(pfo_mom[0],pfo_mom[1],pfo_mom[2], 'PID: ',[id.getPDG() for id in pfo.getParticleIDs()])
            #         print("dist1: ", utils.spacialDistance(pfo.getReferencePoint(), trkSt1.getReferencePoint() ))
            #         print('angle1:',pfo_mom.Angle(momenta[0]))
            #     if utils.spacialDistance(pfo.getReferencePoint(), trkSt2.getReferencePoint() ) < 1.e-3:
            #         # print(pfo_mom[0],pfo_mom[1],pfo_mom[2], 'PID: ',[id.getPDG() for id in pfo.getParticleIDs()])
            #         print("dist2: ", utils.spacialDistance(pfo.getReferencePoint(), trkSt2.getReferencePoint() ))
            #         print('angle2:',pfo_mom.Angle(momenta[1]))
            # print(' ')
                
            # check if vertex was not caused by interaction with detector material
            # and calculate energy inside the cone around the vertex
            secondary_interaction = 0
            coneP = 0
            for trk in trackCollection:
                startpoint = trk.getTrackState(2).getReferencePoint() # 2 = first hit
                endpoint = trk.getTrackState(3).getReferencePoint() # 3 = last hit
                trk_fhit = utils.R(startpoint[0],startpoint[1])
                p_trk = utils.getTrackMomentum(trk)
                p3_trk = ROOT.TVector3()
                p3_trk.SetXYZ(p_trk[0],p_trk[1],p_trk[2])

                # if trk.id() == trksVec[iTrk].id() or trk.id() == trksVec[iTrk+1].id():
                if trk in [trksVec[iTrk],trksVec[iTrk+1]]: # skip track if it's coming out of this vtx
                    # print('found track assigned to vertex ending at: ', [endpoint[edx] for edx in range(3)])
                    continue
                
                if cat == 'trsm':  # reject event if a hard track not from a vtx from IP is found
                    if (abs(trk.getD0()) < 10 or abs(trk.getZ0()) < 10) and trk not in trksVec and p3_trk.Mag() > 2 \
                        and ((trk_fhit < 155 and abs(startpoint[2]) > 215) or trk_fhit < 20):
                        # print(iEvt, [startpoint[x] for x in range(3)], [endpoint[x] for x in range(3)], trk.getD0(), len(trksVec))
                        return False

                distance = []
                time = utils.getDistanceToPoint(trk,pos,distance)
                # if a track passes through (or very close to) the vtx, it is probably the source of this vtx
                if distance[2] < 30 and trk.getD0() < 10 and math.cos(p3_trk.Angle(p_vtx)) > 0.:
                    # print('found a track ending at ', [endpoint[edx] for edx in range(3)], \
                    # 	' passing by vertex ', [pos[ijk] for ijk in range(3)], \
                    # 	' at distance: ', [distance[xyz] for xyz in range(3)])
                    # print('time: ', time)
                    secondary_interaction = 1
                    # break
                # skip track if it starts or ends just next to the vtx (probably comes from the same vtx)
                if utils.spacialDistance(endpoint, pos) < 30 or utils.spacialDistance(startpoint, pos) < 30:
                    continue
                    # print('found a track ending near vertex at: ', [endpoint[edx] for edx in range(3)])
                    # secondary_interaction = 1
                    # break
                # the rest contributes to isolation, if the angular condition is fulfilled
                if math.cos(p3_trk.Angle(p_vtx)) > 0.98:
                    coneP += p3_trk.Mag()

            if secondary_interaction > 0:
                # print(i, "cut on secondary interaction")
                iTrk+=2
                continue
            # if coneP > 1:
            # 	# print(i, "cut on isolation")
            # 	iTrk+=2
            # 	continue
            
            
            if cat == 'trsm': # additional cut on vtx_pT for higgs decays
                if p_vtx.Pt() < 10:
                    # print(i, "cut on isolation")
                    iTrk+=2
                    continue


            # print(p_vtx.Mag(),coneP)
            histograms.vtxIso.Fill(p_vtx.Mag(),coneP)
            histograms.vtxConeP.Fill(coneP)

            # print ''
            if dist0 < 30:
                histograms.vtxP.Fill( p_vtx.Mag() )
                histograms.vtxPt.Fill( p_vtx.Pt() )
                histograms.vtxTheta.Fill( p_vtx.Theta() )
                histograms.vtxPhi.Fill( p_vtx.Phi() )
                # print('True opening angle for a good vtx: ' + str(cosOpenAngleTruth))

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

            # if distRef12 > 50:
            #     print (iEvt, distRef12, pos[0],pos[1],pos[2], utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]))
            #     print (iEvt, vtx_cand.getPosition()[0], vtx_cand.getPosition()[1], vtx_cand.getPosition()[2], vtx_cand.getParameters()[0], distRef12)
            #     print (iEvt, trkSt1.getReferencePoint()[0], trkSt1.getReferencePoint()[1], trkSt1.getReferencePoint()[2],trksVec[iTrk].getNdf())
            #     print (iEvt, trkSt2.getReferencePoint()[0], trkSt2.getReferencePoint()[1], trkSt2.getReferencePoint()[2],trksVec[iTrk+1].getNdf())
            #     print ('')

            cosOpenAngle = math.cos(momenta[0].Angle(momenta[1]))
            histograms.cosOpenAngleVsRefPtDist_true.Fill(cosOpenAngleTruth,distRef12)
            histograms.cosOpenAngleVsRefPtDist_reco.Fill(cosOpenAngle,distRef12)
            

            vtx_refpoint1 = utils.spacialDistance( trkSt1.getReferencePoint(),vtx_cand.getPosition() )
            lh_loc1 = 3 if trkSt1.getLocation() == 2 else 2
            lastHit1 = trksVec[iTrk].getTrackState(lh_loc1)
            vtx_lh1 = utils.spacialDistance( lastHit1.getReferencePoint(), vtx_cand.getPosition() )
            refpoint1_lh1 = utils.spacialDistance( trkSt1.getReferencePoint(), lastHit1.getReferencePoint() )
            rFirstHit1 = utils.R(trkSt1.getReferencePoint()[0],trkSt1.getReferencePoint()[1])
            rLastHit1 = utils.R(lastHit1.getReferencePoint()[0],lastHit1.getReferencePoint()[1])

            vtx_refpoint2 = utils.spacialDistance( trkSt2.getReferencePoint(),vtx_cand.getPosition() )
            lh_loc2 = 3 if trkSt2.getLocation() == 2 else 2
            lastHit2 = trksVec[iTrk+1].getTrackState(lh_loc2)
            vtx_lh2 = utils.spacialDistance( lastHit2.getReferencePoint(), vtx_cand.getPosition() )
            refpoint2_lh2 = utils.spacialDistance( trkSt2.getReferencePoint(), lastHit2.getReferencePoint() )
            rFirstHit2 = utils.R(trkSt2.getReferencePoint()[0],trkSt2.getReferencePoint()[1])
            rLastHit2 = utils.R(lastHit2.getReferencePoint()[0],lastHit2.getReferencePoint()[1])

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
            if distRef12 > 0:
                logCentresRefsDists = 300 * math.log10(distRef12) - centresDist
            else:
                logCentresRefsDists = 0
            

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
                # histograms.cosOpenAngleVtx_reco.Fill(cosOpenAngle)
                histograms.vtxPt1vsPt2.Fill(leadingP.Pt(),secondP.Pt())
                histograms.vtxP1vsP2.Fill(leadingP.Mag(),secondP.Mag())
                histograms.vtxPz1vsPz2.Fill(leadingP.Pz(),secondP.Pz())
                histograms.vtxOmega1vsOmega2.Fill(abs(omega1),abs(omega2))
                histograms.vtxTrksNdf1VsNdf2.Fill(leadingTrack.getNdf(),secondTrack.getNdf())
                histograms.vtxTrksPtVsNdf.Fill(leadingP.Pt(),leadingTrack.getNdf())
                histograms.vtxTrksPtVsNdf.Fill(secondP.Pt(),secondTrack.getNdf())
                histograms.vtxTrksD01VsD02.Fill(math.fabs(leadingTrack.getD0()),math.fabs(secondTrack.getD0()))
                histograms.vtxTrksZ01VsZ02.Fill(math.fabs(leadingTrack.getZ0()),math.fabs(secondTrack.getZ0()))
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
                histograms.corrCentresRefDistsVsPtLog.Fill(centresRefsDists, p_vtx.Pt())
                # if logCentresRefsDists < -2000:
                # 	if p_vtx.Pt() > 1.9:
                histograms.corrLogCentresRefDistsVsPt.Fill(logCentresRefsDists, p_vtx.Pt())
                if utils.R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 329:
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
            histograms.vtxTrksD0AtIP.Fill( math.fabs(trksVec[iTrk].getD0()) )
            histograms.vtxTrksD0AtIP.Fill( math.fabs(trksVec[iTrk+1].getD0()) )
            histograms.vtxTrksZ0AtIP.Fill( math.fabs(trksVec[iTrk].getZ0()) )
            histograms.vtxTrksZ0AtIP.Fill( math.fabs(trksVec[iTrk+1].getZ0()) )

            histograms.vtxTrksD0.Fill( math.fabs(trkSt1.getD0()) )
            histograms.vtxTrksD0.Fill( math.fabs(trkSt2.getD0()) )
            histograms.vtxTrksZ0.Fill( math.fabs(trkSt1.getZ0()) )
            histograms.vtxTrksZ0.Fill( math.fabs(trkSt2.getZ0()) )

            if (rFirstHit1 > rLastHit1 and (math.fabs(trksVec[iTrk].getD0()) < 50 and math.fabs(trksVec[iTrk].getZ0()) < 50)) \
                or (rFirstHit2 > rLastHit2 and (math.fabs(trksVec[iTrk+1].getD0()) < 50 and math.fabs(trksVec[iTrk+1].getZ0()) < 50)):
                    print("cut on one of the tracks coming from IP")
                    iTrk+=2
                    continue

            ee_mass = utils.get_inv_mass( \
                momenta[0],constants.ELECTRON_MASS,momenta[1],constants.ELECTRON_MASS)
            pipi_mass = utils.get_inv_mass( \
                momenta[0],constants.PION_MASS,momenta[1],constants.PION_MASS)
            proton_pion_mass = utils.get_inv_mass( \
                momenta[0],constants.PROTON_MASS,momenta[1],constants.PION_MASS)
            pion_proton_mass = utils.get_inv_mass( \
                momenta[0],constants.PION_MASS,momenta[1],constants.PROTON_MASS)
            
            if ee_mass < constants.GAMMA_MASS + 0.15:
                counter.increment("_eeMassCut", True, 1)
                iTrk+=2
                continue
            if pipi_mass > constants.K0S_MASS - 0.05 and pipi_mass < constants.K0S_MASS + 0.05:
                counter.increment("_pipiMassCut", True, 1)
                iTrk+=2
                continue
            if proton_pion_mass > constants.LAMBDA_MASS - 0.05 and proton_pion_mass < constants.LAMBDA_MASS + 0.05:
                counter.increment("_lambdaMassCut", True, 1)
                iTrk+=2
                continue
            if pion_proton_mass > constants.LAMBDA_MASS - 0.05 and pion_proton_mass < constants.LAMBDA_MASS + 0.05:
                counter.increment("_lambdaMassCut", True, 1)
                iTrk+=2
                continue

            # # # cuts excluding ma = 300 MeV scenario
            # if ee_mass < constants.GAMMA_MASS + 0.7:
            #     counter.increment("_eeMassCut", True, 1)
            #     iTrk+=2
            #     continue
            # if pipi_mass < 0.7:
            #     counter.increment("_pipiMassCut", True, 1)
            #     iTrk+=2
            #     continue

            # if proton_pion_mass > constants.LAMBDA_MASS - 0.02 and proton_pion_mass < constants.LAMBDA_MASS + 0.02:
            #     counter.increment("_lambdaMassCut", True, 1)
            #     iTrk+=2
            #     continue
            # if pion_proton_mass > constants.LAMBDA_MASS - 0.02 and pion_proton_mass < constants.LAMBDA_MASS + 0.02:
            #     counter.increment("_lambdaMassCut", True, 1)
            #     iTrk+=2
            #     continue

            histograms.ee_hypo_mass.Fill( ee_mass )
            histograms.pipi_hypo_mass.Fill( pipi_mass )
            histograms.lambda_mass_diff.Fill( \
							min(math.fabs(proton_pion_mass-constants.LAMBDA_MASS), math.fabs(pion_proton_mass-constants.LAMBDA_MASS)) \
						)
            histograms.lambda_hypo_mass.Fill( proton_pion_mass, pion_proton_mass )

            histograms.ee_mass_vs_R.Fill( ee_mass, utils.R(pos[0],pos[1]) )
            histograms.pipi_mass_vs_R.Fill( pipi_mass, utils.R(pos[0],pos[1]) )

            histograms.vtxEffVsR_reco.Fill( utils.R(pos[0],pos[1]) )

            reco_in_TPC =  (utils.R(pos[0], pos[1]) >= 329) \
                            and (utils.R(pos[0], pos[1]) <= 1808) \
                            and abs(pos[2]) <= 2350
                                
            # FIXME: fake vertices for now don't work for more than 1 true vtx in evt
            if vtx_cand != bestVtx or dist >= 30: 
                counter.increment("_nFakeVertices", True, 1)
                counter.increment("_nFakeVerticesTPC", reco_in_TPC, 1)

            # counter.increment("_nEvSecVerticesTPC", reco_in_TPC and nVtxInEvent < 2, 1)

            if reco_in_TPC and nVtxInEvent < 2:
                nRecoVtxInTPC += 1

            if vtx_cand == bestVtx:
                # count passing vertices
                counter.increment("_nMatchingVertices", dist < 30, 1)
                counter.increment("_matchingVerticesTPC", dist < 30 and reco_in_TPC, 1)

                # # count weighted passing vertices
                # if nMatchedVtxInEvent < 1: # (but only once)
                #     counter.increment("_matchingWeightedEventsTPC", dist < 30 and reco_in_TPC, weights)
                #     counter.increment("_matchingWeightedEventsTPCErr", dist < 30 and reco_in_TPC, weights**2)

                if dist < 30:
                    histograms.vtxEffRVsZ_reco.Fill(trueZ,trueR)
                    histograms.vtxEffVsPt_reco.Fill((mclep1+mclep2).Pt())
                    histograms.vtxEffVsP_reco.Fill((mclep1+mclep2).P())
                    # if (mclep1+mclep2).Pt() < 40:
                    #     print(iEvt, (mclep1+mclep2).Pt())
                    # print(iEvt, llp_endpoint[0], llp_endpoint[1], llp_endpoint[2], 'matched to ', pos[0], pos[1],  pos[2])
                    
                    if reco_in_TPC:
                        nMatchedVtxInEvent += 1

                # return dist < 30 and reco_in_TPC, reco_in_TPC

            # print(iEvt, pos, "nEvVtx: ", vertexCollection.getNumberOfElements(), " is best vtx: ", vtx_cand == bestVtx)

        iTrk+=2

    if nMatchedVtxInEvent > 0:
        return True
    else:
        return False

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

def match_converted_photons(mc_collection, vtx_collection, reco_conversions):
    
    n_vtx_conversions = 0
    vtx_conversions = []
    for particle in mc_collection:
        if particle.getPDG() == 22:
            dec_prod = particle.getDaughters()
            if len(dec_prod) == 2:
                if abs(dec_prod[0].getPDG()) == 11 and abs(dec_prod[1].getPDG()) == 11:
                    vtx_conversions.append(particle.getEndpoint())
                elif abs(dec_prod[0].getPDG()) == 211 and abs(dec_prod[1].getPDG()) == 211:
                    vtx_conversions.append(particle.getEndpoint())

    for vtx in vtx_collection:
        reco_pos = vtx.getPosition()
        for photon_vtx in vtx_conversions:
            if utils.spacialDistance(photon_vtx, reco_pos) < 30:
                n_vtx_conversions+=1
                reco_conversions.append(vtx)

    return n_vtx_conversions

def match_v0s(mc_collection, vtx_collection, reco_v0s):
    
    n_vtx_v0s = 0
    vtx_v0s = []
    for particle in mc_collection:
        if math.fabs(particle.getPDG()) in [310,311,130,3122]:
            dec_prod = particle.getDaughters()
            if len(dec_prod) == 2:
                # print(particle.getPDG(), dec_prod[0].getPDG(), dec_prod[1].getPDG())
                # print(particle.getEndpoint()[0],particle.getEndpoint()[1],particle.getEndpoint()[2])
                # print('\n')
                if abs(dec_prod[0].getPDG()) == 11 and abs(dec_prod[1].getPDG()) == 11:
                    vtx_v0s.append(particle)
                elif abs(dec_prod[0].getPDG()) == 211 and abs(dec_prod[1].getPDG()) == 211:
                    vtx_v0s.append(particle)
                elif abs(dec_prod[0].getPDG()) in [2212,211] and abs(dec_prod[1].getPDG()) in [2212,211]:
                    vtx_v0s.append(particle)

    for vtx in vtx_collection:
        reco_pos = vtx.getPosition()
        for v0_vtx in vtx_v0s:
            if utils.spacialDistance(v0_vtx.getEndpoint(), reco_pos) < 30:
                n_vtx_v0s+=1
                reco_v0s.append(vtx)

    return vtx_v0s

    

