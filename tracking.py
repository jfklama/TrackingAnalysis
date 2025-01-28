import math
from pyLCIO import IOIMPL, EVENT, IMPL, ROOT, IOException 
import utils

def analyse_tracks(mcCollection, trackCollection, \
# def analyse_tracks(mcCollection, trackCollection, \
				llp_pdg, pdg1, pdg2, gen_status_llp, gen_status, \
				track_col, fill_hists, counter, histograms, cat, iEvt):
	vMatchingTracks = []
	matchingTracksInEvent = 0
	matchingTracks = 0
	p_llp = 0
	mcllpP4 = ROOT.TLorentzVector()
	trueR = 0
	mclepP4s = []
	
	for particle in mcCollection:
		# if matchingTracksInEvent > 1:
		# 	break
		
		if(abs(particle.getPDG()) == llp_pdg and particle.getGeneratorStatus() == gen_status_llp):
			e = particle.getEnergy()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p_llp = utils.R3(px,py,pz)
			mcllpP4.SetPxPyPzE(px,py,pz,e)
			llp_endpoint = particle.getEndpoint()
			trueR = utils.R(llp_endpoint[0], llp_endpoint[1])
			trueZ = llp_endpoint[2]
			trueL = utils.spacialDistance(llp_endpoint, [0,0,0])

		mclepP4 = ROOT.TLorentzVector()
		if( particle.getGeneratorStatus() == gen_status):
			if abs(particle.getPDG()) != pdg1 and abs(particle.getPDG()) != pdg2:
				continue
			if utils.R(particle.getVertex()[0],particle.getVertex()[1]) != trueR:
				continue
			vtx_mc = particle.getVertex()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = utils.R3(px,py,pz)
			e = particle.getEnergy()
			mclepP4.SetPxPyPzE(px,py,pz,e)
			cosTheta = abs( math.cos(mclepP4.Theta()) )
			deltaTheta = abs( mcllpP4.Theta() - mclepP4.Theta() )
			deltaPhi = utils.phiDistance( mcllpP4.Phi(), mclepP4.Phi() )
			deltaAlpha = mcllpP4.Angle( mclepP4.Vect() )
			mclepP4s.append(mclepP4)

			if cat == 'sm' and abs(mclepP4.Eta()) > 2.4:
				continue

			# cut = mclepP4.Pt() > 10 and (mclepP4.Theta() > 1.4 and mclepP4.Theta() < 1.75) \
					# and (mclepP4.Phi() > 1.4 and mclepP4.Phi() < 1.75)
			cut = True

			# for trackRel in trackMCLinkCollection:
			for track in trackCollection:

				# track = trackRel.getFrom()

				# ts_IP = track.getTrackState(1) # 1 = atIP, 2 = atFirstHit, 3 = atLastHit
				# ts_IP = utils.fixTrackStateDirection( track, ts_IP.getLocation() )
				# momentum = utils.getTrackMomentum(ts_IP)
				# p_ip = ROOT.TVector3()
				# p_ip.SetXYZ(momentum[0],momentum[1],momentum[2])
				# # dist_fh = utils.angularDistance( p_fh.Theta(), mclepP4.Theta(), p_fh.Phi(), mclepP4.Phi()  )
				# dist_ip = p_ip.Angle(mclepP4.Vect())

				ts_FirstHit = track.getTrackState(2) # 1 = atIP, 2 = atFirstHit, 3 = atLastHit
				if ts_FirstHit.getOmega() == 0:
					continue
				ts_FirstHit = utils.fixTrackStateDirection( track, ts_FirstHit.getLocation() )
				momentum = utils.getTrackMomentum(ts_FirstHit)
				p_fh = ROOT.TVector3()
				p_fh.SetXYZ(momentum[0],momentum[1],momentum[2])
				# dist_fh = utils.angularDistance( p_fh.Theta(), mclepP4.Theta(), p_fh.Phi(), mclepP4.Phi()  )
				dist_fh = p_fh.Angle(mclepP4.Vect())
				vtx_fh = ts_FirstHit.getReferencePoint()
				vtx_dist_fh = utils.spacialDistance(vtx_mc, vtx_fh)

				ts_LastHit = track.getTrackState(3) # 1 = atIP, 2 = atFirstHit, 3 = atLastHit
				if ts_LastHit.getOmega() == 0:
					continue
				ts_LastHit = utils.fixTrackStateDirection( track, ts_LastHit.getLocation() )
				momentum = utils.getTrackMomentum(ts_LastHit)
				p_lh = ROOT.TVector3()
				p_lh.SetXYZ(momentum[0],momentum[1],momentum[2])
				# dist_lh = utils.angularDistance( p_lh.Theta(), mclepP4.Theta(), p_lh.Phi(), mclepP4.Phi()  )
				dist_lh = p_lh.Angle(mclepP4.Vect())
				vtx_lh = ts_LastHit.getReferencePoint()
				vtx_dist_lh = utils.spacialDistance(vtx_mc, vtx_lh)

				ang_dist = dist_fh
				vtx_dist = vtx_dist_fh
				trackState = ts_FirstHit
				if vtx_dist_lh < vtx_dist_fh:
					ang_dist = dist_lh
					vtx_dist = vtx_dist_lh
					trackState = ts_LastHit
				# if trueR < 300:
				# 	ang_dist = dist_ip
				# 	trackState = ts_IP
				# histograms.h_angDist.Fill(ang_dist)

				# hits = track.getTrackerHits()
				# if hits.size() < 4:
				# 	print('')
				# 	print(iEvt, " TRACK CONTAINS < 4 HITS")
				# 	print('')

				if utils.getTrackMomentum(trackState) not in vMatchingTracks:
					vMatchingTracks.append( utils.getTrackMomentum(trackState) )

				# if trackRel.getTo() == particle \
				# and trackRel.getTo().getCharge() * utils.getTrackCharge(trackState) > 0 \
				# and ang_dist < 0.2:
				if particle.getCharge() * utils.getTrackCharge(trackState) > 0 \
				and ang_dist < 0.1:
				#and vtx_dist < 100.:

					if fill_hists:
						if cut:
							histograms.momTrack_reco.Fill(p)
							histograms.ptTrack_reco.Fill(mclepP4.Pt())
							histograms.momLLP_reco.Fill( p_llp )
							histograms.thetaTrack_reco.Fill(mclepP4.Theta())
							histograms.DeltaThetaTrackLLP_reco.Fill( deltaTheta )
							histograms.DeltaPhiTrackLLP_reco.Fill( abs(deltaPhi) )
							histograms.DeltaAlphaTrackLLP_reco.Fill( abs(deltaAlpha) )
							histograms.trackEffiPvsTheta_reco.Fill(mclepP4.Theta(),p)
							histograms.trackEffVsR_reco.Fill(trueR)
							histograms.trackEffVsL_reco.Fill(trueL)
							histograms.trackEffRVsZ_reco.Fill(trueZ,trueR)

							# if trueR > 500:
							# 	print(iEvt, trueR)

							# print(iEvt, "track:", trueR, trueZ)

					counter.increment("_nMatchingTracks", \
										track_col == 'MarlinTrk', 1)

					counter.increment("_nMatchingTracksReffited", \
										track_col == 'Refitted', 1)
					
					matchingTracksInEvent+=1
					break

				# if trueR < 100:
				# 	print(ang_dist)
	# if len(mclepP4s) > 1:
	# 	if matchingTracksInEvent < 2 and \
	# 		((mclepP4s[0].Theta() < 0.6 or mclepP4s[0].Theta() > math.pi - 0.6) or (mclepP4s[1].Theta() < 0.6 or mclepP4s[1].Theta() > math.pi - 0.6)):
	# 		print (iEvt, matchingTracksInEvent, trueR, trueZ, mclepP4s[0].Theta(), mclepP4s[1].Theta())
	# 		print(mclepP4s[0].Px(),mclepP4s[0].Py(),mclepP4s[0].Pz() )
	# 		print(mclepP4s[1].Px(),mclepP4s[1].Py(),mclepP4s[1].Pz() )

	# if len(mclepP4s) > matchingTracksInEvent and trueR > 330:
	# 		print(track_col)
	# 		print (iEvt, matchingTracksInEvent, trueR, trueZ, mclepP4s[0].Theta(), mclepP4s[1].Theta())
	# 		print(mclepP4s[0].Px(),mclepP4s[0].Py(),mclepP4s[0].Pz() )
	# 		print(mclepP4s[1].Px(),mclepP4s[1].Py(),mclepP4s[1].Pz() )
	# 		print('reco tracks:')
	# 		for mom in vMatchingTracks:
	# 			print(mom[0], mom[1], mom[2])
	# 		print('')
