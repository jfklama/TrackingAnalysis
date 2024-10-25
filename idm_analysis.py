import sys
import os

import math
import numpy as np
import re

from scipy.stats import expon


from pyLCIO import IOIMPL, EVENT, IMPL, ROOT, IOException
#from ROOT import Math
import ROOT
ROOT.gROOT.Macro( os.path.expanduser( './rootlogon.C' ) )

import utils
import constants
import counters
import histo_class
import tracking
import vertex_finding

def run_analysis(reader, cat, fs, counter, histograms, decay_lengths):
	for event in reader:
		process_event(event, counter, histograms, cat, fs, decay_lengths)
	reader.close()

def process_event(event,counter, histograms, cat, fs, decay_lengths):

	iEvt = event.getEventNumber()

	# if iEvt > 100:
	# 	return

	mcCollection = event.getCollection('MCParticle')
	#tpcCollection = event.getCollection('TPCCollection')
	if cat.find('Pixel') != -1:
		trackCollection = event.getCollection('ClupatraTracks')
		trackToMCLinkCollection = event.getCollection('ClupatraTracksMCTruthLink')
		# trackCollection = event.getCollection('MarlinTrkTracks')
		# trackToMCLinkCollection = event.getCollection('MarlinTrkTracksMCTruthLink')
		refitTrackCollection = event.getCollection('RefittedMarlinTrkTracks')
		refitTrackToMCLinkCollection = event.getCollection('RefittedMarlinTrkTracksMCTruthLink')
		subsetTrackCollection = event.getCollection('SubsetTracks')
	elif cat.find('Silicon') == -1:
		trackCollection = event.getCollection('MarlinTrkTracks')
		trackToMCLinkCollection = event.getCollection('MarlinTrkTracksMCTruthLink')
		refitTrackCollection = event.getCollection('RefittedMarlinTrkTracks')
		refitTrackToMCLinkCollection = event.getCollection('RefittedMarlinTrkTracksMCTruthLink')
		subsetTrackCollection = event.getCollection('SubsetTracks')
		# siTrackCollection = event.getCollection('SiTracks')
		# TPCHitRelations = event.getCollection('TPCTrackerHitRelations')

	else:
		trackCollection = event.getCollection('SiTracks')
		refitTrackCollection = event.getCollection('SiTracks_Refitted')
		trackToMCLinkCollection = event.getCollection('SiTracksMCTruthLink')
		refitTrackToMCLinkCollection = event.getCollection('SiTracksMCTruthLink')
	pfoCollection = event.getCollection('PandoraPFOs')

	if 'LLPVertices' in event.getCollectionNames():
		vertexCollection = event.getCollection('LLPVertices')
		counter.increment("_nSecVertices", True, int(vertexCollection.getNumberOfElements()))

	if iEvt%100==0:
			print (iEvt, str(trackCollection.getNumberOfElements() ) + ' tracks in event')
	
	if subsetTrackCollection.getNumberOfElements() > 0:
		for trk in subsetTrackCollection:
			histograms.siTrackD0.Fill(trk.getD0())
			histograms.siTrackZ0.Fill(trk.getZ0())
	# 	print(iEvt, "Subset tracks: ", str(subsetTrackCollection.getNumberOfElements()))
	# if siTrackCollection.getNumberOfElements() > 0:
	# 	print(iEvt, "Si tracks: ", str(siTrackCollection.getNumberOfElements()))

	if int(pfoCollection.getNumberOfElements()) > 0:
		# _nPFOs += pfoCollection.getNumberOfElements()
		counter.increment("_nPFOs", True, int(pfoCollection.getNumberOfElements()))

	#check if LLP decayed within tracking acceptance and we have 2 truth tracks
	trueRs = []
	trueZs = []
	trueLs = []
	mcllpP4s = []
	trueR = 0
	trueZ = 0
	trueL = 0
	nFSParticle1 = 0
	nFSParticle2 = 0
	llp_endpoint = 0

	for particle in mcCollection:
		if(particle.getPDG() == llp_pdg and particle.getGeneratorStatus() == gen_status_llp):
			if len(particle.getDaughters()) == 1:
				continue
			llp_endpoint = particle.getEndpoint()
			trueR = utils.R(llp_endpoint[0], llp_endpoint[1])
			trueZ = llp_endpoint[2]
			trueL = utils.spacialDistance(llp_endpoint, [0,0,0])
			trueRs.append(trueR)
			trueZs.append(trueZ)
			trueLs.append(trueL)
			E = particle.getEnergy()
			PX = particle.getMomentum()[0]
			PY = particle.getMomentum()[1]
			PZ = particle.getMomentum()[2]
			mcllpP4 = ROOT.TLorentzVector()
			mcllpP4.SetPxPyPzE(PX,PY,PZ,E)
			mcllpP4s.append(mcllpP4)

	nLLPsInEvent = [0] # using a list since it's mutable object
	nLLPsInTPC = [0] # using a list since it's mutable object
	nMatchedLLPsInEvent = [0] # using a list since it's mutable object
	for particle in mcCollection:
		if(particle.getPDG() == llp_pdg and particle.getGeneratorStatus() == gen_status_llp):
			if len(particle.getDaughters()) == 1:
				continue
			llp_endpoint = particle.getEndpoint()
			trueR = utils.R(llp_endpoint[0], llp_endpoint[1])
			trueZ = llp_endpoint[2]
			trueL = utils.spacialDistance(llp_endpoint, [0,0,0])
			e = particle.getEnergy()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = utils.R3(px,py,pz)
			m_llp = particle.getMass()
			mcllpP4 = ROOT.TLorentzVector()
			mcllpP4.SetPxPyPzE(px,py,pz,e)

			process_llp(iEvt, event, counter, histograms, cat, fs, decay_lengths, \
				mcCollection, trackCollection, trackToMCLinkCollection, refitTrackCollection, refitTrackToMCLinkCollection, \
				pfoCollection, nLLPsInEvent, nMatchedLLPsInEvent, nLLPsInTPC, trueRs, trueZs, trueLs, trueR, trueZ, trueL, \
				mcllpP4s, mcllpP4, llp_endpoint, m_llp, nFSParticle1, nFSParticle2)
			
	nLLPsInEvent = [0]
	nLLPsInTPC = [0]
	nMatchedLLPsInEvent = [0]

def process_llp(iEvt, event, counter, histograms, cat, fs, decay_lengths, \
				mcCollection, trackCollection, trackToMCLinkCollection, refitTrackCollection, refitTrackToMCLinkCollection, \
				pfoCollection, nLLPsInEvent, nMatchedLLPsInEvent, nLLPsInTPC, trueRs, trueZs, trueLs, trueR, trueZ, trueL, \
				mcllpP4s, mcllpP4, llp_endpoint, m_llp, nFSParticle1, nFSParticle2):
	
	nVTX = 0
	if 'LLPVertices' in event.getCollectionNames():
		# vertexCollection = event.getCollection('LLPVertices')
		# vertexTracksRelations = event.getCollection('LLPVtxToTracksLink')
		vertexCollection = event.getCollection('LLPVertices')
		vertexTracksRelations = event.getCollection('LLPVtxToTracksLink')
		nVTX = vertexCollection.getNumberOfElements()
		# if infileName.find("idm") == -1:
		# 	_nSecVertices += int(vertexCollection.getNumberOfElements())
		# 	_nEvSecVertices += 1

		#if vertexCollection.getNumberOfElements() > 1:
		#	print iEvt, str(vertexCollection.getNumberOfElements()) + ' vertices'

	# # if it's a background sample set whatever true LLP values
	# if llp_pdg == 0:
	# 	trueR = 0
	# 	trueZ = 0
	# 	trueL = 0
	# 	e = 0
	# 	px = 0
	# 	py = 0
	# 	pz = 0
	# 	p = 0
	# 	m_llp = 0
	# 	mcllpP4.SetPxPyPzE(px,py,pz,e)


	# if nFSParticle1<1 and nFSParticle2<1:
	# 	#print (iEvt, 'WARNING: no truth tracks in the event!')
	# 	return
	# if nFSParticle1+nFSParticle2 != 2:
	# 	#print (iEvt, nFSParticle1+nFSParticle2, " not 2 truth tracks in the event!")
	# 	return

	gammas = [p4.Gamma() for p4 in mcllpP4s] # gamma factor to get decay length in LLP rest frame
	beta_gamma = [gammas[i] * mcllpP4s[i].Beta() for i in range(len(mcllpP4s))] # factor to get decay length in LLP rest frame
	# print(iEvt, gammas,[mcllpP4s[i].Beta() for i in range(len(mcllpP4s))], beta_gamma)

	# decay lengths used to generate samples have to be hardcoded
	# as this information is not stored in the samples
	if cat.find("idm") != -1:
		tau0 = 1.e3
	elif cat.find("alp") != -1:
		tau0 = m_llp * 10.
	elif cat.find("trsm") != -1:
		if m_llp < 10:
			tau0 = 10.
		else:
			tau0 = 1.e3
	elif cat.find("overlay") != -1 or cat.find("qqbar") != -1:
		tau0 = 1
		beta_gamma = 1
	else:
		print("ERROR: make sure you use correct mean decay length as tau0")
		return

	# calculate weights corresponding to different lifetimes
	weights = np.array([ \
		math.exp( -1. * sum(trueLs[i] / beta_gamma[i] for i in range(len(trueLs))) * (1./tau - 1/tau0) ) \
			* (tau0 / tau)**(len(trueLs)) \
			for tau in decay_lengths ])
	# print(weights)
	
	l_max = 3000
	# w_max = np.ones(weights.shape) \
	w_max = np.array([ \
		expon.cdf( sum(l_max / beta_gamma[i] for i in range(len(beta_gamma))) * (1./tau) )  \
			for tau in decay_lengths ])
	# print (expon.cdf(  sum(l_max * (1./1000) / beta_gamma[i] for i in range(len(beta_gamma))) ))
	# print(w_max)
	# print('')

	if cat.find("overlay") != -1 or cat.find("qqbar") != -1:
		w_max = np.ones(weights.size)
	
	# # if weights[7] > 1.e3 and nLLPsInEvent < 1 and all(l < l_max for l in trueLs):
	# if iEvt == 3919:
		# print(iEvt)
		# print("P1 = (",mcllpP4s[0].Px(),mcllpP4s[0].Py(),mcllpP4s[0].Pz(), mcllpP4s[0].E(), ")")
		# print("P2 = (",mcllpP4s[1].Px(),mcllpP4s[1].Py(),mcllpP4s[1].Pz(), mcllpP4s[1].E(), ")")
		# print("L = " + str(trueL), 'l_max = ' + str(l_max), all(l < l_max for l in trueLs))
		# print('nLLPsInEvent = ', nLLPsInEvent[0])
		# print('nLLPsInTPC = ', nLLPsInTPC[0])
		# print("L = " + str(trueLs), "gamma*beta = " + str(beta_gamma))
		# print('len(mcllpP4s) =',len(mcllpP4s))
		# print(decay_lengths)
		# print(weights,'\n')

	if trueR > 200 and trueR < 300:
		print(iEvt)
		print('trueRs =', trueRs)
		print('secondary vertices: ', nVTX)
		
	# print("Sum of weighted events: ", counter._weightedEvents[4])
	
	# count total number of events corresponding to different lifetimes
	# to get relative number of decays inside TPC
	# counter.increment("_weightedEvents", trueL >= l_max, weights)
	# counter.increment("_weightedEventsErr", trueL >= l_max, (weights)**2)
	# print("_weightedEvents", all(l < l_max for l in trueLs), weights/w_max)
	# print("_weightedEventsErr", all(l < l_max for l in trueLs), (weights)**2,'\n')

	in_acceptance = (trueR <= 1808) and abs(trueZ) <= 2350
	in_TPC = in_acceptance and trueR >= 329
	counter.increment("_nWithinAcceptance", in_acceptance, 2)
	counter.increment("_nWithinTPC", in_TPC, 2)

	# get weighted events in TPC acceptance
	if nLLPsInEvent[0] < 1: # (but only once)
		counter.increment("_weightedEvents", all(l < l_max for l in trueLs), weights/w_max)
		counter.increment("_weightedEventsErr", all(l < l_max for l in trueLs), (weights/w_max)**2)
	if nLLPsInTPC[0] < 1 and all(l < l_max for l in trueLs):
		counter.increment("_weightedEventsInTPC", in_TPC, weights)
	if in_TPC and all(l < l_max for l in trueLs):
		nLLPsInTPC[0] += 1


	if all(l < l_max for l in trueLs):
		nLLPsInEvent[0] += 1

		# histograms.vtxEffRVsZ_true.Fill( trueZ, trueR, (weights)[3] )
		# histograms.vtxEffVsL_true.Fill(trueL,(weights)[3] )
		histograms.vtxEffRVsZ_true.Fill( trueZ, trueR )
		histograms.vtxEffVsL_true.Fill(trueL )
		histograms.vtxEffVsR_true.Fill(trueR )
		histograms.vtxMultiplicity.Fill(nVTX)

		# print(iEvt)
		# print('weights:', weights)
		# print('w_max:', w_max)
		# print('weights/w_max:', weights/w_max)
		# print('(weights/w_max)**2:', (weights/w_max)**2)
		# print('')

	if nLLPsInEvent[0] < 2:

		counter.increment("_nEvSecVertices", True, 1)

		counter.increment("_twoTruthTracks", True, 1)

		#count events with at least 2 reconstructed tracks
		numberOfTracks = trackCollection.getNumberOfElements()
		# if int(numberOfTracks) >= 2:
		# 	_twoTracks+=1
		counter.increment("_twoTracks", int(numberOfTracks) >= 2, 1)
		#if numberOfTracks > 2:
			#print iEvt, numberOfTracks
		histograms.h_nTracks.Fill(numberOfTracks)

		#count reconstructed tracks
		# for track in trackCollection:
			# _allTracks+=1
		counter.increment("_allTracks", True, int(numberOfTracks))

	#count tracks matching to MC particles
	matchingTracks = 0
	for particle in mcCollection:
		if(particle.getGeneratorStatus() == gen_status):
			for track in trackToMCLinkCollection:
				histograms.h_relWeights.Fill(track.getWeight())
				#if track.getTo() == particle and track.getWeight() >= 1.0:
					#_nMatchingTracks+=1
					#print track.getWeight()


	if not in_acceptance and (infileName.find("idm") != -1 or infileName.find("alp") != -1 \
						   or infileName.find("trsm") != -1): # cut on acceptance only if it's BSM sample
		return

	mclepP4 = ROOT.TLorentzVector()
	mclep1 = ROOT.TLorentzVector()
	mclep2 = ROOT.TLorentzVector()
	mcTracks = []
	cosTheta = 0
	deltaTheta = 0
	deltaPhi = 0
	nFSParticle1 = 0
	nFSParticle2 = 0
	for particle in mcCollection:	# loop to fill all quantities in acceptance

		e = particle.getEnergy()
		px = particle.getMomentum()[0]
		py = particle.getMomentum()[1]
		pz = particle.getMomentum()[2]
		p = utils.R3(px,py,pz)

		if(particle.getPDG() == llp_pdg and particle.getGeneratorStatus() == gen_status_llp):
			histograms.momLLP_true.Fill(p,2.)

		if(particle.getGeneratorStatus() == gen_status and (abs( particle.getPDG() ) == pdg1 or abs( particle.getPDG() ) == pdg2)):
			if utils.R(particle.getVertex()[0],particle.getVertex()[1]) != trueR:
				# print(iEvt, trueR, particle.getVertex()[0],particle.getVertex()[1], 'skippin track')
				continue
			mclepP4.SetPxPyPzE(px,py,pz,e)
			if not mclep1.E():
				mclep1.SetPxPyPzE(px,py,pz,e)
			elif not mclep2.E():
				mclep2.SetPxPyPzE(px,py,pz,e)
			cosTheta = abs(math.cos(mclepP4.Theta()))
			deltaTheta = abs( mcllpP4.Theta() - mclepP4.Theta() )
			deltaPhi = utils.phiDistance( mcllpP4.Phi(), mclepP4.Phi() )
			deltaAlpha = mcllpP4.Angle( mclepP4.Vect() )

			if cat == 'sm' and abs(mclepP4.Eta()) > 2.4:
				return

			# cut = mclepP4.Pt() > 10 and (mclepP4.Theta() > 1.4 and mclepP4.Theta() < 1.75) \
			# 	and (mclepP4.Phi() > 1.4 and mclepP4.Phi() < 1.75)
			cut = True
			
			if abs( particle.getPDG() ) == pdg1:
				nFSParticle1+=1
				if cut:
					histograms.trackEffVsR_true.Fill(trueR)
					histograms.trackEffVsL_true.Fill(trueL)
					histograms.trackEffRVsZ_true.Fill(trueZ,trueR)
					histograms.momTrack_true.Fill(p)
					histograms.ptTrack_true.Fill(mclepP4.Pt())
					histograms.thetaTrack_true.Fill(mclepP4.Theta())
					histograms.DeltaThetaTrackLLP_true.Fill( deltaTheta )
					histograms.DeltaPhiTrackLLP_true.Fill( abs(deltaPhi) )
					histograms.DeltaAlphaTrackLLP_true.Fill( abs(deltaAlpha) )
					histograms.trackEffiPvsTheta_true.Fill(mclepP4.Theta(),p)
					histograms.trackEffiPvsR_true.Fill(p,trueR)
					histograms.vtxEffVsPt_true.Fill((mclep1+mclep2).Pt())
					histograms.vtxEffVsP_true.Fill((mclep1+mclep2).P())

					# print (iEvt, trueR, trueZ)

			elif abs( particle.getPDG() ) == pdg2:
				nFSParticle2+=1
				if cut:
					histograms.trackEffVsR_true.Fill(trueR)
					histograms.trackEffVsL_true.Fill(trueL)
					histograms.trackEffRVsZ_true.Fill(trueZ,trueR)
					histograms.momTrack_true.Fill(p)
					histograms.ptTrack_true.Fill(mclepP4.Pt())
					histograms.thetaTrack_true.Fill(mclepP4.Theta())
					histograms.eltaThetaTrackLLP_true.Fill( deltaTheta )
					histograms.DeltaPhiTrackLLP_true.Fill( abs(deltaPhi) )
					histograms.DeltaAlphaTrackLLP_true.Fill( abs(deltaAlpha) )
					histograms.trackEffiPvsTheta_true.Fill(mclepP4.Theta(),p)
					histograms.trackEffiPvsR_true.Fill(p,trueR)

					# print (iEvt, trueR, trueZ)

		# if nFSParticle1 > 1 or nFSParticle2 > 1 or (nFSParticle1>0 and nFSParticle2>0):
		# 	break
		#if nFSParticle1>0 and nFSParticle2>0:
		#	print (iEvt, 'electron and muon!')
	

	histograms.cosOpenAngleVtx_true.Fill( math.cos(mclep1.Angle(mclep2.Vect())) )

	if nLLPsInEvent[0] < 2:
	
		#########################################
		############### MARLINTRK ###############
		#########################################

		tracking.analyse_tracks(mcCollection, trackCollection, \
				llp_pdg, pdg1, pdg2, gen_status_llp, gen_status, \
					"MarlinTrk", useTracks == "MarlinTrk", counter, histograms, cat, iEvt)
		
		#########################################
		############### REFITTED ################
		#########################################

		tracking.analyse_tracks(mcCollection, refitTrackCollection, \
				llp_pdg, pdg1, pdg2, gen_status_llp, gen_status, \
					"Refitted", useTracks == "Refitted", counter, histograms, cat, iEvt)


	#########################################
	############### VERTICES ################
	#########################################

	if 'LLPVertices' in event.getCollectionNames():

		match_success = False

		if useTracks == "Refitted":

			match_success = vertex_finding.vertex_analysis(event, cat, vertexCollection, refitTrackCollection, mcCollection, \
								vertexTracksRelations, refitTrackToMCLinkCollection, \
									llp_endpoint, mclep1, mclep2, weights, \
									counter, histograms, iEvt)
			
		if useTracks == "MarlinTrk":

			match_success = vertex_finding.vertex_analysis(event, cat, vertexCollection, trackCollection, mcCollection, \
								vertexTracksRelations, trackToMCLinkCollection, \
									llp_endpoint, mclep1, mclep2, weights, \
									counter, histograms, iEvt)
			
		# count weighted passing vertices (but only once)
		if all(l < l_max for l in trueLs):
			counter.increment("_matchingWeightedEventsTPC", match_success == True and nMatchedLLPsInEvent[0] < 1, weights)
			counter.increment("_matchingWeightedEventsTPCErr", match_success == True and nMatchedLLPsInEvent[0] < 1, weights**2)
		if match_success:
			nMatchedLLPsInEvent[0] += 1
			# print(iEvt, str(counter._weightedEvents)+ '\n' + str(counter._weightedEventsErr) + "\n" \
		 	# 					 + str(counter._matchingWeightedEventsTPC) +'\n' \
			# 					 + str(counter._matchingWeightedEventsTPCErr) + '\n')

			# print(iEvt)
			# print('weights:', weights)
			# print('(weights)**2:', (weights)**2)
			# print('')


	#########################################
	################## PFOs #################
	#########################################
	
	if int(pfoCollection.getNumberOfElements()) > 0:
		for particle in pfoCollection:
			if particle.getParticleIDs().size() > 0:
				pid = particle.getParticleIDs()[0].getPDG()
				# if abs(pid) == pdg1 or abs(pid) == pdg2:
				# 	_nCorrPID+=1
				counter.increment("_nCorrPID", abs(pid) == pdg1 or abs(pid) == pdg2, 1)
			# print(iEvt)
			# print(trueR, math.sqrt(particle.getMomentum()[0]**2+particle.getMomentum()[1]**2))
			# print([id.getPDG() for id in particle.getParticleIDs()],'\n')
			# 	print ''
			# 	for id in particle.getParticleIDs():
			# 		print id.getAlgorithmType(),
			# 	print ''
			# 	print iEvt, particle.getParticleIDUsed()
			# 	print iEvt, particle.getGoodnessOfPID()
			# print ''

def print_output(counter, histograms):

	centresRefDistsVsPtEv = histograms.corrCentresRefDistsVsPt.Integral()
	# print('linear int: ', centresRefDistsVsPtEv)
	# print("log int: ", histograms.corrLogCentresRefDistsVsPt.Integral())
	print('true vertices in TPC:', histograms.vtxEffVsL_true.Integral( \
		histograms.vtxEffVsL_true.FindBin(329),histograms.vtxEffVsL_true.FindBin(1808) \
		))
	print('true vertices outside:', histograms.vtxEffVsL_true.Integral( \
		histograms.vtxEffVsL_true.FindBin(1808),histograms.vtxEffVsL_true.FindBin(12000) \
																))
	# print('reco vertices with small d0 and first hit further than the last one:', \
	#    histograms.vtxTrksD01VsD02.Integral(),histograms.vtxTrksZ01VsZ02.Integral())

	#print('Decays within tracker acceptance', 'Ev. with min. two tracks', 'Ev. with min. two reco. tracks', 'All reco. tracks', 'Tracks matching to MC', 'Two reco. leptons')
	#print(_nWithinAcceptance, _twoTruthTracks, _twoTracks, _allTracks, _nMatchingTracks, _twoRecoLeptons)

	print('')
	print('Tracks within acceptance: ' + str(counter._nWithinAcceptance))
	print('Ev. with two truth tracks: ' + str(counter._twoTruthTracks))
	print('Decays inside TPC: ' + str(float(counter._nWithinTPC)/2))
	print('Ev. with min. two reco. tracks: ' + str(counter._twoTracks))
	print('All reco. tracks: ' + str(counter._allTracks))
	print('Number of PFOs: ' + str(counter._nPFOs))
	print('Number of PFOs with correct PID: ' + str(counter._nCorrPID))
	print('PFO ID eff.: ', float(counter._nCorrPID) / float(counter._nWithinAcceptance))
	print('Tracks matching to MC (MarlinTrk): ' + str(counter._nMatchingTracks))
	print('Tracks matching to MC (Refit): ' + str(counter._nMatchingTracksReffited))
	print('Tracks matching only with MarlinTrk: ' + str(counter._onlyMarlinTrkTracks))
	print('Tracks matching only with Refit: ' + str(counter._onlyRefittedTracks))
	if counter._nWithinAcceptance > 0:
		print('MarlinTrk reco. eff.: ', float(counter._nMatchingTracks) / float(counter._nWithinAcceptance))
		print('Refitted reco. eff.:  ', float(counter._nMatchingTracksReffited) / float(counter._nWithinAcceptance))
	print('Number of sec. vtx:  ', counter._nSecVertices)
	print('Number of sec. vtx inside TPC:  ', counter._nSecVerticesTPC)
	print('Ev. with >=1 sec. vtx:  ', counter._nEvSecVertices)
	print('Ev. with >=1 sec. vtx inside TPC:  ', counter._nEvSecVerticesTPC)
	print('Matching vertices:  ', counter._nMatchingVertices)
	print('Fake vertices:  ', counter._nFakeVertices)
	print('Vertices from photon conversions: ', counter._nPhotonConversionsTPC)
	print('Removed by ee mass cut: ', counter._eeMassCut)
	print('Removed by pipi mass cut: ', counter._pipiMassCut)
	print('Removed by lambda mass cut: ', counter._lambdaMassCut)
	print('Removed by sec. int. cut: ', counter._secondaryIntCut)
	if counter._nWithinAcceptance > 0:
		print('Vtx finding eff.:  ', float(counter._nMatchingVertices) / (float(counter._nWithinAcceptance)/2))
		if counter._nSecVertices > 0:
			print('Vtx finding purity:  ', float(counter._nMatchingVertices) / float(counter._nMatchingVertices+counter._nFakeVertices))
			if counter._nWithinTPC > 0:
				print('* Eff. (inside TPC):  ', float(counter._matchingVerticesTPC) / (float(counter._nWithinTPC)/2))
				print('* Purity (inside TPC):  ', float(counter._matchingVerticesTPC) / float(counter._matchingVerticesTPC+counter._nFakeVerticesTPC))
				print('Eff. inside TPC after pT, Ndf cuts: ', float(centresRefDistsVsPtEv) / (float(counter._nWithinTPC)/2))
	print('')

def save_to_file(histo_array):
	try:
		outFileName = sys.argv[-1]
		if outFileName.find("root") != -1:
			outHistFile = ROOT.TFile.Open(outFileName ,"RECREATE")
			outHistFile.cd()
			for hist in histo_array:
				hist.Write()
			outHistFile.Close()
			print('*** File ' + outFileName + ' saved! ***')
	except IndexError:
		pass

def print_to_file(counter, scenario, infileName, decay_lengths):
	try:
		fileName = sys.argv[-1]

		if fileName.find(".txt") != -1 or fileName.find(".log") != -1:
			if os.path.exists(fileName):
				logFile = open(fileName,'a')
			else:
				logFile = open(fileName,'a')
				logFile.write("Scenario   DecayLen   Efficiency   Purity   Efficiency (TPC)   Purity (TPC)   N (all)   N (TPC)   N_pass (TPC)   err (rej)   err (pass)   N_w2   Neq_all   Neq_pass" \
				  			+ infileName + "\n")
			for i in range(counter._weightedEventsInTPC.size):
				logFile.write(scenario + " " + str(decay_lengths[i]) \
				  			+ " " + str( float(counter._nMatchingVertices) / (float(counter._nWithinAcceptance)/2) ) \
							+ " " + str( float(counter._nMatchingVertices) / float(counter._nMatchingVertices+counter._nFakeVertices) ) )
				logFile.write( " " + str( float(counter._matchingVerticesTPC) / (float(counter._nWithinTPC)/2) ) \
							+ " " + str( float(counter._matchingVerticesTPC) / float(counter._matchingVerticesTPC+counter._nFakeVerticesTPC) ) )
				logFile.write( " " + str(counter._weightedEvents[i]) \
							 + " " + str(counter._weightedEventsInTPC[i]) \
							 + " " + str(counter._matchingWeightedEventsTPC[i]) )
				logFile.write( " " + str( math.sqrt(counter._weightedEventsErr[i]-counter._matchingWeightedEventsTPCErr[i]) ) \
							 + " " + str( math.sqrt(counter._matchingWeightedEventsTPCErr[i]) ) )
				logFile.write( " " + str( counter._weightedEventsErr[i] ) \
							 + " " + str( counter._weightedEvents[i]*counter._weightedEvents[i]/counter._weightedEventsErr[i] ) \
							 + " " + str( counter._matchingWeightedEventsTPC[i]*counter._matchingWeightedEventsTPC[i]/counter._matchingWeightedEventsTPCErr[i] ) )
				logFile.write('\n')
			logFile.close()
		else:
			print("Logfile not created.")

	except IndexError:
		pass

if __name__ == '__main__':

	try:
		arg1 = sys.argv[1]
	except IndexError:
		print("Usage: " + os.path.basename(__file__) + " <inputFile(s)>" + " <outputHistFile>")
		sys.exit(1)

	print('Start analysis')

	print('Opening input file ' + str(sys.argv[1]) + '...')
	infile = sys.argv[1]
	reader = IOIMPL.LCFactory.getInstance().createLCReader()
	try:
		reader.open(infile)
	except IOException as e:
		print("Caught exception:", e)
	print('input file successfully open')

	infileName = str(infile)

	fs_id = infileName[infileName.find("simulation")+12]
	fs_dict = {"4": ("pi",211), "5": ("e",11), "6": ("mu",13), "7": ("tau",15), "a": ("lambda",2212)}
	fs, fs_pdg = fs_dict[fs_id]

	pdg1 = fs_pdg
	pdg2 = fs_pdg
	gen_status = 1
	gen_status_llp = 2
	llp_pdg = 36

	if infileName.find("lambda") != -1:
		pdg1 = 2212
		pdg2 = 211
		gen_status = 0
		gen_status_llp = 1
		llp_pdg = 3122
		fs = "lambda"

	if infileName.find("trsm") != -1:
		pdg1 = fs_pdg
		pdg2 = fs_pdg
		gen_status = 1
		gen_status_llp = 2
		llp_pdg = 25

	# useTracks = 'Refitted'
	useTracks = 'MarlinTrk'

	cat = 'idm'
	scenario = infileName[infileName.find("dM"):infileName.find("dM")+4]
	decay_lengths = np.array([3.e-2, 1.e-1, 3.e-1, 1., 3., 1.e1, 1.e2, 1.e3, 1.e4]) * 1.e3
	if infileName.find("alp") != -1:
		cat = 'alp'
		scenario = infileName[infileName.find("_M"):infileName.find("_M")+7]
		scenario = re.findall(r'\d+', scenario)[0]
		xfactor = 0.001 # IDM mass splitting is in 10 GeV and ALP masses in MeV
		decay_lengths = np.array([3.e-2, 8.e-2, 0.3, 1., 3., 1.e1, 1.e2, 1.e3, 1.e4]) * xfactor * float(scenario) * 10
		# decay_lengths = np.array([1., 1.e1, 1.e2, 1.e3, 1.e4]) 
	if infileName.find("trsm") != -1:
		cat = 'trsm'
		scenario = infileName[infileName.find("_M"):infileName.find("_M")+4]
		scenario = re.findall(r'\d+', scenario)[0]
		decay_lengths = np.array([6.e-2, 1.e-1, 3.e-1, 1., 3., 1.e1, 1.e2, 1.e3, 1.e4]) * 1.e3
		if int(scenario) < 1:
			decay_lengths = np.array([5.e-4, 1.e-3, 3.e-3, 1.e-2, 4.e-2, 1.e-1, 1., 1.e1, 1.e2]) * 1.e3
		elif int(scenario) < 10:
			decay_lengths = np.array([1.e-3, 1.e-2, 1.e-1, 1., 3., 1.e1, 1.e2, 1.e3, 1.e4]) * 1.e3

	# if infileName.find("sm") != -1:
	# 	cat = 'sm'
	# 	llp_pdg = 0
	if infileName.find("BB") != -1 \
	or infileName.find("BW") != -1 \
	or infileName.find("WB") != -1 \
	or infileName.find("WW") != -1 \
	or infileName.find("eepairs") != -1:
		cat = 'overlay'
		llp_pdg = 0
		decay_lengths = np.array([1.])
	if infileName.find("qqbar") != -1:
		cat = 'qqbar'
		llp_pdg = 0
		decay_lengths = np.array([1.])
	if infileName.find("Silicon") != -1:
		cat = 'SiliconTracker/'+cat
	if infileName.find("Pixel") != -1:
		cat = 'PixelTPC/'+cat

	if len(re.findall(r'\d+', scenario)) > 0:
		scenario = re.findall(r'\d+', scenario)[0]
	else:
		scenario = ''

	print("LLP ("+ scenario +"): " + str(llp_pdg) + ", with final state: " + str(pdg1) + ", " + str(pdg2))

	counter = counters.Counters()
	histograms = histo_class.Histograms()

	run_analysis(reader=reader, cat=cat, fs=fs, \
			  counter=counter, histograms=histograms, \
				decay_lengths=decay_lengths)
	
	print_output(counter, histograms)
	print_to_file(counter, scenario, infileName, decay_lengths)
	histo_array = histograms.print_histos(cat, fs, scenario)
	save_to_file(histo_array)