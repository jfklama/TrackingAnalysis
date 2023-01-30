import sys
import os

import math
import numpy as np

from pyLCIO import IOIMPL, EVENT, IMPL, ROOT, IOException
#from ROOT import Math

import utils
import constants
import counters
import histo_class
import tracking
import vertex_finding

def run_analysis(reader, cat, fs):
	counter = counters.Counters()
	histograms = histo_class.Histograms()
	for event in reader:
		process_event(event, counter, histograms)
	reader.close()
	print_output(counter, histograms)
	histo_array = histograms.print_histos(cat, fs)
	save_to_file(histo_array)

def process_event(event,counter, histograms):

	iEvt = event.getEventNumber()

	mcCollection = event.getCollection('MCParticlesSkimmed')
	#tpcCollection = event.getCollection('TPCCollection')
	trackCollection = event.getCollection('MarlinTrkTracks')
	trackToMCLinkCollection = event.getCollection('MarlinTrkTracksMCTruthLink')
	refitTrackCollection = event.getCollection('RefittedMarlinTrkTracks')
	refitTrackToMCLinkCollection = event.getCollection('RefittedMarlinTrkTracksMCTruthLink')
	pfoCollection = event.getCollection('PandoraPFOs')
	TPCHitRelations = event.getCollection('TPCTrackerHitRelations')

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

	if iEvt%100==0:
			print (iEvt, str(trackCollection.getNumberOfElements() ) + ' tracks in event')

	
	# _nPFOs += pfoCollection.getNumberOfElements()
	counter.increment("_nPFOs", True, int(pfoCollection.getNumberOfElements()))

	#check if LLP decayed within tracking acceptance and we have 2 truth tracks
	trueR = 0
	trueZ = 0
	trueL = 0
	nFSParticle1 = 0
	nFSParticle2 = 0
	llp_endpoint = 0
	for particle in mcCollection:
		if abs(particle.getPDG()) == llp_pdg:
			llp_endpoint = particle.getEndpoint()
			trueR = utils.R(llp_endpoint[0], llp_endpoint[1])
			trueZ = llp_endpoint[2]
			trueL = utils.spacialDistance(llp_endpoint, [0,0,0])
		if(particle.getGeneratorStatus() == gen_status):
			if abs( particle.getPDG() ) == pdg1:
				nFSParticle1+=1
			elif abs( particle.getPDG() ) == pdg2:
				nFSParticle2+=1

	if nFSParticle1<1 and nFSParticle2<1:
		#print (iEvt, 'WARNING: no truth tracks in the event!')
		return
	if nFSParticle1+nFSParticle2 != 2:
		#print (iEvt, nFSParticle1+nFSParticle2, " not 2 truth tracks in the event!")
		return


	in_acceptance = (trueR <= 1808) and abs(trueZ) <= 2350
	in_TPC = in_acceptance and trueR >= 329
	counter.increment("_nWithinAcceptance", in_acceptance, 2)
	counter.increment("_nWithinTPC", in_TPC, 2)
	if not in_acceptance and infileName.find("idm") != -1: # cut on TPC acceptance only if it's IDM sample
		return
	# if (trueR <= 1808) and abs(trueZ) <= 2350:
	# 	_nWithinAcceptance+=2
	# 	if trueR >= 329:
	# 		_nWithinTPC+=2
	# if (trueR < 329) and abs(trueZ) <= 2350:
	# 	_nWithinAcceptance+=2
	# 	_nWithinTPC+=2

	histograms.vtxEffRVsZ_true.Fill( trueZ, trueR )
	histograms.vtxEffVsL_true.Fill(trueL)
	histograms.vtxMultiplicity.Fill(nVTX)

	if 'LLPVertices' in event.getCollectionNames():
		counter.increment("_nSecVertices", True, int(vertexCollection.getNumberOfElements()))
		counter.increment("_nEvSecVertices", True, 1)
		counter.increment("_nSecVerticesTPC", in_TPC, 1)

	# if 'LLPVertices' in event.getCollectionNames():
	# 	nVTX = vertexCollection.getNumberOfElements()
	# 	_nSecVertices += int(vertexCollection.getNumberOfElements())
	# 	_nEvSecVertices += 1
	# 	if (trueR >= 329 and trueR <= 1808) and abs(trueZ) <= 2350:
	# 		_nSecVerticesTPC += 1
	#if nVTX != 3 and nVTX != 6:
	#	return

	# if sm_background, check if tracks are within detector acceptance
	if cat == 'sm':
		for particle in mcCollection:
			if particle.getGeneratorStatus() == gen_status \
			and (abs( particle.getPDG() ) == pdg1 or abs( particle.getPDG() ) == pdg2):
				mctrackP4 = ROOT.TLorentzVector()
				e = particle.getEnergy()
				px = particle.getMomentum()[0]
				py = particle.getMomentum()[1]
				pz = particle.getMomentum()[2]
				mctrackP4.SetPxPyPzE(px,py,pz,e)
				# if abs( mctrackP4.Eta() ) <= 2.4:
					# _nWithinAcceptance += 1
				counter.increment("_nWithinAcceptance", abs( mctrackP4.Eta() ) <= 2.4, 1)

	# _twoTruthTracks+=1
	counter.increment("_twoTruthTracks", True, 1)

	mclepP4 = ROOT.TLorentzVector()
	mclep1 = ROOT.TLorentzVector()
	mclep2 = ROOT.TLorentzVector()
	mcllpP4 = ROOT.TLorentzVector()
	mcTracks = []
	cosTheta = 0
	deltaTheta = 0
	deltaPhi = 0
	nFSParticle1 = 0
	nFSParticle2 = 0
	for particle in mcCollection:	# loop to fill all quantities in acceptance

		if(particle.getPDG() == llp_pdg and particle.getGeneratorStatus() == gen_status_llp):
			e = particle.getEnergy()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = utils.R3(px,py,pz)
			mcllpP4.SetPxPyPzE(px,py,pz,e)
			histograms.momLLP_true.Fill(p,2.)

		if(particle.getGeneratorStatus() == gen_status and (abs( particle.getPDG() ) == pdg1 or abs( particle.getPDG() ) == pdg2)):
			e = particle.getEnergy()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = utils.R3(px,py,pz)
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

		if nFSParticle1 > 1 or nFSParticle2 > 1 or (nFSParticle1>0 and nFSParticle2>0):
			break
		#if nFSParticle1>0 and nFSParticle2>0:
		#	print (iEvt, 'electron and muon!')

	histograms.cosOpenAngleVtx_true.Fill( math.cos(mclep1.Angle(mclep2.Vect())) )

	#count events with at least 2 reconstructed tracks
	numberOfTracks = trackCollection.getNumberOfElements()
	# if int(numberOfTracks) >= 2:
	# 	_twoTracks+=1
	counter.increment("_twoTracks", int(numberOfTracks) >= 2, 2)
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
	
	#########################################
	############### MARLINTRK ###############
	#########################################

	tracking.analyse_tracks(mcCollection, trackToMCLinkCollection, \
			llp_pdg, pdg1, pdg2, gen_status_llp, gen_status, \
				"MarlinTrk", useTracks == "MarlinTrk", counter, histograms, cat)
	
	#########################################
	############### REFITTED ################
	#########################################

	tracking.analyse_tracks(mcCollection, refitTrackToMCLinkCollection, \
			llp_pdg, pdg1, pdg2, gen_status_llp, gen_status, \
				"Refitted", useTracks == "Refitted", counter, histograms, cat)


	#########################################
	############### VERTICES ################
	#########################################

	if 'LLPVertices' in event.getCollectionNames():

		vertex_finding.vertex_analysis(vertexCollection, refitTrackCollection, \
				 vertexTracksRelations, refitTrackToMCLinkCollection, \
					llp_endpoint, counter, histograms)


	#########################################
	################## PFOs #################
	#########################################

	for particle in pfoCollection:
		pid = particle.getParticleIDs()[0].getPDG()
		# if abs(pid) == pdg1 or abs(pid) == pdg2:
		# 	_nCorrPID+=1
		counter.increment("_nCorrPID", abs(pid) == pdg1 or abs(pid) == pdg2, 1)
	# 	print iEvt,
	# 	for id in particle.getParticleIDs():
	# 		print id.getPDG(),
	# 	print ''
	# 	for id in particle.getParticleIDs():
	# 		print id.getAlgorithmType(),
	# 	print ''
	# 	print iEvt, particle.getParticleIDUsed()
	# 	print iEvt, particle.getGoodnessOfPID()
	# print ''

def print_output(counter, histograms):

	centresRefDistsVsPtEv = histograms.corrCentresRefDistsVsPt.Integral()
	print('linear int: ', centresRefDistsVsPtEv)
	print("log int: ", histograms.corrLogCentresRefDistsVsPt.Integral())

	#print('Decays within tracker acceptance', 'Ev. with min. two tracks', 'Ev. with min. two reco. tracks', 'All reco. tracks', 'Tracks matching to MC', 'Two reco. leptons')
	#print(_nWithinAcceptance, _twoTruthTracks, _twoTracks, _allTracks, _nMatchingTracks, _twoRecoLeptons)

	print('')
	print('Tracks within acceptance: ' + str(counter._nWithinAcceptance))
	print('Ev. with two truth tracks: ' + str(counter._twoTruthTracks))
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
	print('Ev. with >=1 sec. vtx:  ', counter._nEvSecVertices)
	print('Matching vertices:  ', counter._nMatchingVertices)
	print('Fake vertices:  ', counter._nFakeVertices)
	if counter._nWithinAcceptance > 0:
		print('Vtx finding eff.:  ', float(counter._nMatchingVertices) / (float(counter._nWithinAcceptance)/2))
		print('Vtx finding purity:  ', float(counter._nMatchingVertices) / float(counter._nSecVertices))
		print('* Eff. (inside TPC):  ', float(counter._matchingVerticesTPC) / (float(counter._nWithinTPC)/2))
		print('* Purity (inside TPC):  ', float(counter._matchingVerticesTPC) / float(counter._nSecVerticesTPC))
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

	useTracks = 'Refitted'
	#useTracks = 'MarlinTrk'

	cat = 'idm'
	if infileName.find("sm") != -1:
		cat = 'sm'
		llp_pdg = 0

	print("LLP: " + str(llp_pdg) + ", with final state: " + str(pdg1) + ", " + str(pdg2))

	run_analysis(reader=reader, cat=cat, fs=fs)