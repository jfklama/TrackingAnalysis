import sys
import os

import math
import numpy as np

from pyLCIO import IOIMPL, ROOT
#from ROOT import Math

import utils

def R(x,y):
	return math.sqrt(x**2 + y**2)
def R3(x,y,z):
	return math.sqrt(x**2 + y**2 + z**2)

try:
    arg1 = sys.argv[1]
except IndexError:
    print "Usage: " + os.path.basename(__file__) + " <inputFile>" + " <outputHistFile>"
    sys.exit(1)

print 'Start analysis'

chargedChannel    		  = 0
withinAcceptance  		  = 0
twoTracks         		  = 0
twoTruthTracks   		  = 0
allTracks		  		  = 0
allMatchingTracks 		  = 0
allMatchingTracksRefitted = 0
twoRecoLeptons    		  = 0
onlyMarlinTrkTracks		  = 0
onlyRefittedTracks		  = 0
nPFOs					  = 0
nSecVertices		  	  = 0

nbins = 12
TPCRmax = 1974.
#TPCRmax = 1808.
TPCZmax = 2350.

trackEffVsR_true = ROOT.TH1F('trackEffVsR_true', 'True vertices', nbins, 0., TPCRmax)
trackEffVsR_reco = ROOT.TH1F('trackEffVsR_reco', 'N_{evt} with 2 reco. tracks', nbins, 0., TPCRmax)
trackEffVsL_true = ROOT.TH1F('trackEffVsL_true', 'True vertices', 30, 0., 3000.)
trackEffVsL_reco = ROOT.TH1F('trackEffVsL_reco', 'N_{evt} with 2 reco. tracks', 30, 0., 3000.)

trackEffRVsZ_true = ROOT.TH2F('trackEffRVsZ_true', 'True vertices', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffRVsZ_reco = ROOT.TH2F('trackEffRVsZ_reco', 'N_{evt} with 2 reco. tracks', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffiPvsTheta_true = ROOT.TH2F('trackEffiPvsTheta_true', '', 3, 0., 3.15, 15, 0., 15.)
trackEffiPvsTheta_reco = ROOT.TH2F('trackEffiPvsTheta_reco', '', 3, 0., 3.15, 15, 0., 15.)
trackEffiPvsR_true = ROOT.TH2F('trackEffiPvsR_true', '', 15, 0., 15., nbins, 0., TPCRmax)

d12vsTheta_reco = ROOT.TH2F('d12vsTheta_reco', '', 32, 0., 3.15, 32, 0., 150.)
d12vsTheta_true = ROOT.TH2F('d12vsTheta_true', '', 32, 0., 3.15, 32, 0., 150.)
d12vsDeltaPhi_true = ROOT.TH2F('d12vsDeltaPhi_true', '', 32, 0., 3.15, 32, 0., 150.)
#d12vsDeltaPhi_reco = ROOT.TH2F('d12vsDeltaPhi_reco', '', 32, 0., 3.15, 32, 0., 100.)


momTrack_true = ROOT.TH1F('momTrack_true', 'True momentum of track', 50, 0., 15.)
momTrack_reco = ROOT.TH1F('momTrack_reco', 'True momentum of reco. track', 50, 0., 15.)
ptTrack_true = ROOT.TH1F('ptTrack_true', 'True pT of track', 50, 0., 15.)
ptTrack_reco = ROOT.TH1F('ptTrack_reco', 'True pT of reco. track', 50, 0., 15.)
thetaTrack_true = ROOT.TH1F('thetaTrack_true', 'True polar angle of track', 42, 0., 3.15)
thetaTrack_reco = ROOT.TH1F('thetaTrack_reco', 'True polar angle of reco. track', 42, 0., 3.15)

DeltaThetaTrackLLP_true = ROOT.TH1F('DeltaThetaTrackLLP_true', 'Dist. between polar angles of track and LLP', 35, 0., 3.5)
DeltaThetaTrackLLP_reco = ROOT.TH1F('DeltaThetaTrackLLP_reco', 'Dist. between polar angles of track and LLP', 35, 0., 3.5)
DeltaPhiTrackLLP_true = ROOT.TH1F('DeltaPhiTrackLLP_true', 'Dist. between azimuthal angles of track and LLP', 35, 0., 3.5)
DeltaPhiTrackLLP_reco = ROOT.TH1F('DeltaPhiTrackLLP_reco', 'Dist. between azimuthal angles of track and LLP', 35, 0., 3.5)
DeltaAlphaTrackLLP_true = ROOT.TH1F('DeltaAlphaTrackLLP_true', 'Angle between track and LLP', 35, 0., 3.5)
DeltaAlphaTrackLLP_reco = ROOT.TH1F('DeltaAlphaTrackLLP_reco', 'Angle between track and LLP', 35, 0., 3.5)

momLLP_true = ROOT.TH1F('momLLP_true', 'True momentum of LLP', 50, 0., 300.)
momLLP_reco = ROOT.TH1F('momLLP_reco', 'True momentum of reco. LLP', 50, 0., 300.)
thetaLLP_true = ROOT.TH1F('thetaLLP_true', 'True polar angle of LLP', 42, 0., 3.15)
thetaLLP_reco = ROOT.TH1F('thetaLLP_reco', 'True polar angle of reco. LLP', 42, 0., 3.15)

h_nTracks = ROOT.TH1F('h_nTracks', 'Number of tracks', 20, 0., 20.)
h_relWeights = ROOT.TH1F('h_relWeights', 'Number of weights', 50, 0., 2.)
h_angDist = ROOT.TH1F('h_angDist', 'Angular separation', 70, 0., 7)

recoVertices = ROOT.TH2F('recoVertices', 'N_{evt} with 2 reco. leptons', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
vtxMisdistances = ROOT.TH1F('vtxMisdistances', 'Distance between true and reco. vertex', 50, 0., 200)
hitsDistances = ROOT.TH1F('hitsDistances', 'Distance between two closest first/last hits', 50, 0., 200)
vtxMultiplicity = ROOT.TH1F('vtxMultiplicity', 'Reco. vertex multiplicity', 10, 0., 10)

hmass = ROOT.TH1F('hmass', '#Lambda mass', 150, 0., 15.)

print 'Opening input file ' + str(sys.argv[1]) + '...'
infile = sys.argv[1]
reader=IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(infile)
print 'input file successfully open'

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

if cat == 'sm':
	llp_pdg = 0

print "LLP: " + str(llp_pdg) + ", with final state: " + str(pdg1) + ", " + str(pdg2)

i=0
for event in reader:
	i+=1

	mcCollection = event.getCollection('MCParticlesSkimmed')
	#tpcCollection = event.getCollection('TPCCollection')
	trackCollection = event.getCollection('MarlinTrkTracks')
	trackToMCLinkCollection = event.getCollection('MarlinTrkTracksMCTruthLink')
	refitTrackCollection = event.getCollection('RefittedMarlinTrkTracks')
	refitTrackToMCLinkCollection = event.getCollection('RefittedMarlinTrkTracksMCTruthLink')
	pfoCollection = event.getCollection('PandoraPFOs')
	TPCHitRelations = event.getCollection('TPCTrackerHitRelations')
	#BuildUpVertexCollection = event.getCollection('BuildUpVertex')
	if 'LLPVerticesTest' in event.getCollectionNames():
		vertexCollection = event.getCollection('LLPVerticesTest')
		nVTX = vertexCollection.getNumberOfElements()
		vtxMultiplicity.Fill(nVTX)
		nSecVertices += 1

	nPFOs += pfoCollection.getNumberOfElements()

	if i%100==0:
		print (i, str(trackCollection.getNumberOfElements() ) + ' tracks in event')

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
			trueR = R(llp_endpoint[0], llp_endpoint[1])
			trueZ = llp_endpoint[2]
			trueL = utils.spacialDistance(llp_endpoint, [0,0,0])
		if(particle.getGeneratorStatus() == gen_status):
			if abs( particle.getPDG() ) == pdg1:
				nFSParticle1+=1
			elif abs( particle.getPDG() ) == pdg2:
				nFSParticle2+=1

	if nFSParticle1<1 and nFSParticle2<1:
		#print (i, 'WARNING: no truth tracks in the event!')
		continue
	if nFSParticle1+nFSParticle2 != 2:
		#print (i, nFSParticle1+nFSParticle2, " not 2 truth tracks in the event!")
		continue

	if (trueR >= 329 and trueR <= 1808) and abs(trueZ) <= 2350:
		withinAcceptance+=2
	elif infileName.find("idm") != -1:
		continue

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
				if abs( mctrackP4.Eta() ) <= 2.4:
					withinAcceptance += 1

	twoTruthTracks+=1

	mclepP4 = ROOT.TLorentzVector()
	mcllpP4 = ROOT.TLorentzVector()
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
			p = R3(px,py,pz)
			mcllpP4.SetPxPyPzE(px,py,pz,e)
			momLLP_true.Fill(p,2.)

		if(particle.getGeneratorStatus() == gen_status):
			e = particle.getEnergy()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = R3(px,py,pz)
			mclepP4.SetPxPyPzE(px,py,pz,e)
			cosTheta = abs(math.cos(mclepP4.Theta()))
			deltaTheta = abs( mcllpP4.Theta() - mclepP4.Theta() )
			deltaPhi = utils.phiDistance( mcllpP4.Phi(), mclepP4.Phi() )
			deltaAlpha = mcllpP4.Angle( mclepP4.Vect() )

			if cat == 'sm' and abs(mclepP4.Eta()) > 2.4:
				continue

			cut = True

			if abs( particle.getPDG() ) == pdg1:
				nFSParticle1+=1
				if cut:
					trackEffVsR_true.Fill(trueR)
					trackEffVsL_true.Fill(trueL)
					trackEffRVsZ_true.Fill(trueZ,trueR)
					momTrack_true.Fill(p)
					ptTrack_true.Fill(mclepP4.Pt())
					thetaTrack_true.Fill(mclepP4.Theta())
					DeltaThetaTrackLLP_true.Fill( deltaTheta )
					DeltaPhiTrackLLP_true.Fill( abs(deltaPhi) )
					DeltaAlphaTrackLLP_true.Fill( abs(deltaAlpha) )
					trackEffiPvsTheta_true.Fill(mclepP4.Theta(),p)
					trackEffiPvsR_true.Fill(p,trueR)

			elif abs( particle.getPDG() ) == pdg2:
				nFSParticle2+=1
				if cut:
					trackEffVsR_true.Fill(trueR)
					trackEffVsL_true.Fill(trueL)
					trackEffRVsZ_true.Fill(trueZ,trueR)
					momTrack_true.Fill(p)
					ptTrack_true.Fill(mclepP4.Pt())
					thetaTrack_true.Fill(mclepP4.Theta())
					DeltaThetaTrackLLP_true.Fill( deltaTheta )
					DeltaPhiTrackLLP_true.Fill( abs(deltaPhi) )
					DeltaAlphaTrackLLP_true.Fill( abs(deltaAlpha) )
					trackEffiPvsTheta_true.Fill(mclepP4.Theta(),p)
					trackEffiPvsR_true.Fill(p,trueR)

		if nFSParticle1 > 1 or nFSParticle2 > 1 or (nFSParticle1>0 and nFSParticle2>0):
			break
		#if nFSParticle1>0 and nFSParticle2>0:
		#	print (i, 'electron and muon!')

	#count events with at least 2 reconstructed tracks
	numberOfTracks = trackCollection.getNumberOfElements()
	if int(numberOfTracks) >= 2:
		twoTracks+=1
	#if numberOfTracks > 2:
		#print i, numberOfTracks
	h_nTracks.Fill(numberOfTracks)

	#count reconstructed tracks
	for track in trackCollection:
		allTracks+=1

	#count tracks matching to MC particles
	matchingTracks = 0
	for particle in mcCollection:
		if(particle.getGeneratorStatus() == gen_status):
			for track in trackToMCLinkCollection:
				h_relWeights.Fill(track.getWeight())
				#if track.getTo() == particle and track.getWeight() >= 1.0:
					#allMatchingTracks+=1
					#print track.getWeight()

	#########################################
	############### MARLINTRK ###############
	#########################################

	vMatchingTracks = []
	matchingTracksInEvent = 0
	matchingTracks = 0
	p_llp = 0
	for particle in mcCollection:

		if matchingTracksInEvent > 1:
			break

		if(particle.getPDG() == llp_pdg and particle.getGeneratorStatus() == gen_status_llp):
			e = particle.getEnergy()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p_llp = R3(px,py,pz)
			mcllpP4.SetPxPyPzE(px,py,pz,e)

		if( particle.getGeneratorStatus() == gen_status):
			if abs(particle.getPDG()) != pdg1 and abs(particle.getPDG()) != pdg2:
				continue
			vtx_mc = particle.getVertex()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = R3(px,py,pz)
			e = particle.getEnergy()
			mclepP4.SetPxPyPzE(px,py,pz,e)
			cosTheta = abs( math.cos(mclepP4.Theta()) )
			deltaTheta = abs( mcllpP4.Theta() - mclepP4.Theta() )
			deltaPhi = utils.phiDistance( mcllpP4.Phi(), mclepP4.Phi() )
			deltaAlpha = mcllpP4.Angle( mclepP4.Vect() )

			if cat == 'sm' and abs(mclepP4.Eta()) > 2.4:
				continue

			cut = True

			for trackRel in trackToMCLinkCollection:

				track = trackRel.getFrom()

				ts_FirstHit = track.getTrackState(2) # 1 = atIP, 2 = atFirstHit, 3 = atLastHit
				ts_FirstHit = utils.fixTrackStateDirection( track, ts_FirstHit.getLocation() )
				momentum = utils.getTrackMomentum(ts_FirstHit)
				p_fh = ROOT.TVector3()
				p_fh.SetXYZ(momentum[0],momentum[1],momentum[2])
				dist_fh = utils.angularDistance( p_fh.Theta(), mclepP4.Theta(), p_fh.Phi(), mclepP4.Phi()  )
				vtx_fh = ts_FirstHit.getReferencePoint()
				vtx_dist_fh = utils.spacialDistance(vtx_mc, vtx_fh)

				ts_LastHit = track.getTrackState(3) # 1 = atIP, 2 = atFirstHit, 3 = atLastHit
				ts_LastHit = utils.fixTrackStateDirection( track, ts_LastHit.getLocation() )
				momentum = utils.getTrackMomentum(ts_LastHit)
				p_lh = ROOT.TVector3()
				p_lh.SetXYZ(momentum[0],momentum[1],momentum[2])
				dist_lh = utils.angularDistance( p_lh.Theta(), mclepP4.Theta(), p_lh.Phi(), mclepP4.Phi()  )
				vtx_lh = ts_LastHit.getReferencePoint()
				vtx_dist_lh = utils.spacialDistance(vtx_mc, vtx_lh)

				ang_dist = dist_fh
				vtx_dist = vtx_dist_fh
				trackState = ts_FirstHit
				hitsInTrack = track.getTrackerHits()
				if vtx_dist_lh < vtx_dist_fh:
					ang_dist = dist_lh
					vtx_dist = vtx_dist_lh
					trackState = ts_LastHit
					hitsInTrack = utils.invert(hitsInTrack)

				# distance between 1st and 2nd hits
				if trackRel.getTo() == particle and trackRel.getWeight() >= 1.0 \
				and trackRel.getTo().getCharge() * utils.getTrackCharge(trackState) > 0:

					firstSimHit = 0
					secondSimHit = 0
					idx = 0
					while firstSimHit == 0:
						for rel in TPCHitRelations:
							if hitsInTrack[idx] == rel.getFrom():
								firstSimHit = rel.getTo()
							if hitsInTrack[idx+1] == rel.getFrom():
								secondSimHit = rel.getTo()
							if firstSimHit != 0 and secondSimHit != 0:
								break
						idx += 1
						if idx+1 == len(hitsInTrack):
						#if idx+1 == hitsInTrack.size():
							break
					if firstSimHit and secondSimHit:
						d12 = utils.spacialDistance(firstSimHit.getPosition(), secondSimHit.getPosition())
					else:
						d12 = 0
						print i, ' event: one of the SimTrackerHits not found!'
						print 'firstSimHit:', firstSimHit
						print 'secondSimHit:',secondSimHit
						print ''

					#if d12 > 100:
					#	print i, "d12:", d12

					if firstSimHit != 0 and secondSimHit != 0:
						d12vsTheta_true.Fill(mclepP4.Theta(),d12)
						d12vsDeltaPhi_true.Fill(deltaPhi,d12)

				'''
				if ang_dist > 0.4:
					print i, " MC track vertex:", particle.getVertex()[0],particle.getVertex()[1],particle.getVertex()[2]
					print i, " MarlinTrkTrack at FirstHit:", ts_FirstHit.getReferencePoint()[0],ts_FirstHit.getReferencePoint()[1],ts_FirstHit.getReferencePoint()[2]
					print i, " MarlinTrkTrack at LastHit:", ts_LastHit.getReferencePoint()[0],ts_LastHit.getReferencePoint()[1],ts_LastHit.getReferencePoint()[2]
					print i, mclepP4.Px(), mclepP4.Py(), mclepP4.Pz()
					print i, p_fh.Px(), p_fh.Py(), p_fh.Pz()
					print i, p_lh.Px(), p_lh.Py(), p_lh.Pz()
					print ''
				'''
				if trackRel.getTo() == particle \
				and trackRel.getTo().getCharge() * utils.getTrackCharge(trackState) > 0 \
				and ang_dist < 0.2:
				#and vtx_dist < 100.:
					'''
					print i, " MC track vertex:", particle.getVertex()[0],particle.getVertex()[1],particle.getVertex()[2]
					#print i, " MarlinTrkTrack at IP:", track.getReferencePoint()[0],track.getReferencePoint()[1],track.getReferencePoint()[2], track.getD0(), track.getZ0(), track.getPhi()
					print i, " MarlinTrkTrack at FirstHit:", ts_FirstHit.getReferencePoint()[0],ts_FirstHit.getReferencePoint()[1],ts_FirstHit.getReferencePoint()[2], ts_FirstHit.getD0(), ts_FirstHit.getZ0(), ts_FirstHit.getPhi()
					print i, " MarlinTrkTrack at LastHit:", ts_LastHit.getReferencePoint()[0],ts_LastHit.getReferencePoint()[1],ts_LastHit.getReferencePoint()[2], ts_LastHit.getD0(), ts_LastHit.getZ0(), ts_LastHit.getPhi()
					print i, mclepP4.Px(), mclepP4.Py(), mclepP4.Pz()
					print i, p_fh.Px(), p_fh.Py(), p_fh.Pz()
					print i, p_lh.Px(), p_lh.Py(), p_lh.Pz()
					'''

					if useTracks == 'MarlinTrk':
						if cut:
							momTrack_reco.Fill(p)
							ptTrack_reco.Fill(mclepP4.Pt())
							momLLP_reco.Fill( p_llp )
							thetaTrack_reco.Fill(mclepP4.Theta())
							DeltaThetaTrackLLP_reco.Fill( deltaTheta )
							DeltaPhiTrackLLP_reco.Fill( abs(deltaPhi) )
							DeltaAlphaTrackLLP_reco.Fill( abs(deltaAlpha) )
							trackEffiPvsTheta_reco.Fill(mclepP4.Theta(),p)
							trackEffVsR_reco.Fill(trueR)
							trackEffVsL_reco.Fill(trueL)
							trackEffRVsZ_reco.Fill(trueZ,trueR)



					vMatchingTracks.append( (trackState, ang_dist) )

					allMatchingTracks+=1
					matchingTracksInEvent+=1
					# print i, particle.getPDG(), matchingTracksInEvent
					break

	for track in trackCollection:
		momentum = utils.getTrackMomentum(track)
		trackP3 = ROOT.TVector3()
		trackP3.SetXYZ(momentum[0],momentum[1],momentum[2])

		firstSimHit  = 0
		secondSimHit = 0
		idx = 0
		hitsInTrack = track.getTrackerHits()
		while firstSimHit == 0:
			for rel in TPCHitRelations:
				if hitsInTrack[idx] == rel.getFrom():
					firstSimHit = rel.getTo()
				if hitsInTrack[idx+1] == rel.getFrom():
					secondSimHit = rel.getTo()
				if firstSimHit != 0 and secondSimHit != 0:
					break
			idx += 1
			if idx+1 == hitsInTrack.size():
				break
		if firstSimHit and secondSimHit:
			d12 = utils.spacialDistance(firstSimHit.getPosition(), secondSimHit.getPosition())
		'''
		else:
			d12 = 0
			print i, ' event: at least one of the SimTrackerHits in track not found!'
			print 'firstSimHit:', firstSimHit
			print 'secondSimHit:',secondSimHit
			print ''
		'''
		if firstSimHit !=0 and secondSimHit !=0:
			d12vsTheta_reco.Fill(trackP3.Theta(),d12)

	#########################################
	############### REFITTED ################
	#########################################

	vMatchingRefittedTracks = []
	matchingRefittedTracksInEvent = 0
	matchingTracks = 0
	p_llp = 0
	for particle in mcCollection:

		if matchingRefittedTracksInEvent > 1:
			break

		if(particle.getPDG() == llp_pdg and particle.getGeneratorStatus() == gen_status_llp):
			e = particle.getEnergy()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p_llp = R3(px,py,pz)
			mcllpP4.SetPxPyPzE(px,py,pz,e)

		if( particle.getGeneratorStatus() == gen_status):
			if abs(particle.getPDG()) != pdg1 and abs(particle.getPDG()) != pdg2:
				continue
			vtx_mc = particle.getVertex()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = R3(px,py,pz)
			e = particle.getEnergy()
			mclepP4.SetPxPyPzE(px,py,pz,e)
			cosTheta = abs( math.cos(mclepP4.Theta()) )
			deltaTheta = abs( mcllpP4.Theta() - mclepP4.Theta() )
			deltaPhi = utils.phiDistance( mcllpP4.Phi(), mclepP4.Phi() )
			deltaAlpha = mcllpP4.Angle( mclepP4.Vect() )

			if cat == 'sm' and abs(mclepP4.Eta()) > 2.4:
				continue

			cut = True

			for trackRel in refitTrackToMCLinkCollection:

				track = trackRel.getFrom()

				ts_FirstHit = track.getTrackState(2) # 1 = atIP, 2 = atFirstHit, 3 = atLastHit
				ts_FirstHit = utils.fixTrackStateDirection( track, ts_FirstHit.getLocation() )
				momentum = utils.getTrackMomentum(ts_FirstHit)
				p_fh = ROOT.TVector3()
				p_fh.SetXYZ(momentum[0],momentum[1],momentum[2])
				dist_fh = utils.angularDistance( p_fh.Theta(), mclepP4.Theta(), p_fh.Phi(), mclepP4.Phi()  )
				vtx_fh = ts_FirstHit.getReferencePoint()
				vtx_dist_fh = utils.spacialDistance(vtx_mc, vtx_fh)

				ts_LastHit = track.getTrackState(3) # 1 = atIP, 2 = atFirstHit, 3 = atLastHit
				ts_LastHit = utils.fixTrackStateDirection( track, ts_LastHit.getLocation() )
				momentum = utils.getTrackMomentum(ts_LastHit)
				p_lh = ROOT.TVector3()
				p_lh.SetXYZ(momentum[0],momentum[1],momentum[2])
				dist_lh = utils.angularDistance( p_lh.Theta(), mclepP4.Theta(), p_lh.Phi(), mclepP4.Phi()  )
				vtx_lh = ts_LastHit.getReferencePoint()
				vtx_dist_lh = utils.spacialDistance(vtx_mc, vtx_lh)

				ang_dist = dist_fh
				vtx_dist = vtx_dist_fh
				trackState = ts_FirstHit
				if vtx_dist_lh < vtx_dist_fh:
					ang_dist = dist_lh
					vtx_dist = vtx_dist_lh
					trackState = ts_LastHit
				h_angDist.Fill(ang_dist)

				hits = track.getTrackerHits()
				if hits.size() < 4:
					print ''
					print i, "REFITTED TRACK CONTAINS < 4 HITS"
					print ''

				vMatchingRefittedTracks.append( (trackState, ang_dist) )

				if trackRel.getTo() == particle \
				and trackRel.getTo().getCharge() * utils.getTrackCharge(trackState) > 0 \
				and ang_dist < 0.2:
				#and vtx_dist < 100.:


					if useTracks == 'Refitted':
						if cut:
							momTrack_reco.Fill(p)
							ptTrack_reco.Fill(mclepP4.Pt())
							momLLP_reco.Fill( p_llp )
							thetaTrack_reco.Fill(mclepP4.Theta())
							DeltaThetaTrackLLP_reco.Fill( deltaTheta )
							DeltaPhiTrackLLP_reco.Fill( abs(deltaPhi) )
							DeltaAlphaTrackLLP_reco.Fill( abs(deltaAlpha) )
							trackEffiPvsTheta_reco.Fill(mclepP4.Theta(),p)
							trackEffVsR_reco.Fill(trueR)
							trackEffVsL_reco.Fill(trueL)
							trackEffRVsZ_reco.Fill(trueZ,trueR)
							#print i, deltaPhi

					allMatchingTracksRefitted+=1
					matchingRefittedTracksInEvent+=1
					break
				'''
				elif trackRel.getTo() == particle \
				and trackRel.getTo().getCharge() * utils.getTrackCharge(trackState) > 0 \
				and deltaPhi > 1.5 and deltaPhi < 2.:

					print i, " deltaPhi:", deltaPhi, "angDist:", ang_dist
					print i, " MC track vertex:	    ", particle.getVertex()[0],particle.getVertex()[1],particle.getVertex()[2]
					print i, " RefitTrack at FirstHit:", ts_FirstHit.getReferencePoint()[0],ts_FirstHit.getReferencePoint()[1],ts_FirstHit.getReferencePoint()[2]
					print i, " RefitTrack at LastHit: ", ts_LastHit.getReferencePoint()[0],ts_LastHit.getReferencePoint()[1],ts_LastHit.getReferencePoint()[2]
					print i, " MC track P_i:			", mclepP4.Px(), mclepP4.Py(), mclepP4.Pz()
					print i, " RefitTrack track P_i (FirstHit):", p_fh.Px(), p_fh.Py(), p_fh.Pz()
					print i, " RefitTrack track P_i (LastHit): ", p_lh.Px(), p_lh.Py(), p_lh.Pz()
					hits = track.getTrackerHits()
					#for hit in hits:
					#	print (hit.getPosition()[0],hit.getPosition()[1],hit.getPosition()[2]),
					print ''
				'''
	closest_dists = []
	for idx in range( refitTrackCollection.getNumberOfElements() - 1 ):
		firstTrack = refitTrackCollection.getElementAt(idx)
		for jdx in range(idx+1, refitTrackCollection.getNumberOfElements() ):
			secondTrack = refitTrackCollection.getElementAt(jdx)

			ftFH = firstTrack.getTrackState(2).getReferencePoint()
			ftLH = firstTrack.getTrackState(3).getReferencePoint()
			stFH = secondTrack.getTrackState(2).getReferencePoint()
			stLH = secondTrack.getTrackState(3).getReferencePoint()

			dists = []
			dists.append( utils.spacialDistance(ftFH,stFH) )
			dists.append( utils.spacialDistance(ftFH,stLH) )
			dists.append( utils.spacialDistance(ftLH,stFH) )
			dists.append( utils.spacialDistance(ftLH,stFH) )
			closest_dists.append( min(dists) )

	for dist in closest_dists:
		hitsDistances.Fill( dist )
	#if mclepP4.Theta() > 1.2 and mclepP4.Theta() > 1.2 and matchingRefittedTracksInEvent < 2:
	#	print i
	'''
	if matchingRefittedTracksInEvent < 2:
		print matchingRefittedTracksInEvent, " matching Refitted Tracks in event ", i
		print ''

	if matchingTracksInEvent < matchingRefittedTracksInEvent:
		onlyRefittedTracks += matchingRefittedTracksInEvent - matchingTracksInEvent
		print i, str(matchingTracksInEvent) + " matching MarlinTrkTracks and " + str(matchingRefittedTracksInEvent) + " RefittedTracks in event " + str(i)
		print ''
		for particle in mcCollection:
			if ( particle.getGeneratorStatus() == gen_status):
				if abs(particle.getPDG()) != pdg1 and abs(particle.getPDG()) != pdg2:
					continue
				print i, 'MC particle momentum:    ', particle.getMomentum()[0],particle.getMomentum()[1],particle.getMomentum()[2], ", vertex: ", \
				particle.getVertex()[0],particle.getVertex()[1],particle.getVertex()[2]
		for (trackState, ang_dist) in vMatchingTracks:
			momentum = utils.getTrackMomentum(trackState)
			print i, 'MarlinTrk momentum:      ', momentum[0],momentum[1],momentum[2], ", vertex: ", \
			trackState.getReferencePoint()[0], trackState.getReferencePoint()[1], trackState.getReferencePoint()[2], ", ang_dist: ", ang_dist
		for (trackState, ang_dist) in vMatchingRefittedTracks:
			momentum = utils.getTrackMomentum(trackState)
			print i, 'Refitted Track momentum: ', momentum[0],momentum[1],momentum[2], ", vertex: ", \
			trackState.getReferencePoint()[0], trackState.getReferencePoint()[1], trackState.getReferencePoint()[2], ", ang_dist: ", ang_dist
		print ''
	'''
	'''
	if matchingTracksInEvent > matchingRefittedTracksInEvent or (matchingRefittedTracksInEvent < 2):
		onlyMarlinTrkTracks += matchingTracksInEvent - matchingRefittedTracksInEvent
		print i, str(matchingTracksInEvent) + " matching MarlinTrkTracks and " + str(matchingRefittedTracksInEvent) + " RefittedTracks in event " + str(i)
		print ''
		for particle in mcCollection:
			if ( particle.getGeneratorStatus() == gen_status):
				if abs(particle.getPDG()) != pdg1 and abs(particle.getPDG()) != pdg2:
					continue
				print i, 'MC particle momentum:    ', particle.getMomentum()[0],particle.getMomentum()[1],particle.getMomentum()[2], ", vertex: ", \
				particle.getVertex()[0],particle.getVertex()[1],particle.getVertex()[2]
		for (trackState, ang_dist) in vMatchingTracks:
			momentum = utils.getTrackMomentum(trackState)
			print i, 'MarlinTrk momentum:      ', momentum[0],momentum[1],momentum[2], ", vertex: ", \
			trackState.getReferencePoint()[0], trackState.getReferencePoint()[1], trackState.getReferencePoint()[2], ", ang_dist: ", ang_dist
		for (trackState, ang_dist) in vMatchingRefittedTracks:
			momentum = utils.getTrackMomentum(trackState)
			print i, 'Refitted Track momentum: ', momentum[0],momentum[1],momentum[2], ", vertex: ", \
			trackState.getReferencePoint()[0], trackState.getReferencePoint()[1], trackState.getReferencePoint()[2], ", ang_dist: ", ang_dist
		print ''
	'''

	#########################################
	############### VERTICES ################
	#########################################
	if 'LLPVerticesTest' in event.getCollectionNames():
		for vtx in vertexCollection:
			dist = utils.spacialDistance(llp_endpoint, vtx.getPosition() )
			vtxMisdistances.Fill(dist)

reader.close()

#print('Decays within tracker acceptance', 'Ev. with min. two tracks', 'Ev. with min. two reco. tracks', 'All reco. tracks', 'Tracks matching to MC', 'Two reco. leptons')
#print(withinAcceptance, twoTruthTracks, twoTracks, allTracks, allMatchingTracks, twoRecoLeptons)

print ''
print 'Tracks within acceptance: ' + str(withinAcceptance)
print 'Ev. with two truth tracks: ' + str(twoTruthTracks)
print 'Ev. with min. two reco. tracks: ' + str(twoTracks)
print 'All reco. tracks: ' + str(allTracks)
print 'Number of PFOs: ' + str(nPFOs)
print 'Tracks matching to MC (MarlinTrk): ' + str(allMatchingTracks)
print 'Tracks matching to MC (Refit): ' + str(allMatchingTracksRefitted)
print 'Tracks matching only with MarlinTrk: ' + str(onlyMarlinTrkTracks)
print 'Tracks matching only with Refit: ' + str(onlyRefittedTracks)
print 'MarlinTrk reco. eff.: ', float(allMatchingTracks) / float(withinAcceptance)
print 'Refitted reco. eff.:  ', float(allMatchingTracksRefitted) / float(withinAcceptance)
print 'Ev. with >=1 sec. vtx:  ', nSecVertices
print ''

cutText = ' (True)'

ROOT.gStyle.SetOptStat(0)

c0=ROOT.TCanvas('c0', 'c0', 600, 400)

trackEffVsR = ROOT.TH1F('trackEffVsR', 'Track reconstruction efficiency', nbins, 0., TPCRmax)
trackEffVsR.SetXTitle('True vertex distance from beam axis')
trackEffVsR.SetYTitle('Track reco. efficiency')
trackEffVsR.Divide(trackEffVsR_reco, trackEffVsR_true,1,1,"cl=0.683 b(1,1) mode")
trackEffVsR.SetMarkerStyle(20)
trackEffVsR.SetMarkerSize(0.5)
trackEffVsR.Draw('e')
c0.SaveAs('plots/'+cat+'/'+fs+'/trackEffVsR_'+fs+'.pdf')

trackEffVsL = ROOT.TH1F('trackEffVsL', 'Track reconstruction efficiency', 30, 0., 3000.)
trackEffVsL.SetXTitle('#beta c t_{LLP} [mm]')
trackEffVsL.SetYTitle('Track reco. efficiency')
trackEffVsL.Divide(trackEffVsL_reco, trackEffVsL_true,1,1,"cl=0.683 b(1,1) mode")
trackEffVsL.SetMarkerStyle(20)
trackEffVsL.SetMarkerSize(0.5)
trackEffVsL.Draw('e')
c0.SaveAs('plots/'+cat+'/'+fs+'/trackEffVsL_'+fs+'.pdf')

momTrack = ROOT.TEfficiency(momTrack_reco,momTrack_true)
momTrack.SetTitle('Track reconstruction efficiency;True momentum [GeV];Efficiency')
momTrack.SetMarkerStyle(20)
momTrack.SetMarkerSize(0.5)
momTrack.Draw('ap')
c0.SaveAs('plots/'+cat+'/'+fs+'/momTrack_'+fs+'.pdf')

ptTrack = ROOT.TEfficiency(ptTrack_reco,ptTrack_true)
ptTrack.SetTitle('Track reconstruction efficiency;True p_{T} [GeV];Efficiency')
ptTrack.SetMarkerStyle(20)
ptTrack.SetMarkerSize(0.5)
ptTrack.Draw('ap')
c0.SaveAs('plots/'+cat+'/'+fs+'/ptTrack_'+fs+'.pdf')

momLLP = ROOT.TEfficiency(momLLP_reco,momLLP_true)
momLLP.SetTitle('Track reconstruction efficiency;True LLP momentum [GeV];Efficiency')
momLLP.SetMarkerStyle(20)
momLLP.SetMarkerSize(0.5)
momLLP.Draw('ap')
c0.SaveAs('plots/'+cat+'/'+fs+'/momLLP_'+fs+'.pdf')

thetaTrack = ROOT.TEfficiency(thetaTrack_reco,thetaTrack_true)
thetaTrack.SetTitle('Track reconstruction efficiency;True polar angle;Efficiency')
thetaTrack.SetMarkerStyle(20)
thetaTrack.SetMarkerSize(0.5)
thetaTrack.Draw('ap')
c0.SaveAs('plots/'+cat+'/'+fs+'/thetaTrack_'+fs+'.pdf')
'''
print 'DeltaTheta: ', DeltaThetaTrackLLP_reco.GetEntries(), DeltaThetaTrackLLP_true.GetEntries()
DeltaThetaTrackLLP = ROOT.TEfficiency(DeltaThetaTrackLLP_reco,DeltaThetaTrackLLP_true)
DeltaThetaTrackLLP.SetTitle('Track reconstruction efficiency;|#theta_{LLP} - #theta_{track}|;Efficiency')
DeltaThetaTrackLLP.SetMarkerStyle(20)
DeltaThetaTrackLLP.SetMarkerSize(0.5)
DeltaThetaTrackLLP.Draw('ap')
c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaThetaTrackLLP_'+fs+'.pdf')

print 'DeltaPhi: ', DeltaPhiTrackLLP_reco.GetEntries(), DeltaPhiTrackLLP_true.GetEntries()
DeltaPhiTrackLLP = ROOT.TEfficiency(DeltaPhiTrackLLP_reco,DeltaPhiTrackLLP_true)
DeltaPhiTrackLLP.SetTitle('Track reconstruction efficiency;|#phi_{LLP} - #phi_{track}|;Efficiency')
DeltaPhiTrackLLP.SetMarkerStyle(20)
DeltaPhiTrackLLP.SetMarkerSize(0.5)
DeltaPhiTrackLLP.Draw('ap')
c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaPhiTrackLLP_'+fs+'.pdf')

print 'DeltaAlpha: ', DeltaAlphaTrackLLP_reco.GetEntries(), DeltaAlphaTrackLLP_true.GetEntries()
DeltaAlphaTrackLLP = ROOT.TEfficiency(DeltaAlphaTrackLLP_reco,DeltaAlphaTrackLLP_true)
DeltaAlphaTrackLLP.SetTitle('Track reconstruction efficiency;|#Delta#alpha_{LLP-track}|;Efficiency')
DeltaAlphaTrackLLP.SetMarkerStyle(20)
DeltaAlphaTrackLLP.SetMarkerSize(0.5)
DeltaAlphaTrackLLP.Draw('ap')
c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaAlphaTrackLLP_'+fs+'.pdf')
'''
trackEffVsR_true.Scale( 1./trackEffVsR_true.Integral() )
trackEffVsR_true.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/R_true_'+fs+'.pdf')
trackEffVsR_reco.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/R_reco_'+fs+'.pdf')
momTrack_true.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/momTrack_true_'+fs+'.pdf')
thetaTrack_true.Scale( 1./thetaTrack_true.Integral() )
thetaTrack_true.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/thetaTrack_true_'+fs+'.pdf')
DeltaThetaTrackLLP_true.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaThetaTrackLLP_true_'+fs+'.pdf')
DeltaThetaTrackLLP_reco.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaThetaTrackLLP_reco_'+fs+'.pdf')
DeltaPhiTrackLLP_true.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaPhiTrackLLP_true_'+fs+'.pdf')
DeltaAlphaTrackLLP_true.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaAlphaTrackLLP_true_'+fs+'.pdf')
h_nTracks.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/nTracks_'+fs+'.pdf')
h_relWeights.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/relWeights_'+fs+'.pdf')
h_angDist.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/angularDistance_'+fs+'.pdf')
vtxMisdistances.Draw('hist')
vtxMisdistances.SetXTitle('#Delta R_{vtx} [mm]')
c0.SetLogy()
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxMisdistance_'+fs+'.pdf')
vtxMultiplicity.Draw('hist')
vtxMultiplicity.SetXTitle('N_{vtx}')
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxMultiplicity_'+fs+'.pdf')
hitsDistances.Draw('hist')
hitsDistances.SetXTitle('l_{ij} [mm]')
c0.SaveAs('plots/'+cat+'/'+fs+'/hitsDistances_'+fs+'.pdf')

c1=ROOT.TCanvas('c1', 'c1', 800, 400)
c1.SetRightMargin(0.15)
'''
trackEffRVsZ = ROOT.TH2F('trackEffRVsZ', 'Track reconstruction efficiency', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffRVsZ.SetYTitle('R [mm]')
trackEffRVsZ.SetXTitle('z [mm]')
trackEffRVsZ.SetZTitle('Track reco. efficiency')
trackEffRVsZ.Divide(trackEffRVsZ_reco, trackEffRVsZ_true,1,1,"cl=0.683 b(1,1) mode")
trackEffRVsZ.SetStats(0)
'''
trackEffRVsZ = ROOT.TEfficiency(trackEffRVsZ_reco,trackEffRVsZ_true)
trackEffRVsZ.SetTitle('Track reconstruction efficiency;z [mm];R [mm];Efficiency')
trackEffRVsZ.Draw('colz')
line1 = ROOT.TLine(-2350.,329.0,2350.0,329.0);
line2 = ROOT.TLine(-2350.,1808.0,2350.0,1808.0);
line3 = ROOT.TLine(-2350.,329.0,-2350.0,1808.0);
line4 = ROOT.TLine(2350.,329.0,2350.0,1808.0);
line1.SetLineColor(2);
line2.SetLineColor(2);
line3.SetLineColor(2);
line4.SetLineColor(2);
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/trackEffRVsZ_'+fs+'.pdf')

trackEffiPvsTheta = ROOT.TEfficiency(trackEffiPvsTheta_reco,trackEffiPvsTheta_true)
trackEffiPvsTheta.SetTitle('Track reconstruction efficiency;#theta_{track};Momentum;Efficiency')
trackEffiPvsTheta.Draw('colz text')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/trackEffiPvsTheta_'+fs+'.pdf')

trackEffRVsZ_true.SetYTitle('R [mm]')
trackEffRVsZ_true.SetXTitle('z [mm]')
trackEffRVsZ_true.SetZTitle('Number of vertices')
trackEffRVsZ_true.SetStats(0)
trackEffRVsZ_true.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/trueVerticesRVsZ_'+fs+'.pdf')

recoVertices.SetYTitle('R [mm]')
recoVertices.SetXTitle('z [mm]')
recoVertices.SetZTitle('Number of vertices')
recoVertices.SetStats(0)
recoVertices.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/recoVerticesRVsZ_'+fs+'.pdf')

trackEffiPvsTheta_true.SetYTitle('p [GeV]')
trackEffiPvsTheta_true.SetXTitle('#theta_{track}')
trackEffiPvsTheta_true.SetZTitle('Number of tracks')
trackEffiPvsTheta_true.SetStats(0)
#trackEffiPvsTheta_true.Scale( 1./trackEffiPvsTheta_true.Integral() )
trackEffiPvsTheta_true.Draw('colz text')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/trackEffiPvsTheta_true_'+fs+'.pdf')

trackEffiPvsTheta_reco.SetYTitle('p [GeV]')
trackEffiPvsTheta_reco.SetXTitle('#theta_{track}')
trackEffiPvsTheta_reco.SetZTitle('Number of tracks')
trackEffiPvsTheta_reco.SetStats(0)
#trackEffiPvsTheta_reco.Scale( 1./trackEffiPvsTheta_reco.Integral() )
trackEffiPvsTheta_reco.Draw('colz text')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/trackEffiPvsTheta_reco_'+fs+'.pdf')

trackEffiPvsR_true.SetYTitle('R [mm]')
trackEffiPvsR_true.SetXTitle('p [GeV]')
trackEffiPvsR_true.SetZTitle('Number of tracks')
trackEffiPvsR_true.SetStats(0)
trackEffiPvsR_true.Draw('colz text')
c1.SaveAs('plots/'+cat+'/'+fs+'/trackEffiPvsR_true_'+fs+'.pdf')

d12vsTheta_true.SetYTitle('d_{12} [mm]')
d12vsTheta_true.SetXTitle('Track polar angle (truth)')
d12vsTheta_true.SetZTitle('Number of tracks')
d12vsTheta_true.SetStats(0)
d12vsTheta_true.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/d12vsTheta_true_'+fs+'.pdf')

d12vsTheta_reco.SetYTitle('d_{12} [mm]')
d12vsTheta_reco.SetXTitle('Track polar angle')
d12vsTheta_reco.SetZTitle('Number of tracks')
d12vsTheta_reco.SetStats(0)
d12vsTheta_reco.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/d12vsTheta_reco_'+fs+'.pdf')

d12vsDeltaPhi_true.SetYTitle('d_{12} [mm]')
d12vsDeltaPhi_true.SetXTitle('|#phi_{LLP} - #phi_{track}|')
d12vsDeltaPhi_true.SetZTitle('Number of tracks')
d12vsDeltaPhi_true.SetStats(0)
d12vsDeltaPhi_true.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/d12vsDeltaPhi_true_'+fs+'.pdf')

histo_array = []
histo_array.append(trackEffVsL)
histo_array.append(ptTrack)
histo_array.append(thetaTrack)
#histo_array.append(DeltaAlphaTrackLLP)
#histo_array.append(DeltaThetaTrackLLP)
#histo_array.append(DeltaPhiTrackLLP)

try:
	arg2 = sys.argv[2]
	outFileName = arg2
	outHistFile = ROOT.TFile.Open(outFileName ,"RECREATE")
	outHistFile.cd()
	for hist in histo_array:
		hist.Write()
	outHistFile.Close()
except IndexError:
    pass

'''
# save histograms to output file
if sys.argv[2]:
	outFileName = sys.argv[2]
	outHistFile = ROOT.TFile.Open(outFileName ,"RECREATE")
	outHistFile.cd()

	for hist in histo_array:
		hist.Write()

	outHistFile.Close()
else:
	reader.close()
'''
reader.close()
