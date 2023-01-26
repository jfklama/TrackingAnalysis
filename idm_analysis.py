import sys
import os

import math
import numpy as np

from pyLCIO import IOIMPL, EVENT, IMPL, ROOT
#from ROOT import Math

import utils

def R(x,y):
	return math.sqrt(x**2 + y**2)
def R3(x,y,z):
	return math.sqrt(x**2 + y**2 + z**2)

try:
    arg1 = sys.argv[1]
except IndexError:
    print "Usage: " + os.path.basename(__file__) + " <inputFile(s)>" + " <outputHistFile>"
    sys.exit(1)

print 'Start analysis'

chargedChannel    		  = 0
withinAcceptance  		  = 0
withinTPC		  		  = 0
twoTracks         		  = 0
twoTruthTracks   		  = 0
allTracks		  		  = 0
allMatchingTracks 		  = 0
allMatchingTracksRefitted = 0
allMatchingVertices		  = 0
matchingVerticesTPC		  = 0
twoRecoLeptons    		  = 0
onlyMarlinTrkTracks		  = 0
onlyRefittedTracks		  = 0
nPFOs					  = 0
nCorrPID				  = 0
nEvSecVertices		  	  = 0
nSecVertices		  	  = 0
nSecVerticesTPC		  	  = 0
nFakeVertices		  	  = 0

nbins = 30
TPCRmax = 1974.
TPCZmax = 2350.
# TPCRmax = 100.
# TPCZmax = 200.

# track efficiency
trackEffVsR_true = ROOT.TH1F('trackEffVsR_true', 'True vertices', nbins, 0., TPCRmax)
trackEffVsR_reco = ROOT.TH1F('trackEffVsR_reco', 'N matching reco. tracks', nbins, 0., TPCRmax)
trackEffVsL_true = ROOT.TH1F('trackEffVsL_true', 'True vertices', 30, 0., 3000.)
trackEffVsL_reco = ROOT.TH1F('trackEffVsL_reco', 'N matching reco. tracks', 30, 0., 3000.)

# 2-D track efficiency
trackEffRVsZ_true = ROOT.TH2F('trackEffRVsZ_true', 'True vertices', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffRVsZ_reco = ROOT.TH2F('trackEffRVsZ_reco', 'N matching reco. tracks', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffiPvsTheta_true = ROOT.TH2F('trackEffiPvsTheta_true', '', 3, 0., 3.15, 15, 0., 15.)
trackEffiPvsTheta_reco = ROOT.TH2F('trackEffiPvsTheta_reco', '', 3, 0., 3.15, 15, 0., 15.)
trackEffiPvsR_true = ROOT.TH2F('trackEffiPvsR_true', '', 15, 0., 15., nbins, 0., TPCRmax)

# first two hits distances
d12vsTheta_reco = ROOT.TH2F('d12vsTheta_reco', '', 32, 0., 3.15, 32, 0., 150.)
d12vsTheta_true = ROOT.TH2F('d12vsTheta_true', '', 32, 0., 3.15, 32, 0., 150.)
d12vsDeltaPhi_true = ROOT.TH2F('d12vsDeltaPhi_true', '', 32, 0., 3.15, 32, 0., 150.)
#d12vsDeltaPhi_reco = ROOT.TH2F('d12vsDeltaPhi_reco', '', 32, 0., 3.15, 32, 0., 100.)

# track eff. vs p, pT, theta
momTrack_true = ROOT.TH1F('momTrack_true', 'True momentum of track', 50, 0., 15.)
momTrack_reco = ROOT.TH1F('momTrack_reco', 'True momentum of reco. track', 50, 0., 15.)
ptTrack_true = ROOT.TH1F('ptTrack_true', 'True pT of track', 50, 0., 15.)
ptTrack_reco = ROOT.TH1F('ptTrack_reco', 'True pT of reco. track', 50, 0., 15.)
thetaTrack_true = ROOT.TH1F('thetaTrack_true', 'True polar angle of track', 42, 0., 3.15)
thetaTrack_reco = ROOT.TH1F('thetaTrack_reco', 'True polar angle of reco. track', 42, 0., 3.15)

# track eff. vs decay angle
DeltaThetaTrackLLP_true = ROOT.TH1F('DeltaThetaTrackLLP_true', 'Dist. between polar angles of track and LLP', 35, 0., 3.5)
DeltaThetaTrackLLP_reco = ROOT.TH1F('DeltaThetaTrackLLP_reco', 'Dist. between polar angles of track and LLP', 35, 0., 3.5)
DeltaPhiTrackLLP_true = ROOT.TH1F('DeltaPhiTrackLLP_true', 'Dist. between azimuthal angles of track and LLP', 35, 0., 3.5)
DeltaPhiTrackLLP_reco = ROOT.TH1F('DeltaPhiTrackLLP_reco', 'Dist. between azimuthal angles of track and LLP', 35, 0., 3.5)
DeltaAlphaTrackLLP_true = ROOT.TH1F('DeltaAlphaTrackLLP_true', 'Angle between track and LLP', 35, 0., 3.5)
DeltaAlphaTrackLLP_reco = ROOT.TH1F('DeltaAlphaTrackLLP_reco', 'Angle between track and LLP', 35, 0., 3.5)

# LLP momentum and theta
momLLP_true = ROOT.TH1F('momLLP_true', 'True momentum of LLP', 50, 0., 300.)
momLLP_reco = ROOT.TH1F('momLLP_reco', 'True momentum of reco. LLP', 50, 0., 300.)
thetaLLP_true = ROOT.TH1F('thetaLLP_true', 'True polar angle of LLP', 42, 0., 3.15)
thetaLLP_reco = ROOT.TH1F('thetaLLP_reco', 'True polar angle of reco. LLP', 42, 0., 3.15)

# track reco. features
h_nTracks = ROOT.TH1F('h_nTracks', 'Number of tracks', 20, 0., 20.)
h_relWeights = ROOT.TH1F('h_relWeights', 'Number of weights', 50, 0., 2.)
h_angDist = ROOT.TH1F('h_angDist', 'Angular separation', 70, 0., 7)

# vertices reconstruction
vtxMisdistances = ROOT.TH1F('vtxMisdistances', 'Distance between true and reco. vertex', 50, 0., 200)
hitsDistances = ROOT.TH1F('hitsDistances', 'Distance between two closest first/last hits', 200, 0., 800)
helixDistances = ROOT.TH1F('helixDistances', 'Distance between helices for a vertex', 50, 0., 25)
vtxMultiplicity = ROOT.TH1F('vtxMultiplicity', 'Reco. vertex multiplicity', 10, 0., 10)

# vertex reco. eff.
vtxEffVsL_true = ROOT.TH1F('vtxEffVsL_true', 'True vertices', 30, 0., 3000.)
vtxEffVsL_reco = ROOT.TH1F('vtxEffVsL_reco', 'N matching reco. vertices', 30, 0., 3000.)
vtxVsR_reco = ROOT.TH1F('vtxVsR_reco', 'Reco. vertices', 100, 0., 2000.)
vtxVsL_reco = ROOT.TH1F('vtxVsL_reco', 'Reco. vertices', 75, 0., 300.)

logx = np.logspace(-1.0, 3.7, num=51)
logx2 = np.logspace(-1.0, 3.5, num=51)

# 2-D vertices plots
recoVertices = ROOT.TH2F('recoVertices', 'Number of found vertices', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
fakeVertices = ROOT.TH2F('fakeVertices', 'Vertices not matched to MC', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
vtxEffRVsZ_true = ROOT.TH2F('vtxEffRVsZ_true', 'True vertices', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
vtxEffRVsZ_reco = ROOT.TH2F('vtxEffRVsZ_reco', 'True vtx corresponding to reco. vtx', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
vtxPt1vsPt2 = ROOT.TH2F('vtxPt1vsPt2', 'pT of tracks at reco. vtx', 40, 0.0, 4.0, 40, 0.0, 4.0)
vtxP1vsP2 = ROOT.TH2F('vtxP1vsP2', 'p of tracks at reco. vtx', 40, 0.0, 4.0, 40, 0.0, 4.0)
vtxPz1vsPz2 = ROOT.TH2F('vtxPz1vsPz2', 'p_{z} of tracks at reco. vtx', 40, -4.0, 4.0, 40, -4.0, 4.0)
vtxOmega1vsOmega2 = ROOT.TH2F('vtxOmega1vsOmega2', '#Omega of tracks at reco. vtx', 40, 0.0, 0.02, 40, 0.0, 0.02)
vtxRefPointDist1VsDist2 = ROOT.TH2F('vtxRefPointDist1VsDist2', 'Distance of ref. points from reco. vtx', 50, 0., 2000., 50, 0., 2000.)
vtxTrksNdf1VsNdf2 = ROOT.TH2F('vtxTrksNdf1VsNdf2', 'N_{df} for tracks at reco. vtx', 50, 0, 500, 50, 0, 500)
vtxTrksPtVsNdf = ROOT.TH2F('vtxTrksPtVsNdf', 'p_{T} and N_{df} for tracks at reco. vtx', 50, 0, 15, 50, 0, 500)
tracksVtxPvsTheta = ROOT.TH2F('tracksVtxPvsTheta', '', 16, 0., 3.15, 50, 0., 15.)
tracksVtxPtVsTheta = ROOT.TH2F('tracksVtxPtVsTheta', '', 16, 0., 3.15, 50, 0., 15.)
vtxTrksD01VsD02 = ROOT.TH2F('vtxTrksD01VsD02', 'd_{0} (w.r.t. IP) of tracks in vtx', 50, 0, 2000, 50, 0, 2000)
vtxTrksZ01VsZ02 = ROOT.TH2F('vtxTrksZ01VsZ02', 'z_{0} (w.r.t. IP) of tracks in vtx', 50, 0, 2500, 50, 0, 2500)
vtxTrksCentresRefPointsDist = ROOT.TH2F('vtxTrksCentresRefPointsDist', '', 50, logx, 50, 0, 5000)
vtxTrksCentresRefPointsDist_scaled = ROOT.TH2F('vtxTrksCentresRefPointsDist_scaled', '', 50, logx, 50, 0, 5000)
vtxRefPointDistsVsCentres = ROOT.TH2F('vtxRefPointDistsVsCentres', 'Distance of ref. points from reco. vtx vs. distance between circle centres', 50, logx2, 50, 0., 5000.)
corrCentresRefDistsVsPt = ROOT.TH2F('corrCentresRefDistsVsPt', '', 50, -15000, 3000, 50, 0, 15)
corrLogCentresRefDistsVsPt = ROOT.TH2F('corrLogCentresRefDistsVsPt', '', 50, -15000, 3000, 50, 0, 15)
corrCentresRefDistsVsR = ROOT.TH2F('corrCentresRefDistsVsR', '', 50, -3000, 3000, 50, 0, 2000)
cosOpenAngleVsCurvRatio = ROOT.TH2F('cosOpenAngleVsCurvRatio', '', 50, -1, 1, 50, 0., 1.)
vtxRefPhiVsLhPhi = ROOT.TH2F('vtxRefPhiVsLhPhi', 'Ref. point and last hit phi w.r.t. vtx phi', 50, -0.5, 0.5, 50, -0.5, 0.5)
vtxRefPhiVsArcPhi = ROOT.TH2F('vtxRefPhiVsArcPhi', 'Phi from vtx to ref. pt. vs. arc phi', 50, -6.3, 6.3, 50, 0, 7)
vtxRefZVsLhZ = ROOT.TH2F('vtxRefZVsLhZ', 'Vtx dist. in z from ref. pt. and last hit', 50, -0.5, 1, 50, -1, 1)
vtxRefZVsLenZ = ROOT.TH2F('vtxRefZVsLenZ', 'Vtx dist. from ref. pt. and trk len. in z', 50, -0.5, 1, 50, logx)

# tracks assoc. to vertex
vtxTracksTheta = ROOT.TH1F('vtxTracksTheta', 'Theta of the tracks assoc. to good vertex', 42, 0., 3.15)
vtxTracksPhi = ROOT.TH1F('vtxTracksPhi', 'Phi of the tracks assoc. to good vertex', 84, -3.15, 3.15)
vtxTheta = ROOT.TH1F('vtxTheta', 'Reco. theta of the dilepton system', 42, 0., 3.15)
vtxPhi = ROOT.TH1F('vtxPhi', 'Reco. phi of the dilepton system', 84, -3.15, 3.15)
vtxTracksPt = ROOT.TH1F('vtxTracksPt', 'pT of tracks assoc. to good vertex', 50, 0., 15.)
vtxTracksPtSum = ROOT.TH1F('vtxTracksPtSum', 'pT sum of tracks assoc. to good vertex', 30, 0., 6.)
vtxTracksP = ROOT.TH1F('vtxTracksP', 'Momentum of tracks assoc. to good vertex', 50, 0., 15.)
vtxPt = ROOT.TH1F('vtxPt', 'pT of the dilepton system', 50, 0., 15.)
vtxP = ROOT.TH1F('vtxP', 'Momentum of the dilepton system', 50, 0., 15.)
vtxRefPointDist = ROOT.TH1F('vtxRefPointDist', 'Distance between vtx and track ref. point', 60, 0., 300.)
vtxTrksRefPointsDist = ROOT.TH1F('vtxTrksRefPointsDist', 'Distance between track ref. points in vtx', 80, 0., 2000.)
vtxTrksCurvRatio = ROOT.TH1F('vtxTrksCurvRatio', 'Ratio of the two track curvature', 50, 0., 1.)
vtxTrksOmega = ROOT.TH1F('vtxTrksOmega', '#Omega param. of tracks in vtx', 50, 0, 0.04)
vtxTrksRefPointZ = ROOT.TH1F('vtxTrksRefPointZ', 'Distance in z for tracks in vtx', 50, logx)
vtxTrksChi2 = ROOT.TH1F('vtxTrksChi2', '#chi^{2} of tracks in vtx', 20, 0, 40)
vtxTrksChi2Ndf = ROOT.TH1F('vtxTrksChi2Ndf', '#chi^{2}/N_{df} of tracks in vtx', 50, 0, 20)
vtxTrksNdf = ROOT.TH1F('vtxTrksNdf', 'N_{df} of tracks in vtx', 50, 0, 500)
vtxTrksD0 = ROOT.TH1F('vtxTrksD0', 'd_{0} of tracks in vtx', 50, 0, 50)
vtxTrksZ0 = ROOT.TH1F('vtxTrksZ0', 'z_{0} of tracks in vtx', 50, 0, 50)
vtxTrksD0AtIP = ROOT.TH1F('vtxTrksD0AtIP', 'd_{0} (w.r.t. IP) of tracks in vtx', 50, 0, 25)
vtxTrksZ0AtIP = ROOT.TH1F('vtxTrksZ0AtIP', 'z_{0} (w.r.t. IP) of tracks in vtx', 50, 0, 25)
vtxTrksTruePDG = ROOT.TH1F('vtxTrksTruePDG', 'PDG of tracks in vtx', 212, 0, 212)
vtxTrksAreBothCurlers = ROOT.TH1F('vtxTrksAreBothCurlers', '1 if both tracks in vtx are curlers', 2, 0, 2)
vtxTrksCentresDist = ROOT.TH1F('vtxTrksCentresDist', 'Distance between helix circle centres', 50, 0, 5000)
vtxTrksLogCentresRefDists = ROOT.TH1F('vtxTrksLogCentresRefDists', '', 50, -4000, 3000)

# vertex angular analysis
cosOpenAngleVtx_true = ROOT.TH1F('cosOpenAngleVtx_true', 'MC track opening angle at vtx', 40, -1, 1)
cosOpenAngleVtx_reco = ROOT.TH1F('cosOpenAngle_recoVtx', 'Reco. track opening angle at a good vtx', 40, -1, 1)


print 'Opening input file ' + str(sys.argv[1]) + '...'
infile = sys.argv[1]
reader = IOIMPL.LCFactory.getInstance().createLCReader()
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
	llp_pdg = 0

print "LLP: " + str(llp_pdg) + ", with final state: " + str(pdg1) + ", " + str(pdg2)

i=0
for event in reader:
	i+=1
	#if i>999:
	#	break

	mcCollection = event.getCollection('MCParticlesSkimmed')
	#tpcCollection = event.getCollection('TPCCollection')
	trackCollection = event.getCollection('MarlinTrkTracks')
	trackToMCLinkCollection = event.getCollection('MarlinTrkTracksMCTruthLink')
	refitTrackCollection = event.getCollection('RefittedMarlinTrkTracks')
	refitTrackToMCLinkCollection = event.getCollection('RefittedMarlinTrkTracksMCTruthLink')
	pfoCollection = event.getCollection('PandoraPFOs')
	TPCHitRelations = event.getCollection('TPCTrackerHitRelations')

	if 'LLPVertices' in event.getCollectionNames():
		# vertexCollection = event.getCollection('LLPVertices')
		# vertexTracksRelations = event.getCollection('LLPVtxToTracksLink')
		vertexCollection = event.getCollection('LLPVertices')
		vertexTracksRelations = event.getCollection('LLPVtxToTracksLink')
		# if infileName.find("idm") == -1:
		# 	nSecVertices += int(vertexCollection.getNumberOfElements())
		# 	nEvSecVertices += 1

		#if vertexCollection.getNumberOfElements() > 1:
		#	print i, str(vertexCollection.getNumberOfElements()) + ' vertices'

	nVTX = 0
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

	if (trueR <= 1808) and abs(trueZ) <= 2350:
		withinAcceptance+=2
		if trueR >= 329:
			withinTPC+=2
	# if (trueR < 329) and abs(trueZ) <= 2350:
	# 	withinAcceptance+=2
	# 	withinTPC+=2
	elif infileName.find("idm") != -1: # cut on TPC acceptance only if it's IDM sample
		continue

	vtxEffRVsZ_true.Fill( trueZ, trueR )
	vtxEffVsL_true.Fill(trueL)
	if 'LLPVertices' in event.getCollectionNames():
		nVTX = vertexCollection.getNumberOfElements()
		nSecVertices += int(vertexCollection.getNumberOfElements())
		nEvSecVertices += 1
		if (trueR >= 329 and trueR <= 1808) and abs(trueZ) <= 2350:
			nSecVerticesTPC += 1
	#if nVTX != 3 and nVTX != 6:
	#	continue
	vtxMultiplicity.Fill(nVTX)

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
			p = R3(px,py,pz)
			mcllpP4.SetPxPyPzE(px,py,pz,e)
			momLLP_true.Fill(p,2.)

		if(particle.getGeneratorStatus() == gen_status and (abs( particle.getPDG() ) == pdg1 or abs( particle.getPDG() ) == pdg2)):
			e = particle.getEnergy()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = R3(px,py,pz)
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

	cosOpenAngleVtx_true.Fill( math.cos(mclep1.Angle(mclep2.Vect())) )

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
				'''
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
				if trackRel.getTo() == particle \
				and trackRel.getTo().getCharge() * utils.getTrackCharge(trackState) > 0 \
				and ang_dist < 0.2:
				#and vtx_dist < 100.:

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
	'''
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

		else:
			d12 = 0
			print i, ' event: at least one of the SimTrackerHits in track not found!'
			print 'firstSimHit:', firstSimHit
			print 'secondSimHit:',secondSimHit
			print ''

		if firstSimHit !=0 and secondSimHit !=0:
			d12vsTheta_reco.Fill(trackP3.Theta(),d12)
	'''

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


	#########################################
	############### VERTICES ################
	#########################################

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


	#if closest_dist < 50:
	#	continue

	if 'LLPVertices' in event.getCollectionNames():

		dist0 = 50000
		helixDist0 = 50000
		bestVtx = IMPL.VertexImpl()

		for vtx in vertexCollection:

			helixDist = vtx.getParameters()[0]
			helixDistances.Fill(helixDist )
			if len(closest_dists) > 0:
				closest_dist = min(closest_dists)
				hitsDistances.Fill( closest_dist )
			pos = vtx.getPosition()
			vtxVsR_reco.Fill( R(pos[0],pos[1]) )
			vtxVsL_reco.Fill( R3(pos[0],pos[1],pos[2]) )
			dist = utils.spacialDistance(llp_endpoint, pos )
			#vtxMisdistances.Fill(dist)
			# if dist < 30:
			# 	print i, dist0, dist, helixDist

			if dist < dist0:
				dist0 = dist
				pos0 = pos
				helixDist0 = helixDist
				bestVtx = vtx
			if dist > 30:
				fakeVertices.Fill( pos[2], R(pos[0],pos[1]) )

			'''
			if helixDist < helixDist0:
				dist0 = dist
				pos0 = pos
				helixDist0 = helixDist
			'''
		vtxMisdistances.Fill(dist0)
		if dist0 < 30:
			recoVertices.Fill( pos0[2], R(pos0[0],pos0[1]) )
			vtxEffRVsZ_reco.Fill(trueZ,trueR)
			vtxEffVsL_reco.Fill(trueL)
			allMatchingVertices+=1
			if (trueR >= 329 and trueR <= 1808) and abs(trueZ) <= 2350:
				matchingVerticesTPC+=1
			nFakeVertices += vertexCollection.getNumberOfElements() - 1
		else:
			#print 'no good vtx in ev. ' + str(i) + ' with ' + str(vertexCollection.getNumberOfElements()) + ' vertices'
			nFakeVertices += vertexCollection.getNumberOfElements()


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
					vtxTrksTruePDG.Fill(abs(particle.getPDG()))

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
			vtxTrksAreBothCurlers.Fill( abs(omega1) > curlerOmega and abs(omega2) > curlerOmega )

			trkStates = [trkSt1, trkSt2]
			momenta = []

			if  R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 330 and R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) < 1800:
			# if True:

				# print i, llp_endpoint[0],llp_endpoint[1], llp_endpoint[2]
				# print i, vtx_cand.getPosition()[0],vtx_cand.getPosition()[1], vtx_cand.getPosition()[2], vtx_cand.getParameters()[0]
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

					# print i, mom[0],mom[1],mom[2]

					vtxTracksP.Fill(p.Mag())
					vtxTracksPt.Fill(p.Pt())
					vtxTracksTheta.Fill(p.Theta())
					vtxTracksPhi.Fill(p.Phi())
					vtx_refpoint = utils.spacialDistance( trkSt.getReferencePoint(),vtx_cand.getPosition() )
					vtxRefPointDist.Fill(vtx_refpoint)

					tracksVtxPtVsTheta.Fill(p.Theta(),p.Pt())
					tracksVtxPvsTheta.Fill(p.Theta(),p.Mag())
				# print ''
				if dist0 < 30:
					vtxP.Fill( p_vtx.Mag() )
					vtxPt.Fill( p_vtx.Pt() )
					vtxTheta.Fill( p_vtx.Theta() )
					vtxPhi.Fill( p_vtx.Phi() )

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

				vtxTracksPtSum.Fill(leadingP.Pt()+secondP.Pt())

				vtx_refpoint1 = utils.spacialDistance( trkSt1.getReferencePoint(),vtx_cand.getPosition() )
				vtx_refpoint2 = utils.spacialDistance( trkSt2.getReferencePoint(),vtx_cand.getPosition() )
				vtxRefPointDist1VsDist2.Fill(vtx_refpoint1,vtx_refpoint2)

				# if vtx_refpoint1 > 200 and vtx_refpoint2 > 200:
				# 	print i, vtx_cand.getPosition()[0],vtx_cand.getPosition()[1], vtx_cand.getPosition()[2]

				vtxTrksRefPointZ.Fill( abs(trkSt1.getReferencePoint()[2]-trkSt2.getReferencePoint()[2]) )
				distRef12 = utils.spacialDistance( trkSt1.getReferencePoint(), trkSt2.getReferencePoint() )
				vtxTrksRefPointsDist.Fill( distRef12 )
				# print i, vtx_cand.getPosition()[0], vtx_cand.getPosition()[1], vtx_cand.getPosition()[2], vtx_cand.getParameters()[0], distRef12
				# print i, trkSt1.getReferencePoint()[0], trkSt1.getReferencePoint()[1], trkSt1.getReferencePoint()[2],trksVec[iTrk].getNdf()
				# print i, trkSt2.getReferencePoint()[0], trkSt2.getReferencePoint()[1], trkSt2.getReferencePoint()[2],trksVec[iTrk+1].getNdf()
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
					print arc1, arc2
	
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
				# 	continue
				# if ( (refpoint_lhZ1 < 100 and cutPhi1) \
				# or (refpoint_lhZ2 < 100 and cutPhi2) ):
				# 	iTrk += 2
				# 	continue

				vtxRefZVsLhZ.Fill(sign_pz1*vtx_refpointZ1/refpoint_lhZ1, sign_pz2*vtx_refpointZ2/refpoint_lhZ2)
				vtxRefZVsLenZ.Fill(sign_pz1*vtx_refpointZ1/refpoint_lhZ1, refpoint_lhZ1)
				vtxRefZVsLenZ.Fill(sign_pz2*vtx_refpointZ2/refpoint_lhZ2, refpoint_lhZ2)

				if q1 < 0:
					if refpoint_lhZ1 < 100000 and refpoint_lhZ2 < 100000:
						vtxRefPhiVsLhPhi.Fill(q1*phiRef1/arc1, q2*(2*math.pi-phiRef2)/arc2)
						vtxRefPhiVsArcPhi.Fill(q1*phiRef1/arc1, arc1)
					if refpoint_lhZ2 < 100000:
						# vtxRefPhiVsLhPhi.Fill(q2*(2*math.pi-phiRef2)/arc2, q2*(2*math.pi-phiLH2)/arc2)
						vtxRefPhiVsArcPhi.Fill(q2*(2*math.pi-phiRef2)/arc2,arc2)
				else:
					if refpoint_lhZ1 < 100000 and refpoint_lhZ2 < 100000:
						vtxRefPhiVsLhPhi.Fill(q1*(2*math.pi-phiRef1)/arc1, q2*phiRef2/arc2)
						vtxRefPhiVsArcPhi.Fill(q1*(2*math.pi-phiRef1)/arc1, arc1)
					if refpoint_lhZ2 < 100000:
						# vtxRefPhiVsLhPhi.Fill(q2*phiRef2/arc2, q2*phiLH2/arc2)
						vtxRefPhiVsArcPhi.Fill(q2*phiRef2/arc2, arc2)
				# if not (leadingP.Pt() > 1.5 and leadingTrack.getNdf() < 70):
				if True:
					vtxTrksCentresDist.Fill(centresDist)
					cosOpenAngleVtx_reco.Fill(cosOpenAngle)
					vtxPt1vsPt2.Fill(leadingP.Pt(),secondP.Pt())
					vtxP1vsP2.Fill(leadingP.Mag(),secondP.Mag())
					vtxPz1vsPz2.Fill(leadingP.Pz(),secondP.Pz())
					vtxOmega1vsOmega2.Fill(abs(omega1),abs(omega2))
					vtxTrksNdf1VsNdf2.Fill(leadingTrack.getNdf(),secondTrack.getNdf())
					vtxTrksPtVsNdf.Fill(leadingP.Pt(),leadingTrack.getNdf())
					vtxTrksPtVsNdf.Fill(secondP.Pt(),secondTrack.getNdf())
					vtxTrksD01VsD02.Fill(leadingTrack.getD0(),secondTrack.getD0())
					vtxTrksZ01VsZ02.Fill(leadingTrack.getZ0(),secondTrack.getZ0())
					# if 550 * math.log10(distRef12) + 550 < centresDist:
					# if 2.2 * distRef12 + 900 < centresDist:
					vtxTrksCentresRefPointsDist.Fill(distRef12, centresDist)
					vtxTrksCentresRefPointsDist_scaled.Fill(distRef12, centresDistScaled)
					vtxRefPointDistsVsCentres.Fill(vtx_refpoint1,centresDist)
					vtxRefPointDistsVsCentres.Fill(vtx_refpoint2,centresDist)
					# if R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 330 and R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) < 1800:
					# if p_vtx.Pt() > 1.9:
						# if centresRefsDists < -2000:
					corrCentresRefDistsVsPt.Fill(centresRefsDists, p_vtx.Pt())
					# if logCentresRefsDists < -2000:
					# 	if p_vtx.Pt() > 1.9:
					corrLogCentresRefDistsVsPt.Fill(logCentresRefsDists, p_vtx.Pt())
					if R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]) > 350:
						corrCentresRefDistsVsR.Fill(centresRefsDists, R(vtx_cand.getPosition()[0],vtx_cand.getPosition()[1]))
					vtxTrksOmega.Fill(abs(omega1))
					vtxTrksOmega.Fill(abs(omega2))
					vtxTrksCurvRatio.Fill(ratio)
					cosOpenAngleVsCurvRatio.Fill(cosOpenAngle, ratio)
					vtxTrksChi2.Fill( trksVec[iTrk].getChi2() )
					vtxTrksChi2.Fill( trksVec[iTrk+1].getChi2() )
					vtxTrksNdf.Fill( trksVec[iTrk].getNdf() )
					vtxTrksNdf.Fill( trksVec[iTrk+1].getNdf() )
					if trksVec[iTrk].getNdf() > 0:
						vtxTrksChi2Ndf.Fill( trksVec[iTrk].getChi2() / trksVec[iTrk].getNdf() )
					if trksVec[iTrk+1].getNdf() > 0:
						vtxTrksChi2Ndf.Fill( trksVec[iTrk+1].getChi2() / trksVec[iTrk+1].getNdf() )
				vtxTrksD0AtIP.Fill( trksVec[iTrk].getD0() )
				vtxTrksD0AtIP.Fill( trksVec[iTrk+1].getD0() )
				vtxTrksZ0AtIP.Fill( trksVec[iTrk].getZ0() )
				vtxTrksZ0AtIP.Fill( trksVec[iTrk+1].getZ0() )
				vtxTrksD0.Fill( trkSt1.getD0() )
				vtxTrksD0.Fill( trkSt2.getD0() )
				vtxTrksZ0.Fill( trkSt1.getZ0() )
				vtxTrksZ0.Fill( trkSt2.getZ0() )

			iTrk+=2


		for particle in pfoCollection:
			pid = particle.getParticleIDs()[0].getPDG()
			if abs(pid) == pdg1 or abs(pid) == pdg2:
				nCorrPID+=1
		# 	print i,
		# 	for id in particle.getParticleIDs():
		# 		print id.getPDG(),
		# 	print ''
		# 	for id in particle.getParticleIDs():
		# 		print id.getAlgorithmType(),
		# 	print ''
		# 	print i, particle.getParticleIDUsed()
		# 	print i, particle.getGoodnessOfPID()
		# print ''


reader.close()

centresRefDistsVsPtEv = corrCentresRefDistsVsPt.Integral()
print 'linear int: ', centresRefDistsVsPtEv
print "log int: ", corrLogCentresRefDistsVsPt.Integral()

#print('Decays within tracker acceptance', 'Ev. with min. two tracks', 'Ev. with min. two reco. tracks', 'All reco. tracks', 'Tracks matching to MC', 'Two reco. leptons')
#print(withinAcceptance, twoTruthTracks, twoTracks, allTracks, allMatchingTracks, twoRecoLeptons)

print ''
print 'Tracks within acceptance: ' + str(withinAcceptance)
print 'Ev. with two truth tracks: ' + str(twoTruthTracks)
print 'Ev. with min. two reco. tracks: ' + str(twoTracks)
print 'All reco. tracks: ' + str(allTracks)
print 'Number of PFOs: ' + str(nPFOs)
print 'Number of PFOs with correct PID: ' + str(nCorrPID)
print 'PFO ID eff.: ', float(nCorrPID) / float(withinAcceptance)
print 'Tracks matching to MC (MarlinTrk): ' + str(allMatchingTracks)
print 'Tracks matching to MC (Refit): ' + str(allMatchingTracksRefitted)
print 'Tracks matching only with MarlinTrk: ' + str(onlyMarlinTrkTracks)
print 'Tracks matching only with Refit: ' + str(onlyRefittedTracks)
if withinAcceptance > 0:
	print 'MarlinTrk reco. eff.: ', float(allMatchingTracks) / float(withinAcceptance)
	print 'Refitted reco. eff.:  ', float(allMatchingTracksRefitted) / float(withinAcceptance)
print 'Number of sec. vtx:  ', nSecVertices
print 'Ev. with >=1 sec. vtx:  ', nEvSecVertices
print 'Matching vertices:  ', allMatchingVertices
print 'Fake vertices:  ', nFakeVertices
if withinAcceptance > 0:
	print 'Vtx finding eff.:  ', float(allMatchingVertices) / (float(withinAcceptance)/2)
	print 'Vtx finding purity:  ', float(allMatchingVertices) / float(nSecVertices)
	print '* Eff. (inside TPC):  ', float(matchingVerticesTPC) / (float(withinTPC)/2)
	print '* Purity (inside TPC):  ', float(matchingVerticesTPC) / float(nSecVerticesTPC)
	print 'Eff. inside TPC after pT, Ndf cuts: ', float(centresRefDistsVsPtEv) / (float(withinTPC)/2)
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

vtxEffVsL = ROOT.TH1F('vtxEffVsL', 'Vertex reconstruction efficiency', 30, 0., 3000.)
vtxEffVsL.SetXTitle('#beta c t_{LLP} [mm]')
vtxEffVsL.SetYTitle('Vertex reco. efficiency')
vtxEffVsL.Divide(vtxEffVsL_reco, vtxEffVsL_true,1,1,"cl=0.683 b(1,1) mode")
vtxEffVsL.SetMarkerStyle(20)
vtxEffVsL.SetMarkerSize(0.5)
vtxEffVsL.Draw('e')
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxEffVsL_'+fs+'.pdf')

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
print '1'
DeltaThetaTrackLLP.SetTitle('Track reconstruction efficiency;|#theta_{LLP} - #theta_{track}|;Efficiency')
print '2'
DeltaThetaTrackLLP.SetMarkerStyle(20)
print '3'
DeltaThetaTrackLLP.SetMarkerSize(0.5)
print '4'
#DeltaThetaTrackLLP.Draw('ap')
print '5'
#c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaThetaTrackLLP_'+fs+'.pdf')

print 'DeltaPhi: ', DeltaPhiTrackLLP_reco.GetEntries(), DeltaPhiTrackLLP_true.GetEntries()
DeltaPhiTrackLLP = ROOT.TEfficiency(DeltaPhiTrackLLP_reco,DeltaPhiTrackLLP_true)
DeltaPhiTrackLLP.SetTitle('Track reconstruction efficiency;|#phi_{LLP} - #phi_{track}|;Efficiency')
DeltaPhiTrackLLP.SetMarkerStyle(20)
DeltaPhiTrackLLP.SetMarkerSize(0.5)
#DeltaPhiTrackLLP.Draw('ap')
#c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaPhiTrackLLP_'+fs+'.pdf')

print 'DeltaAlpha: ', DeltaAlphaTrackLLP_reco.GetEntries(), DeltaAlphaTrackLLP_true.GetEntries()
DeltaAlphaTrackLLP = ROOT.TEfficiency(DeltaAlphaTrackLLP_reco,DeltaAlphaTrackLLP_true)
DeltaAlphaTrackLLP.SetTitle('Track reconstruction efficiency;|#Delta#alpha_{LLP-track}|;Efficiency')
DeltaAlphaTrackLLP.SetMarkerStyle(20)
DeltaAlphaTrackLLP.SetMarkerSize(0.5)
#DeltaAlphaTrackLLP.Draw('ap')
#c0.SaveAs('plots/'+cat+'/'+fs+'/DeltaAlphaTrackLLP_'+fs+'.pdf')
'''
trackEffVsR_true.Scale( 1./trackEffVsR_true.Integral() )
trackEffVsR_true.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/R_true_'+fs+'.pdf')
trackEffVsR_reco.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/R_reco_'+fs+'.pdf')
momTrack_true.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/momTrack_true_'+fs+'.pdf')
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

c0.SetLogy()

#thetaTrack_true.Scale( 1./thetaTrack_true.Integral() )
thetaTrack_true.Draw('hist')
c0.SaveAs('plots/'+cat+'/'+fs+'/thetaTrack_true_'+fs+'.pdf')
cosOpenAngleVtx_true.Draw('hist')
cosOpenAngleVtx_true.SetXTitle('cos(#Delta#alpha)')
c0.SaveAs('plots/'+cat+'/'+fs+'/cosOpenAngleVtx_true_'+fs+'.pdf')
cosOpenAngleVtx_reco.Draw('hist')
cosOpenAngleVtx_reco.SetXTitle('cos(#Delta#alpha)')
c0.SaveAs('plots/'+cat+'/'+fs+'/cosOpenAngleVtx_reco_'+fs+'.pdf')
vtxMisdistances.Draw('hist')
vtxMisdistances.SetXTitle('#Delta R_{vtx} [mm]')
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxMisdistance_'+fs+'.pdf')
vtxMultiplicity.Draw('hist')
vtxMultiplicity.SetXTitle('N_{vtx}')
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxMultiplicity_'+fs+'.pdf')
hitsDistances.Draw('hist')
hitsDistances.SetXTitle('l_{ij} [mm]')
c0.SaveAs('plots/'+cat+'/'+fs+'/hitsDistances_'+fs+'.pdf')
helixDistances.Draw('hist')
helixDistances.SetXTitle('helix distance [mm]')
c0.SaveAs('plots/'+cat+'/'+fs+'/helixDistances_'+fs+'.pdf')
vtxEffVsL_true.SetMinimum(0.1)
vtxEffVsL_true.Draw('hist')
vtxEffVsL_true.SetXTitle('#beta c t_{LLP} [mm]')
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxEffVsL_true_'+fs+'.pdf')
vtxEffVsL_reco.Draw('hist')
vtxEffVsL_reco.SetXTitle('#beta c t_{LLP} [mm]')
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxEffVsL_reco_'+fs+'.pdf')
vtxVsR_reco.Draw('hist')
vtxVsR_reco.SetXTitle('R [mm]')
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxVsR_reco_'+fs+'.pdf')
vtxVsL_reco.Draw('hist')
vtxVsL_reco.SetXTitle('L [mm]')
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxVsL_reco_'+fs+'.pdf')
vtxP.Draw('hist')
vtxP.SetXTitle("p [GeV]")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxP_'+fs+'.pdf')
vtxPt.Draw('hist')
vtxPt.SetXTitle("p_{T} [GeV]")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxPt_'+fs+'.pdf')
vtxTheta.Draw('hist')
vtxTheta.SetXTitle("Polar angle")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTheta_'+fs+'.pdf')
vtxPhi.Draw('hist')
vtxPhi.SetXTitle("Azimuthal angle")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxPhi_'+fs+'.pdf')
vtxTracksPt.Draw('hist')
vtxTracksPt.SetXTitle("p_{T} [GeV]")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTracksPt_'+fs+'.pdf')
vtxTracksPtSum.Draw('hist')
vtxTracksPtSum.SetXTitle("p_{T}^{(1)}+p_{T}^{(2)} [GeV]")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTracksPtSum_'+fs+'.pdf')
vtxTracksP.Draw('hist')
vtxTracksP.SetXTitle("p [GeV]")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTracksP_'+fs+'.pdf')
vtxTracksTheta.Draw('hist')
vtxTracksTheta.SetXTitle("Polar angle")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTracksTheta_'+fs+'.pdf')
vtxTracksPhi.Draw('hist')
vtxTracksPhi.SetXTitle("Azimuthal angle")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTracksPhi_'+fs+'.pdf')
vtxRefPointDist.Draw('hist')
vtxRefPointDist.SetXTitle("#DeltaR(vtx-trackState) [mm]")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxRefPointDist_'+fs+'.pdf')
vtxTrksRefPointsDist.Draw('hist')
vtxTrksRefPointsDist.SetXTitle("#DeltaR(vtx-trackState) [mm]")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksRefPointsDist_'+fs+'.pdf')
vtxTrksCurvRatio.Draw('hist')
vtxTrksCurvRatio.SetXTitle("|#Omega_{1}|")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksCurvRatio'+fs+'.pdf')
vtxTrksOmega.Draw('hist')
vtxTrksOmega.SetXTitle("|#Omega|")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksOmega'+fs+'.pdf')
c0.SetLogx()
vtxTrksRefPointZ.Draw('hist')
vtxTrksRefPointZ.SetXTitle("#Deltaz_{trk} [mm]")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksRefPointZ'+fs+'.pdf')
c0.SetLogx(0)
vtxTrksChi2.Draw('hist')
vtxTrksChi2.SetXTitle("#chi^{2}")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksChi2_'+fs+'.pdf')
vtxTrksChi2Ndf.Draw('hist')
vtxTrksChi2Ndf.SetXTitle("#chi^{2} / N_{df}")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksChi2Ndf_'+fs+'.pdf')
vtxTrksNdf.Draw('hist')
vtxTrksNdf.SetXTitle("N_{df}")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksNdf_'+fs+'.pdf')
vtxTrksD0.Draw('hist')
vtxTrksD0.SetXTitle("d_{0}")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksD0_'+fs+'.pdf')
vtxTrksZ0.Draw('hist')
vtxTrksZ0.SetXTitle("z_{0}")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksZ0_'+fs+'.pdf')
vtxTrksD0AtIP.Draw('hist')
vtxTrksD0AtIP.SetXTitle("d_{0}")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksD0AtIP_'+fs+'.pdf')
vtxTrksZ0AtIP.Draw('hist')
vtxTrksZ0AtIP.SetXTitle("z_{0}")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksZ0AtIP_'+fs+'.pdf')
vtxTrksTruePDG.Draw('hist')
vtxTrksTruePDG.SetXTitle("|PDG|")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksTruePDG_'+fs+'.pdf')
vtxTrksAreBothCurlers.Draw('hist')
vtxTrksAreBothCurlers.SetXTitle("|PDG|")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksAreBothCurlers_'+fs+'.pdf')
vtxTrksCentresDist.Draw('hist')
vtxTrksCentresDist.SetXTitle("d_{C}")
c0.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksCentresDist_'+fs+'.pdf')

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
trackEffiPvsTheta.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/trackEffiPvsTheta_'+fs+'.pdf')

#vtxEffRVsZ = ROOT.TEfficiency(vtxEffRVsZ_reco,vtxEffRVsZ_true)
vtxEffRVsZ = vtxEffRVsZ_reco.Clone()
vtxEffRVsZ.Divide(vtxEffRVsZ_true)
vtxEffRVsZ.SetTitle('Vertex reconstruction efficiency;z [mm];R [mm];Efficiency')
vtxEffRVsZ.SetMinimum(0.)
vtxEffRVsZ.SetMaximum(1.)
vtxEffRVsZ.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
'''
ROOT.gPad.Update()
graph = vtxEffRVsZ.GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1)
ROOT.gPad.Update()
'''
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxEffRVsZ_'+fs+'.pdf')

vtxEffRVsZ_true.SetYTitle('R [mm]')
vtxEffRVsZ_true.SetXTitle('z [mm]')
vtxEffRVsZ_true.SetZTitle('Number of vertices')
vtxEffRVsZ_true.SetStats(0)
vtxEffRVsZ_true.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxEffRVsZ_true_'+fs+'.pdf')

vtxPt1vsPt2.SetYTitle('p_{T}^{2} [GeV]')
vtxPt1vsPt2.SetXTitle('p_{T}^{1} [GeV]')
vtxPt1vsPt2.SetZTitle('Number of vertices')
vtxPt1vsPt2.SetStats(0)
vtxPt1vsPt2.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxPt1vsPt2_'+fs+'.pdf')

vtxP1vsP2.SetYTitle('p^{(2)} [GeV]')
vtxP1vsP2.SetXTitle('p^{(1)} [GeV]')
vtxP1vsP2.SetZTitle('Number of vertices')
vtxP1vsP2.SetStats(0)
vtxP1vsP2.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxP1vsP2_'+fs+'.pdf')

vtxPz1vsPz2.SetYTitle('p_{z}^{(2)} [GeV]')
vtxPz1vsPz2.SetXTitle('p_{z}^{(1)} [GeV]')
vtxPz1vsPz2.SetZTitle('Number of vertices')
vtxPz1vsPz2.SetStats(0)
vtxPz1vsPz2.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxPz1vsPz2_'+fs+'.pdf')

vtxOmega1vsOmega2.SetYTitle('#Omega_{2}')
vtxOmega1vsOmega2.SetXTitle('#Omega_{1}')
vtxOmega1vsOmega2.SetZTitle('Number of vertices')
vtxOmega1vsOmega2.SetStats(0)
vtxOmega1vsOmega2.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxOmega1vsOmega2_'+fs+'.pdf')

vtxRefPointDist1VsDist2.SetYTitle('d_{2}')
vtxRefPointDist1VsDist2.SetXTitle('d_{1}')
vtxRefPointDist1VsDist2.SetZTitle('Number of vertices')
vtxRefPointDist1VsDist2.SetStats(0)
vtxRefPointDist1VsDist2.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxRefPointDist1VsDist2_'+fs+'.pdf')

vtxTrksNdf1VsNdf2.SetYTitle('N_{df}^{(2)}')
vtxTrksNdf1VsNdf2.SetXTitle('N_{df}^{(1)}')
vtxTrksNdf1VsNdf2.SetZTitle('Number of vertices')
vtxTrksNdf1VsNdf2.SetStats(0)
vtxTrksNdf1VsNdf2.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksNdf1VsNdf2_'+fs+'.pdf')

vtxTrksPtVsNdf.SetYTitle('N_{df}')
vtxTrksPtVsNdf.SetXTitle('p_{T} [GeV]')
vtxTrksPtVsNdf.SetZTitle('Number of tracks')
vtxTrksPtVsNdf.SetStats(0)
vtxTrksPtVsNdf.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksPtVsNdf_'+fs+'.pdf')

vtxTrksD01VsD02.SetYTitle('d_{0}^{(2)}')
vtxTrksD01VsD02.SetXTitle('d_{0}^{(1)}')
vtxTrksD01VsD02.SetZTitle('Number of vertices')
vtxTrksD01VsD02.SetStats(0)
vtxTrksD01VsD02.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksD01VsD02_'+fs+'.pdf')

vtxTrksZ01VsZ02.SetYTitle('z_{0}^{(2)}')
vtxTrksZ01VsZ02.SetXTitle('z_{0}^{(1)}')
vtxTrksZ01VsZ02.SetZTitle('Number of vertices')
vtxTrksZ01VsZ02.SetStats(0)
vtxTrksZ01VsZ02.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksZ01VsZ02_'+fs+'.pdf')

c1.SetLogx()
vtxTrksCentresRefPointsDist_scaled.SetYTitle('d_{C}/#sqrt{p_{T}^{1} p_{T}^{2}}')
vtxTrksCentresRefPointsDist_scaled.SetXTitle('d_{ref}')
vtxTrksCentresRefPointsDist_scaled.SetZTitle('Number of vertices')
vtxTrksCentresRefPointsDist_scaled.SetStats(0)
vtxTrksCentresRefPointsDist_scaled.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksCentresRefPointsDist_scaled_'+fs+'.pdf')

vtxTrksCentresRefPointsDist.SetYTitle('d_{C}')
vtxTrksCentresRefPointsDist.SetXTitle('d_{ref}')
vtxTrksCentresRefPointsDist.SetZTitle('Number of vertices')
vtxTrksCentresRefPointsDist.SetStats(0)
vtxTrksCentresRefPointsDist.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxTrksCentresRefPointsDist_'+fs+'.pdf')

vtxRefPointDistsVsCentres.SetYTitle('d_{C}')
vtxRefPointDistsVsCentres.SetXTitle('d_{ref-vtx}')
vtxRefPointDistsVsCentres.SetZTitle('Number of tracks')
vtxRefPointDistsVsCentres.SetStats(0)
vtxRefPointDistsVsCentres.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxRefPointDistsVsCentres_'+fs+'.pdf')
c1.SetLogx(0)
# print vtxTrksCentresRefPointsDist.Integral()

corrCentresRefDistsVsPt.SetYTitle('p_{T}^{vtx}')
corrCentresRefDistsVsPt.SetXTitle('2.2d_{ref} - d_{C}')
corrCentresRefDistsVsPt.SetZTitle('Number of vertices')
corrCentresRefDistsVsPt.SetStats(0)
corrCentresRefDistsVsPt.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/corrCentresRefDistsVsPt_'+fs+'.pdf')
# print 'pT vs. dd integral:', centresRefDistsVsPtEv

corrLogCentresRefDistsVsPt.SetYTitle('p_{T}^{vtx}')
corrLogCentresRefDistsVsPt.SetXTitle('300log_{10}(d_{ref}) - d_{C}')
corrLogCentresRefDistsVsPt.SetZTitle('Number of vertices')
corrLogCentresRefDistsVsPt.SetStats(0)
corrLogCentresRefDistsVsPt.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/corrLogCentresRefDistsVsPt_'+fs+'.pdf')

corrCentresRefDistsVsR.SetYTitle('R_{vtx}')
corrCentresRefDistsVsR.SetXTitle('2.2d_{ref} - d_{C}')
corrCentresRefDistsVsR.SetZTitle('Number of vertices')
corrCentresRefDistsVsR.SetStats(0)
corrCentresRefDistsVsR.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/corrCentresRefDistsVsR_'+fs+'.pdf')
# print 'R vs. dd integral:', corrCentresRefDistsVsR.Integral()

cosOpenAngleVsCurvRatio.SetXTitle('cos(#Delta#alpha)')
cosOpenAngleVsCurvRatio.SetYTitle('|#Omega_{1}/#Omega_{2}|')
cosOpenAngleVsCurvRatio.SetZTitle('Number of vertices')
cosOpenAngleVsCurvRatio.SetStats(0)
cosOpenAngleVsCurvRatio.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/cosOpenAngleVsCurvRatio_'+fs+'.pdf')

vtxRefPhiVsLhPhi.SetXTitle('q_{1} #phi^{(1)}_{ref}/#phi^{(1)}_{arc}')
vtxRefPhiVsLhPhi.SetYTitle('q_{2} #phi^{(2)}_{ref}/#phi^{(2)}_{arc}')
vtxRefPhiVsLhPhi.SetZTitle('Number of tracks')
vtxRefPhiVsLhPhi.SetStats(0)
vtxRefPhiVsLhPhi.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxRefPhiVsLhPhi_'+fs+'.pdf')

vtxRefPhiVsArcPhi.SetXTitle('#phi_{ref}/#phi_{arc}')
vtxRefPhiVsArcPhi.SetYTitle('#phi_{arc}')
vtxRefPhiVsArcPhi.SetZTitle('Number of tracks')
vtxRefPhiVsArcPhi.SetStats(0)
vtxRefPhiVsArcPhi.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxRefPhiVsArcPhi_'+fs+'.pdf')

vtxRefZVsLhZ.SetXTitle('sgn(p_{z})(z_{ref}-z_{vtx})/|z_{ref}-z_{lh}|')
vtxRefZVsLhZ.SetYTitle('sgn(p_{z})(z_{lh}-z_{vtx})/|z_{ref}-z_{lh}|')
vtxRefZVsLhZ.SetZTitle('Number of tracks')
vtxRefZVsLhZ.SetStats(0)
vtxRefZVsLhZ.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxRefZVsLhZ_'+fs+'.pdf')

c1.SetLogy()
vtxRefZVsLenZ.SetXTitle('sgn(p_{z})(z_{ref}-z_{vtx})/|z_{ref}-z_{lh}|')
vtxRefZVsLenZ.SetYTitle('|z_{ref}-z_{lh}|')
vtxRefZVsLenZ.SetZTitle('Number of tracks')
vtxRefZVsLenZ.SetStats(0)
vtxRefZVsLenZ.Draw('colz')
c1.SaveAs('plots/'+cat+'/'+fs+'/vtxRefZVsLenZ_'+fs+'.pdf')
c1.SetLogy(0)

tracksVtxPvsTheta.SetYTitle('p [GeV]')
tracksVtxPvsTheta.SetXTitle('#theta_{track}')
tracksVtxPvsTheta.SetZTitle('Number of tracks')
tracksVtxPvsTheta.SetStats(0)
#tracksVtxPvsTheta.Scale( 1./tracksVtxPvsTheta.Integral() )
tracksVtxPvsTheta.Draw('colz text')
c1.SaveAs('plots/'+cat+'/'+fs+'/tracksVtxPvsTheta_'+fs+'.pdf')

tracksVtxPtVsTheta.SetYTitle('p [GeV]')
tracksVtxPtVsTheta.SetXTitle('#theta_{track}')
tracksVtxPtVsTheta.SetZTitle('Number of tracks')
tracksVtxPtVsTheta.SetStats(0)
tracksVtxPtVsTheta.Draw('colz text')
c1.SaveAs('plots/'+cat+'/'+fs+'/tracksVtxPtVsTheta_'+fs+'.pdf')


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

fakeVertices.SetYTitle('R [mm]')
fakeVertices.SetXTitle('z [mm]')
fakeVertices.SetZTitle('Number of vertices')
fakeVertices.SetStats(0)
fakeVertices.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/fakeVerticesRVsZ_'+fs+'.pdf')

trackEffiPvsTheta_true.SetYTitle('p [GeV]')
trackEffiPvsTheta_true.SetXTitle('#theta_{track}')
trackEffiPvsTheta_true.SetZTitle('Number of tracks')
trackEffiPvsTheta_true.SetStats(0)
#trackEffiPvsTheta_true.Scale( 1./trackEffiPvsTheta_true.Integral() )
trackEffiPvsTheta_true.Draw('colz')
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
trackEffiPvsTheta_reco.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/'+cat+'/'+fs+'/trackEffiPvsTheta_reco_'+fs+'.pdf')

trackEffiPvsR_true.SetYTitle('R [mm]')
trackEffiPvsR_true.SetXTitle('p [GeV]')
trackEffiPvsR_true.SetZTitle('Number of tracks')
trackEffiPvsR_true.SetStats(0)
trackEffiPvsR_true.Draw('colz')
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
histo_array.append(vtxEffVsL)
#histo_array.append(DeltaAlphaTrackLLP)
#histo_array.append(DeltaThetaTrackLLP)
#histo_array.append(DeltaPhiTrackLLP)
histo_array.append(vtxMisdistances)
histo_array.append(vtxVsL_reco)
histo_array.append(vtxVsR_reco)
histo_array.append(recoVertices)
histo_array.append(vtxP)
histo_array.append(vtxPt)
histo_array.append(vtxTheta)
histo_array.append(vtxPhi)
histo_array.append(vtxTracksPt)
histo_array.append(vtxTracksP)
histo_array.append(vtxTracksTheta)
histo_array.append(vtxTracksPhi)
histo_array.append(vtxRefPointDist)
histo_array.append(vtxTrksCentresRefPointsDist)
histo_array.append(vtxTrksRefPointsDist)
histo_array.append(vtxTrksCurvRatio)
histo_array.append(cosOpenAngleVtx_true)
histo_array.append(cosOpenAngleVtx_reco)
histo_array.append(vtxTrksChi2)
histo_array.append(vtxTrksChi2Ndf)
histo_array.append(vtxTrksNdf)
histo_array.append(vtxTrksD0AtIP)
histo_array.append(vtxTrksZ0AtIP)

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
