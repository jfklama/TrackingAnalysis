import sys
import os

import math
import numpy as np

from pyLCIO import IOIMPL, ROOT

import utils

def R(x,y):
	return math.sqrt(x**2 + y**2)
def R3(x,y,z):
	return math.sqrt(x**2 + y**2 + z**2)

try:
    arg1 = sys.argv[1]
except IndexError:
    print "Usage: " + os.path.basename(__file__) + " <inputFile>"
    sys.exit(1)



chargedChannel    = 0
withinAcceptance  = 0
twoTracks         = 0
twoTruthTracks    = 0
allTracks		  = 0
allMatchingTracks = 0
twoRecoLeptons    = 0

nbins = 12
TPCRmax = 1974.
#TPCRmax = 1770.
TPCZmax = 2350.

trackEffVsR_true = ROOT.TH1F('trackEffVsR_true', 'True vertices, p_{track} < 2 GeV', nbins, 0., TPCRmax)
trackEffVsR_reco = ROOT.TH1F('trackEffVsR_reco', 'N_{evt} with 2 reco. tracks', nbins, 0., TPCRmax)

trackEffRVsZ_true = ROOT.TH2F('trackEffRVsZ_true', 'True vertices', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffRVsZ_reco = ROOT.TH2F('trackEffRVsZ_reco', 'N_{evt} with 2 reco. tracks', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffiPvsTheta_true = ROOT.TH2F('trackEffiPvsTheta_true', '', 3, 0., 3.15, 15, 0., 15.)
trackEffiPvsTheta_reco = ROOT.TH2F('trackEffiPvsTheta_reco', '', 3, 0., 3.15, 15, 0., 15.)
trackEffiPvsR_true = ROOT.TH2F('trackEffiPvsR_true', '', 15, 0., 15., nbins, 0., TPCRmax)


momTrack_true = ROOT.TH1F('momTrack_true', 'True momentum of track', 50, 0., 15.)
momTrack_reco = ROOT.TH1F('momTrack_reco', 'True momentum of reco. track', 50, 0., 15.)
ptTrack_true = ROOT.TH1F('ptTrack_true', 'True pT of track', 50, 0., 15.)
ptTrack_reco = ROOT.TH1F('ptTrack_reco', 'True pT of reco. track', 50, 0., 15.)
thetaTrack_true = ROOT.TH1F('thetaTrack_true', 'True polar angle of track', 32, 0., 3.15)
thetaTrack_reco = ROOT.TH1F('thetaTrack_reco', 'True polar angle of reco. track', 32, 0., 3.15)

DeltaThetaTrackLLP_true = ROOT.TH1F('DeltaThetaTrackLLP_true', 'Dist. between polar angles of track and LLP', 35, 0., 3.5)
DeltaThetaTrackLLP_reco = ROOT.TH1F('DeltaThetaTrackLLP_reco', 'Dist. between polar angles of track and LLP', 35, 0., 3.5)
DeltaPhiTrackLLP_true = ROOT.TH1F('DeltaPhiTrackLLP_true', 'Dist. between azimuthal angles of track and LLP', 35, 0., 3.5)
DeltaPhiTrackLLP_reco = ROOT.TH1F('DeltaPhiTrackLLP_reco', 'Dist. between azimuthal angles of track and LLP', 35, 0., 3.5)

momLLP_true = ROOT.TH1F('momLLP_true', 'True momentum of LLP', 50, 0., 300.)
momLLP_reco = ROOT.TH1F('momLLP_reco', 'True momentum of reco. LLP', 50, 0., 300.)
thetaLLP_true = ROOT.TH1F('thetaLLP_true', 'True polar angle of LLP', 32, 0., 3.15)
thetaLLP_reco = ROOT.TH1F('thetaLLP_reco', 'True polar angle of reco. LLP', 32, 0., 3.15)

h_nTracks = ROOT.TH1F('h_nTracks', 'Number of tracks', 20, 0., 20.)
h_relWeights = ROOT.TH1F('h_relWeights', 'Number of weights', 50, 0., 2.)
h_angDist = ROOT.TH1F('h_angDist', 'Angular separation', 70, 0., 7)

recoVertices = ROOT.TH2F('recoVertices', 'N_{evt} with 2 reco. leptons', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)

hmass = ROOT.TH1F('hmass', '#Lambda mass', 150, 0., 15.)

infile = sys.argv[1]
reader=IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(infile)

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

print "LLP: " + str(llp_pdg) + ", with final state: " + str(pdg1) + ", " + str(pdg2)

i=0
for event in reader:
	i+=1

	mcCollection = event.getCollection('MCParticlesSkimmed')
	#tpcCollection = event.getCollection('TPCCollection')
	trackCollection = event.getCollection('MarlinTrkTracks')
	trackToMCLinkCollection = event.getCollection('MarlinTrkTracksMCTruthLink')
	refitTrackCollection = event.getCollection('RefittedTracks')
	refitTrackToMCLinkCollection = event.getCollection('RefittedTracksMCP')
	pfoCollection = event.getCollection('PandoraPFOs')
	#vertexCollection = event.getCollection('PandoraPFANewStartVertices')

	if i%100==0:
		print (i, str(trackCollection.getNumberOfElements() ) + ' tracks in event')

	#check if LLP decayed within tracking acceptance and we have 2 truth tracks
	trueR = 0
	trueZ = 0
	nFSParticle1 = 0
	nFSParticle2 = 0
	for particle in mcCollection:
		if abs(particle.getPDG()) == llp_pdg:
			trueR = R(particle.getEndpoint()[0], particle.getEndpoint()[1])
			trueZ = particle.getEndpoint()[2]
		if(particle.getGeneratorStatus() == gen_status):
			if abs( particle.getPDG() ) == pdg1:
				nFSParticle1+=1
			elif abs( particle.getPDG() ) == pdg2:
				nFSParticle2+=1
	if (trueR >= 329 and trueR <= 1770) and abs(trueZ) <= 2350:
		withinAcceptance+=1
	else:
		continue

	if nFSParticle1<1 and nFSParticle2<1:
		#print (i, 'WARNING: no truth tracks in the event!')
		continue
	if nFSParticle1+nFSParticle2 != 2:
		#print (i, nFSParticle1+nFSParticle2, " not 2 truth tracks in the event!")
		continue
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

			if abs( particle.getPDG() ) == pdg1:
				nFSParticle1+=1
				#if trueR < 450 and abs(trueZ) < 400 and (mclepP4.Theta() > 1 and mclepP4.Theta() < 2.14):
				#if mclepP4.Theta() > 1 and mclepP4.Theta() < 2.14:
				trackEffVsR_true.Fill(trueR)
				trackEffRVsZ_true.Fill(trueZ,trueR)
				momTrack_true.Fill(p)
				ptTrack_true.Fill(mclepP4.Pt())
				thetaTrack_true.Fill(mclepP4.Theta())
				DeltaThetaTrackLLP_true.Fill( deltaTheta )
				DeltaPhiTrackLLP_true.Fill( abs(deltaPhi) )
				trackEffiPvsTheta_true.Fill(mclepP4.Theta(),p)
				trackEffiPvsR_true.Fill(p,trueR)
				# print i, 'found: ' + str(pdg1)

			elif abs( particle.getPDG() ) == pdg2:
				nFSParticle2+=1
				#if trueR < 450 and abs(trueZ) < 400 and (mclepP4.Theta() > 1 and mclepP4.Theta() < 2.14):
				#if mclepP4.Theta() > 1 and mclepP4.Theta() < 2.14:
				trackEffVsR_true.Fill(trueR)
				trackEffRVsZ_true.Fill(trueZ,trueR)
				momTrack_true.Fill(p)
				ptTrack_true.Fill(mclepP4.Pt())
				#if trueR < 600 and abs(trueZ) < 400:
				thetaTrack_true.Fill(mclepP4.Theta())
				DeltaThetaTrackLLP_true.Fill( deltaTheta )
				DeltaPhiTrackLLP_true.Fill( abs(deltaPhi) )
				trackEffiPvsTheta_true.Fill(mclepP4.Theta(),p)
				trackEffiPvsR_true.Fill(p,trueR)
				# print i, 'found: ' + str(pdg2)

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


	# fill histos if tracks are matching

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
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = R3(px,py,pz)
			e = particle.getEnergy()
			mclepP4.SetPxPyPzE(px,py,pz,e)
			cosTheta = abs( math.cos(mclepP4.Theta()) )
			deltaTheta = abs( mcllpP4.Theta() - mclepP4.Theta() )
			deltaPhi = utils.phiDistance( mcllpP4.Phi(), mclepP4.Phi() )

			for trackRel in trackToMCLinkCollection:

				trackState = trackRel.getFrom().getTrackState(2) # 1 = atIP, 2 = atFirstHit
				momentum = utils.getTrackMomentum(trackState)
				p_tr = ROOT.TVector3()
				p_tr.SetXYZ(momentum[0],momentum[1],momentum[2])
				dist = utils.angularDistance( p_tr.Theta(), mclepP4.Theta(), p_tr.Phi(), mclepP4.Phi()  )
				h_angDist.Fill(dist)

				if trackRel.getTo() == particle \
				and trackRel.getTo().getCharge() * utils.getTrackCharge(trackState) > 0 \
				and dist < 0.2:

					momTrack_reco.Fill(p)
					ptTrack_reco.Fill(mclepP4.Pt())
					momLLP_reco.Fill( p_llp )
					thetaTrack_reco.Fill(mclepP4.Theta())
					DeltaThetaTrackLLP_reco.Fill( deltaTheta )
					DeltaPhiTrackLLP_reco.Fill( abs(deltaPhi) )
					trackEffiPvsTheta_reco.Fill(mclepP4.Theta(),p)
					trackEffVsR_reco.Fill(trueR)
					trackEffRVsZ_reco.Fill(trueZ,trueR)

					allMatchingTracks+=1
					matchingTracksInEvent+=1
					# print i, particle.getPDG(), matchingTracksInEvent
					break


reader.close()

#print('Decays within tracker acceptance', 'Ev. with min. two tracks', 'Ev. with min. two reco. tracks', 'All reco. tracks', 'Tracks matching to MC', 'Two reco. leptons')
#print(withinAcceptance, twoTruthTracks, twoTracks, allTracks, allMatchingTracks, twoRecoLeptons)

print ''
print 'Decays within TPC: ' + str(withinAcceptance)
print 'Ev. with two truth tracks: ' + str(twoTruthTracks)
print 'Ev. with min. two reco. tracks: ' + str(twoTracks)
print 'All reco. tracks: ' + str(allTracks)
print 'Tracks matching to MC: ' + str(allMatchingTracks)
print 'Two reco. leptons: ' + str(twoRecoLeptons)
print ''


c0=ROOT.TCanvas('c0', 'c0', 600, 400)

trackEffVsR = ROOT.TH1F('trackEffVsR', 'Track reconstruction efficiency', nbins, 0., TPCRmax)
trackEffVsR.SetXTitle('True vertex distance from beam axis')
trackEffVsR.SetYTitle('Track reco. efficiency')
trackEffVsR.Divide(trackEffVsR_reco, trackEffVsR_true,1,1,"cl=0.683 b(1,1) mode")
trackEffVsR.SetMarkerStyle(20)
trackEffVsR.SetMarkerSize(0.5)
trackEffVsR.Draw('e')
c0.SaveAs('plots/idm/'+fs+'/trackEffVsR_'+fs+'.pdf')

momTrack = ROOT.TEfficiency(momTrack_reco,momTrack_true)
momTrack.SetTitle('Track reconstruction efficiency;True momentum [GeV];Efficiency')
momTrack.SetMarkerStyle(20)
momTrack.SetMarkerSize(0.5)
momTrack.Draw('ap')
c0.SaveAs('plots/idm/'+fs+'/momTrack_'+fs+'.pdf')

ptTrack = ROOT.TEfficiency(ptTrack_reco,ptTrack_true)
ptTrack.SetTitle('Track reconstruction efficiency;True p_{T} [GeV];Efficiency')
ptTrack.SetMarkerStyle(20)
ptTrack.SetMarkerSize(0.5)
ptTrack.Draw('ap')
c0.SaveAs('plots/idm/'+fs+'/ptTrack_'+fs+'.pdf')

momLLP = ROOT.TEfficiency(momLLP_reco,momLLP_true)
momLLP.SetTitle('Track reconstruction efficiency;True LLP momentum [GeV];Efficiency')
momLLP.SetMarkerStyle(20)
momLLP.SetMarkerSize(0.5)
momLLP.Draw('ap')
c0.SaveAs('plots/idm/'+fs+'/momLLP_'+fs+'.pdf')

thetaTrack = ROOT.TEfficiency(thetaTrack_reco,thetaTrack_true)
thetaTrack.SetTitle('Track reconstruction efficiency;True polar angle;Efficiency')
thetaTrack.SetMarkerStyle(20)
thetaTrack.SetMarkerSize(0.5)
thetaTrack.Draw('ap')
c0.SaveAs('plots/idm/'+fs+'/thetaTrack_'+fs+'.pdf')

print 'DeltaTheta: ', DeltaThetaTrackLLP_reco.GetEntries(), DeltaThetaTrackLLP_true.GetEntries()
DeltaThetaTrackLLP = ROOT.TEfficiency(DeltaThetaTrackLLP_reco,DeltaThetaTrackLLP_true)
DeltaThetaTrackLLP.SetTitle('Track reconstruction efficiency;|#theta_{LLP} - #theta_{track}|;Efficiency')
DeltaThetaTrackLLP.SetMarkerStyle(20)
DeltaThetaTrackLLP.SetMarkerSize(0.5)
DeltaThetaTrackLLP.Draw('ap')
c0.SaveAs('plots/idm/'+fs+'/DeltaThetaTrackLLP_'+fs+'.pdf')

print 'DeltaPhi: ', DeltaPhiTrackLLP_reco.GetEntries(), DeltaPhiTrackLLP_true.GetEntries()
DeltaPhiTrackLLP = ROOT.TEfficiency(DeltaPhiTrackLLP_reco,DeltaPhiTrackLLP_true)
DeltaPhiTrackLLP.SetTitle('Track reconstruction efficiency;|#phi_{LLP} - #phi_{track}|;Efficiency')
DeltaPhiTrackLLP.SetMarkerStyle(20)
DeltaPhiTrackLLP.SetMarkerSize(0.5)
DeltaPhiTrackLLP.Draw('ap')
c0.SaveAs('plots/idm/'+fs+'/DeltaPhiTrackLLP_'+fs+'.pdf')

trackEffVsR_true.Scale( 1./trackEffVsR_true.Integral() )
trackEffVsR_true.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/R_true_'+fs+'.pdf')
trackEffVsR_reco.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/R_reco_'+fs+'.pdf')
momTrack_true.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/momTrack_true_'+fs+'.pdf')
thetaTrack_true.Scale( 1./thetaTrack_true.Integral() )
thetaTrack_true.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/thetaTrack_true_'+fs+'.pdf')
DeltaThetaTrackLLP_true.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/DeltaThetaTrackLLP_true_'+fs+'.pdf')
DeltaThetaTrackLLP_reco.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/DeltaThetaTrackLLP_reco_'+fs+'.pdf')
DeltaPhiTrackLLP_true.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/DeltaPhiTrackLLP_true_'+fs+'.pdf')
h_nTracks.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/nTracks_'+fs+'.pdf')
h_relWeights.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/relWeights_'+fs+'.pdf')
h_angDist.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/angularDistance_'+fs+'.pdf')
hmass.Draw('hist')
c0.SaveAs('plots/idm/'+fs+'/mass_'+fs+'.pdf')

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
line2 = ROOT.TLine(-2350.,1770.0,2350.0,1770.0);
line3 = ROOT.TLine(-2350.,329.0,-2350.0,1770.0);
line4 = ROOT.TLine(2350.,329.0,2350.0,1770.0);
line1.SetLineColor(2);
line2.SetLineColor(2);
line3.SetLineColor(2);
line4.SetLineColor(2);
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/idm/'+fs+'/trackEffRVsZ_'+fs+'.pdf')

trackEffiPvsTheta = ROOT.TEfficiency(trackEffiPvsTheta_reco,trackEffiPvsTheta_true)
trackEffiPvsTheta.SetTitle('Track reconstruction efficiency;#theta_{track};Momentum;Efficiency')
trackEffiPvsTheta.Draw('colz text')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/idm/'+fs+'/trackEffiPvsTheta_'+fs+'.pdf')

trackEffRVsZ_true.SetYTitle('R [mm]')
trackEffRVsZ_true.SetXTitle('z [mm]')
trackEffRVsZ_true.SetZTitle('Number of vertices')
trackEffRVsZ_true.SetStats(0)
trackEffRVsZ_true.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/idm/'+fs+'/trueVerticesRVsZ_'+fs+'.pdf')

recoVertices.SetYTitle('R [mm]')
recoVertices.SetXTitle('z [mm]')
recoVertices.SetZTitle('Number of vertices')
recoVertices.SetStats(0)
recoVertices.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/idm/'+fs+'/recoVerticesRVsZ_'+fs+'.pdf')

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
c1.SaveAs('plots/idm/'+fs+'/trackEffiPvsTheta_true_'+fs+'.pdf')

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
c1.SaveAs('plots/idm/'+fs+'/trackEffiPvsTheta_reco_'+fs+'.pdf')

trackEffiPvsR_true.SetYTitle('R [mm]')
trackEffiPvsR_true.SetXTitle('p [GeV]')
trackEffiPvsR_true.SetZTitle('Number of tracks')
trackEffiPvsR_true.SetStats(0)
#trackEffiPvsR_true.Scale( 1./trackEffiPvsR_true.Integral() )
trackEffiPvsR_true.Draw('colz text')
#line1.Draw('same')
#line2.Draw('same')
#line3.Draw('same')
#line4.Draw('same')
c1.SaveAs('plots/idm/'+fs+'/trackEffiPvsR_true_'+fs+'.pdf')

reader.close()
