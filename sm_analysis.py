import sys
import os

import math
import numpy as np

from pyLCIO import IOIMPL, ROOT

def R(x,y):
	return math.sqrt(x**2 + y**2)
def R3(x,y,z):
	return math.sqrt(x**2 + y**2 + z**2)

try:
    arg1 = sys.argv[1]
except IndexError:
    print "Usage: " + os.path.basename(__file__) + " <inputFile> <considered momenta> (optional)"
    sys.exit(1)



chargedChannel    = 0
withinAcceptance  = 0
oneTrack         = 0
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
trackEffiPvsTheta_true = ROOT.TH2F('trackEffiPvsTheta_true', '', 3, 0., 3.15, 22, 0., 22.)
trackEffiPvsTheta_reco = ROOT.TH2F('trackEffiPvsTheta_reco', '', 3, 0., 3.15, 22, 0., 22.)
trackEffiPvsR_true = ROOT.TH2F('trackEffiPvsR_true', '', 22, 0.5, 22.5, nbins, 0., TPCRmax)
trackEffiPvsR_reco = ROOT.TH2F('trackEffiPvsR_reco', '', 22, 0.5, 22.5, nbins, 0., TPCRmax)


momTrack_true = ROOT.TH1F('momTrack_true', 'True momentum of track', 50, 0., 25.)
momTrack_reco = ROOT.TH1F('momTrack_reco', 'True momentum of reco. track', 50, 0., 25.)
thetaTrack_true = ROOT.TH1F('thetaTrack_true', 'True polar angle of track', 32, 0., 3.15)
thetaTrack_reco = ROOT.TH1F('thetaTrack_reco', 'True polar angle of reco. track', 32, 0., 3.15)

h_nTracks = ROOT.TH1F('h_nTracks', 'Number of tracks', 20, 0., 20.)
h_relWeights = ROOT.TH1F('h_relWeights', 'Number of weights', 50, 0., 2.)

recoVertices = ROOT.TH2F('recoVertices', 'N_{evt} with 2 reco. leptons', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)

hmass = ROOT.TH1F('hmass', '#Lambda mass', 150, 0., 15.)

infile = sys.argv[1]
reader=IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(infile)

infileName = str(infile)
particleType = infileName[infileName.find("sm_")+3:infileName.find("sm_")+7]

print 'Analysed particle: ' + particleType

sys_argv = sys.argv[2:]
momenta = [float(p) for p in sys_argv]

i=0
for event in reader:
	i+=1

	mcCollection = event.getCollection('MCParticlesSkimmed')
	#tpcCollection = event.getCollection('TPCCollection')
	trackCollection = event.getCollection('MarlinTrkTracks')
	#trackCollection = event.getCollection('ClupatraTracks')
	trackToMCLinkCollection = event.getCollection('MarlinTrkTracksMCTruthLink')
	pfoCollection = event.getCollection('PandoraPFOs')
	#vertexCollection = event.getCollection('PandoraPFANewStartVertices')

	if i%1000==0:
		print (i, str(trackCollection.getNumberOfElements() ) + ' tracks in event')

	#trackEffVsR_true.Fill(trueR, 2.)
	#trackEffRVsZ_true.Fill(trueZ,trueR, 2.)

	trueR = 0
	trueZ = 0
	p = 0
	for particle in mcCollection:
		if particle.getGeneratorStatus() == 1:
			trueR = R(particle.getVertex()[0], particle.getVertex()[1])
			trueZ = particle.getVertex()[2]
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = R3(px,py,pz)
			#print p
			break

	# Use only events with selected particle momenta
	ok = True
	for mom in momenta:
		if abs(mom - p) < 0.01:
			ok = True
			break
		else:
			ok = False
	if not ok:
		continue

	nMuons = 0
	nElectrons = 0
	lepP4 = ROOT.TLorentzVector()
	for particle in mcCollection:
		if(particle.getGeneratorStatus() == 1):
			e = particle.getEnergy()
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = R3(px,py,pz)
			lepP4.SetPxPyPzE(px,py,pz,e)
			cosTheta = abs(math.cos(lepP4.Theta()))
			if abs( particle.getPDG() ) == 11:
				nElectrons+=1
				#if lepP4.Theta() > 0.785 and lepP4.Theta() < math.pi - 0.785:
				trackEffVsR_true.Fill(trueR)
				trackEffRVsZ_true.Fill(trueZ,trueR)
				momTrack_true.Fill(p)
				thetaTrack_true.Fill(lepP4.Theta())
				trackEffiPvsTheta_true.Fill(lepP4.Theta(),p)
				trackEffiPvsR_true.Fill(p,trueR)
			if abs( particle.getPDG() ) == 13:
				nMuons+=1
				#if lepP4.Theta() > 0.785 and lepP4.Theta() < math.pi - 0.785:
				trackEffVsR_true.Fill(trueR)
				trackEffRVsZ_true.Fill(trueZ,trueR)
				momTrack_true.Fill(p)
				thetaTrack_true.Fill(lepP4.Theta())
				trackEffiPvsTheta_true.Fill(lepP4.Theta(),p)
				trackEffiPvsR_true.Fill(p,trueR)
		if nElectrons > 0 or nMuons > 0 or (nElectrons>0 and nMuons>0):
			break
		if nElectrons>0 and nMuons>0:
			print (i, 'electron and muon!')
	if nElectrons<1 and nMuons<1:
		print (i, 'WARNING: no MC electrons and muons!')
		continue

	#count events with at least 2 reconstructed tracks
	numberOfTracks = trackCollection.getNumberOfElements()
	if int(numberOfTracks) >= 1:
		oneTrack+=1
	#if numberOfTracks > 2:
		#print i, numberOfTracks
	h_nTracks.Fill(numberOfTracks)

	#count reconstructed tracks
	for track in trackCollection:
		allTracks+=1

	#count tracks matching to MC particles
	matchingTracks = 0
	for particle in mcCollection:
		if(particle.getGeneratorStatus() == 1):
			for track in trackToMCLinkCollection:
				h_relWeights.Fill(track.getWeight())
				if track.getTo() == particle and track.getWeight() >= 1.0:
					allMatchingTracks+=1
					#print track.getWeight()


	# fill histos if tracks are matching
	#if matchingTracks > -1:
	matchingTracks = 0
	for particle in mcCollection:
		if(particle.getGeneratorStatus() == 1):
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p = R3(px,py,pz)
			e = particle.getEnergy()
			lepP4.SetPxPyPzE(px,py,pz,e)
			cosTheta = abs( math.cos(lepP4.Theta()) )
			#print (i, "genStatus = 1")
			for track in trackToMCLinkCollection:
				if track.getTo() == particle and track.getWeight() >= 1.0:
					#if lepP4.Theta() > 0.785 and lepP4.Theta() < math.pi - 0.785:
					momTrack_reco.Fill(p)
					thetaTrack_reco.Fill(lepP4.Theta())
					trackEffiPvsTheta_reco.Fill(lepP4.Theta(),p)
					trackEffVsR_reco.Fill(trueR)
					trackEffRVsZ_reco.Fill(trueZ,trueR)
					trackEffiPvsR_reco.Fill(p,trueR)
					matchingTracks+=1
					break
				#if track.getWeight() > 0.9 and track.getWeight() < 1.0:
					#print (i, matchingTracks, numberOfTracks)
			#if p >=2 and p <= 3:
				#print (i, numberOfTracks, matchingTracks, p)
		#if matchingTracks == 2:
			#print (i, numberOfTracks, matchingTracks, p)
		'''
		for track in trackToMCLinkCollection:
			if track.getWeight() > 0.9 and track.getWeight() < 1.2:
				particle = track.getTo()
				e = particle.getEnergy()
				px = particle.getMomentum()[0]
				py = particle.getMomentum()[1]
				pz = particle.getMomentum()[2]
				p = R3(px,py,pz)
				lepP4.SetPxPyPzE(px,py,pz,e)
				momTrack_reco.Fill(p)
				thetaTrack_reco.Fill(lepP4.Theta())
				trackEffiPvsTheta_reco.Fill(abs(math.cos(lepP4.Theta())),p)
				trackEffVsR_reco.Fill(trueR)
				trackEffRVsZ_reco.Fill(trueZ,trueR)
		'''
	# check if electron and muon PFOs are reconstructed
	nMuons = 0
	nElectrons = 0
	for particle in pfoCollection:
			if abs( particle.getType() ) == 11:
				nElectrons+=1
			if abs( particle.getType() ) == 13:
				nMuons+=1
	'''
	recoR = 0
	recoZ = 0
	if nMuons >= 2:
		#print (i, nMuons, nElectrons, trueR, trueZ)
		for particle in pfoCollection:
			if abs( particle.getType() ) == 13:
				for vertex in vertexCollection:
					if vertex.id() == particle.id():
						recoR = R(vertex.getPosition()[0], vertex.getPosition()[1])
						recoZ = vertex.getPosition()[2]
	elif nElectrons >= 2:
		#print (i, nMuons, nElectrons, trueR, trueZ)
		for particle in pfoCollection:
			if abs( particle.getType() ) == 11:
				for vertex in vertexCollection:
					if vertex.id() == particle.id():
						recoR = R(vertex.getPosition()[0], vertex.getPosition()[1])
						recoZ = vertex.getPosition()[2]

	recoVertices.Fill(recoZ,recoR)
	'''
	if nMuons > 1 or nElectrons > 1:
		twoRecoLeptons+=1
		#continue
	#calculate lambda mass
	ZP4 = ROOT.TLorentzVector()
	ZP4.SetPxPyPzE(0,0,0,0)
	for particle in pfoCollection:
		if abs( particle.getType() ) == 11 or abs( particle.getType() ) == 13:
			p4 = ROOT.TLorentzVector()
			p4.SetPxPyPzE( particle.getMomentum()[0],particle.getMomentum()[1],particle.getMomentum()[2],particle.getEnergy() )
			ZP4 += p4
	mass = ZP4.M()
	hmass.Fill(mass)
	'''
	if mass < 0.7 or mass > 3.5:
		for particle in pfoCollection:
			print (i, particle.getMomentum()[0], particle.getMomentum()[1], particle.getMomentum()[2], particle.getEnergy(), particle.getMass())
	'''
	#if i%100==0:
	#	print ('Reco. Lambda mass: ' + str(mass))

reader.close()

print('Ev. with min. one track', 'All reco. tracks', 'Tracks matching to MC', 'Two reco. leptons')
print(oneTrack, allTracks, allMatchingTracks, twoRecoLeptons)

c0=ROOT.TCanvas('c0', 'c0', 600, 400)

trackEffVsR = ROOT.TH1F('trackEffVsR', 'Track reconstruction efficiency', nbins, 0., TPCRmax)
trackEffVsR.SetXTitle('True vertex distance from beam axis')
trackEffVsR.SetYTitle('Track reco. efficiency')
trackEffVsR.Divide(trackEffVsR_reco, trackEffVsR_true,1,1,"cl=0.683 b(1,1) mode")
trackEffVsR.SetMarkerStyle(20)
trackEffVsR.SetMarkerSize(0.5)
trackEffVsR.Draw('e')
c0.SaveAs('plots/sm/'+particleType+'/trackEffVsR_'+particleType+'.pdf')

#momTrack = ROOT.TH1F('momTrack', 'Track reconstruction efficiency', 75, 0., 15.)
#momTrack = ROOT.TEfficiency("momTrack", "Track reconstruction efficiency;p_truth;efficiency",75, 0., 15.)
#momTrack.SetTotalHistogram(momTrack_true,'')
#momTrack.SetPassedHistogram(momTrack_reco,'')
#momTrack.SetXTitle('True momentum')
#momTrack.SetYTitle('Track reco. efficiency')
#momTrack.Divide(momTrack_reco, momTrack_true,1,1,"cl=0.683 b(1,1) mode")
momTrack = ROOT.TEfficiency(momTrack_reco,momTrack_true)
momTrack.SetTitle('Track reconstruction efficiency;True momentum [GeV];Efficiency')
momTrack.SetMarkerStyle(20)
momTrack.SetMarkerSize(0.5)
momTrack.Draw('ap')
c0.SaveAs('plots/sm/'+particleType+'/momTrack_'+particleType+'.pdf')

thetaTrack = ROOT.TEfficiency(thetaTrack_reco,thetaTrack_true)
thetaTrack.SetTitle('Track reconstruction efficiency;True polar angle;Efficiency')
thetaTrack.SetMarkerStyle(20)
thetaTrack.SetMarkerSize(0.5)
thetaTrack.Draw('ap')
c0.SaveAs('plots/sm/'+particleType+'/thetaTrack_'+particleType+'.pdf')

trackEffVsR_true.Scale( 1./trackEffVsR_true.Integral() )
trackEffVsR_true.Draw('hist')
c0.SaveAs('plots/sm/'+particleType+'/R_true_'+particleType+'.pdf')
trackEffVsR_reco.Draw('hist')
c0.SaveAs('plots/sm/'+particleType+'/R_reco_'+particleType+'.pdf')
momTrack_true.Draw('hist')
c0.SaveAs('plots/sm/'+particleType+'/momTrack_true_'+particleType+'.pdf')
thetaTrack_true.Scale( 1./thetaTrack_true.Integral() )
thetaTrack_true.Draw('hist')
c0.SaveAs('plots/sm/'+particleType+'/thetaTrack_true_'+particleType+'.pdf')
h_nTracks.Draw('hist')
c0.SaveAs('plots/sm/'+particleType+'/nTracks_'+particleType+'.pdf')
h_relWeights.Draw('hist')
c0.SaveAs('plots/sm/'+particleType+'/relWeights_'+particleType+'.pdf')
hmass.Draw('hist')
c0.SaveAs('plots/sm/'+particleType+'/mass_'+particleType+'.pdf')

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
c1.SaveAs('plots/sm/'+particleType+'/trackEffRVsZ_'+particleType+'.pdf')

trackEffiPvsTheta = ROOT.TEfficiency(trackEffiPvsTheta_reco,trackEffiPvsTheta_true)
trackEffiPvsTheta.SetTitle('Track reconstruction efficiency;#theta_{track};Momentum;Efficiency')
trackEffiPvsTheta.Draw('colz text')
c1.SaveAs('plots/sm/'+particleType+'/trackEffiPvsTheta_'+particleType+'.pdf')

trackEffiPvsR = ROOT.TEfficiency(trackEffiPvsR_reco,trackEffiPvsR_true)
trackEffiPvsR.SetTitle('Track reconstruction efficiency;p [GeV];R [mm];Efficiency')
trackEffiPvsR.Draw('colz text')
c1.SaveAs('plots/sm/'+particleType+'/trackEffiPvsR_'+particleType+'.pdf')

trackEffRVsZ_true.SetYTitle('R [mm]')
trackEffRVsZ_true.SetXTitle('z [mm]')
trackEffRVsZ_true.SetZTitle('Number of vertices')
trackEffRVsZ_true.SetStats(0)
trackEffRVsZ_true.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/sm/'+particleType+'/trueVerticesRVsZ_'+particleType+'.pdf')

recoVertices.SetYTitle('R [mm]')
recoVertices.SetXTitle('z [mm]')
recoVertices.SetZTitle('Number of vertices')
recoVertices.SetStats(0)
recoVertices.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/sm/'+particleType+'/recoVerticesRVsZ_'+particleType+'.pdf')

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
c1.SaveAs('plots/sm/'+particleType+'/trackEffiPvsTheta_true_'+particleType+'.pdf')

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
c1.SaveAs('plots/sm/'+particleType+'/trackEffiPvsTheta_reco_'+particleType+'.pdf')

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
c1.SaveAs('plots/sm/'+particleType+'/trackEffiPvsR_true_'+particleType+'.pdf')

reader.close()
