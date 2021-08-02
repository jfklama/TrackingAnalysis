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
    print "Usage: " + os.path.basename(__file__) + " <inputFile>"
    sys.exit(1)



chargedChannel    = 0
withinAcceptance  = 0
twoTracks         = 0
allTracks		  = 0
matchingTracks    = 0
evWithV0		  = 0
pandoraPID		  = 0
V0FinderPID		  = 0

nbins = 12
TPCRmax = 1974.
#TPCRmax = 1770.
TPCZmax = 2350.

trackEffVsR_true = ROOT.TH1F('trackEffVsR_true', 'True vertices', nbins, 0., TPCRmax)
trackEffVsR_reco = ROOT.TH1F('trackEffVsR_reco', 'N_{evt} with 2 reco. tracks', nbins, 0., TPCRmax)
R_v0_TrueReco = ROOT.TH1F('R_v0_TrueReco', 'Distance of reco. V0 vtx. from true vtx.', 100, 0., 20.)
R_v0_true_IP = ROOT.TH1F('R_v0_true_IP', 'Distance of reco. V0 and true vtx. from IP', 150, 0., 3000.)
R_v0_reco_IP = ROOT.TH1F('R_v0_reco_IP', 'Distance of reco. V0 and true vtx. from IP', 150, 0., 3000.)

momV0_true = ROOT.TH1F('momV0_true', 'True momentum of V0', 11, 0., 55.)
momV0_reco = ROOT.TH1F('momV0_reco', 'True momentum of reco. V0', 11, 0., 55.)

V0TrueVtxRVsZ_true = ROOT.TH2F('V0TrueVtxRVsZ_true', 'True vertices', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
V0TrueVtxRVsZ_reco = ROOT.TH2F('V0TrueVtxRVsZ_reco', 'True vertices', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffRVsZ_true = ROOT.TH2F('trackEffRVsZ_true', 'N_{evt} with 2 reco. tracks', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffRVsZ_reco = ROOT.TH2F('trackEffRVsZ_reco', 'N_{evt} with 2 reco. tracks', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)

V0vtxRVsZ = ROOT.TH2F('V0vtxRVsZ', 'N_{evt} with one V0', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)

hmass = ROOT.TH1F('hmass', '#V0 mass', 100, 0., 5.)
henergy = ROOT.TH1F('henergy', '#V0 truth energy', 10, 0., 50.)

infile = sys.argv[1:]

if infile[0].find('lambda') != -1:
	v0PID = 3122
	daughter1PID = 2212
	daughter2PID = -211
elif infile[0].find('kaon') != -1:
	v0PID = 310
	daughter1PID = 211
	daughter2PID = -211
else:
	print 'Provide correct input file!'
	sys.exit()

reader=IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(infile)

i=0
for event in reader:
	i+=1

	mcCollection = event.getCollection('MCParticlesSkimmed')
	#tpcCollection = event.getCollection('TPCCollection')
	trackCollection = event.getCollection('MarlinTrkTracks')
	trackToMCLinkCollection = event.getCollection('MarlinTrkTracksMCTruthLink')
	pfoCollection = event.getCollection('PandoraPFOs')
	if 'V0Vertices' in event.getCollectionNames():
		v0VtxCollection = event.getCollection('V0Vertices')
		v0RecoPartCollection = event.getCollection('V0RecoParticles')
	else:
		v0VtxCollection = event.getCollection('BuildUpVertex_V0')
		for i in range(v0VtxCollection.getNumberOfElements()):
			v0VtxCollection.removeElementAt(i)
		v0RecoPartCollection = event.getCollection('BuildUpVertex_V0_RP')
		for i in range(v0RecoPartCollection.getNumberOfElements()):
			v0RecoPartCollection.removeElementAt(i)

	if i%100==0:
		print (i, str(trackCollection.getNumberOfElements()) + ' tracks in event')


	### MC TRUTH ###

	for particle in mcCollection:
		if particle.getPDG() == v0PID:
			henergy.Fill(particle.getEnergy())
			break

	#check if event contains proton and pion
	isDaughter1 = False
	isDaughter2 = False
	for particle in mcCollection:
		if particle.getPDG() == daughter1PID:
			isDaughter1 = True
		if particle.getPDG() == daughter2PID:
			isDaughter2 = True

	if isDaughter1 != True or isDaughter2 != True:
		#print 'Suppressed because of decay mode'
		continue
	else:
		chargedChannel+=1

	#check if V0 decayed within tracking acceptance
	trueR = 0
	trueZ = 0
	V0X = 0
	V0Y = 0
	V0Z = 0
	px = 0
	py = 0
	pz = 0
	isV0MC = False
	for particle in mcCollection:
		if abs(particle.getPDG()) == v0PID:
			isV0MC = True
			V0X = particle.getEndpoint()[0]
			V0Y = particle.getEndpoint()[1]
			V0Z = particle.getEndpoint()[2]
			trueR = R(V0X, V0Y)
			trueZ = V0Z
			trueR3 = R3(V0X, V0Y, V0Z)
			px = particle.getMomentum()[0]
			py = particle.getMomentum()[1]
			pz = particle.getMomentum()[2]
			p_truth = R3(px,py,pz)
			break
	if not isV0MC:
		#print i, 'Suppressed because no MC V0'
		continue
	if trueR > 1770 or abs(trueZ) > 2350:
		#print 'Suppressed because out of det. accept.'
		continue
	else:
		withinAcceptance+=1

	trackEffVsR_true.Fill(trueR)
	trackEffRVsZ_true.Fill(trueZ,trueR)

	# events with at least 2 reconstructed tracks
	numberOfTracks = trackCollection.getNumberOfElements()
	if int(numberOfTracks) >= 2:
		twoTracks+=1
		trackEffVsR_reco.Fill(trueR)
		trackEffRVsZ_reco.Fill(trueZ,trueR) # for track reco. efficiency

		V0TrueVtxRVsZ_true.Fill(trueZ,trueR) # for v0 reco. efficiency
		R_v0_true_IP.Fill(trueR3)
		momV0_true.Fill(p_truth)
	else:
		continue

	if v0VtxCollection.getNumberOfElements() >= 1:
		for vertex in v0VtxCollection:
			v0R = R(vertex.getPosition()[0], vertex.getPosition()[1])
			v0Z = vertex.getPosition()[2]
			#if abs( v0R-trueR ) < 3 and abs( v0Z-trueZ ) < 3:
			R_trueReco = R3( vertex.getPosition()[0]-V0X, vertex.getPosition()[1]-V0Y, vertex.getPosition()[2]-V0Z )
			R_v0_TrueReco.Fill(R_trueReco) # distance between true and reco. vtx
			R_reco = R3( vertex.getPosition()[0], vertex.getPosition()[1], vertex.getPosition()[2] ) # vtx dist. from IP
			if R_trueReco < 5:
				V0vtxRVsZ.Fill(v0Z,v0R)
				V0TrueVtxRVsZ_reco.Fill(trueZ,trueR)
				R_v0_reco_IP.Fill( R_reco )
				evWithV0+=1
				break

	isV0 = False
	if v0RecoPartCollection.getNumberOfElements() >= 1:
		for particle in v0RecoPartCollection:
			if particle.getType() == v0PID:
				V0FinderPID+=1
				isV0 = True
				momV0_reco.Fill(p_truth)
				break
			#else:
				#print i

	#count reconstructed tracks
	for track in trackCollection:
		allTracks+=1

	#count tracks matching to MC particles
	for track in trackToMCLinkCollection:
		#allTracks+=1
		for particle in mcCollection:
			if track.getTo() == particle and track.getWeight() > 0.9:
				matchingTracks+=1
			#if trackCollection.getNumberOfElements() == 1:
				#print track.getTo().getPDG()

	#if trackCollection.getNumberOfElements() < trackToMCLinkCollection.getNumberOfElements():
	#	print (i, trackCollection.getNumberOfElements(), trackToMCLinkCollection.getNumberOfElements() )

	#calculate V0 mass
	V0P4 = ROOT.TLorentzVector()
	V0P4.SetPxPyPzE(0,0,0,0)
	for particle in pfoCollection:
		if particle.getType() == v0PID:
			#continue
			#print trackCollection.getNumberOfElements()
			#print i
			pandoraPID+=1
			p4 = ROOT.TLorentzVector()
			p4.SetPxPyPzE( particle.getMomentum()[0],particle.getMomentum()[1],particle.getMomentum()[2],particle.getEnergy() )
			V0P4 = p4
			break
		else:
			continue
			#print trackCollection.getNumberOfElements()
			#p4 = ROOT.TLorentzVector()
			#p4.SetPxPyPzE( particle.getMomentum()[0],particle.getMomentum()[1],particle.getMomentum()[2],particle.getEnergy() )
			#V0P4 += p4
	mass = V0P4.M()
	if mass != 0:
		hmass.Fill(mass)
	else:
		if isV0:
			print i

	#if mass < 0.7 or mass > 3.5:
	#	for particle in pfoCollection:
	#		print (i, particle.getMomentum()[0], particle.getMomentum()[1], particle.getMomentum()[2], particle.getEnergy(), particle.getMass())

	#if i%100==0:
	#	print ('Reco. V0 mass: ' + str(mass))

reader.close()

c0=ROOT.TCanvas('c0', 'c0', 200, 400)

trackEffVsR = ROOT.TH1F('trackEffVsR', 'Track reconstruction efficiency', nbins, 0., TPCRmax)
trackEffVsR.SetXTitle('True vertex distance from beam axis')
trackEffVsR.SetYTitle('Track reco. efficiency')
trackEffVsR.Divide(trackEffVsR_reco, trackEffVsR_true,1,1,"cl=0.683 b(1,1) mode")
trackEffVsR.SetMarkerStyle(20)
trackEffVsR.SetMarkerSize(0.5)
trackEffVsR.Draw('e')
c0.SaveAs('plots/trackEffVsR_' + str(nbins) + 'bins.pdf')


v0EffVsMom = ROOT.TH1F('v0EffVsMom', '', 11, 0., 55.)
v0EffVsMom.SetXTitle('True V0 momentum')
v0EffVsMom.SetYTitle('V0 finding efficiency')
v0EffVsMom.Divide(momV0_reco, momV0_true,1,1,"cl=0.683 b(1,1) mode")
v0EffVsMom.SetMarkerStyle(20)
v0EffVsMom.SetMarkerSize(0.5)
v0EffVsMom.SetMinimum(0.0)
v0EffVsMom.SetMaximum(1.1)
v0EffVsMom.Draw('e')
c0.SaveAs('plots/v0EffVsMom.pdf')

trackEffVsR_true.Draw('hist')
c0.SaveAs('plots/R_true_' + str(nbins) + 'bins.pdf')
trackEffVsR_reco.Draw('hist')
c0.SaveAs('plots/R_reco_' + str(nbins) + 'bins.pdf')
hmass.Draw('hist')
c0.SaveAs('plots/mass.pdf')
henergy.Draw('hist')
c0.SaveAs('plots/energy.pdf')
c0.SetLogy()
R_v0_TrueReco.Draw('hist')
c0.SaveAs('plots/R_v0_TrueReco.pdf')

R_v0_true_IP.Draw('hist')
R_v0_reco_IP.SetLineColor(2)
R_v0_reco_IP.Draw('hist same')
c0.SaveAs('plots/R_v0_TrueReco_IP.pdf')

c1=ROOT.TCanvas('c1', 'c1', 800, 400)
c1.SetRightMargin(0.15)

trackEffRVsZ = ROOT.TH2F('trackEffRVsZ', 'Track reconstruction efficiency', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
trackEffRVsZ.SetYTitle('R [mm]')
trackEffRVsZ.SetXTitle('z [mm]')
trackEffRVsZ.SetZTitle('Track reco. efficiency')
trackEffRVsZ.Divide(trackEffRVsZ_reco, trackEffRVsZ_true,1,1,"cl=0.683 b(1,1) mode")
trackEffRVsZ.SetStats(0)
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
c1.SaveAs('plots/trackEffRVsZ_' + str(nbins) + 'bins.pdf')

V0TrueVtxRVsZ_true.SetYTitle('R [mm]')
V0TrueVtxRVsZ_true.SetXTitle('z [mm]')
V0TrueVtxRVsZ_true.SetZTitle('Number of vertices')
V0TrueVtxRVsZ_true.SetStats(0)
V0TrueVtxRVsZ_true.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/trueVerticesRVsZ_' + str(nbins) + 'bins.pdf')

V0vtxRVsZ.SetYTitle('R [mm]')
V0vtxRVsZ.SetXTitle('z [mm]')
V0vtxRVsZ.SetZTitle('V0 vertices')
V0vtxRVsZ.SetStats(0)
V0vtxRVsZ.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/V0vtxReco_' + str(nbins) + 'bins.pdf')

v0EffRVsZ = ROOT.TH2F('v0EffRVsZ', 'V0 reconstruction efficiency', nbins, -TPCZmax, TPCZmax, nbins, 0., TPCRmax)
v0EffRVsZ.SetYTitle('R [mm]')
v0EffRVsZ.SetXTitle('z [mm]')
v0EffRVsZ.SetZTitle('V0 reco. efficiency')
v0EffRVsZ.Divide(V0TrueVtxRVsZ_reco, V0TrueVtxRVsZ_true,1,1,"cl=0.683 b(1,1) mode")
v0EffRVsZ.SetStats(0)
v0EffRVsZ.Draw('colz')
line1.Draw('same')
line2.Draw('same')
line3.Draw('same')
line4.Draw('same')
c1.SaveAs('plots/v0VerticesEffRVsZ_' + str(nbins) + 'bins.pdf')


print('Ev. with p and pi', 'Decays within tracker acceptance', 'Ev. with two tracks', 'All reco. tracks', 'Tracks matching to MC')
print(chargedChannel, withinAcceptance, twoTracks, allTracks, matchingTracks)
print 'Ev. with correct V0 vtx: ' + str(evWithV0)
print 'Ev. with correct V0 by PandoraPID: ' + str(pandoraPID)
print 'Ev. with correct V0 by V0FinderPID: ' + str(V0FinderPID)

reader.close()
