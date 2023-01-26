import math
import numpy as np

from pyLCIO import IOIMPL, IMPL, EVENT, ROOT
#import dd4hep

'''
description = dd4hep.Detector.getInstance()
bfieldV = [0.,0.,0.]
#print bfieldV
description.field().magneticField( [ 0., 0., 0. ]  , bfieldV  )
bField = bfieldV[2]/dd4hep.tesla
'''

bField = 3.5 # HARDCODED!!!
FCT = 2.99792458E-4
c = 2.99792458E+14

def R(x,y):
	return math.sqrt(x**2 + y**2)

def R3(x,y,z):
	return math.sqrt(x**2 + y**2 + z**2)

def logx(nbins, max):
	xbins = np.empty(nbins-1)
	dx = max/nbins
	l10 = math.log(10)
	for i in range(nbins-1):
		xbins[i] = math.exp(l10*i*dx)
		print xbins[i]
	return xbins

def getTrackMomentum(track):

	momentum = [0.,0.,0.]

  	#d0 = track.getD0();
	#z0 = track.getZ0();
	phi = track.getPhi()
	tanLambda = track.getTanLambda()
	omega = track.getOmega()

	radius = 1./abs(omega)
	pxy = FCT * bField * radius
	momentum[0] = pxy * math.cos(phi)
	momentum[1] = pxy * math.sin(phi)
	momentum[2] = tanLambda * pxy

	return momentum

def getTrackCharge(track):
	omega = track.getOmega()
	charge = omega/abs(omega)
	return charge

def phiDistance(phi1, phi2):

	dPhi = phi1 - phi2
	if abs(dPhi) > math.pi:
		dPhi = 2 * math.pi - abs(dPhi)

	return abs(dPhi)

def getArcLength(phi1, phi2, q):

	if q < 0:
		dPhi = phi2 - phi1
		if dPhi < 0:
			dPhi += 2 * math.pi
	else:
		dPhi = phi1 - phi2
		if dPhi < 0:
			dPhi += 2 * math.pi

	return abs(dPhi)

def angularDistance(theta1, theta2, phi1, phi2):

	dPhi = phi1 - phi2
	if abs(dPhi) > math.pi:
		dPhi = 2 * math.pi - abs(dPhi)
	'''
	theta1 -= math.pi / 2
	theta2 -= math.pi / 2
	phi1   += math.pi
	'''

	dist = math.sqrt( (theta1-theta2)**2 * (math.sin((theta1+theta2)/2))**2 + (dPhi)**2 )
	#dist = math.acos( \
	#	   math.cos(theta1) * math.cos(phi1) * math.cos(theta2) * math.cos(phi2) \
	# 	 + math.cos(theta1) * math.sin(phi1) * math.cos(theta2) * math.sin(phi2) \
	#	 + math.sin(theta1) * math.sin(theta2) \
	#	 )

	return dist

def spacialDistance(v1, v2):

	dist = (v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2
	return math.sqrt(dist)

def getNextTrackState(trackState0, hit, loc):

	x0 = trackState0.getReferencePoint()[0]
	y0 = trackState0.getReferencePoint()[1]
	z0 = trackState0.getReferencePoint()[2]

	d0 = trackState0.getD0()
	dz0 = trackState0.getZ0()
	phi0 = trackState0.getPhi()
	tanLambda = trackState0.getTanLambda()
	omega = trackState0.getOmega()
	alfa_kappa = 1. / (omega * c)
	r = 1. / omega

	# x-y coordinates of the circle of helix
	Xc = x0 + (d0+alfa_kappa) * math.cos(phi0)
	Yc = y0 + (d0+alfa_kappa) * math.sin(phi0)
	#if phi0 > 0:
	#	Xc = x0 - (d0+alfa_kappa) * math.cos(phi0)
	#	Yc = y0 - (d0+alfa_kappa) * math.sin(phi0)
	#	#Xc = x0 - (r-d0) * math.sin(phi0)
	#	#Yc = y0 - (r-d0) * math.cos(phi0)

	# new point in which we want to get a track state
	x = hit.getPosition()[0]
	y = hit.getPosition()[1]
	z = hit.getPosition()[2]

	referencePoint = trackState0.getReferencePoint()
	referencePoint[0] = x
	referencePoint[1] = y
	referencePoint[2] = z

	phi = 0.
	if omega > 0:
		phi = math.atan( (Yc-y)/(Xc-x) )
	else:
		phi = math.atan( (y-Yc)/(x-Xc) )
	dz = z0 - z + dz0 - alfa_kappa * (phi-phi0) * tanLambda
	d = (Xc-x) * math.cos(phi) + (Yc-y) * math.sin(phi) - alfa_kappa
	# d = math.sqrt( (Xc-x)**2 + (Yc-y)**2 ) - alfa_kappa

	trackState = IMPL.TrackStateImpl()

	trackState.setLocation(loc)
	trackState.setD0(d)
	trackState.setPhi(phi)
	trackState.setOmega(omega)
	trackState.setZ0(dz)
	trackState.setTanLambda(tanLambda)
	trackState.setReferencePoint(referencePoint)

	return trackState

def getAnyTrackState(track, myHit, loc):

	x = myHit.getPosition()[0]
	y = myHit.getPosition()[1]
	z = myHit.getPosition()[2]

	hits = track.getTrackerHits()
	ts = track.getTrackState(2)
	if (not ts):
		print 'ERROR: No TrackState in First Hit!'
		return

	myTrackState = IMPL.TrackStateImpl()

	for i in range(len(hits)-1):

		x0 = hit[i+1].getPosition()[0]
		y0 = hit[i+1].getPosition()[1]
		z0 = hit[i+1].getPosition()[2]

		ts = getNextTrackState(ts,hit[i+1],loc)
		#if hit[i+1]==myHit:

def fixTrackStateDirection(track,loc):

	oppos_loc = 3
	if loc == 3:
		oppos_loc = 2

	init_ts = track.getTrackState(loc)
	oppos_ts = track.getTrackState(oppos_loc)

	direction = 1
	if init_ts.getReferencePoint()[2] > oppos_ts.getReferencePoint()[2]:
		direction = -1

	fixed_ts = IMPL.TrackStateImpl()

	if direction * getTrackMomentum(init_ts)[2] > 0:
		fixed_ts = init_ts
	elif direction * getTrackMomentum(init_ts)[2] < 0:
		fixed_ts = flipTrackState(init_ts)
	else:
		print 'ERROR: dir * p_z = ' + str(direction * init_ts.getMomentum()[2])
		return

	return fixed_ts

def flipTrackState(ts):

	loc = ts.getLocation()
	d0 = ts.getD0()
	z0 = ts.getZ0()
	phi0 = ts.getPhi()
	tanLambda = ts.getTanLambda()
	omega = ts.getOmega()
	covMat = ts.getCovMatrix()
	refPoint = ts.getReferencePoint()

	phi = phi0
	if phi0 > 0:
		phi = -1*(math.pi - phi)
	else:
		phi = math.pi + phi

	flipped_ts = IMPL.TrackStateImpl()

	flipped_ts.setLocation(loc)
	flipped_ts.setD0(-d0)
	flipped_ts.setPhi(phi)
	flipped_ts.setOmega(-omega)
	flipped_ts.setZ0(z0)
	flipped_ts.setTanLambda(-tanLambda)
	flipped_ts.setReferencePoint(refPoint)
	flipped_ts.setCovMatrix(covMat)

	return flipped_ts

def invert(list):
	n = list.size()
	out_list = []
	for i in range(n):
		out_list.append( list[n-i-1] )
	return out_list

def findClosestTrackState(track,vtx_mc):

	ts_FirstHit = track.getTrackState(2) # 1 = atIP, 2 = atFirstHit, 3 = atLastHit
	ts_FirstHit = fixTrackStateDirection( track, ts_FirstHit.getLocation() )
	momentum = getTrackMomentum(ts_FirstHit)
	p_fh = ROOT.TVector3()
	p_fh.SetXYZ(momentum[0],momentum[1],momentum[2])
	vtx_fh = ts_FirstHit.getReferencePoint()
	vtx_dist_fh = spacialDistance(vtx_mc, vtx_fh)

	ts_LastHit = track.getTrackState(3) # 1 = atIP, 2 = atFirstHit, 3 = atLastHit
	ts_LastHit = fixTrackStateDirection( track, ts_LastHit.getLocation() )
	momentum = getTrackMomentum(ts_LastHit)
	p_lh = ROOT.TVector3()
	p_lh.SetXYZ(momentum[0],momentum[1],momentum[2])
	vtx_lh = ts_LastHit.getReferencePoint()
	vtx_dist_lh = spacialDistance(vtx_mc, vtx_lh)

	vtx_dist = vtx_dist_fh
	trackState = ts_FirstHit
	if vtx_dist_lh < vtx_dist_fh:
		vtx_dist = vtx_dist_lh
		trackState = ts_LastHit

	return trackState

def findClosestLegalTrackStates(trk1, trk2):

	trkStates = []
	trkStates.append( fixTrackStateDirection(trk1,2) )
	trkStates.append( fixTrackStateDirection(trk1,3) )
	trkStates.append( fixTrackStateDirection(trk2,2) )
	trkStates.append( fixTrackStateDirection(trk2,3) )

	goodTS1 = IMPL.TrackStateImpl()
	goodTS2 = IMPL.TrackStateImpl()
	minDist = 50000
	for iTS in range(2):
		for jTS in range(2,4):
			refPoint1 = trkStates[iTS].getReferencePoint()
			refPoint2 = trkStates[jTS].getReferencePoint()
			dist = spacialDistance( refPoint1,refPoint2 )
			# print dist,

			charge1 = trkStates[iTS].getOmega()/abs(trkStates[iTS].getOmega())
			charge2 = trkStates[jTS].getOmega()/abs(trkStates[jTS].getOmega())

			if dist < minDist and charge1*charge2 < 0:
				minDist = dist
				goodTS1 = trkStates[iTS]
				goodTS2 = trkStates[jTS]
	# print 'minDist: ', minDist
	# print ''

	return goodTS1, goodTS2

def getHelixXYC(trackState):

	d0 = trackState.getD0()
	phi0 = trackState.getPhi()
	referencePoint = trackState.getReferencePoint()
	omega = trackState.getOmega()
	charge = omega / abs(omega)
	radius = 1. / abs(omega)

	xAtPCA = -d0 * math.sin(phi0) + referencePoint[0]
	yAtPCA = d0 * math.cos(phi0) + referencePoint[1]

	xCentre = xAtPCA + radius * math.cos(phi0 - 0.5*math.pi*charge)
	yCentre = yAtPCA + radius * math.sin(phi0 - 0.5*math.pi*charge)

	return xCentre, yCentre

def rotateAngle(angle, rotation):
	rotAngle = angle + rotation
	if rotAngle < 0:
		rotAngle += 2*math.pi
	return rotAngle
