import math
import numpy as np

from pyLCIO import IOIMPL, ROOT
import dd4hep

'''
description = dd4hep.Detector.getInstance()
bfieldV = [0.,0.,0.]
#print bfieldV
description.field().magneticField( [ 0., 0., 0. ]  , bfieldV  )
bField = bfieldV[2]/dd4hep.tesla
'''

bField = 3.5 # HARDCODED!!!
FCT = 2.99792458E-4

def R(x,y):
	return math.sqrt(x**2 + y**2)

def R3(x,y,z):
	return math.sqrt(x**2 + y**2 + z**2)

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

def angularDistance(theta1, theta2, phi1, phi2):

	dPhi = phi1 - phi2
	if abs(dPhi) > math.pi:
		dPhi = 2 * math.pi - abs(dPhi)
	'''
	theta1 -= math.pi / 2
	theta2 -= math.pi / 2
	phi1   += math.pi
	'''

	dist = math.sqrt( (theta1-theta2)**2 * (math.sin(theta1-theta2))**2 + (dPhi)**2 )
	#dist = math.acos( \
	#	   math.cos(theta1) * math.cos(phi1) * math.cos(theta2) * math.cos(phi2) \
	# 	 + math.cos(theta1) * math.sin(phi1) * math.cos(theta2) * math.sin(phi2) \
	#	 + math.sin(theta1) * math.sin(theta2) \
	#	 )

	return dist
