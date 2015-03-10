"""
	The parameters for gaussian sound wave packet
	k0 - wavenumber
	zc - center of the gauss function
	W width of the gauss function

"""

import numpy as np
from constants import z0, zf
from math import pi,sqrt

k0 = 60.0
zc = z0 + 0.2 * (zf - z0)
W = 0.05

from soundwave_medium_params import mediumType
if mediumType == "inhomog":
	from soundwave_medium_params import inhomogSubtype
	if inhomogSubtype == 2:
		#change 	
		k0 = 15.0#second exp of inhom
		#zc = z0 + 0.2 * (zf - z0) homog
		#zc = z0 + 0.3 * (zf - z0)
		zc = z0 + (3.0/20.0)*(zf - z0)#second exp of inhom and first new
		#W = 0.01
		W = 0.25 #second exp of inhom

def getSoundWaveGaussFunction(zc, W):
	"""
	returns the gauss function (of the envelope) - this is a gaussian wave packet!	
		Parameters
		----------
		zc : real
			center of the gauss function    
		W : real
			width of the gauss function
	returns g so that
	g(z) = -(z-zc)**2/W**2

	"""
	def gaussFunction(z):
		t2 = np.subtract(z,zc) ** 2
		return np.exp(-np.divide(t2, W**2))	 
	return gaussFunction


def getSoundWaveFunction(k0, zc, W):
	"""
	returns the resulting wave function obtained by multiplying the gauss function with a cos function:
		Parameters
		----------
		k0 : int  
			wavenumber
		zc : real
			center of the gauss function    
		W : real
			width of the gauss function
	returns f so that:
	f(z) =  g(z) / cos(2 pi k0 (z-z0)/(zf-z0) )
	with g the gauss function defined above

	"""
	def gaussPacketFunction(z):
		return np.multiply(getSoundWaveGaussFunction(zc, W)(z),  np.cos(2.0 * pi * k0 * (z-z0)/ (zf - z0) ) )
	return gaussPacketFunction




def getSoundWaveFFTAnalytical(k0, zc, W):
	"""
		Parameters
		----------
		k0 : int  
			wavenumber
		zc : real
			center of the gauss function    
		W : real
			width of the gauss function
	
		returns analytical expression of the sound wave function defined above: [-(z-zc)**2/W**2] *  cos(2 pi k0 (z-z0)/(zf-z0) )
		calculated with mathematica 
				
	
	"""
	from constants import z0, zf
	def analyticFFT(k):
		return ((np.exp(-((pi*(k0**2*pi*W**2 - 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)) +  np.exp(-((pi*(k0**2*pi*W**2 + 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)))*sqrt(pi)*W)/2
	return analyticFFT
