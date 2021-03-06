import numpy as np
from soundwave_medium_params import mediumType
from initcond_soundwave import getWFunction, fromCurvesToVals
from soundwave_boundary_conditions import getPeriodicXArray

def getVals(z, curves):
	from soundwave_medium_params import p00, v00
	if mediumType == "homog":
		from soundwave_medium_params import rho00
		rhoIni = rho00
	elif mediumType == "inhomog":
		from soundwave_medium_params import densFunc
		rhoIni = densFunc(z)
	presCurve = curves['pres']
	rhoCurve = curves['rho']
	velCurve = curves['vel']
	return fromCurvesToVals(p00, rhoIni, v00, presCurve, rhoCurve, velCurve)


if mediumType == "homog":
	def getWFunctionVals(functiontype, csSign, z, t):
		from soundwave_medium_params import cs00		
		newZ = getPeriodicXArray(z - csSign * cs00 * t)
		return getWFunction(functiontype)(newZ)


elif mediumType == "inhomog":
	from common import getZArray	
	newZ = getZArray()

	print("creating newZ for analytic sol and  " )

	
	#method = 1  #NEVER USE!
	#method = 2
	method = 3
	#method = 4  #NEVER USE!

	#useC = True
	useC = False #is weave faster than numpy operations?

	from constants import loopType
	if useC and loopType == "cython":
		import pyximport
		pyximport.install(setup_args={'include_dirs':[np.get_include()]})

	def getWFunctionVals(functiontype, csSign, z, time):
		if(functiontype != "wavepacket"):
			print("analytical method not implemented")
			import sys
			sys.exit(0)
		from constants import z0, zf
		from soundwave_medium_params import cs00
		from sound_wave_packet_params import getSoundWaveGaussFunction, getSoundWaveFunction, zc, W, k0
		from math import pi

		if method == 1:
			curve = np.zeros(z.shape)
			def getZ0Index(zval):
				from common import getDz
				delta = getDz()
				for index in range(newZ.shape[0]):
					if (abs(newZ[index] - zval)<delta):
						return index
			indexZ0 = np.zeros(newZ.shape)		
			for index in range(newZ.shape[0]):
				indexZ0[index] = getZ0Index(z[index])
				#print(indexZ0[index])
			nans, x= np.isnan(indexZ0), lambda z: z.nonzero()[0]
			indexZ0[nans]= np.interp(x(nans), x(~nans), indexZ0[~nans])
			indexZ0 = indexZ0.astype(int)
			firstZ = z[indexZ0]
			csZ0 = cs00(firstZ)
			csZt = cs00(z)
			#ampIni = getSoundWaveGaussFunction(zc, W)(firstZ)
			ampIni = getSoundWaveFunction(k0, zc, W)(firstZ)
			#print("AMP INI")
			#print(ampIni)
			for index in range(firstZ.shape[0]):
				#curve[index] =  ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) * np.cos( 2 * pi * k0 * (csZ0[index]/csZt[index]) / (zf - z0) * ( z[index] - csZt[index] * time) - 2 * pi * k0 * z0 / (zf - z0))
				curve[index] =  ampIni[index] * (csZ0[index] / csZt[index]) ** (0.5)
	
		elif method == 2:
			from common import getZIndex
			curve = np.full(z.shape, np.nan)
			csZ0 = cs00(z)
			csZt = cs00(newZ)
			#ampIni = getSoundWaveGaussFunction(zc, W)(z)
			ampIni = getSoundWaveFunction(k0,zc, W)(z)
			from constants import nint
			if  useC and (loopType == "weave" or loopType == "cython"):
				if loopType == "weave":
					from scipy.weave import inline, converters
					from constants import z0, zf
					code = """
					for(int index = 0;index<nint+2; index++) {
						int zIndexNewZ =  int(float(nint)*(newZ(index) - z0)/(zf - z0) + 0.5 );
						curve(zIndexNewZ) = ampIni(index) * pow((csZ0(index) / csZt(index)), 0.5);
					}
					"""
					inline(code, ['curve', 'ampIni', 'csZ0', 'csZt', 'newZ', 'nint', 'z0', 'zf'],type_converters=converters.blitz)
				elif loopType == "cython":
					from cython_as import as2
					as2(curve, ampIni, csZ0, csZt, newZ, nint, z0, zf)


					
			else:
				for index in range(z.shape[0]):
				#curve[getZIndex(newZ[index])] = ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) *  np.cos( 2 * pi * k0 * (csZ0[index]/csZt[index]) / (zf - z0) * ( newZ[index] - csZt[index] * time) - 2 * pi * k0 * z0 / (zf - z0))
					curve[getZIndex(newZ[index])] = ampIni[index] * (csZ0[index]  / csZt[index] )** (0.5)

			#print("curve before interp")
			#print(curve)
			nans, x= np.isnan(curve), lambda z: z.nonzero()[0]
			curve[nans]= np.interp(x(nans), x(~nans), curve[~nans])

		elif method == 3:
			from common import getZIndex
			curve = np.full(z.shape, np.nan)
			csZ0 = cs00(z)
			csZt = cs00(newZ)
			#ampIni = getSoundWaveGaussFunction(zc, W)(z)
			ampIni = getSoundWaveFunction(k0, zc, W)(z)
			from constants import loopType,nint
			if  useC and (loopType == "weave" or loopType == "cython"):
				if loopType == "weave":
					from scipy.weave import inline, converters
					code = """
					for(int index = 0;index<nint+2; index++) {
						curve(index) = ampIni(index) * pow((csZ0(index) / csZt(index)), 0.5);
					}
					"""
					inline(code, ['curve', 'ampIni', 'csZ0', 'csZt', 'nint'],type_converters=converters.blitz)
				elif loopType == "cython":
					from cython_as import as3
					as3(curve, ampIni, csZ0, csZt, nint)
			else:
#				for index in range(z.shape[0]):
#				#curve[index] = ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) *  np.cos( 2 * pi * k0 * (csZ0[index]/csZt[index]) / (zf - z0) * ( newZ[index] - csZt[index] * time) - 2 * pi * k0 * z0 / (zf - z0))
#					curve[index] = ampIni[index] * (csZ0[index] / csZt[index]) ** (0.5)
				#numpy array operations
				curve = ampIni * (csZ0 / csZt) ** (0.5)
				 
		elif method == 4:
			from common import getZIndex
			curve = np.full(z.shape, np.nan)
			csZ0 = cs00(z)
			csZt = cs00(newZ)
			ampIni = getSoundWaveGaussFunction(zc, W)(z)
			#ampIni = getSoundWaveFunction(k0,zc, W)(z)
			for index in range(z.shape[0]):
				curve[getZIndex(newZ[index])] = ampIni[index] * (csZt[index] / csZ0[index] )** (0.5) *  np.cos( 2 * pi * k0 * (csZ0[index]/csZt[index]) / (zf - z0) * ( newZ[index] - csZt[index] * time) - 2 * pi * k0 * z0 / (zf - z0))
				#curve[getZIndex(newZ[index])] = ampIni[index] * (csZ0[index]  * csZt[index] )** (0.5)
			#print("curve before interp")
			#print(curve)
			nans, x= np.isnan(curve), lambda z: z.nonzero()[0]
			curve[nans]= np.interp(x(nans), x(~nans), curve[~nans])

		#print("curve")

		#print(curve)
		#print("MAX CURVE = %e , IND MAX CURVE =%d, min CURVE = %e , ind MIN curve = %d" % (np.max(curve), np.argmax(curve), np.min(curve), np.argmin(curve) ) )
		
		

		return curve
	

	def updateNewZ(modelObj, dt, csNumerical):
		for index in range(newZ.shape[0]):
			newZ[index] = modelObj.getNewPoint(newZ[index],dt, csNumerical)





def getCurves(z, t):
	"""
		calculates values of generic perturbations A * h(z)
		if there is a superposition pres and rho perturbation will be different from velocity perturbation	
	"""
	from soundwave_perturbation_params import perturbationType
	if perturbationType == "superposition":
		from soundwave_perturbation_params import  init_functions_generation
		velPert = np.zeros(z.shape)
		presRhoPert = np.zeros(z.shape)
		for fdef in init_functions_generation:
			wFuncVals =  getWFunctionVals(fdef['functiontype'],fdef['csSign'] , z, t)
			velPert += fdef['csSign'] * fdef['A'] * wFuncVals
			presRhoPert += fdef['A'] * wFuncVals
		return {'pres':presRhoPert, 'rho':presRhoPert, 'vel': velPert}
	elif perturbationType == "one":
		from soundwave_perturbation_params import  A, functiontype
		velPresRhoPert = A * getWFunctionVals(functiontype ,1 , z, t)    
		return {'pres':velPresRhoPert, 'rho':velPresRhoPert, 'vel': velPresRhoPert}


	


