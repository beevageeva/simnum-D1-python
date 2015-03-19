from constants import problemType

#useLimits = False
useLimits = True


if problemType == "riemann" or not useLimits :

	def getXLimits(title):
		return None

	def getYLimits(title):
		return None

else:

	def getXLimits(title):
		"""
			get axis x limits (if I want fixed)
			if returns None by default the autoscale is done 
		"""
		if(title == "velFFT" or title == "presFFT"):
			return {"minX":-80, "maxX":80}
		from constants import z0, zf
		return {"minX" : z0, "maxX" : zf}
	
	def getYLimits(title):
		"""
			get axis y limits (if I want fixed)
			if returns None by default the autoscale is done 
		"""
		from soundwave_medium_params import mediumType
		if mediumType == "homog":
			medType = "homog"
		else:	
			from soundwave_medium_params import inhomogSubtype
			if inhomogSubtype == 1:
				medType = "inhomog1"
			else:
				medType = "inhomog2"
		if(title == "pres"):
			if medType == "homog":
				return { "maxY": 1.0006, "minY": 0.9994} #homog
			elif medType == "inhomog1":
				return { "maxY": 1.0006, "minY": 0.9995} #inhomog1
			elif medType == "inhomog2":
				return {  "maxY": 1.002, "minY": 0.998} #inhomog2
		elif(title == "vel"):
			if medType == "homog":
				return	{ "maxY": 0.0004, "minY": -0.00035} 
			elif medType == "inhomog1":
				return	{ "maxY": 0.0015, "minY": -0.0015} 
			elif medType == "inhomog2":
				return	 { "maxY": 0.0035, "minY": -0.004}
		elif(title == "rho"):
			if medType == "homog":
				return { "maxY": 1.0004, "minY": 0.9996} 
			elif medType == "inhomog1":
				return { "maxY": 1.0004, "minY": 0}	 
			elif medType == "inhomog2":
				return { "maxY": 1.3, "minY": 0} 	 
		elif(title == "rhoCurve"):
			if medType == "homog":
				return { "maxY": 0.00025, "minY": -0.00025}		
			elif medType == "inhomog1":
				return	{ "maxY": 0.00035, "minY": -0.00035}	
			elif medType == "inhomog2":
				return	{ "maxY": 0.0020, "minY": -0.0020}
		#I should always relim y for pres fft beacuse of the central value in k = 0 (= 1 ,=  p00, = mean value of p)
		#there is no need for vel fft beacuse mean vel = 0
		elif(title == "presFFT"):
			if medType == "inhomog1":
				return	{ "maxY": 0.00015, "minY": 0}
			if medType == "inhomog2":
				return	{ "maxY": 0.00015, "minY": 0}
			elif medType == "homog":
				return	{ "maxY": 0.000025, "minY": 0}
	
		return None
