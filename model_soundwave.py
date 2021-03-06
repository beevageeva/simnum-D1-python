import numpy as np
import sys
from math import pi
from constants import gamma
from base_model import BaseModel
from notifier_params import plotPresCurve, plotVelCurve, plotRhoCurve, markPoints, plotPresAn, plotRhoAn, plotVelAn, plotVelFFT, plotVelFFTAnal,  plotPresFFT
from soundwave_medium_params import mediumType, cs00

"""
Parameters:
	addMarkPoints a hash with the points I want to follow: at each iteration I call getNewPoint to get their new value

	calcWidth if True it calculates the width of the packet as it is evolving : it calculates it by difference of
		self.addMarkPoint["right"] - self.addMarkPoint["left"] so I have to have points labelled with "left" 
		and "right" in addMarkPoints
		tries to get a relationship between fourier width and packet width (it will append  time maxAmp/widthPacket to ww.txt)
		it also appends time widthPachet/variation length of density to 	wdcl.txt
		

	calcKc calculates kc and appends time kc * cs(self.addMarkPoints[0]) to kc.txt
	(I refer here to self.addMarkPoints[0] as the point which is first in the hash - and calculated in each step as it is displaced with cs: 
	I put zc)
	plotWidthOnGraph plots the self.addMarkPoint["left"] and self.addMarkPoint["right"] on the graph

	showErr = True it will append at each iteration time error to preserr.txt, rhoerr.txt or velerr.txt

	all the files *.txt generated here can be plot with 
		python plotValuesFromText.py file.txt 
 	

"""

showErr = False
#calcWidth = True
calcWidth = True
#plotWidthOnGraph = False
plotWidthOnGraph = True
calcKc = True
#calcKc = False
from soundwave_perturbation_params import perturbationType
if(perturbationType != "one"):
	calcKc = False
else:		
	from soundwave_perturbation_params import functiontype
	if(functiontype != "wavepacket" or mediumType != "inhomog"):
		calcKc = False


#addMarkPoint = None
#uncomment this to add a new mark point
#the following for the wave packet
from sound_wave_packet_params import zc,W
addMarkPoints = {'zc':zc, 'left':zc-1.666*W, 'right':zc+1.666*W} 


plotCsMaxMin = True
#plotCsMaxMin = False




class Model(BaseModel):
	

	def __init__(self):
		BaseModel.__init__(self)
		if addMarkPoints:
			self.addMarkPoints = addMarkPoints


	def getNewPoint(self, zval, dt, csNumerical = True):
		"""
			calculates the new z point 
			Parameters:
				zval - the  value to calculate
				dt 
				csNumerical if True it will calculate cs from current values of p, rho and v, otherwise it will take cs from soundwave_medium_params
				if the perturbation is a superposition csNumerical=False won't work, but otherwise 
				taking it from soundwave_medium_params may be faster
		"""
		from common import getZIndex
		from math import sqrt
		from soundwave_medium_params import v00, p00
		from soundwave_boundary_conditions import getPeriodicX
		if(csNumerical):
			zIndex = getZIndex(zval)	
			#from sound_wave_params import csSign 
			#I should not import from here: this should be used for generating initial conditions only
			# I have to calculate it from actual values
			#Imagine that it should work for a superposition of wave travelling right and left
			if (self.pres[zIndex]< p00 and self.vel[zIndex] > v00 ) or (self.pres[zIndex]> p00 and self.vel[zIndex] < v00 ):
				csSign = -1
			else:
				csSign = 1
			cs = csSign * sqrt(gamma * self.pres[zIndex] / self.rho[zIndex])
		else:
			from soundwave_medium_params import cs00
			if(mediumType == "homog"):
				cs = cs00
			else:
				cs = cs00(zval)	
		v = v00 +  cs	
		newz = getPeriodicX(zval + v * dt)
		return newz


	def updateValuesModel(self, dt, time):
		"""
			called in each iteration of mainLoop methos of parent class base_model.py
		"""
		#markPoints used to plot maximum and minumum for rho, pres and vel , in hom medium they travel at cs speed (phase vel), but in inhom at group speed
		if(markPoints):
			self.maxPresZ = self.getNewPoint(self.maxPresZ,dt)
		if(calcKc):
			self.presFFT = self.getPresFFTVals(False)
			F = self.presFFT[1]
			Y = self.presFFT[0]
			indCenter = np.argmax(Y[1:])+1
			kc = abs(F[indCenter])
			from common import testKeyInDict
			if(calcWidth and addMarkPoints is not None and testKeyInDict("left",addMarkPoints) and testKeyInDict("right",addMarkPoints)):
				#TODO test wavepcket
				import os	
				maxAmp = Y[indCenter]
#				from analyze_functions import getFirstIndexDifferentLeft,getFirstIndexDifferentRight, getGaussianLimitRight
#				delta = 0.005
#				rightIndex = getGaussianLimitRight(Y, indCenter, delta)
#				#leftIndex = 2 * indCenter - rightIndex
#				leftIndex = 0
#				widthFourier =  F[rightIndex] - F[leftIndex]
#				delta = 0.00005
#				widthFunction =  self.z[getFirstIndexDifferentRight(self.pres, delta)] -  self.z[getFirstIndexDifferentLeft(self.pres, delta)]
#				#print("widthFunction = %e == %e?, widthFourier = %e, w1*w2 = %e, w1*w2gr=%e" % (widthFunction, widthFunc2, widthFourier,widthFourier * widthFunc2, widthFourier * widthFunction ))
#				os.system("echo %e %e >> %s" % (time, widthFourier * widthFunction , "ww.txt"))				
#				os.system("echo %e %e >> %s" % (time, maxAmp/widthFunction, "adW.txt"))			
				widthFunction = self.addMarkPoints["right"] - self.addMarkPoints["left"]	
				#VARLENGTH
				os.system("echo %e %e >> %s" % (time, maxAmp/widthFunction, "ww.txt"))				
				from soundwave_medium_params import getDensVarLength	
				cl =	getDensVarLength(self.z)
				os.system("echo %e %e >> %s" % (time, widthFunction/cl , "wdcl.txt"))				
				#print("1/kc= %e, wp = %e, cl = %e" % (1.0/kc, widthFunction, cl))

		if hasattr(self, "addMarkPoints"):
			for k in self.addMarkPoints:
				self.addMarkPoints[k] = self.getNewPoint(self.addMarkPoints[k],dt)
		from soundwave_perturbation_params import perturbationType
		if(perturbationType == "one"):
			from soundwave_perturbation_params import functiontype
			if(functiontype == "wavepacket" and mediumType == "inhomog"):
				if(plotPresAn or plotVelAn or plotRhoAn):
					from analytic_solution import updateNewZ
					updateNewZ(self, dt, False)
				if(calcKc and not addMarkPoints is None):
					localMarkPoint = self.addMarkPoints[list(addMarkPoints.keys())[0]] 
					from sound_wave_packet_params import k0, getSoundWaveGaussFunction, zc, W
					from common import getDz, getZIndex
					from constants import z0, zf
					from initcond_soundwave import fromValsToCurvePres
					from soundwave_perturbation_params import A	
					from soundwave_medium_params import cs00
					#print("dz: %e, dt*minCs=%e" % (getDz(), dt * np.min(cs00(self.z))))
					cszp0 = cs00(localMarkPoint)
					cszp = 	cs00(localMarkPoint)
					#kt = k0 *  cszp0/cszp
					#from soundwave_medium_params import csderAnal
					#ktnew =  k0 * np.exp(-csderAnal(localMarkPoint) * time)
					#print("kc=%e,kt=%e,ktnew=%e,kc/kt = %e, kc*cs(zp(t))=%e,ktnew*cszp=%e,k0*csz0=%e" % (kc, kt,ktnew,kc / kt, kc * cszp, ktnew*cszp, k0*cszp0))
					import os
					os.system("echo %e %e >> %s" % (time, kc * cszp, "kc.txt"))				

						
	
	
#					gaussFunc = getSoundWaveGaussFunction( zc, W)
#					#gaussFunc = getSoundWaveFunction(k0, zc, W)
#					ampzp0 = A * gaussFunc(addMarkPoint)
#					#phaseZp =  2 * pi * kc / (zf - z0) * (self.addMarkPoint - cszp * time) - 2 * pi * k0 * z0 / (zf - z0) 
#					phaseZp =  2 * pi * k0 / (zf - z0) * cszp0 * (self.addMarkPoint/cszp - time) - 2 * pi * k0 * z0 / (zf - z0) 
#					#ampzp = fromValsToCurvePres(self.pres[getZIndex(self.addMarkPoint)]) / np.cos(phaseZp) 
#					ampzp = fromValsToCurvePres(self.pres[getZIndex(self.addMarkPoint)])  
#					#print("ampz0**2*csz0 = %e == %e = ampzp**2*cszp" % (ampzp0**2*cszp0,  ampzp**2*cszp))
#					print("ampz0*csz0 = %e == %e = ampzp*cszp" % (ampzp0*cszp0,  ampzp*cszp))
#					os.system("echo %e %e >> %s" % (time, ampzp**2*cszp, "amp.txt"))					
#					os.system("echo %e %e >> %s" % (time, kc *self.addMarkPoint, "kz.txt"))					
				


	def updateValuesNotifier(self, dt, time):
		"""
			called in every nStepsPlot(defined in notifier_params.py) iterations of mainLoop methos of parent class base_model.py
		"""
		#TODO simpl
		if(plotPresCurve):
			from initcond_soundwave import fromValsToCurvePres
		if(plotVelCurve):
			from initcond_soundwave import fromValsToCurveVel
		if(plotRhoCurve):
			from initcond_soundwave import fromValsToCurveRho
		#TODO
		newZNA = None
		if(plotPresAn or plotVelAn or plotRhoAn):
			from analytic_solution import getCurves, getVals
			curves = getCurves(self.z, time)
			anValuesCalc = False
			if(mediumType == "inhomog"):
				from analytic_solution import method
				if method==3:
					from analytic_solution import newZ
					newZNA = ((self.z, False),(newZ, True))
					anVals = getVals(newZ, curves)
					anValuesCalc = True
			if not anValuesCalc:
				anVals = getVals(self.z, curves)
		
		if(plotPresAn):
			presc = curves['pres']
			anPres =  anVals['pres']
			#print("in update vgals notif max pres an ind = %d  max = %e" % (np.argmax(anPres),np.max(anPres) ))
			presNewVals = [self.pres, anPres]
			if(showErr):
				err = np.max(np.absolute(np.subtract(self.pres, anPres)))
				#print("time=%E,max abs (self.pres - anPres) = %E" % (time,err))
				os.system("echo %e %e >> %s" % (time, err , "preserr.txt"))				
			if(plotPresCurve):
				presCurveNewVals = [fromValsToCurvePres(self.pres), presc]
		else:
			presNewVals = self.pres
			if(plotPresCurve):
				presCurveNewVals = fromValsToCurvePres(self.pres)
		if(plotRhoAn):
			rhoc = curves["rho"]
			anRho =  anVals["rho"]
			rhoNewVals = [self.rho, anRho]
			if(showErr):
				err = np.max(np.absolute(np.subtract(self.rho, anRho)))
				os.system("echo %e %e >> %s" % (time, err , "rhoerr.txt"))				
				#print("max abs (self.rho - rhoAn) = %E" % err)
			if(plotRhoCurve):
				rhoCurveNewVals = [fromValsToCurveRho(self.rho, self.z), rhoc]
		else:
			rhoNewVals = self.rho
			if(plotRhoCurve):
				rhoCurveNewVals = fromValsToCurveRho(self.rho, self.z)
		if(plotVelAn):
			velc = curves["vel"]
			anVel =  anVals["vel"]
			velNewVals = [self.vel, anVel]
			if(showErr):
				err = np.max(np.absolute(np.subtract(self.vel, anVel)))
				#print("max abs (self.vel - velAn) = %E" % err)
				os.system("echo %e %e >> %s" % (time, err , "velerr.txt"))				
			if(plotVelCurve):
				velCurveNewVals = [fromValsToCurveVel(self.vel), velc]
		else:
			velNewVals = self.vel
			if(plotVelCurve):
				velCurveNewVals = fromValsToCurveVel(self.vel)

		

		self.notifier.updateValues("rho", rhoNewVals, time, newZNA)
		self.notifier.updateValues("pres", presNewVals, time, newZNA)
		self.notifier.updateValues("vel", velNewVals,  time, newZNA)
		if(plotPresCurve):
			self.notifier.updateValues("presCurve", presCurveNewVals, time, newZNA)
		if(plotVelCurve):
			self.notifier.updateValues("velCurve", velCurveNewVals,  time, newZNA)
		if(plotRhoCurve):
			self.notifier.updateValues("rhoCurve", rhoCurveNewVals, time, newZNA)
		if(plotVelFFT):
			self.notifier.updateValues("velFFT", self.getVelFFTVals(True)[0:-1], time)
		if(plotPresFFT):
			#if ('presFFT' in vars()):
			if (hasattr(self, 'presFFT')):
				presFFT = self.presFFT
			else:
				print("Calculate presFFT ")
				presFFT = self.getPresFFTVals(False)
			self.notifier.updateValues("presFFT", presFFT[0], time)
#			if(calcWidth and plotWidthOnGraph):
#				Y = presFFT[0]
#				F = presFFT[1]
#				indCenter = np.argmax(Y[1:])+1
#				delta = 0.005
#				from analyze_functions import getFirstIndexDifferentLeft,getFirstIndexDifferentRight, getGaussianLimitRight
#				rightIndex = getGaussianLimitRight(Y, indCenter, delta)
#				#leftIndex = 2 * indCenter - rightIndex
#				leftIndex = 0
#				widthFourier =  F[rightIndex] - F[leftIndex]
#				self.notifier.markPoint("presFFT", "presFFTLeft", F[leftIndex])
#				self.notifier.markPoint("presFFT", "presFFTRight", F[rightIndex])

		if(markPoints):
			#only for pres	
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
		if hasattr(self, "addMarkPoints"):
			for k in self.addMarkPoints:
				self.notifier.markPoint("pres", "mark_%s" % k, self.addMarkPoints[k])
#		if(calcWidth and plotWidthOnGraph):
#			from analyze_functions import getFirstIndexDifferentLeft, getFirstIndexDifferentRight
#			delta = 0.00005
#			self.notifier.markPoint("pres", "presLeft", self.z[getFirstIndexDifferentLeft(self.pres, delta)])
#			self.notifier.markPoint("pres", "presRight", self.z[getFirstIndexDifferentRight(self.pres, delta)]) 
#			#print("LZ = %e , RZ = %e" % (self.z[getFirstIndexDifferentLeft(self.pres, delta)], self.z[getFirstIndexDifferentRight(self.pres, delta)]))		

	def getInitialValues(self):
		"""
			the initial values of pres, rho and vel as an array (used in order to plot)
		"""	
		if plotPresAn:
			presIniVal = [self.pres, self.pres]
		else:
			presIniVal = self.pres
		if plotRhoAn:
			rhoIniVal = [self.rho, self.rho]
		else:
			rhoIniVal = self.rho
		if plotVelAn:
			velIniVal = [self.vel, self.vel]
		else:
			velIniVal = self.vel
		#return [[self.pres, self.pres], [self.rho, self.rho], [self.vel, self.vel]]
		#If I don't want analitical solution plotted for velocity:
		#return  [[self.pres, self.pres], [self.rho, self.rho], self.vel]
		return  [presIniVal, rhoIniVal, velIniVal]

		
	def getVelFFTVals(self, middlePoints):
		"""
			calculates f1(k) as defined in pdf of velocity
			Parameters:
				middlePoints = True it will use values (v(i) + v(i+1))/2 instead of original values
		"""	
		if(middlePoints):
			vals = (self.vel[1:] + self.vel[:-1]) / 2
		else:
			vals = self.vel	
		from constants import z0, zf
		intlen = zf - z0
		from scipy.fftpack import fft,fftfreq
		Y=fft(vals)/len(vals)
		F=fftfreq(len(vals), self.z[1] - self.z[0]) 
		#print("max freq in plotVelFFT %e" % intlen * F[np.argmax(Y[1:]) + 1 ])
		if(plotVelFFTAnal):
			from initcond_soundwave import getVelFFTAn
			#when using first form
			#vals = [intlen * abs(Y), abs(getVelFFTAn(F)), F * intlen]
			vals = [intlen * (1.0 / np.sqrt(2 * pi)) * abs(Y), abs(getVelFFTAn(F*(-2 * pi))), F * 2 * pi]
		else:
			#when using first form
			#vals = [intlen * abs(Y), F * intlen]	
			vals = [intlen * (1.0 / np.sqrt(2 * pi)) *abs(Y), F * 2 * pi]	
		return vals

	def getPresFFTVals(self, middlePoints):
		"""
			calculates f1(k) as defined in pdf of pres
			Parameters:
				middlePoints = True it will use values (p(i) + p(i+1))/2 instead of original values
		"""	
		if(middlePoints):
			vals = (self.pres[1:] + self.pres[:-1]) / 2
		else:
			vals = self.pres	
		from constants import z0, zf
		intlen = zf - z0
		from scipy.fftpack import fft,fftfreq
		Y=fft(vals)/len(vals)
		F=fftfreq(len(vals), self.z[1] - self.z[0])
		#print("getPresFFT")
		#print(abs(Y)) 
		#print("max freq in plotPresFFT %e " % intlen * F[np.argmax(Y[1:]) + 1 ] )
		#when using the first form 
		#return  [intlen * abs(Y), F * intlen]	
		return [intlen * (1.0 / np.sqrt(2 * pi)) *abs(Y), F * 2 * pi]	
	
	

	def additionalInit(self):
		"""
			called in constructor of parent class base_model.py
		"""
		#TODO all initial values of curves are calcutaed 2 times: 2 function calls for each
		#plot Curves of pression , vel, density
		if(plotPresCurve):
			from initcond_soundwave import fromValsToCurvePres
			prescvals = fromValsToCurvePres(self.pres)
			self.notifier.addGraph(self.z, "presCurve",[prescvals, prescvals] if plotPresAn else prescvals)
		if(plotVelCurve):
			from initcond_soundwave import fromValsToCurveVel
			velcvals = fromValsToCurveVel(self.vel, self.z)
			self.notifier.addGraph(self.z, "velCurve", [velcvals, velcvals] if plotVelAn else velcvals )
		if(plotRhoCurve):
			from initcond_soundwave import fromValsToCurveRho
			rcvals = fromValsToCurveRho(self.rho, self.z)
			self.notifier.addGraph(self.z, "rhoCurve", [rcvals,rcvals] if plotRhoAn else rcvals)
		if(plotVelFFT):
			vals = self.getVelFFTVals(True)
			self.notifier.addGraph(vals[-1], "velFFT", vals[0:-1] , "k" , 'None')
		if(plotPresFFT):
			#use the same boolean parameter as in update	
			vals = self.getPresFFTVals(False)
			self.notifier.addGraph(vals[-1], "presFFT", vals[0] , "k", 'None')
	
		if(markPoints):
			self.maxPresZ = self.z[np.argmax(self.pres)]	
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
		if(addMarkPoints):
			for k in addMarkPoints:
				self.notifier.markPoint("pres", "mark_%s" % k, addMarkPoints[k])
		if(mediumType=="inhomog"):
			from soundwave_medium_params import cs00
			self.notifier.plotAxisTwin(self.z, "vel",cs00(self.z) , "cs00")
		if(plotCsMaxMin):
			if(mediumType=="inhomog"):

				maxZIndex = np.argmax(self.pres)
				minZIndex = np.argmin(self.pres)

				from soundwave_medium_params import cs00 as cs00Func
				from soundwave_medium_params import p00, densFunc

				cs00 = cs00Func(self.z)
				#for alpha in np.arange(-2.5, 2.5, 0.5):
				for alpha in [-0.5]:
					vals = np.power(cs00, alpha)
					self.notifier.plotAxis(self.z,"pres",  p00 + ((self.pres[minZIndex] - p00)/ vals[minZIndex]) * vals, "(min)%1.1f" % alpha)
					self.notifier.plotAxis(self.z,"pres",  p00 + ((self.pres[maxZIndex] - p00)/ vals[maxZIndex]) * vals, "(max)%1.1f" % alpha)
				for alpha in [ 0.5]:
					vals = np.power(cs00, alpha)
					self.notifier.plotAxis(self.z,"vel", np.multiply(self.vel[minZIndex] / vals[minZIndex],vals), "(min)%1.1f" % alpha)
					self.notifier.plotAxis(self.z,"vel", np.multiply(self.vel[maxZIndex] / vals[maxZIndex],vals), "(max)%1.1f" % alpha)
					if(plotRhoCurve):
						for alpha in [-0.5]:
							vals = np.power(cs00, alpha)
							self.notifier.plotAxis(self.z,"rhoCurve", ((self.rho[minZIndex] - densFunc(self.z[minZIndex])) / (vals[minZIndex] * densFunc(self.z[minZIndex])))* vals , "(min)%1.1f" % alpha)
							self.notifier.plotAxis(self.z,"rhoCurve", ((self.rho[maxZIndex] - densFunc(self.z[maxZIndex])) / (vals[maxZIndex] * densFunc(self.z[maxZIndex])))* vals , "(max)%1.1f" % alpha)

			else:
				print("csmin max only inhom")
		



