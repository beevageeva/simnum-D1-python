import numpy as np
import math
from constants import gamma

rho00 = 1.0
#rho00 = 0.3  #second exp of inhom

#mediumType = "homog"
mediumType = "inhomog"  #variable density rho00 to test with wave packet
if(mediumType=="inhomog"):
	from constants import z0, zf
	rho01 = 0.01
	#rho01 = 1.2#second exp of inhom
	#ze = 0.5*(z0 + zf)
	ze = z0 + 0.2*(zf - z0) #first exp of inhom new
	#ze = z0 + 0.7*(zf - z0) #second exp of inhom
	#we = 0.4
	we = 0.2
	densFunc = lambda z: rho00 + 0.5 * (rho01-rho00) * (1 + np.tanh((z-ze)/we))



#functiontype = 'sine'
functiontype = 'defined'
#periodicType = "repeat" 
periodicType = "refl" 

p00 = 1.0


v00 = 0.0

#cs00 = math.sqrt(gamma * p00 / rho00)
#v00 = - cs00 /  5.5
#v00 =  0.5 * cs00
#v00 = - cs00
#v00 =  cs00

A = 3.0 * 10.0 ** (-4)
#A = 5.0 * 10.0 ** (-2)
#A = 10.0 ** (-2) #this will work fine as well with fg scheme

#init_functions_generation = [{'csSign':1, 'A': A}, {'csSign':-1, 'A': 0.5*A}] #SUPERPOSITION wave travelling right with amp A and left with amp 0.5 * A
#init_functions_generation = [{'csSign':-1, 'A': A}, {'csSign':1, 'A': 0.5*A}] #SUPERPOSITION wave travelling left with amp A and right with amp 0.5 * A
#init_functions_generation = [{'csSign':-1, 'A': A}, {'csSign':1, 'A': A}] #SUPERPOSITION wave travelling left with amp A and right with amp  A
#init_functions_generation = [{'csSign':-1, 'A': A}] #travelling left
init_functions_generation = [{'csSign':1, 'A': A}] #travelling right



#The following are used for the notifier(taken from notifier_params as they were wavesound specific)
plotPresCurve = False
#plotRhoCurve = False
plotRhoCurve = True
plotVelCurve = False
#plotPresAn = False
#plotRhoAn = False
#plotVelAn = False
plotPresAn = True
plotRhoAn = True
plotVelAn = True
markPoints = False
#plotVelFFT = True
plotVelFFT = False

if(mediumType == "inhomog"):
	plotPresCurve = False
	plotVelCurve = False


if periodicType == "refl" or mediumType == "inhomog":
	#analytical function does not make sense
	plotPresAn = False
	plotRhoAn = False
	plotVelAn = False
	markPoints = False




















