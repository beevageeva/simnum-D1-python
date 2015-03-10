import numpy as np
from constants import z0, zf


def getPeriodicX(xval, a=z0, b=zf):
	"""
		returns x in [z0,zf]  in concordance with boundary conditions repeat
	"""
	p = float(b - a)
	k = int((xval-a)/p)
	res = xval - k * p
	if(res < a):
		res+=p
	if(res > b):
		res-=p
	return res


def  lrBoundaryConditions(array, skip=0):
	n = array.shape[0] - 1
	array = np.insert(array, 0, array[0])
	array = np.insert(array, n+2, array[-1])
	return array
