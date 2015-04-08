import numpy as np
cimport numpy as np

import cython
cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t




@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef as2(np.ndarray[DTYPE_t, ndim=1] curve, np.ndarray[DTYPE_t, ndim=1] ampIni, np.ndarray[DTYPE_t, ndim=1] csZ0, np.ndarray[DTYPE_t, ndim=1] csZt, np.ndarray[DTYPE_t, ndim=1] newZ, int nint, float z0, float zf):
	cdef int zIndexNewZ

	for index in range(nint + 2):
		zIndexNewZ =  int(float(nint)*(newZ(index) - z0)/(zf - z0) + 0.5 );
		curve[zIndexNewZ] = ampIni[index] * (csZ0[index]  * csZt[index] )** (0.5)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef as3(np.ndarray[DTYPE_t, ndim=1] curve, np.ndarray[DTYPE_t, ndim=1] ampIni, np.ndarray[DTYPE_t, ndim=1] csZ0, np.ndarray[DTYPE_t, ndim=1] csZt, int nint):

	for index in range(nint + 2):
		curve[index] = ampIni[index] * (csZ0[index]  / csZt[index] )** (0.5)
