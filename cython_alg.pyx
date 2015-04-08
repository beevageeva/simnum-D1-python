import numpy as np
cimport numpy as np

import cython
cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t




@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_interm_u_array(np.ndarray[DTYPE_t, ndim=1] res, np.ndarray[DTYPE_t, ndim=1] u, np.ndarray[DTYPE_t, ndim=1] f, int nint, float lambdaParam):
	cdef int i
	for i in range(1, nint+2):
		res[i-1] = 0.5 * (u[i] + u[i-1]) - 0.5 * lambdaParam  * (f[i] - f[i-1])


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_final_u_array(np.ndarray[DTYPE_t, ndim=1] res, np.ndarray[DTYPE_t, ndim=1] u, np.ndarray[DTYPE_t, ndim=1] intermF, int n, float lambdaParam, int skip):
	cdef int i
	for i in range(0, n):
		res[i] =  u[i+skip] - lambdaParam  * (intermF[i+1] - intermF[i]) 


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_singlestep_u_array(np.ndarray[DTYPE_t, ndim=1] res, np.ndarray[DTYPE_t, ndim=1] u, np.ndarray[DTYPE_t, ndim=1] f, int nint, float lambdaParam):
	cdef int i
	for i in range(1, nint+1):
		res[i-1] =  0.5 * (u[i-1] + u[i+1]) - 0.5 * lambdaParam  * (f[i+1] - f[i-1]) 



