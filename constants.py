
"""
Parameters:
	gamma = constant (5/3)
	problemType = "soundwave" | "riemann"
	[z0, zf] domain
	nint - number of intervals
	timeEnd - the default timeEnd if not specified as parameter to main.py
	schemeType = "lf" | "fg"
	loopType = "python" | "weave"  (if set to weave the loops are run calling C from python using weave). 
This is used in alg.py (when calculating time-derivatives) and 
analytic_solution.py when calculating analytic pressure in inhomogeneous medium

	fcfl - is different for each scheme type
	In the case of schemeType = "fg" we can define bcStep parameter: "interm" or "final" (the step when to apply the boundary conditions)
	

"""


gamma = 5.0/3

problemType = "soundwave"
#problemType = "riemann"
if problemType == "soundwave":
	z0 = 3.1
	zf = 7.4
elif problemType == "riemann":
	z0 = -5
	zf = 10
	from riemann_params import riemann_problemType
	if riemann_problemType == "complete":
		gamma = 1.4

#nint = 32
#nint = 64
#nint=128
#nint =  256
#nint = 1024
nint = 2048
#nint = 4096
#nint = 32

timeEnd = 1.0 #if not set as program argument it's taken from here

#schemeType = "lf"  # scheme type may be lf(Lax - Fr) or fg (first generation)
schemeType = "fg" 
#loopType = "python" 
loopType = "weave" 
if schemeType == "lf":
	fcfl = 0.99 #use this for lax - fr scheme type
elif schemeType == "fg":
	fcfl = 0.97#use this for first generation scheme
	bcStep = "interm" 
	#bcStep = "final" 


