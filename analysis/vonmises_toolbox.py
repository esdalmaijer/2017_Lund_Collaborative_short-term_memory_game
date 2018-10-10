import copy

import numpy
import scipy
from matplotlib import pyplot
import scipy.special
import scipy.stats


# FITTING FUNCTION

def fullspacestimate(responses, targetoris, nontargetoris=None, fituniform=True):
	
	"""Does a parameter estimate for the combination of two Von Mises
	distributions and a uniform distribution that fits the responses to
	targets and non-targets best, by sampling a large portion of the parameter
	space (SD from 0.01 to 3, pU and pNT from 0 to 1, all in steps of 0.01)

	Reference:
	Bays PM, Catalao RFG & Husain M. (2009). The precision of visual working
		memory is set by allocation of a shared resource. Journal of Vision,
		9(10):7, p. 1-11
	
	arguments
	
	reponses		-	a participant's recall responses (NumPy array or list)
	targetoris	-	target orientations (NumPy array or list)
	
	keyword arguments
	
	nontargetoris	-	orientations of all non-targets (NumPy array or list;
					each index being a number, or another array/list); or
					None	when non-targets were not present or should be
					ignored (default = None)
	fituniform	-	Boolean indicating whether a uniform distribution
					should be fitted to the data
	
	returns
	
	[kT, kNT, pT, pNT, pU], LL
	k 			-	kappa, concentration parameter dictating the
					distribution's spread (equivalent to a normal
					distribution's standard deviation)
					NOTE: kNT will be None if no nontargetoris were passed
	p 			-	proportion of responses to the Target, NonTargets, or
					to a random orientation (Uniform)
					NOTE: pNT will be None if no nontargetoris were passed
	LL 			-	log likelihood value of the estimate
	"""

	# number of nontargets (updated as required)
	nn = 1

	# check input
	if type(responses) in [list, tuple]:
		responses = numpy.array(responses)
	if type(targetoris) in [list, tuple]:
		targetoris = numpy.array(targetoris)
	if nontargetoris != None:
		if type(nontargetoris) in [tuple, list, numpy.array, numpy.ndarray]:
			for i in range(len(nontargetoris)):
				if type(nontargetoris[i]) in [list, tuple, numpy.array, numpy.ndarray]:
					nontargetoris[i] = numpy.array([nontargetoris[i]])
					nn = len(nontargetoris[i])
			nontargetoris = numpy.array(nontargetoris)
	
	# number of trials
	if len(responses) == len(targetoris):
		n = len(responses)
	else:
		raise Exception("ERROR in baystimate: unequal number of responses (N=%d) and target orientations (N=%d)" % (len(responses), len(targetoris)))

	# check if we need to take non-target responses into account
	if nontargetoris is None:
		nn = 0
	
	# calculate errors to target orientations
	err = calc_angular_dist(targetoris, responses)
	# calculate errors to nontarget orientations
	if nn == 1:
		nerr = calc_angular_dist(nontargetoris, responses)
	elif nn > 1:
		nerr = []
		for i in range(nn):
			nerr.append(calc_angular_dist(nontargetoris[:,i], responses))
		nerr = numpy.array(nerr)
	else:
		nerr = None
	
	# starting values
	B = [numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan]
	LL = -numpy.inf

	# parameter estimate ranges
	sdr = numpy.arange(0.01, 1.01, 0.025)
	pur = numpy.arange(0, 1.01, 0.025)
	if nn > 0:
		ptr = numpy.arange(0, 1.01, 0.025)
	else:
		ptr = numpy.zeros(1)
	if not fituniform:
		pur = numpy.zeros(1)
		sdr = numpy.arange(0.01, 10.01, 0.025)

	# explore parameter space
	for sd in sdr:
		kappa = sd2kappa(sd)
		for pU in pur:
			for pNT in ptr:
				
				# skip if the proportions don't add up
				if pU + pNT > 1:
					continue

				# create target error distribution
				pT = 1 - (pU + pNT)
				tardist = pT * vonmisespdf(err, 0, kappa)
				# create uniform error distribution
				unidist = pU * numpy.ones(n) / (2 * numpy.pi)
				# create single nontarget error distribution
				if nn == 1:
					nondist = pNT * vonmisespdf(nerr, 0, kappa)
				# create nontarget error distribution for every non-target
				elif nn > 1:
					nondist = []
					for i in range(nn):
						nondist.append( (pNT/nn) * vonmisespdf(nerr[i], 0, kappa))
					# sum all non-target distributions
					nondist = numpy.sum(nondist, axis=0)
				# create dummy nontarget error distribution
				else:
					nondist = numpy.zeros(n)
				
				# sum distributions
				distsum = tardist + unidist + nondist
				
				# calculate log likelihood
				ll = numpy.sum(numpy.log(distsum))
				
				# store parameters if the log likelihood if the current
				# estimate is better than the stored one
				if ll > LL:
					if nn > 0:
						B = [0+kappa, 0+kappa, 0+pT, 0+pNT, 0+pU]
					else:
						B = [0+kappa, None, 0+pT, 0+pNT, 0+pU]
					LL = 0 + ll
	
	# return values	
	return B, LL


def baysfit(responses, targetoris, nontargetoris=None):
	
	"""

	arguments
	
	reponses		-	a participant's recall responses (NumPy array or list)
	targetoris	-	target orientations (NumPy array or list)
	
	keyword arguments
	
	nontargetoris	-	orientations of all non-targets (NumPy array or list;
					each index being a number, or another array/list); or
					None	when non-targets were not present or should be
					ignored (default = None)
	
	returns
	
	[kT, kNT, pT, pNT, pU], LL
	k 			-	kappa, concentration parameter dictating the
					distribution's spread (equivalent to a normal
					distribution's standard deviation)
					NOTE: kNT will be None if no nontargetoris were passed
	p 			-	proportion of responses to the Target, NonTargets, or
					to a random orientation (Uniform)
					NOTE: pNT will be None if no nontargetoris were passed
	LL 			-	log likelihood value of the estimate
	"""

	# make sure there is an equal number of responses and targets
	if len(responses) != len(targetoris):
		raise Exception("ERROR in baysfit: unequal number of responses (N=%d) and target orientations (N=%d)" % (len(responses), len(targetoris)))
	
	# starting parameter values
	kstarts = [1, 10, 100]
	nPstarts = [0.01, 0.1, 0.4]
	uPstarts = [0.01, 0.1, 0.4]
	
	# check if we need to take non-target responses into account
	if nontargetoris is None:
		nPstarts = [0]
	
	# starting values for log likelihood and parameter estimates
	LL = -numpy.inf
	B = [numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan]

	# loop through all starting values for all parameters
	for ks in kstarts:
		for ns in nPstarts:
			for us in uPstarts:
				# parameter estimation
				b, ll = baystimate(responses, targetoris, nontargetoris, ks, 1-(ns+us), ns, us)
				# save parameters if the log likelihood is higher
				if ll > LL and not numpy.isnan(ll):
					LL = copy.deepcopy(ll)
					B = copy.deepcopy(b)
	
	return B, LL


def baystimate(responses, targetoris, nontargetoris, kTstart, pTstart, pNTstart, pUstart):
	
	"""Does a parameter estimate for the combination of two Von Mises
	distributions and a uniform distribution that fits the responses to
	targets and non-targets best, provided the passed starting points.

	Reference:
	Bays PM, Catalao RFG & Husain M. (2009). The precision of visual working
		memory is set by allocation of a shared resource. Journal of Vision,
		9(10):7, p. 1-11
	
	arguments
	
	reponses		-	a participant's recall responses (NumPy array or list)
	targetoris	-	target orientations (NumPy array or list)
	nontargetoris	-	orientations of all non-targets (NumPy array or list;
					each index being a number, or another array/list); or
					None	when non-targets were not present or should be
					ignored
	kTstart 		-	starting value for the kappa estimate
	pTstart 		-	starting value for the proportion of target responses
	pNTstart 		-	starting value for the proportion of non-target
					responses
	pUstart 		-	starting value for the proportion of random responses
	
	returns
	
	[kT, kNT, pT, pNT, pU], LL
	k 			-	kappa, concentration parameter dictating the
					distribution's spread (equivalent to a normal
					distribution's standard deviation)
					NOTE: kNT will be None if no nontargetoris were passed
	p 			-	proportion of responses to the Target, NonTargets, or
					to a random orientation (Uniform)
					NOTE: pNT will be None if no nontargetoris were passed
	LL 			-	log likelihood value of the estimate
	"""

	# number of nontargets (updated as required)
	nn = 1

	# check input
	if type(responses) in [list, tuple]:
		responses = numpy.array(responses)
	if type(targetoris) in [list, tuple]:
		targetoris = numpy.array(targetoris)
	if nontargetoris != None:
		if type(nontargetoris) in [tuple, list]:
			for i in range(len(nontargetoris)):
				if type(nontargetoris[i]) in [list, tuple]:
					nontargetoris[i] = numpy.array([nontargetoris[i]])
					nn = len(nontargetoris[i])
			nontargetoris = numpy.array(nontargetoris)
	
	# number of trials
	if len(responses) == len(targetoris):
		n = len(responses)
	else:
		raise Exception("ERROR in baystimate: unequal number of responses (N=%d) and target orientations (N=%d)" % (len(responses), len(targetoris)))

	# check if we need to take non-target responses into account
	if nontargetoris is None:
		nn = 0
	
	# calculate errors to target orientations
	err = calc_angular_dist(targetoris, responses)
	# calculate errors to nontarget orientations
	if nn == 1:
		nerr = calc_angular_dist(nontargetoris, responses)
	elif nn > 1:
		nerr = []
		for i in range(len(nn)):
			nerr.append(calc_angular_dist(nontargetoris[i], responses))
		nerr = numpy.array(nerr)
	else:
		nerr = None
	
	# maximum values
	maxdLL = 10**-4
	maxi = 10**3

	# starting values
	kappa = 0 + kTstart
	pT = 0 + pTstart
	pNT = 0 + pNTstart
	pU = 0 + pUstart
	LL = numpy.nan
	dLL = numpy.nan
	
	# starting log likelihood at a complete uniform distribution
	#LL = numpy.sum(numpy.log(numpy.ones(n) / (2 * numpy.pi)))
	
	# run until we reach a good level of certainty, or the maximal amount of
	# iterations
	i = 0
	while True:
		
		#print("i=%4.0f, pT=%.2f, pU=%.2f" % (i,pT,pU))

		# update iteration number
		i += 1

		# create target error distribution
		tardist = pT * vonmisespdf(err, 0, kappa)
		# create uniform error distribution
		unidist = pU * numpy.ones(n) / (2 * numpy.pi)
		# create single nontarget error distribution
		if nn == 1:
			nondist = pNT * vonmisespdf(nerr, 0, kappa)
		# create nontarget error distribution for every non-target
		elif nn > 1:
			nondist = []
			for i in range(nn):
				nondist.append( (pNT/nn) * vonmisespdf(nerr[i], 0, kappa))
			# sum all non-target distributions
			nondist = numpy.sum(nondist, axis=0)
		# create dummy nontarget error distribution
		else:
			nondist = numpy.zeros(n)
		
		# sum distributions
		distsum = tardist + unidist + nondist
		
		# calculate log-likelihood difference
		dLL = LL - numpy.sum(numpy.log(distsum))
		# calculate log likelihood
		LL = numpy.sum(numpy.log(distsum))
		
		# stop looping if we reached a local minimum in log likelihood space,
		# or if we have reached the maximal amount of iterations
		if numpy.abs(dLL) > maxdLL or i > maxi:
			break
		if numpy.isnan(LL):
			i = maxi
			break
		
		# update distribution proportions
		pT = numpy.sum(tardist / distsum) / n
		pNT = numpy.sum(nondist / distsum) / n
		pU = numpy.sum(unidist / distsum) / n
		
		# do some magic (fmin search)
		rw = [distsum / tardist]
		s = [numpy.sin(err)]
		c = [numpy.cos(err)]
		if nn == 1:
			rw.append(distsum / nondist)
			s.append(numpy.sin(nerr))
			c.append(numpy.cos(nerr))
		elif nn > 1:
			for i in range(nn):
				rw.append(distsum / nondist[i])
				s.append(numpy.sin(nerr[i]))
				c.append(numpy.cos(nerr[i]))
		rw = numpy.array(rw)
		s = numpy.array(s)
		c = numpy.array(c)

		r = [numpy.sum(numpy.sum(s * rw, axis=0)), numpy.sum(numpy.sum(c * rw, axis=0))]
		
		# update kappa
		if numpy.sum(rw) == 0:
			kappa = 0
		else:
			R = (r[0]**2 + r[1]**2)**0.5 / numpy.sum(rw)
			if 0 <= R < 0.53:
				kappa = 2 * R + R**3 + (5 * R**5)/6
			elif R < 0.85:
				kappa = -0.4 + 1.39 * R + 0.43/(1 - R)
			else:
				kappa = 1/(R**3 - 4 * R**2 + 3 * R)
		
		# kappa correction for small sample sizes
		if n < 15:
			if kappa < 2:
				kappa = numpy.max(kappa - 2/(n*kappa), 0)
			else:
				kappa = kappa * (n-1)**3/(n**3+n)

	# return values
	if nontargetoris is None:
		kappaNT = None
	else:
		kappaNT = 0 + kappa
	if i < maxi:
		B = [kappa, kappaNT, pT, pNT, pU]
	else:
		B = [numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan]
		LL = numpy.nan
	
	return B, LL


# CORE FUNCTIONS

def vonmisespdf(x, mu, kappa):
	
	"""Creates a probability density function for a Von Mises distribution.
	NOTE: x and mu should be in radians, x between -pi and +pi
	
	arguments

	x		-	x values to be used for the pdf
	mu		-	distribution centre (average)
	kappa	-	distribution spread (equivalent to a normal distribution's
				standard deviation)
	
	keyword arguments

	mode		-	input mode, either 'degrees' or 'radians'
	
	returns

	pdf		-	numpy.array of the probability density of the Von Mises
				distribution with the passed parameters
	"""

	return numpy.exp(kappa * numpy.cos(x-mu)) / (2*numpy.pi * scipy.special.iv(0,kappa))

def entropy(kappa, mode='natural'):
	
	"""Calculates the entropy of a Von Mises pdf with given kappa.
	
	arguments

	kappa		-	kappa value (spread parameter) of the Von Mises pdf
	
	keyword arguments
	
	mode		-	string indicating the type of logarithm to use to
				calculate the entropy; 'natural' for ln, or 'binary'
				for a base of 2 (Shannon entropy)
	"""
	
	if mode == 'natural':
		H = numpy.log(2*numpy.pi * scipy.special.iv(0,kappa)) - kappa * (scipy.special.iv(1,kappa) / scipy.special.iv(0,kappa))
	elif mode == 'binary':
		H = numpy.log2(2*numpy.pi * scipy.special.iv(0,kappa)) - kappa * (scipy.special.iv(1,kappa) / scipy.special.iv(0,kappa))
	
	return H

def kappa2sd(kappa):
	
	"""Converts a kappa value of a Von Mises distribution to the equivalent
	standard deviation of the corresponding normal distribution
	"""

	if type(kappa) not in [numpy.array, numpy.ndarray]:
		kappa = numpy.array([kappa])
	
	zero = kappa == 0
	infinite = kappa == numpy.inf
	other =  numpy.invert(zero | infinite)
	
	sd = numpy.zeros(kappa.shape)

	sd[zero] = numpy.inf
	sd[infinite] = 0
	sd[other] = numpy.sqrt(-2 * numpy.log(scipy.special.iv(1,kappa[other]) / scipy.special.iv(0,kappa[other])))

	if sd.shape == ():
		return sd
	elif len(sd) == 1:
		return sd[0]
	else:
		return sd

def sd2kappa(sd):
	
	"""Converts a standard deviation of a normal distribution to the
	equivalent kappa value of the corresponding Von Mises distribution
	"""
	
	if type(sd) not in [numpy.array, numpy.ndarray]:
		sd = numpy.array([sd])
	
	r = numpy.exp((-1*(sd)**2) / 2)
	k = 1 / (r**3 - 4*r**2 + 3*r)
	selone = numpy.array(r < 0.85)
	k[selone] = -0.4 + 1.39 * r[selone] + 0.43/(1 - r[selone])
	seltwo = numpy.array(r < 0.53)
	k[seltwo] = 2 * r[seltwo] + r[seltwo]**3 + (5 * r[seltwo]**5)/6;
	
	if k.shape == ():
		return k
	elif len(k) == 1:
		return k[0]
	else:
		return k

def calc_angular_dist(ori1, ori2, zero=360):
	
	"""Calculates the angular distance between two angles (in DEGREES!),
	each angle representing a point on an imaginary circle (default). Note
	that the origin of the round distribution can be determined by setting the
	'zero' keyword.
	
	arguments

	ori1	-	first orientation(s) (the target orientation), in degrees
	ori2	-	second orientation (the response), in degrees
	
	keyword arguments
	
	zero	-	origin of the distribution (360 for a full circle, 180 for lines
			along the diameter of a circle)
	"""
	
	if type(ori1) in [int,float, numpy.int, numpy.int16, numpy.int32, numpy.int64, numpy.float, numpy.float16, numpy.float32, numpy.float64]:
		ori1 = [ori1]
	elif type(ori1) == numpy.ndarray and ori1.shape == ():
		ori1 = [ori1]
	ori1 = numpy.array(ori1)
	if type(ori2) in [int,float, numpy.int, numpy.int16, numpy.int32, numpy.int64, numpy.float, numpy.float16, numpy.float32, numpy.float64]:
		ori2 = [ori2]
	elif type(ori2) == numpy.ndarray and ori2.shape == ():
		ori2 = [ori2]
	ori2 = numpy.array(ori2)
	
	# moduli after division
	ori1 = ori1 % zero
	ori2 = ori2 % zero
	
	# calculate the difference
	error = ori2 - ori1
	# where the difference is larger then a clockwise 90 degree rotation, it
	# should be counted as a counter-clockwise (negative) rotation
	error[error>zero/2] = -1 * (zero - error[error>zero/2])
	error[error<=-zero/2] = (zero + error[error<=-zero/2])
	
	return error

