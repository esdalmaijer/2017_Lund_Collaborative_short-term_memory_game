import copy
import numpy
from scipy import integrate, stats

def create_cdf(X):
	
	"""Create the cummulative density function of a continuous random
	variable, e.g. observed data.
	
	arguments

	X	-	Observed data from which a cdf should be constructed; please
			provide the data as a vector (Nx1 NumPy array or list).
	
	returns
	
	(bins, cdf)

	bins	-	All unique values in X, used as bins for the cdf.

	cdf	-	The cummulative density for each value in bins.
	"""
	
	# convert the data to a NumPy array.
	data = numpy.array(X)
	# the bins are all unique values in the data
	bins = copy.deepcopy(data)
	bins.sort()
	# calculate the cummulative frequency for each unique value in the data
	cumfreq, lowerlim, binsize, extra = stats.cumfreq(data, numbins=len(bins))
	# transform the cummulative frequencies to a cdf
	cdf = cumfreq / numpy.max(cumfreq)

	return bins, cdf

def create_pdf(X, mode='hist'):
	
	"""Create a probability density function of a continuous random variable,
	e.g. observed data.
	
	arguments

	X	-	Observed data from which a cdf should be constructed; please
			provide the data as a vector (Nx1 NumPy array or list).

	keyword arguments

	mode	-	String indicating the mode: either 'hist' to create a (rather
			pointy) pdf, or 'kde' to use a Gaussian kernel density
			estimate (for a smoother result).
	
	returns
	
	(bins, pdf)

	bins	-	All unique values in X, used as bins for the cdf.

	pdf	-	The cummulative density for each value in bins.
	"""
	
	# covert the data to a NumPy array
	data = numpy.array(X)
	# the bins are all unique values
	bins = copy.deepcopy(data)
	bins.sort()
	# create a pdf-like histogram
	if mode == 'hist':
		pdf, binedges = numpy.histogram(data, bins=len(bins), density=True)
	elif mode == 'kde':
		# create a probability density funtion function
		kde = stats.gaussian_kde(data)
		# create the values of the pdf at each value in the data
		pdf = kde(bins)
	
	return bins, pdf

def pdf2cdf(bins, pdf, mode='x'):
	
	"""Create a cummulative density function from a probability density
	function.
	
	arguments

	bins	-	The values (bins) of which the pdf gives the probability
			density.

	pdf	-	The probability density for each value in bins.

	keyword arguments

	mode	-	String indicating the mode: either 'x' to use bins as actual x
			values when integrating, or 'dx' to use the average distance
			between bins.
	
	returns
	bins	-	Same bins.

	cdf	-	The cummulative density for each value in x.
	"""
	
	# create an empty list to contain values
	cdf = []
	# loop through all bins
	for i in range(len(bins)):
		# integrate the pdf from the first to the current bin
		if mode== 'x':
			cdf.append(integrate.trapz(pdf[0:i],x=bins[0:i]))
		elif mode == 'dx':
			cdf.append(integrate.trapz(pdf[0:i],dx=numpy.mean(numpy.diff(bins[0:i]))))
	
	return bins, numpy.array(cdf)

def cdf2pdf(bins, cdf):
	
	"""Creates a probability density function out of a cummulative density
	function by using the slope between each bin in the cdf.
	
	arguments

	bins	-	All unique values in X, used as bins for the cdf (and thus
			ordered from lowest to highest!)

	cdf	-	The cummulative density for each value in bins.

	returns
	
	(bins, pdf)

	bins	-	All unique values in X, used as bins for the pdf.

	pdf	-	The probability density for each value in bins.
	"""
	
	# calculate average bin size
	dx = numpy.mean(numpy.diff(bins))
	# calculate the gradients between each point in the cdf, assuming equal
	# bin size
	pdf = numpy.gradient(cdf, dx)
	
	return bins, pdf

#	# calculate the distances between bins
#	db = numpy.diff(bins)
#	# calculate the distances between frequency densities of the values in
#	# each bin
#	dc = numpy.diff(cdf)
#	# the slope of the cdf between each bin value is 
#	pdf = dc / db
#	
#	return bins, pdf

def entropy_from_pdf(bins, pdf, mode='nats'):
	
	"""Calculates the entropy of a probability density function, defined as
	sum(p*log(p)), where p is the integral of each point in the pdf.
	
	arguments
	
	bins	-	The values (bins) of which the pdf gives the probability
			density.

	pdf	-	The probability density for each value in bins.

	keyword arguments

	mode	-	String indicating the mode: either 'nats' to use a natural
			logarithm (with base Euler), or 'bits' to use a binary
			logarithm (with base 2)
	
	returns
	
	(e, fish)

	e	-	Entropy.

	fish	-	Fisher self-information.
	"""
	
	# choose log function
	if mode == 'nats':
		log = numpy.log
	elif mode == 'bits':
		log = numpy.log2
	
	# calculate P for every value in the bins
	p = numpy.zeros(len(bins))
	for i in range(1,len(bins)):
		p[i] = integrate.trapz(pdf[0:i], dx=numpy.mean(numpy.diff(bins[0:i]))) - integrate.trapz(pdf[0:i-1], dx=numpy.mean(numpy.diff(bins[0:i-1])))
	
	# calculate entropy
	non0 = p > 0
	e = numpy.sum(p[non0] * log(p[non0]))
	
	# calculate the integral of log(p)
	gr = numpy.diff(log(p[non0]))
	# calculate Fisher self-information
	fish = numpy.sum(p[non0][1:] * gr * gr);
	
	return e, fish

def entropy_from_cdf(bins, cdf, mode='nats'):
	
	"""Calculates the entropy of a cummulative density function, defined as
	sum(p*log(p)), where p is the difference between each point in the cdf.
	
	arguments
	
	bins	-	The values (bins) of which the pdf gives the probability
			density.

	cdf	-	The cummulative density for each value in bins.

	keyword arguments

	mode	-	String indicating the mode: either 'nats' to use a natural
			logarithm (with base Euler), or 'bits' to use a binary
			logarithm (with base 2)
	
	returns
	
	(e, fish)

	e	-	Entropy.

	fish	-	Fisher self-information.
	"""
	
	# choose log function
	if mode == 'nats':
		log = numpy.log
	elif mode == 'bits':
		log = numpy.log2
	
	# calculate the difference between each point and the previous point in
	# the cdf (this is the probability of the value at that point)
	p = cdf[1:] - cdf[:-1]
	bins = bins[1:]

	# calculate entropy
	non0 = p > 0
	e = numpy.sum(p[non0] * log(p[non0]))
	
	# calculate the integral of log(p)
	gr = numpy.diff(log(p[non0]))
	# calculate Fisher self-information
	fish = numpy.sum(p[non0][1:] * gr * gr);
	
	return e, fish

def KL_divergence_continuous(pdfP, pdfQ, dx, mode='nats'):
	
	"""Calculates the Kullback-Leibler divergence between two probability
	density distributions (P and Q, both assumed to be continuous random
	variables). The KL divergence is the amount of extra information that is
	required to identify an event from distribution P when distribution Q is
	used to decode said event. Note that the divergence of P and Q is not
	necessarily the same as the divergence of Q and P!
	
	arguments
	
	pdfP	-	The probability density distribution of P (the 'true'
			distribution).

	pdfQ	-	The probability density distribution of Q (the 'wrong'
			distribution).

	dx	-	The step size on the x-axis for each pdf.

	keyword arguments

	mode	-	String indicating the mode: either 'nats' to use a natural
			logarithm (with base Euler), or 'bits' to use a binary
			logarithm (with base 2). (default = 'nats')
	
	returns

	Dkl	-	Kullback-Leibler divergence (in nats or bits, depending on what
			mode was selected).
	"""
	
	# The KL divergence is the integral of (P * log(P/Q) * dx), where the base
	# of the log is either Eulers number (mode=='nats') or 2 (mode=='bits')
	if mode == 'nats':
		y = pdfP * numpy.log(pdfP/pdfQ)
	elif mode == 'bits':
		y = pdfP * numpy.log2(pdfP/pdfQ)
    
    # When there are 0s in pdfP, they will result in log being called on 0,
    # which will result in -inf. Then, when multiplying 0 with -inf, a NaN
    # will be the result. Hence, we need to replace NaN values by 0.
	y[numpy.isnan(y)] = 0

	return integrate.trapz(y, dx=dx)

	
