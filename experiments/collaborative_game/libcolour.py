# -*- coding: utf-8 -*-

import numpy
from matplotlib import pyplot

def pol2cart(theta, radius, units='deg'):
	"""Convert from polar to cartesian coordinates
	
	**usage**:
		x,y = pol2cart(theta, radius, units='deg')
	"""
	if units in ['deg', 'degs']:
		theta = theta*numpy.pi/180.0
	xx = radius*numpy.sin(theta)
	yy = radius*numpy.cos(theta)
	
	return xx,yy


def lab2xyz(lab):
	
	"""Converts a CIE Lab triplet to a CIE XYZ triplet, using a 2 degree
	observer and a D65 illuminant.
	More information:
	CIE Lab: http://en.wikipedia.org/wiki/Lab_color_space
	CIE XYZ: http://en.wikipedia.org/wiki/CIE_1931_color_space
	http://www.easyrgb.com/index.php?X=MATH&H=08#text8
	reference values: http://www.easyrgb.com/index.php?X=MATH&H=15#text15
	"""
	
	var = {'Y': (lab[0] + 16.0) / 116.0}
	var['X'] = (lab[1] / 500.0) + var['Y']
	var['Z'] = var['Y'] - (lab[2] / 200.0)
	
	for g in ['X','Y','Z']:
		if var[g] > 0.008856:
			var[g] = var[g]**3
		else:
			var[g] = ((var[g] - 16.0) / 116.0) / 7.787
	
	# references for observer=2 degrees, illuminant=D65
	ref = {'X':95.047}
	ref['Y'] = 100.000
	ref['Z'] = 108.883
	
	out = {}
	for g in var.keys():
		out[g] = ref[g] * var[g]
	
	return (out['X'], out['Y'], out['Z'])


def xyz2rgb(xyz, mode='255'):
	
	"""Converts a CIE XYZ triplet to a standard RGB triplet, using a 2 degree
	observer and a D65 illuminant.
	More information:
	CIE XYZ: http://en.wikipedia.org/wiki/CIE_1931_color_space
	sRGB: http://en.wikipedia.org/wiki/SRGB
	conversion: http://www.easyrgb.com/index.php?X=MATH&H=01#text1
	reference values: http://www.easyrgb.com/index.php?X=MATH&H=15#text15
	"""
	
	# variables divided by 100, asuming the XYZ values are in the range 0-100
	# (which should be the case if the lab2xyz function was used)
	var = {'X': xyz[0] / 100.0}
	var['Y'] = xyz[1] / 100.0
	var['Z'] = xyz[2] / 100.0

	var['R'] = var['X'] *  3.2406 + var['Y'] * -1.5372 + var['Z'] * -0.4986
	var['G'] = var['X'] * -0.9689 + var['Y'] *  1.8758 + var['Z'] *  0.0415
	var['B'] = var['X'] *  0.0557 + var['Y'] * -0.2040 + var['Z'] *  1.0570
	
	for g in ['R','G','B']:
		if var[g] > 0.0031308:
			var[g] = 1.055 * ( var[g] ** ( 1 / 2.4 ) ) - 0.055
		else:
			var[g] = 12.92 * var[g]
	
	if mode == '1':
		return (var['R'], var['G'], var['B'])
	else:
		return (int(var['R']*255), int(var['G']*255), int(var['B']*255))


def create_colourwheel(L, r, savefile=None):
	
	"""Produces a colourwheel, and returns a NumPy array with the RGB colours
	(shape = (360,3)); if a location for a save file is provided, it will
	also store an image of the colour wheel
	"""
	
	theta = numpy.arange(0,360,1)
	figsize = (10.0,10.0)
	dpi = 300.0
	a = numpy.zeros(theta.shape)
	b = numpy.zeros(theta.shape)
	cols = numpy.zeros((10, len(theta), 3))
	for i in range(len(theta)):
		a[i], b[i] = pol2cart(theta[i], r, units='deg')
		(R,G,B) = xyz2rgb(lab2xyz((L,a[i],b[i])))
		if 0 <= R <= 255:
			if 0 <= G <= 255:
				if 0 <= B <= 255:
					cols[:,i,0] = R
					cols[:,i,1] = G
					cols[:,i,2] = B
	
	if savefile:
		fig, ax = pyplot.subplots(figsize=figsize, dpi=dpi)
		ax.imshow(cols)
		ax.axes.set_yticks([])
		ax.axes.set_yticklabels([])
		ax.set_title("colour wheel at L=%d, r=%d" % (L, r))
		ax.set_xlabel("angle (degs)")
		fig.savefig(savefile)
		pyplot.close(fig)
	
	return cols[0,:,:]

if __name__ == "__main__":
	cw = []
	Ls = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]
	rs = [11, 13, 15, 18, 20, 22, 24, 27, 29, 31, 33, 36, 38, 39, 30, 22, 14]
	for i in range(len(Ls)):
		cw.append(create_colourwheel(Ls[i], rs[i], savefile="colourwheel_L%d_r%d.png" % (Ls[i],rs[i])))
		