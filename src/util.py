# Copyright (C) 2016 Peloton
#########################
# Contain classes and routines for
# 	* binning
# 	* file management
# author: julien@sussex
# See 1611.01446
#########################
import glob,sys,os
import cPickle as pickle
import errno
import numpy as np

########################################################################
# Class to manipulate binned quantities

class binner4cov():
	def __init__(self,lmin,lmax,step=20,method='CMB',bins=None):
		"""
		Binning routines for covariance matrices.
		The binning is pre-defined, unless specified by the user (not yet avail)
		Input:
			* lmin: int, minimum multipole
			* lmax: int, maximum multipole
			* step: int, bin size. Default is 20.
			* method: string, the type of the data (CMB. phi, or phiCMB).
			* bins: 1D array, flat array containing bin boundaries (e.g. [2,20,40,60,...,lmax]).
				Not yet available.
		"""
		if bins is None:
			if step > lmin:
				self.bin_boundaries_all = np.array([lmin] + range(step,lmax,step) + [lmax])
			else:
				self.bin_boundaries_all = np.array([lmin] + range(step+lmin,lmax,step) + [lmax])
			self.bin_boundaries = np.array([[a,b] for a,b in zip(self.bin_boundaries_all[0:-1],self.bin_boundaries_all[1:])])
			self.bin_centers = np.mean(self.bin_boundaries,axis=1)
			self.bin_size = np.diff(self.bin_boundaries_all)
			self.nbins = len(self.bin_centers)
			self.lmax = lmax
			self.lmin = lmin

			if method == 'CMB':
				self.alpha_x = 1.0
				self.alpha_y = 1.0
			elif method == 'phi':
				self.alpha_x = 2.0
				self.alpha_y = 2.0
			elif method == 'phiCMB':
				self.alpha_x = 1.0
				self.alpha_y = 2.0
		else:
			print 'Not yet available'
			sys.exit()

	def P_operator(self,alpha):
		'''
		Binning operator (size n bins x n modes).
		Takes alpha in input (the power of the ell).
		'''
		P = np.zeros((self.nbins,self.lmax+1))
		for index_b,[b_l,b_r] in enumerate(self.bin_boundaries):
			for ell in range(b_l,b_r):
				norm = (self.bin_centers[index_b]*(self.bin_centers[index_b]+1.))**alpha * self.bin_size[index_b]
				P[index_b][ell] = (ell*(ell+1))**alpha / norm
		return P

	def bin_this_cov(self,cov):
		'''
		Take unbinned covariance and bin it.
		'''
		cov_bin = np.zeros((self.nbins,self.nbins))
		Px = self.P_operator(alpha=self.alpha_x)
		Py = self.P_operator(alpha=self.alpha_y)
		cov_bin = dot(Py,cov[0:self.lmax+1,0:self.lmax+1],np.transpose(Px))
		return cov_bin

	def bin_this_spec(self,spec):
		'''
		Take unbinned power spectrum and bin it.
		'''
		spec_bin = np.zeros(self.nbins)
		P = self.P_operator(alpha=self.alpha_x)
		spec_bin = dot(P,spec[0:self.lmax+1])
		return spec_bin

class binner():
	def __init__(self, bins_l, bins_r):
		"""
		Binning routines. Left and right inclusive.
		For most general situation
		:param bins_l: left edges (inclusive)
		:param bins_r: right edges (inclusive)
		"""
		assert (len(bins_l) == len(bins_r)), "inconsistent inputs"
		assert (np.all(bins_r - bins_l > 0.)), "inconsistent input"
		self.bins_l = np.array(bins_l)
		self.bins_r = np.array(bins_r)

	def Nbins(self):
		return len(self.bins_l)

	def bin_centers(self):
		return 0.5 * self.bins_l + 0.5 * self.bins_r

	def bin_that(self, x, y, weights=None):
		ret = np.zeros(self.Nbins())
		if weights is None: weights = np.ones(len(x))
		assert (len(x) == len(y) and len(x) == len(weights)), "inconsistent inputs"
		for i, bin_l, bin_r in zip(xrange(self.Nbins()), self.bins_l, self.bins_r):
			idc = np.array(np.where((x >= bin_l) & (x <= bin_r)))
			if idc.size > 0.: ret[i] = np.sum(y[idc] * weights[idc]) / idc.size
		return ret

########################################################################
# Functions to manipulate binned spectra

def binvec(vec, ell, lbins, lmax,t=lambda l : 1.):
	""" rebins vec with non-uniform binning
		 * lbins		= list definining the bin edges [lbins[0], lbins[1]], [lbins[1], lbins[2]], ...
		 * (optional) w = l-dependent scaling to apply when accumulating into bins (in addition to number of modes in each bin).
	"""
	l = np.arange(0, lmax+1, 1)
	l = 0.5*(l[:-1] + l[1:]) # get bin centers
	t = t(l)

	nm, bins = np.histogram(ell, bins=np.arange(0, lmax+1, 1))
	modes = np.nan_to_num(nm)

	norm, bins = np.histogram(l, bins=lbins, weights=modes) # get number of modes in each l-bin.
	spec, bins = np.histogram(l, bins=lbins, weights=t*modes*np.nan_to_num(vec)) # bin the spectrum.

	# normalize the spectrum
	spec[np.nonzero(norm)] /= norm[np.nonzero(norm)]

	return spec

def rebin(cl,old_bin_centers,new_bin_boundaries,lmax,weight=1.0):
	interp_cl = np.interp(range(0,lmax+1),old_bin_centers,cl)
	return bin_spectra_theory(interp_cl,new_bin_boundaries,weight=weight)

def bin_spectra(sl,cl,lbins_boundaries,cl_transf=None,t=lambda l : 1.):
	if cl_transf is None:
		cl_transf = np.ones_like(cl)

	## Initialize
	ell = np.arange(len(cl))
	ellmax = len(cl[sl])
	bins_l = np.transpose(lbins_boundaries)[0]
	bins_u = np.transpose(lbins_boundaries)[1]-1
	binner_func = binner(bins_l, bins_u)

	bin_func = lambda x, y: binner_func.bin_that(x, y)
	xord = lambda ell: binner_func.bin_centers()

	## Bin spectrum
	bins = xord(ell[sl])
	cb = bin_func(ell[sl], cl[sl] * t(ell[sl]) )

	## debias
	if cl_transf is not None:
		## Bin beam function
		ell = np.arange(len(cl_transf))
		beamb = bin_func(ell, cl_transf**2 )
		cb = cb / beamb

	return cb, bins

def bin_spectra_theory(cl,lbins,weight=1.0):
	lbins_boundaries = np.array([[int(lbins[i]),int(lbins[i+1])] for i in range(len(lbins)-1)])
	ell = np.arange(len(cl))
	bins_l = np.transpose(lbins_boundaries)[0]
	bins_u = np.transpose(lbins_boundaries)[1]-1

	binner_func = binner(bins_l, bins_u)

	bin_func = lambda x, y: binner_func.bin_that(x, y)
	xord = lambda ell: binner_func.bin_centers()

	bins = xord(ell)
	cb = bin_func(ell, cl * (ell*(ell+1))**weight )
	cb = cb / (bins * (bins+1))**weight

	return cb

def linear_binning(lmin, lmax, nbins):
	lbins	  = np.linspace(lmin, lmax, nbins)	   # bin boundaries
	lbins_center	  = 0.5 * (lbins[:-1] + lbins[1:]) # bin centers
	lbins_boundaries = np.array([[int(lbins[i]),int(lbins[i+1])] for i in range(len(lbins)-1)])
	return lbins, lbins_center, lbins_boundaries

def Pbinning(bins,ell):
	P=np.zeros((len(bins)-1,len(ell)))
	shift=ell[0]
	minmax=zip(bins[0:-1],bins[1:])
	for i in range(len(bins)-1):
		for l in range(int(minmax[i][0]),int(minmax[i][1])):
			P[i][l-shift]=1./(2.*np.pi)*float(l)*(float(l)+1)/float(minmax[i][1]-minmax[i][0])
	return P

def Qbinning(bins,ell):
	Q=np.zeros((len(ell),len(bins)-1))
	shift=ell[0]
	minmax=zip(bins[0:-1],bins[1:])
	for i in range(len(bins)-1):
		for l in range(int(minmax[i][0]),int(minmax[i][1])):
			Q[l-shift][i]=2.*np.pi/float(l*(l+1))
	return Q

########################################################################
# Class and functions to manipulate files

def pickle_save(d,fn):
	'''
	Save data d into fn (.pkl) file
	Input:
		* d: dictionary, the data to be saved.
		* fn: string, the name of the file where data will be stored.
	'''
	with open(fn,'wb') as f:
		pickle.dump(d,f,protocol=2)

def pickle_load(fn):
	'''
	Load fn (.pkl) file
	Input:
		* fn: string, the name of the file (.pkl) containing data
	'''
	with open(fn,'rb') as f:
		x = pickle.load(f)
	return x

def safe_mkdir(path):
	'''
	Create a path and catch the race condition between path exists and mkdir.
	Input:
		* path: string, name of the folder to be created.
	'''
	path = os.path.abspath(path)
	if not os.path.exists(path):
		try:
			os.makedirs(path)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise

class file_manager(object):
	""" class which manages the cache. """
	def __init__(self, mode, exp, spec, lmax, force_recomputation=False, folder='cache',foldersim='sims',rank=0):
		"""
		Look into folder/ if the file you are trying to compute already exists on the disk.
		Input
			 * mode: string, quantity you want to compute.
			 * exp: string, name of the experimental configuration.
			 * spec: string, name of CMB spectrum you are looking at.
			 * lmax: int, maximum multipole to load (all multipoles in file will be loaded by default).
			 * force_recomputation: boolean, if True, recompute no matter the file already exists on the disk
			 * folder: string, name of the folder which contains the files
			 * rank: int, the rank of the processor calling this class (by default 0 = master).
		"""
		## Make a safe mkdir
		self.folder = folder
		safe_mkdir(self.folder)

		self.rank = rank
		self.exp = exp

		self.basefn = '%s_%s_%d_%s'%(mode,spec,lmax,exp)
		self.fn = os.path.join(self.folder,self.basefn)

		# self.foldersim = os.path.join(foldersim,self.basefn)
		# safe_mkdir(self.foldersim)

		if force_recomputation is False:
			self.FileExist = self.check_if_file_exist()
		else:
			self.FileExist = False

		if self.FileExist is True and force_recomputation is False:
			## Need to define the format of the data: dic, or just array? dic would be better
			self.data = self.load_data_from_disk()
		else:
			self.data = {}

	def check_if_file_exist(self):
		'''
		Check if the file corresponding to the experimental configuration exists on the disk.
		Output
			* boolean
		'''
		fn = self.fn + '.pkl'
		if os.path.isfile(fn):
			if self.rank==0: print 'Data stored on the disk'
			return True
		else:
			if self.rank==0: print 'Data not stored on the disk'
			return False

	def load_data_from_disk(self,fn=None):
		'''
		Load the file corresponding to the experimental configuration from the disk.
		Input:
			* (optional) fn: string, external file to load. Should be of the same format that save_data_on_disk() returns.
		Output:
			* array: 1D or 2D-array, the data stored in the file (vector or matrix)
		'''
		if fn is not None:
			if self.rank==0: print 'Loading external data',fn
			data = pickle_load(fn)
			array = data['data']
		elif self.FileExist is False:
			if self.rank==0: print 'Cannot load the data'
			array = []
		else:
			fn = self.fn + '.pkl'
			data = pickle_load(fn)
			array = data['data']
		return array

	def save_data_on_disk(self,array,name=None):
		'''
		Save the data into a file(.pkl) on the disk for later re-use.
		The name of the file is given by self.fn to have quick knowlegde of what's inside.
		For the moment, the routine is very basic and store only the array.
		Further informations could be stored.
		Input:
			* array: 1D or 2D-array, the data stored in the file (vector or matrix)
		'''
		data = {}

		## Ideally put more than just the array
		data['data'] = array

		if name is not None:
			fn = name
		else:
			fn = self.fn
			fn += '.pkl'
		pickle_save(data, fn)
		if self.rank==0: print 'data saved on disk at', fn

	def save_figure_on_disk(self,plot_folder,extension='.pdf',name=None):
		'''
		Save a figure into a file on the disk, and close the figure to avoid overplot.
		Input:
			* plot_folder: string, the folder which contain the plots
			* extension: string, the extension for the file (default = .pdf)
		'''
		if name is not None:
			fn = os.path.join(plot_folder,name)
		else:
			fn = os.path.join(plot_folder,self.basefn)
		fn += extension
		import pylab as pl
		pl.savefig(fn)
		pl.clf()
		if self.rank==0: print 'figure saved on disk at', fn

########################################################################
# Misc

def dot(*l):
	''' np.dot wrapper for dotting several linear operators together '''
	return reduce(np.dot,l)

def split_block_n0(block):
	'''
	Return name of spectra involve in N0 computation.
	'''
	x1 = block[0].upper()
	x2 = block[1].upper()
	x3 = block[2].upper()
	x4 = block[3].upper()
	flavor1 = 'cl%s%s'%(x1.lower(),x2.lower())
	flavor2 = 'cl%s%s'%(x3.lower(),x4.lower())
	flavor_n0 = 'cl%s'%block.lower()
	return flavor1, flavor2, flavor_n0
