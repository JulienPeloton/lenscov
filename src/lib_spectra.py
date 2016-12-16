# Copyright (C) 2016 Peloton
#########################
# Library for spectra
# author: julien@sussex
# Date of creation: 10/2015
# See 1611.01446
#########################
import glob, copy, sys, os
import numpy as np
import pylab as pl

import lib_covariances
from loop_lensing import loop_lensing
from util import binner, split_block_n0

from matplotlib import cm, gridspec
from scipy import interpolate

########################################################################
# Noise & beam

def bl(fwhm_arcmin, lmax):
	""" returns the map-level transfer function for a symmetric Gaussian beam.
	Input
		 * fwhm_arcmin: float, beam full-width-at-half-maximum (fwhm) in arcmin.
		 * lmax: int, maximum multipole.
	Output
		* bell: 1D array, the beam function
	"""
	ls = np.arange(0, lmax+1)
	return np.exp( -(fwhm_arcmin * np.pi/180./60.)**2 / (16.*np.log(2.)) * ls*(ls+1.) )

def nl(noise_uK_arcmin, fwhm_arcmin, lmax):
	""" returns the beam-deconvolved (white) noise power spectrum in units of uK^2
	Input
		 * noise_uK_arcmin: float, map noise level in uK.arcmin
		 * fwhm_arcmin: float, beam full-width-at-half-maximum (fwhm) in arcmin.
		 * lmax: int, maximum multipole.
	Output
		 * nell: 1D array, the noise power spectrum
	"""
	return (noise_uK_arcmin * np.pi/180./60.)**2 / bl(fwhm_arcmin, lmax)**2

def nlcorr(noise_uK_arcmin, fwhm_arcmin, lmax, TTcorrelatednoise=100.):
	""" returns the beam-deconvolved noise power spectrum in units of uK^2
	with a correlated component (low-ell contamination).
	Input
		 * noise_uK_arcmin: float, map noise level in uK.arcmin
		 * fwhm_arcmin: float, beam full-width-at-half-maximum (fwhm) in arcmin.
		 * lmax: int, maximum multipole.
		 * TTcorrelatednoise: cutoff for the correlations in ell space
	Output
		 * nell: 1D array, the noise power spectrum
	"""
	ls = np.arange(0, lmax+1)
	correlation = ( 1. + 1./(ls/float(TTcorrelatednoise))**.5 )
	correlation[0] = correlation[1]
	return correlation * (noise_uK_arcmin * np.pi/180./60.)**2 / bl(fwhm_arcmin, lmax)**2

def add_noise_to_cl(cl,noise_uK_arcmin,fwhm_arcmin,lmax,spec='',TTcutoff=False,TTcorr=False):
	""" returns the signal plus the beam-deconvolved noise power spectrum in units of uK^2
	Input
		 * cl: 1D array, signal in uK^2 (length lmax+1)
		 * noise_uK_arcmin: float, map noise level in uK.arcmin
		 * fwhm_arcmin: float, beam full-width-at-half-maximum (fwhm) in arcmin.
		 * lmax: int, maximum multipole.
		 * observable: string, T or P, or TP.
	Output
		 * cl_noisy: 1D array, the sum of signal and noise (length lmax+1)
	"""
	assert len(cl) == lmax+1
	if spec == '':
		print 'you need to specify a spec for noise'
		sys.exit()
	if spec == 'cltt':
		if TTcutoff:
			## Trick to throw away all contributions for ell>3000
			print 'You apply a cutoff at ell_cut=3000'
			return cl + nl(noise_uK_arcmin, fwhm_arcmin, lmax) + nl(1.e-6, 15, lmax)
		elif TTcorr:
			print 'You load 1/f noise'
			return cl + nlcorr(noise_uK_arcmin, fwhm_arcmin, lmax, TTcorrelatednoise=100.)
		else:
			return cl + nl(noise_uK_arcmin, fwhm_arcmin, lmax)
	elif spec in ['clee', 'clbb']:
		return cl + nl(noise_uK_arcmin*np.sqrt(2), fwhm_arcmin, lmax)
	elif spec == 'clte':
		return cl

########################################################################
# Spin values used for wigner

def load_spin_values_wigner(flavor):
	'''
	Load spin values used for Wigner computation (see loop_lensing.f90).
	'''
	xfield = flavor[-2]
	yfield = flavor[-1]

	if xfield == 't':
		spinl2_x=0; spinl3_x=0
	elif xfield in ['e','b']:
		spinl2_x=0; spinl3_x=2
	if yfield == 't':
		spinl2_y=0; spinl3_y=0
	elif yfield in ['e','b']:
		spinl2_y=0; spinl3_y=2

	return spinl2_x, spinl3_x, spinl2_y, spinl3_y

########################################################################
# Weights for the quadratic estimator

def load_weights(cls_class,flavor, noise_uK_arcmin, fwhm_arcmin, lmax, extra=''):
	'''
	Load the weights used for N0 computation.
	Input:
		* cls_class: camb_clfile object, containing lensed power spectra.
		* flavor: string, name of the spectrum as defined in camb_clfile.
		* noise_uK_arcmin: float, white noise level to be used (in uK.arcmin).
		* fwhm_arcmin: float, beam FWHM to be used.
		* lmax: int, the maximum multipole. Currently accepted: lmax, 2*lmax (extra=_long)
	Output
		* cl_XX: 2D-array, contain cl_xx and cl_xx_noisy.
	'''
	xfield = flavor[-2]
	yfield = flavor[-1]

	cl_lensed_XX = getattr(cls_class, 'cl%s%s%s'%(xfield,xfield,extra))
	cl_lensed_YY = getattr(cls_class, 'cl%s%s%s'%(yfield,yfield,extra))
	if flavor == 'cltb':
		cl_lensed_XY = getattr(cls_class, 'clte%s'%extra) ## TB gets TE for weight Xspec
	elif flavor == 'cleb':
		cl_lensed_XY = np.zeros_like(cl_lensed_XX)
	else:
		cl_lensed_XY = getattr(cls_class, 'cl%s%s%s'%(xfield,yfield,extra))

	if flavor =='cltb':
		cl_lensed_XX_noisy = add_noise_to_cl(cl_lensed_XX,noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,lmax=lmax,spec='cltt')
		cl_lensed_YY_noisy = add_noise_to_cl(cl_lensed_YY,noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,lmax=lmax,spec='clbb')
		cl_lensed_XY_noisy = np.zeros_like(cl_lensed_XY)
	elif flavor =='cleb':
		cl_lensed_XX_noisy = add_noise_to_cl(cl_lensed_XX,noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,lmax=lmax,spec='clee')
		cl_lensed_YY_noisy = add_noise_to_cl(cl_lensed_YY,noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,lmax=lmax,spec='clbb')
		cl_lensed_XY_noisy = np.zeros_like(cl_lensed_XY)
	elif flavor=='clte':
		cl_lensed_XY_noisy = cl_lensed_XY
		cl_lensed_XX_noisy = add_noise_to_cl(cl_lensed_XX,noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,lmax=lmax,spec='cltt')
		cl_lensed_YY_noisy = add_noise_to_cl(cl_lensed_YY,noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,lmax=lmax,spec='clee')
	else:
		cl_lensed_XY_noisy = add_noise_to_cl(cl_lensed_XY,noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,lmax=lmax,spec=flavor)
		cl_lensed_XX_noisy = add_noise_to_cl(cl_lensed_XX,noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,lmax=lmax,spec=flavor)
		cl_lensed_YY_noisy = add_noise_to_cl(cl_lensed_YY,noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,lmax=lmax,spec=flavor)

	cl_XX = np.zeros((2,lmax+1),dtype=float)
	cl_YY = np.zeros((2,lmax+1),dtype=float)
	cl_XY = np.zeros((2,lmax+1),dtype=float)

	cl_XX[0] = cl_lensed_XX
	cl_XX[1] = cl_lensed_XX_noisy

	cl_YY[0] = cl_lensed_YY
	cl_YY[1] = cl_lensed_YY_noisy

	cl_XY[0] = cl_lensed_XY
	cl_XY[1] = cl_lensed_XY_noisy

	return cl_XX, cl_YY, cl_XY

########################################################################
# Half the mean-squared deflection

def compute_R(cls_unlensed):
	'''
	Half the mean-squared deflection (Eq 63 from astro-ph/0001303)
	Input:
		cls_unlensed: camb_clfile object
	Output
		R: float, half the mean-squared deflection
	'''
	ells = cls_unlensed.ls
	R = 1./(8*np.pi)*np.sum(cls_unlensed.clpp*(2*ells+1)*(ells+1)*ells)
	return R

########################################################################
# Compute biases and minimum variance quantities

def compute_N0_XYWZ(cls_lensed,lmin=2,blocks=['TTTT'],noise_uK_arcmin=0.0,fwhm_arcmin=0.0,MPI=None):
	'''
	Compute the gaussian bias (N0) for a given experimental setup.
	Input:
		* cls_lensed: camb_clfile object, containing lensed power spectra.
		* lmin: int, minimum multipole.
		* block: string, 4-pt name (e.g. TTTT, or EBEB).
		* noise_uK_arcmin: float, white noise level to be used (in uK.arcmin).
		* fwhm_arcmin: float, beam FWHM to be used.
		* MPI: module, module to pass if you want to parallelize the computation. If None, the computation is done serially.
	Output
		* N0: 1D array (length lmax+1)
	'''
	try:
		comm = MPI.COMM_WORLD
		rank = comm.rank
	except:
		## No MPI
		comm = None
		rank = 0

	lmax = int(cls_lensed.lmax)
	N0_tot = np.array([np.zeros(lmax+1) for i in range(len(blocks))])

	#############################################################
	# Loop over blocks
	#############################################################
	for position_block,block in enumerate(blocks):
		flavor1, flavor2,flavor_n0 = split_block_n0(block)

		if rank == 0: print 'Doing block %s (%s, %s, %s)\n'%(block,flavor1,flavor2, flavor_n0)

		## Load the weights and spin values used in N0 computation
		cl_XX, cl_YY, cl_XY = load_weights(cls_lensed, flavor1, noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		spinl2_x, spinl3_x, spinl2_y, spinl3_y = load_spin_values_wigner(flavor1)
		if block == 'TTEE':
			cl_XX, junk1, junk2 = load_weights(cls_lensed, 'cltt', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
			cl_junk1, cl_YY, junk2 = load_weights(cls_lensed, 'clee', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
			junk1, junk2, cl_XY = load_weights(cls_lensed, 'clte', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
			spinl2_x, spinl3_x, spinl2_y, spinl3_y = load_spin_values_wigner('clte')
		if block in ['TTTE','EETE']:
			cl_XX, cl_YY, cl_XY = load_weights(cls_lensed, 'clte', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
			spinl2_x, spinl3_x, spinl2_y, spinl3_y = load_spin_values_wigner('clte')
		if block == 'EBTB':
			## clxx is EE, clyy is BB, clxy is clte**2/(clee*clbb)
			cl_XX, cl_YY, cl_XY = load_weights(cls_lensed, 'cleb', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
			TT, EE, TE = load_weights(cls_lensed, 'clte', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
			cl_XY[0] = TE[0]**2 / (TT[0]*EE[0])
			cl_XY[1] = TE[1]**2 / (TT[1]*EE[1])
			spinl2_x, spinl3_x, spinl2_y, spinl3_y = load_spin_values_wigner('cleb')

		if comm is None:
			l2range = range(lmin, lmax+1, 1)
			l2dim = lmax - lmin

			## /!\ compute_bias_zero_xyxy_mpi just returns sum(g_ell f_ell)
			tmp = loop_lensing.compute_bias_zero_xyxy_mpi(cl_XX,cl_YY,cl_XY,l2range,flavor_n0,lmin,spinl2_x,spinl3_x,spinl2_y,spinl3_y,l2dim,lmax)
			N0_tot[position_block] = (2*cls_lensed.ls + 1.0)/tmp
		else:
			n_tot = comm.size
			l2range = range(lmin+rank,lmax+1,n_tot)
			l2dim = len(l2range)-1

			## /!\ compute_bias_zero_xyxy_mpi returns sum(g_ell f_ell)
			tmp = loop_lensing.compute_bias_zero_xyxy_mpi(cl_XX,cl_YY,cl_XY,l2range,flavor_n0,lmin,spinl2_x,spinl3_x,spinl2_y,spinl3_y,l2dim,lmax)
			comm.Reduce([tmp, MPI.DOUBLE],[N0_tot[position_block], MPI.DOUBLE],op = MPI.SUM,root = 0)
			comm.barrier()

			if rank==0:
				if block == 'TTEE':
					index_TT = blocks.index('TTTT')
					index_EE = blocks.index('EEEE')
					N0_tot[position_block] = N0_tot[index_TT] * N0_tot[index_EE] / (2*cls_lensed.ls + 1.0) * N0_tot[position_block]
				elif block == 'TTTE':
					index_TT = blocks.index('TTTT')
					index_TE = blocks.index('TETE')
					N0_tot[position_block] = N0_tot[index_TT] * N0_tot[index_TE] / (2*cls_lensed.ls + 1.0) * N0_tot[position_block]
				elif block == 'EETE':
					index_EE = blocks.index('EEEE')
					index_TE = blocks.index('TETE')
					N0_tot[position_block] = N0_tot[index_EE] * N0_tot[index_TE] / (2*cls_lensed.ls + 1.0) * N0_tot[position_block]
				elif block == 'EBTB':
					index_EB = blocks.index('EBEB')
					index_TB = blocks.index('TBTB')
					N0_tot[position_block] = N0_tot[index_EB] * N0_tot[index_TB] / (2*cls_lensed.ls + 1.0) * N0_tot[position_block]
				else:
					N0_tot[position_block] = (2*cls_lensed.ls + 1.0)/N0_tot[position_block]

	return N0_tot, blocks

def compute_minimum_variance_weights(N0_array,N0_names,checkit=False):
	'''
	Compute the variance of the minimum variance estimator and the associated weights.
	Input:
		* N0_array: ndarray, contain the N0s to combine
		* N0_names: ndarray of string, contain the name of the N0s to combine (['TTTT', 'EEEE', etc.])
	Output:
		* minimum_variance_n0: 1D array, the MV N0
		* weights*minimum_variance_n0: ndarray, the weights for each spectrum (TT, EE, etc.)
		* N0_names_ordered: 1D array, contain the name of the spectra (TT, EE, etc.)
	'''
	n_elements = len(N0_array[0])
	n_estimator = len([i for i in N0_names if i[0:2]==i[2:]]) ## Take the diag
	N0_array = np.nan_to_num(N0_array)

	## Fill matrix
	N0_names_ordered = ['TT','EE','EB','TE','TB','BB']
	sub_vec = [[name,pos] for pos,name in enumerate(N0_names_ordered)]
	dic_mat = {'%s%s'%(XY,ZW):[i,j] for XY,i in sub_vec for ZW,j in sub_vec}

	## Build the inverse matrix for each ell
	def build_inv_submatrix(vector_ordered,names_ordered,dic_mat):
		mat = np.zeros((6,6))
		for pos,name in enumerate(names_ordered):
			mat[dic_mat[name][0]][dic_mat[name][1]] = mat[dic_mat[name][1]][dic_mat[name][0]] = vector_ordered[pos]
		return np.linalg.pinv(mat)

	inv_submat_array = np.array([build_inv_submatrix(vector_ordered,N0_names,dic_mat) for vector_ordered in np.transpose(N0_array)])
	inv_N0_array = np.array([ np.sum(submat) for submat in inv_submat_array ])
	minimum_variance_n0 = 1./inv_N0_array

	weights = np.array([ [np.sum(submat[i]) for submat in inv_submat_array] for i in range(6) ])

	if checkit:
		print np.sum(weights*minimum_variance_n0)/len(minimum_variance_n0) == 1.
		print np.sum(weights*minimum_variance_n0)/len(minimum_variance_n0)

	return minimum_variance_n0, weights*minimum_variance_n0, N0_names_ordered

def compute_minimum_variance_N1(bin_centers,bin_function,N1,N0_array,N0_names):
	'''
	Takes all N1 and form the mimimum variance estimator.
	Assumes N1 structure is coming from Biases_n1mat.f90
	Input:
		* bin_centers: 1D array, contain the bin centers for the output
		* bin_function: function, function to bin the input into the output
		* N1: ndarray, contain the N1 (output of Biases_n1mat.f90)
		* N0_array: ndarray, contain the N0s
		* N0_names: ndarray of string, contain the name of the N0s (['TTTT', 'EEEE', etc.])
	'''
	n0_array = [N0_array[N0_names.index(i)] for i in N0_names]
	minimum_variance_n0, weights_for_MV, names_for_weights = compute_minimum_variance_weights(n0_array,N0_names)

	step = 6
	n1_ell = N1[0]
	n1_mat = np.reshape(N1[1:],(step,step,len(n1_ell)))

	## Ordering: i_TT=0,i_EE=1,i_EB=2,i_TE=3,i_TB=4, i_BB=5 (from Biases.f90)
	names_N1 = [  'TTTT','TTEE', 'TTTE', 'TTBB', 'TTEB','TTTB',
				  'EETT','EEEE', 'EETE', 'EEBB', 'EEEB','EETB',
				  'TETT','TEEE', 'TETE', 'TEBB', 'TEEB','TETB',
				  'BBTT','BBEE', 'BBTE', 'BBBB', 'BBEB','BBTB',
				  'EBTT','EBEE', 'EBTE', 'EBBB', 'EBEB','EBTB',
				  'TBTT','TBEE', 'TBTE', 'TBBB', 'TBEB','TBTB']
	indices = {'TT':0,'EE':1,'EB':2,'TE':3,'TB':4,'BB':5}

	n1_tot = np.zeros_like(bin_centers)

	for estimator_name in names_N1:
		## Indices for arrays
		index_x = indices[estimator_name[0:2]]
		index_y = indices[estimator_name[2:]]

		## Interpolate N1
		n1_not_interp = n1_mat[index_x][index_y]
		n1_interp = np.interp(bin_centers,n1_ell,n1_not_interp)

		## Weights
		wXY_index = names_for_weights.index(estimator_name[0:2])
		wZW_index = names_for_weights.index(estimator_name[2:4])

		## Update N1
		if bin_function is not None:
			n1_tot += bin_function(weights_for_MV[wXY_index]) * bin_function(weights_for_MV[wZW_index]) * n1_interp
		else:
			n1_tot += weights_for_MV[wXY_index] * weights_for_MV[wZW_index] * n1_interp
	return n1_tot

def compute_minimum_variance_N1_mat(bin_object,path_to_N1matrix,N0_array,N0_names):
	'''
	Compute the N1 matrix to perform N1 deconvolution (see App. C in the paper).
	Assumes N1 structure is coming from Biases.f90
	Input:
		* bin_object: object, contain binning routines, bin centers, etc.
		* path_to_N1matrix: strig, where the output of Biases_n1mat are stored
		* N0_array: ndarray, contain the N0s
		* N0_names: ndarray of string, contain the name of the N0s (['TTTT', 'EEEE', etc.])
	Output:
		* dN1_over_dCpp_tot * bin_object.bin_size: 2D array, the N1 MV matrix
	'''
	n0_array = [N0_array[N0_names.index(i)] for i in N0_names]
	minimum_variance_n0, weights_for_MV, names_for_weights = compute_minimum_variance_weights(n0_array,N0_names)

	## Ordering: i_TT=0,i_EE=1,i_EB=2,i_TE=3,i_TB=4, i_BB=5 (from Biases.f90)
	names_N1 = [  'TTTT','TTEE', 'TTTE', 'TTBB', 'TTEB','TTTB',
				  'EETT','EEEE', 'EETE', 'EEBB', 'EEEB','EETB',
				  'TETT','TEEE', 'TETE', 'TEBB', 'TEEB','TETB',
				  'BBTT','BBEE', 'BBTE', 'BBBB', 'BBEB','BBTB',
				  'EBTT','EBEE', 'EBTE', 'EBBB', 'EBEB','EBTB',
				  'TBTT','TBEE', 'TBTE', 'TBBB', 'TBEB','TBTB']
	indices = {'TT':0,'EE':1,'EB':2,'TE':3,'TB':4,'BB':5}

	nbins = len(bin_object.bin_centers)
	dN1_over_dCpp_tot = np.zeros((nbins,nbins))

	for estimator_name in names_N1:
		## Load corresponding dN1_over_dCpp matrix
		tag = '_'+path_to_N1matrix.split('/')[1]

		# N1 matrix from code and samplings
		try:
			Lout, Ls, dLs, dN1_over_dCpp_binned = lib_covariances.readN1_matrix(ntag='All',atag=tag, matrixtag = estimator_name, path_to_N1matrix=path_to_N1matrix)
		except:
			Lout, Ls, dLs, dN1_over_dCpp_binned = lib_covariances.readN1_matrix(ntag='All',atag=tag, matrixtag = estimator_name[2:]+estimator_name[:2], path_to_N1matrix=path_to_N1matrix)

		# dN1_over_dCpp_binned = np.array([[dN1_over_dCpp_binned[list(Lout).index(j)][list(Ls).index(i)] * (i*(i+1.))**2 / (j*(j+1.))**2 for i in Ls] for j in Lout])
		dN1_over_dCpp_binned = np.nan_to_num(dN1_over_dCpp_binned)

		for i,_ in enumerate(Lout):
			dN1_over_dCpp_binned[i,:] /= dLs

		## Interpolate
		matsp = interpolate.RectBivariateSpline(Lout, Ls, dN1_over_dCpp_binned,kx=5,ky=5)
		dN1_over_dCpp = matsp(bin_object.bin_centers,bin_object.bin_centers, grid = True)

		## Weights
		wXY_index = names_for_weights.index(estimator_name[0:2])
		wZW_index = names_for_weights.index(estimator_name[2:4])

		## Update N1
		dN1_over_dCpp_tot += lib_covariances.vecmat(bin_object.bin_this_spec(weights_for_MV[wXY_index]) * bin_object.bin_this_spec(weights_for_MV[wZW_index]), dN1_over_dCpp)
		# dN1_over_dCpp_tot += lib_covariances.vecmat(weights_for_MV[wXY_index] * weights_for_MV[wZW_index], dN1_over_dCpp)

	## N1 matrix operate on and return (L(L+1)]^2C^pp /2pi normalized quantities, so we need to undo this.
	dN1_over_dCpp_tot = np.array([[dN1_over_dCpp_tot[list(bin_object.bin_centers).index(j)][list(bin_object.bin_centers).index(i)] * (i*(i+1.))**2 / (j*(j+1.))**2 for i in bin_object.bin_centers] for j in bin_object.bin_centers])
	return dN1_over_dCpp_tot * bin_object.bin_size

def compute_N1(cls_unlensed,cls_lensed,lmin=2,flavor='cltt',noise_uK_arcmin=0.0,fwhm_arcmin=0.0,MPI=None):
	'''
	!!!!!!!!!!!!!!!!!!!!!
	!!!! Not working !!!!
	!!!!!!!!!!!!!!!!!!!!!
	Compute the gaussian bias (N0) for a given experimental setup.
	The routine is parallelized in a dump way (and within python), by splitting evenly the ells of the first sum
	within processors. For ell > 3000, better to use many procs (few minutes on 1 core).
	/!\ spin conventions are hardcoded in order to fit with fortran routines (see loop_lensing.f90).
	If you modify the fortran routines, the values of spin can change. /!\
	/!\ Due to parallelization the fortran routine returns only sum(g_ell f_ell), then the coaddition
	and the inversion are done in python. /!\
	Input:
		* cls_lensed: camb_clfile object, containing lensed power spectra.
		* lmin: int, minimum multipole.
		* flavor: string, name of the spectrum as defined in camb_clfile.
		* noise_uK_arcmin: float, white noise level to be used (in uK.arcmin).
		* fwhm_arcmin: float, beam FWHM to be used.
		* MPI: module, module to pass if you want to parallelize the computation. If None, the computation is done serially.
	Output
		* N0: 1D array (length lmax+1)
	'''
	pass

########################################################################
# Functions and class to manipulate power spectra

def get_camb_cls(fname=None, prefix=None, lmax=None):
	""" loads and returns a "scalar Cls" file produced by CAMB (camb.info).
	Input
		 * fname: string, file name to load.
		 * prefix: string, directory in quicklens/data/cl directory to pull the *_scalCls.dat from. defaults to 'planck_wp_highL' (only used if fname==None).
		 * lmax: int, maximum multipole to load (all multipoles in file will be loaded by default).
	Output
		 * cls: camb_clfile object, object containing everything about input cls.
	"""
	tf = glob.glob( fname )
	assert(len(tf) == 1)

	return camb_clfile( tf[0], lmax=lmax )

def is_camb_clfile(obj):
	""" check (by ducktyping) if obj is a Cls file produced by CAMB """
	if not hasattr(obj, 'lmax'):
		return False
	if not hasattr(obj, 'ls'):
		return False
	return set(object.__dict__.keys()).issubset( set( ['lmax', 'ls', 'cltt', 'clee', 'clte', 'clpp', 'cltp', 'clbb', 'clep', 'cleb' ] ) )

class camb_clfile(object):
	""" class to load and store Cls from the output files produced by CAMB. """
	def __init__(self, tfname, lmax=None):
		"""
		load Cls.
		Input
			 * tfname: string, file name to load from.
			 * lmax: int, maximum multipole to load (all multipoles in file will be loaded by default).
		"""

		tarray = np.loadtxt(tfname)
		lmin   = int(tarray[0, 0])
		assert(lmin in [1,2])

		if lmax == None:
			lmax = np.shape(tarray)[0]-lmin+1
			if lmax > 10000:
				lmax = 10000
			else:
				assert(tarray[-1, 0] == lmax)

		ncol = np.shape(tarray)[1]
		ell  = np.arange(lmin, lmax+1, dtype=np.float)

		## interpolate values beyond 10000
		if 3*lmax > 10000:
			print 'Interpolate high ell values for spectra'
			lmax_file = tarray[-1, 0]
			tarray_tmp = np.zeros((ncol,lmax_file-lmin+1))
			full_lrange = range(int(lmin),int(lmax_file+1))
			tarray = tarray.T ## transpose to ease the computation
			for i in range(ncol):
				tarray_tmp[i] = np.interp(full_lrange, tarray[0],tarray[i])
			tarray = tarray_tmp
			tarray = tarray.T ## restore convention

		assert( (np.shape(tarray)[0]+1) >= lmax )


		## Double the lmax in order to be able to have all wigner computation correct (l1<=l2+l3)
		lmax_spec = 2*lmax
		ell_spec  = np.arange(lmin, lmax_spec+1, dtype=np.float)

		self.lmax = lmax
		self.ls   = np.concatenate( [ np.arange(0, lmin), ell ] )
		self.ls_long   = np.concatenate( [ np.arange(0, lmin), ell_spec ] )
		if ncol == 5:																			# _lensedCls
			self.cltt = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax-lmin+1),1]*2.*np.pi/ell/(ell+1.) ] )
			self.clee = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax-lmin+1),2]*2.*np.pi/ell/(ell+1.) ] )
			self.clbb = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax-lmin+1),3]*2.*np.pi/ell/(ell+1.) ] )
			self.clte = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax-lmin+1),4]*2.*np.pi/ell/(ell+1.) ] )

			self.cltt_long = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax_spec-lmin+1),1]*2.*np.pi/ell_spec/(ell_spec+1.) ] )
			self.clee_long = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax_spec-lmin+1),2]*2.*np.pi/ell_spec/(ell_spec+1.) ] )
			self.clbb_long = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax_spec-lmin+1),3]*2.*np.pi/ell_spec/(ell_spec+1.) ] )
			self.clte_long = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax_spec-lmin+1),4]*2.*np.pi/ell_spec/(ell_spec+1.) ] )

		elif ncol == 6:																		  # _scalCls
			tcmb  = 2.726*1e6 #uK
			self.tcmb  = 2.726*1e6 #uK

			self.cltt = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax-lmin+1),1]*2.*np.pi/ell/(ell+1.) ] )
			self.clee = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax-lmin+1),2]*2.*np.pi/ell/(ell+1.) ] )
			self.clte = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax-lmin+1),3]*2.*np.pi/ell/(ell+1.) ] )
			self.clbb = np.concatenate( [ np.zeros(lmin), np.zeros_like(tarray[0:(lmax-lmin+1),3]) ] )
			self.clpp = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax-lmin+1),4]/ell**4/tcmb**2 ] )
			self.cltp = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax-lmin+1),5]/ell**3/tcmb ] )

			self.cltt_long = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax_spec-lmin+1),1]*2.*np.pi/ell_spec/(ell_spec+1.) ] )
			self.clee_long = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax_spec-lmin+1),2]*2.*np.pi/ell_spec/(ell_spec+1.) ] )
			self.clte_long = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax_spec-lmin+1),3]*2.*np.pi/ell_spec/(ell_spec+1.) ] )
			self.clbb_long = np.concatenate( [ np.zeros(lmin), np.zeros_like(tarray[0:(lmax_spec-lmin+1),3]) ] )
			self.clpp_long = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax_spec-lmin+1),4]/ell_spec**4/tcmb**2 ] )
			self.cltp_long = np.concatenate( [ np.zeros(lmin), tarray[0:(lmax_spec-lmin+1),5]/ell_spec**3/tcmb ] )


	def copy(self, lmax=None, lmin=None):
		""" clone this object.
		Input
			 * lmax: int, restrict copy to L<=lmax.
			 * lmin: int, set spectra in copy to zero for L<lmin.
		Output
			 * ret: camb_clfile object, the copy
		"""
		if (lmax == None):
			return copy.deepcopy(self)
		else:
			assert( lmax <= self.lmax )
			ret	  = copy.deepcopy(self)
			ret.lmax = lmax
			ret.ls   = np.arange(0, lmax+1)
			for k, v in self.__dict__.items():
				if k[0:2] == 'cl':
					setattr( ret, k, copy.deepcopy(v[0:lmax+1]) )

			if lmin != None:
				assert( lmin <= lmax )
				for k in self.__dict__.keys():
					if k[0:2] == 'cl':
						getattr( ret, k )[0:lmin] = 0.0
			return ret

	def hashdict(self):
		""" return a dictionary uniquely associated with the contents of this clfile. """
		ret = {}
		for attr in ['lmax', 'cltt', 'clee', 'clte', 'clpp', 'cltp', 'clbb', 'clep', 'cleb' ]:
			if hasattr(self, attr):
				ret[attr] = getattr(self, attr)
		return ret

	def __add__(self, other):
		"""
		Sum two clfile objects.
		lmax must be the same for both.
		Spectra which are not common between the two are ignored.
		Input
			* other: camb_clfile object, the cls that you want to add
		"""
		if is_camb_clfile(other):
			assert( self.lmax == other.lmax )
			ret = self.copy()
			zs  = np.zeros(self.lmax+1)
			for attr in ['cltt', 'clee', 'clte', 'clpp', 'cltp', 'clbb', 'clep', 'cleb' ]:
				if (hasattr(self, attr) or hasattr(other, attr)):
					setattr(ret, attr, getattr(self, attr, zs) + getattr(other, attr, zs) )
			return ret
		else:
			assert(0)

	def __eq__(self, other):
		""" compare two clfile objects. """
		try:
			for key in self.__dict__.keys()+other.__dict__.keys():
				if type(self.__dict__[key]) == np.ndarray:
					assert( np.all( self.__dict__[key] == other.__dict__[key] ) )
				else:
					assert( self.__dict__[key] == other.__dict__[key] )
		except:
			return False

		return True

	def plot(self, spec='cltt', p=pl.plot, t=lambda l:1., **kwargs):
		"""
		Plot the spectrum
		Input
			 * spec: string, name of the spectrum to display (e.g. cltt, clee, clte, etc.)
			 * p: function, plotting function to use p(x,y,**kwargs)
			 * t: function, scaling to apply to the plotted C_ell -> t(ell)*C_ell
		"""
		p( self.ls, t(self.ls) * getattr(self, spec), **kwargs )

########################################################################
# Old and ugly (probably not working)

def analytic_lensed_spectra_XX_fortran(cls_unlensed,lmin=2,flavor='cltt',MPI=None):
	'''
	Compute (perturbatively) lensed spectrum given an unlensed spectrum (in the form XX, where X= T, E, or B).
	Use Eq. 62 & 76 from astro-ph/0001303 (Hu)
	The routine returns separately the two leading terms.
	The routine is parallelized in a dump way (and within python), by splitting evenly the ells of the first sum
	within processors. For ell > 3000, better to use many procs (few minutes on 1 core).
	/!\ spin conventions are hardcoded in order to fit with fortran routines (see loop_lensing.f90).
	If you modify the fortran routines, the values of spin can change. /!\
	Input:
		* cls_unlensed: camb_clfile object, containing unlensed power spectra.
		* lmin: int, minimum multipole.
		* flavor: string, name of the spectrum as defined in camb_clfile.
		* MPI: module, module to pass if you want to parallelize the computation. If None, the computation is done serially.
	Output
		* C_l1_analytic_part1: 1D array (length lmax+1)
		* C_l2_analytic_part1: 1D array (length lmax+1)
	'''
	try:
		comm = MPI.COMM_WORLD
		rank = comm.rank
	except:
		## No MPI
		comm = None
		rank = 0

	R = compute_R(cls_unlensed)

	if flavor == 'cltt':
		C_l1_analytic_part1 = np.array( [ cls_unlensed.cltt[i] * (1 - i*(i+1)*R) for i in cls_unlensed.ls] )

		cl = np.zeros((2,2*cls_unlensed.lmax+1),dtype=float)
		cl[0] = cls_unlensed.cltt_long
		cl[1] = np.zeros_like(cls_unlensed.cltt_long)
		spinl2=0; spinl3=0

	elif flavor == 'clee':
		try:
			cls_unlensed.clbb
		except:
			if rank==0:
				print 'Assume no scalar B-modes'
			cls_unlensed.clbb=np.zeros_like(cls_unlensed.clee)
			cls_unlensed.clbb_long=np.zeros_like(cls_unlensed.clee_long)

		C_l1_analytic_part1 = np.array( [ cls_unlensed.clee[i] * (1 - (i**2 + i - 4)*R) for i in cls_unlensed.ls] )

		cl = np.zeros((2,2*cls_unlensed.lmax+1),dtype=float)
		cl[0] = (cls_unlensed.clee_long + cls_unlensed.clbb_long)/2.
		cl[1] = (cls_unlensed.clee_long - cls_unlensed.clbb_long)/2.
		spinl2=2; spinl3=0

	elif flavor == 'clte':
		C_l1_analytic_part1 = np.array( [ cls_unlensed.clte[i] * (1 - (i**2 + i - 2)*R) for i in cls_unlensed.ls] )

		cl = np.zeros((2,2*cls_unlensed.lmax+1),dtype=float)
		cl[0] = cls_unlensed.clte_long
		cl[1] = np.zeros_like(cls_unlensed.clte_long)
		spinl2_0=0; spinl3_0=0
		spinl2_2=2; spinl3_2=0

	elif flavor == 'clbb':
		try:
			cls_unlensed.clbb
		except:
			print 'Assume no scalar B-modes'
			cls_unlensed.clbb_long=np.zeros_like(cls_unlensed.clee_long)
			cls_unlensed.clbb=np.zeros_like(cls_unlensed.clee)

		C_l1_analytic_part1 = np.array( [ cls_unlensed.clbb[i] * (1 - (i**2 + i - 4)*R) for i in cls_unlensed.ls] )

		cl = np.zeros((2,2*cls_unlensed.lmax+1),dtype=float)
		cl[0] = (cls_unlensed.clee_long + cls_unlensed.clbb_long)/2.
		cl[1] = -(cls_unlensed.clee_long - cls_unlensed.clbb_long)/2.
		spinl2=2; spinl3=0

	clpp = cls_unlensed.clpp_long
	lmax = int(cls_unlensed.lmax)
	C_l1_analytic_part2_tot = np.zeros_like(C_l1_analytic_part1)

	if comm is None:
		l2range = range(lmin, lmax+1, 1)
		l2dim = lmax - lmin
		if flavor != 'clte':
			C_l1_analytic_part2=loop_lensing.compute_lensed_spectra_mpi(cl,clpp,l2range,flavor,lmin,spinl2,spinl3,l2dim,lmax)
		else:
			C_l1_analytic_part2=loop_lensing.compute_lensed_spectra_xy_mpi(cl,clpp,l2range,flavor,lmin,spinl2_0,spinl3_0,spinl2_2,spinl3_2,l2dim,lmax)
		return C_l1_analytic_part1 , C_l1_analytic_part2
	else:
		n_tot = comm.size
		l2range = range(lmin+rank,lmax+1,n_tot)
		l2dim = len(l2range)-1
		if flavor != 'clte':
			C_l1_analytic_part2=loop_lensing.compute_lensed_spectra_mpi(cl,clpp,l2range,flavor,lmin,spinl2,spinl3,l2dim,lmax)
		else:
			C_l1_analytic_part2=loop_lensing.compute_lensed_spectra_xy_mpi(cl,clpp,l2range,flavor,lmin,spinl2_0,spinl3_0,spinl2_2,spinl3_2,l2dim,lmax)
		comm.Reduce([C_l1_analytic_part2, MPI.DOUBLE],[C_l1_analytic_part2_tot, MPI.DOUBLE],op = MPI.SUM,root = 0)
		if rank==0:
			return C_l1_analytic_part1 , C_l1_analytic_part2_tot
		else:
			return 0,0

def compute_fisher_errorbars(clbinned,lbins,binsize,noise_uK_arcmin,fwhm_arcmin,fsky):
	'''
	Return simple analytic errorbars.
	'''
	prefac = np.sqrt( 2. / (2.*lbins + 1.) / binsize / fsky)
	return prefac * clbinned
