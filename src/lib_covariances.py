# Copyright (C) 2016 Peloton
#########################
# Library for covariances
# author: julien@sussex
# Date of creation: 10/2015
# See 1611.01446
#########################
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pylab as pl
pl.ioff()
import sys, copy, os, glob

from scipy.weave import inline, converters
from scipy import interpolate

from loop_lensing import loop_lensing
from correlation_functions import correlation_functions
import lib_spectra, util

########################################################################
# Spit out combinations X,Y,Z,W,U,V for different terms

def combination_trispectrumB(block):
	'''
	Define the combination of fields X,Y,Z,W and U,V to be used
	when computing the trispectrum type B contribution to the
	cross-covariance. This is based on Eq. 37 in 1611.01446
	Input:
		* block: string, the block XYZW_UV
	Output:
		* list_to_return, list of lists of strings, contain names for
			- C_ell^{XU}
			- Amplitude A_ell^{YZ}
			- product of weights g^{YZ}*f^{WV}
	'''
	X,Y,Z,W = [i for i in block.split('_')[0]]
	U,V = [i for i in block.split('_')[1]]

	lists = [ [X+U, X+Y, X+Y+V+Y],
			  [Y+U, Y+X, Y+X+V+X],
			  [Z+U, Z+W, Z+W+V+W],
			  [W+U, W+Z, W+Z+V+Z],
			  [X+V, X+Y, X+Y+U+Y],
			  [Y+V, Y+X, Y+X+U+X],
			  [Z+V, Z+W, Z+W+U+W],
			  [W+V, W+Z, W+Z+U+Z],]

	nonzero_spectra = ['TT', 'EE', 'BB', 'TE']
	transposed_spectra = ['ET', 'BE', 'BT']
	list_to_return = []
	for sublist in lists:
		cl_name, amplitude_name, mat_name = sublist

		## find the correct ordering
		if cl_name in transposed_spectra:
			cl_name = cl_name[1]+cl_name[0]
		if amplitude_name in transposed_spectra:
			amplitude_name = amplitude_name[1]+amplitude_name[0]

		## if EB or TB
		if cl_name not in nonzero_spectra:
			continue

		n0_name = amplitude_name+amplitude_name
		list_to_return.append([cl_name, n0_name, mat_name])

	return list_to_return

def combination_noise_signal_trispA(block):
	'''
	Define the combination of fields X,Y,Z,W and U,V to be used
	when computing the noise, signal and trispectrum type A
	contributions to the cross-covariance.
	This is based on Eqs. A11, 36, A12 in 1611.01446
	Input:
		* block: string, the block XYZW_UV
	Output:
		* noise, list of strings, contain names for the noise term
		* trispA, list of strings, contain names for the trispA term
		* signal, string, contain names for the signal term
	'''
	CMB = ['TT', 'EE', 'TE', 'BB']
	noise = [ [block[0:5]+i, i+block[5:]] for i in CMB]
	trispA = [ [block[0:5]+i, i+block[5:]] for i in CMB]
	signal = block[5:]+block[5:]
	return noise, trispA, signal

def combination_disconnected_8pt(block):
	'''
	Define the combination to use for the disconnected 8pt function
	contribution to the lensing potential power spectrum auto-covariance.
	see Eq. A19 in 1611.01446
	Input:
		* block: string, the block XYZW_UV
	Output:
		* list_to_return, list of lists of strings
	'''
	X,Y,Z,W = [i for i in block.split('_')[0]]
	Xp,Yp,Zp,Wp = [i for i in block.split('_')[1]]
	derivativeXYZW = '%s%s%s%s'%(X,Y,Z,W)
	derivativeXpYpZpWp = '%s%s%s%s'%(Xp,Yp,Zp,Wp)
	list_tot = []
	for x,y in [[X,Y],[Y,X]]:
		for z,w in [[Z,W],[W,Z]]:
			for xp,yp in [[Xp,Yp],[Yp,Xp]]:
				for zp,wp in [[Zp,Wp],[Wp,Zp]]:
					for xf,yf,zf,wf,xpf,ypf,zpf,wpf in [[x,y,z,w,xp,yp,zp,wp],[z,w,x,y,zp,wp,xp,yp]]:
						list_tot.append('%s%s_%s%s_%s%s_%s%s'%(xf,zf,yf,ypf,wf,wpf,xpf,zpf))
	list_tot = np.unique(list_tot)
	list_elements = []
	for sublist in list_tot:
		ele1, ele2, ele3, ele4 = sublist.split('_')
		if ele1 == 'ET' or ele1 == 'BE' or ele1 == 'BT':
			ele1 = ele1[1]+ele1[0]
		if ele2 == 'ET' or ele2 == 'BE' or ele2 == 'BT':
			ele2 = ele2[1]+ele2[0]
		if ele3 == 'ET' or ele3 == 'BE' or ele3 == 'BT':
			ele3 = ele3[1]+ele3[0]
		if ele4 == 'ET' or ele4 == 'BE' or ele4 == 'BT':
			ele4 = ele4[1]+ele4[0]

		## Exception with EBEB derivatives
		if ele1 == 'BB' and block.split('_')[0] == 'EBEB': ele1 = 'EE'
		elif ele1 == 'EE' and block.split('_')[0] == 'EBEB': ele1 = 'BB'
		if ele4 == 'BB'  and block.split('_')[1] == 'EBEB': ele4 = 'EE'
		elif ele4 == 'EE'  and block.split('_')[1] == 'EBEB': ele4 = 'BB'

		## Exception with EETE derivatives
		if ele1 == 'TE' and block.split('_')[0] == 'EETE': ele1 = 'EE'
		elif ele1 == 'EE' and block.split('_')[0] == 'EETE': ele1 = 'TE'
		if ele4 == 'TE' and block.split('_')[1] == 'EETE': ele4 = 'EE'
		elif ele4 == 'EE' and block.split('_')[1] == 'EETE': ele4 = 'TE'

		## Exception with TTTE derivatives
		if ele1 == 'TE'  and block.split('_')[0] == 'TTTE': ele1 = 'TT'
		elif ele1 == 'TT'  and block.split('_')[0] == 'TTTE': ele1 = 'TE'
		if ele4 == 'TE'  and block.split('_')[1] == 'TTTE': ele4 = 'TT'
		elif ele4 == 'TT'  and block.split('_')[1] == 'TTTE': ele4 = 'TE'

		## Exception with TTTE derivatives
		if ele1 == 'EE'  and block.split('_')[0] == 'TETE': ele1 = 'TT'
		elif ele1 == 'TT'  and block.split('_')[0] == 'TETE': ele1 = 'EE'
		if ele4 == 'EE'  and block.split('_')[1] == 'TETE': ele4 = 'TT'
		elif ele4 == 'TT'  and block.split('_')[1] == 'TETE': ele4 = 'EE'

		list_elements.append(['%s_%s'%(derivativeXYZW,ele1), ele2, ele3, '%s_%s'%(derivativeXpYpZpWp,ele4)])
	return list_elements

########################################################################
# Generate simple Gaussian variance

def gaussian_variance(cl11,cl22=None,cl12=None,ells=[],nzcount=None):
	'''
	Compute a simple gaussian variance from a given C_ell
	Input
		* cl11: 1D array, auto-spectrum C_ell for the flavor 11.
		* cl11: 1D array, auto-spectrum C_ell for the flavor 22.
		* cl11: 1D array, cross-spectrum C_ell for the flavor 12.
		* ells: 1D array, the range of ells.
		* nzcount: 1D array, if specified it replaces the number of modes used to normalized.
	Output
		* diag(gauss_var): 2D array, variance of the C_ell (diagonal matrix)
	'''
	if nzcount is not None:
		gauss_var = np.array( [(cl11[i]*cl22[i] + cl12[i]**2) / nzcount[i] for i in ells] ,dtype=float)
		return np.diag(gauss_var)

	if cl22 is None and cl12 is None:
		# print 'Gaussian variance for auto-spectrum'
		gauss_var = np.array( [2.0*cl11[i]**2 / (2.0*i + 1.0) for i in ells] ,dtype=float)
	else:
		# print 'Gaussian variance for cross-spectrum'
		gauss_var = np.array( [(cl11[i]*cl22[i] + cl12[i]**2) / (2.0*i + 1.0) for i in ells] ,dtype=float)
	return np.diag(gauss_var)

def cross_gaussian_variance(cl1,cl2,ells):
	'''
	Compute a simple gaussian variance from two C_ells
	Input
		* cl1: 1D array, the first C_ell
		* cl2: 1D array, the second C_ell
		* ells: 1D array, the range of ells
	Output
		* diag(gauss_var): 2D array, variance of the C_ells (diagonal matrix)
	'''
	gauss_var = np.array( [2.0*cl1[i]*cl2[i] / (2.0*i + 1.0) for i in ells] ,dtype=float)
	return np.diag(gauss_var)

########################################################################
# Main functions to compute the covariances

def precompute_trispB(cls_unlensed,cls_lensed,combination='TTTT',lmin=2,noise_uK_arcmin=0.0,
				fwhm_arcmin=0.0,MPI=None,exp='CMB-S4',folder_cache='cache'):
	'''
	This routine returns the terms sum_l g_XY f_ZW necessary for the Type B trispectrum contribution for the cross-covariance.
	See Eq. 37 in 1611.01446
	Those terms are then combined later in the phixCMB routine.
	Input:
		* cls_unlensed: camb_clfile object, containing unlensed power spectra.
		* cls_lensed: camb_clfile object, containing lensed power spectra.
		* combination: string, names XYZW for the weight functions (e.g. TTTT, EEEE, BBBB, EEBB, etc).
		* lmin: int, minimum multipole.
		* noise_uK_arcmin: float, white noise level to be used (in uK.arcmin).
		* fwhm_arcmin: float, beam FWHM to be used.
		* MPI: module, module to pass if you want to parallelize the computation. If None, the computation is done serially.
		* exp: string, name of the experiment you chose according to misc.py
		* folder_cache: string, where to save the data.
	Output
		* trispB_tot: 2D array, sum_l g_XY f_ZW
	'''
	try:
		comm = MPI.COMM_WORLD
		rank = comm.rank
		n_tot = comm.size
		Barrier = comm.Barrier
	except:
		## No MPI
		comm = None
		rank = 0
		n_tot = 1
		Barrier = lambda: -1

	lmax = int(cls_unlensed.lmax)

	trispB_tot = np.zeros((lmax+1,lmax+1))
	l2range = range(lmin+rank,lmax+1,n_tot)
	l2dim = len(l2range)-1

	flavor1 = combination[0:2].lower()
	flavor2 = combination[2:4].lower()

	## If ET or BE, load as if it was TE and EB, and inside the code, everything is done correctly.
	if flavor1 in ['et', 'be']:
		cl_len_XX1, cl_len_YY1, cl_len_XY1 = lib_spectra.load_weights(cls_lensed, 'cl'+flavor1[1]+flavor1[0], noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		spinl2_x1, spinl3_x1, spinl2_y1, spinl3_y1 = lib_spectra.load_spin_values_wigner('cl'+flavor1[1]+flavor1[0])
	elif flavor1 == 'tb':
		cl_len_XX1, cl_len_YY1, junk = lib_spectra.load_weights(cls_lensed, 'cl'+flavor1, noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		junk1, junk2, cl_len_XY1 = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		spinl2_x1, spinl3_x1, spinl2_y1, spinl3_y1 = lib_spectra.load_spin_values_wigner('cl'+flavor1)
	elif flavor1 == 'bt':
		cl_len_XX1, cl_len_YY1, junk = lib_spectra.load_weights(cls_lensed, 'cl'+flavor1[1]+flavor1[0], noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		junk1, junk2, cl_len_XY1 = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		spinl2_x1, spinl3_x1, spinl2_y1, spinl3_y1 = lib_spectra.load_spin_values_wigner('cl'+flavor1[1]+flavor1[0])
	else:
		cl_len_XX1, cl_len_YY1, cl_len_XY1 = lib_spectra.load_weights(cls_lensed, 'cl'+flavor1, noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		spinl2_x1, spinl3_x1, spinl2_y1, spinl3_y1 = lib_spectra.load_spin_values_wigner('cl'+flavor1)

	if flavor2 in ['et', 'be']:
		cl_len_XX2, cl_len_YY2, cl_len_XY2 = lib_spectra.load_weights(cls_lensed, 'cl'+flavor2[1]+flavor2[0], noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		spinl2_x2, spinl3_x2, spinl2_y2, spinl3_y2 = lib_spectra.load_spin_values_wigner('cl'+flavor2[1]+flavor2[0])
	elif flavor2 == 'tb':
		cl_len_XX2, cl_len_YY2, junk = lib_spectra.load_weights(cls_lensed, 'cl'+flavor2, noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		junk1, junk2, cl_len_XY2 = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		spinl2_x2, spinl3_x2, spinl2_y2, spinl3_y2 = lib_spectra.load_spin_values_wigner('cl'+flavor2)
	elif flavor2 == 'bt':
		cl_len_XX2, cl_len_YY2, junk = lib_spectra.load_weights(cls_lensed, 'cl'+flavor2[1]+flavor2[0], noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		junk1, junk2, cl_len_XY2 = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		spinl2_x2, spinl3_x2, spinl2_y2, spinl3_y2 = lib_spectra.load_spin_values_wigner('cl'+flavor2[1]+flavor2[0])
	else:
		cl_len_XX2, cl_len_YY2, cl_len_XY2 = lib_spectra.load_weights(cls_lensed, 'cl'+flavor2, noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		spinl2_x2, spinl3_x2, spinl2_y2, spinl3_y2 = lib_spectra.load_spin_values_wigner('cl'+flavor2)

	## Computation of the derivatives
	Bl2l3 = loop_lensing.precompute_trispb_fortran(cl_len_XX1, cl_len_YY1, cl_len_XY1,
													cl_len_XX2, cl_len_YY2, cl_len_XY2,
													l2range,flavor1,flavor2,lmin,
													spinl2_x1,spinl3_x1,spinl2_y1,spinl3_y1,
													spinl2_x2, spinl3_x2, spinl2_y2, spinl3_y2,
													l2dim,lmax)
	comm.Reduce([Bl2l3, MPI.DOUBLE],[trispB_tot, MPI.DOUBLE],op = MPI.SUM,root = 0)
	if rank == 0:
		## Because fortran routine returns the transposed
		trispB_tot = trispB_tot.T

	Barrier()

	return trispB_tot

def analytic_covariances_CMBxCMB(cls_unlensed,cls_lensed,lmin=2,blocks=['TTTT'],noise_uK_arcmin=0.0,
		fwhm_arcmin=0.0,MPI=None,use_corrfunc=False,exp='CMB-S4',folder_cache='cache'):
	'''
	This routine computes the autocovariance of lensed CMB spectra at second order in C^{phi phi}.
	This corresponds to equations 27, 29, 30 of 1611.01446
	We pre-compute first the derivatives, and then update the covariance.
	The routine is parallelized in a dump way (and within python), by splitting evenly the ells of the first sum
	within processors. For ell > 3000, better to use many procs (few minutes on 1 core).
	/!\ spin conventions are hardcoded in order to fit with fortran routines (see loop_lensing.f90).
	If you modify the fortran routines, the values of spin can change. /!\
	Input:
		* cls_unlensed: camb_clfile object, containing unlensed power spectra.
		* cls_lensed: camb_clfile object, containing lensed power spectra.
		* lmin: int, minimum multipole.
		* blocks: list of strings, names for the two 2-pt functions (e.g. [TTTT, EEEE, BBBB, EEBB]).
		* noise_uK_arcmin: float, white noise level to be used (in uK.arcmin).
		* fwhm_arcmin: float, beam FWHM to be used.
		* MPI: module, module to pass if you want to parallelize the computation. If None,
						the computation is done serially.
		* use_corrfunc: boolean, if True use correlation function method to compute derivatives.
									Use series-expansion otherwise.
		* exp: string, name of the experiment you chose according to misc.py
		* folder_cache: string, where to save the data.
	Output
		* cov_order0_tot, cov_order1_tot, cov_order2_tot: 2D arrays, covariance matrices at zeroth,
														first and second order in C^{phi phi}.
			Note that for covariance involving only B-modes,
						it contains only terms at second order (written in cov1 and cov2).
	'''
	try:
		comm = MPI.COMM_WORLD
		rank = comm.rank
	except:
		## No MPI
		comm = None
		rank = 0

	## Maximum multipole
	lmax = int(cls_unlensed.lmax)

	## Input lensing potential power spectrum
	clpp = cls_unlensed.clpp_long

	## Choice of method for the computation of derivatives
	## By default, it uses the correlation function method (see App. B in the paper)
	## Note that if the derivatives are already computed, the code won't recompute them.
	if use_corrfunc:
		ACCURACY_BOOST 		 = 4 # 4 is conservative. oversampling of the theta grid
		accurate_lensing 	 = True # Gauss-Legendre (true) vs simple integration (False)
		EXTRA_MULTIPOLES = 0 # If you want to add more multipole to compute dC^{CMB}/dC^{\phi \phi} (not necessary)
	else:
		EXTRA_MULTIPOLES = 0 # If you want to add more multipole to compute dC^{CMB}/dC^{\phi \phi} (not necessary)

	#############################################################
	# Initialization of containers
	#############################################################
	cov_order0_tot = np.array([covmat(0,lmax) for i in range(len(blocks))]) ## Gaussian variance
	cov_order1_tot = np.array([covmat(0,lmax) for i in range(len(blocks))]) ## O(clpp)
	cov_order2_tot = np.array([covmat(0,lmax) for i in range(len(blocks))]) ## O(clpp^2)

	## We have 4 derivatives dC^{CMB}/dC^{phiphi} to compute: TTTT, EEEE, BBBB, TETE
	names_dCMB_over_dcpp_tot = ['TTTT','EEEE','BBBB','TETE']
	dCMB_over_dcpp_tot = np.array([np.zeros((lmax+1,lmax+1)) for i in range(len(names_dCMB_over_dcpp_tot))])


	## We approximate other derivatives (e.g. dlensedEE/dEE) by ones.
	dBB_over_dEE_tot = np.zeros((lmax+1,lmax+1))

	#############################################################
	# Gaussian variance
	#############################################################
	for position_block,block in enumerate(blocks):
		index_block = blocks.index(block)
		flavor = 'cl%s%s'%(block[0].lower(),block[2].lower())

		if block in ['EEBB','TTBB','TEBB']:
			## There is no Gaussian variance contribution for those terms
			continue

		if rank==0: print 'Gaussian Variance: doing block %s (%s)\n'%(block,flavor)

		## Load weights (lensed spectra, and their noisy version)
		if block in ['TETE', 'TTTE', 'EETE']:
			cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, 'clte',
										noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		else:
			cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, flavor,
										noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')

		## Gaussian variance
		if block == 'TTEE':
			cov_order0_tot[index_block].data = cross_gaussian_variance(cl1=cl_len_XY[1],cl2=cl_len_XY[1],ells=cls_lensed.ls)
		elif block == 'TTTE':
			cov_order0_tot[index_block].data = cross_gaussian_variance(cl1=cl_len_XX[1],cl2=cl_len_XY[1],ells=cls_lensed.ls)
		elif block == 'EETE':
			cov_order0_tot[index_block].data = cross_gaussian_variance(cl1=cl_len_YY[1],cl2=cl_len_XY[1],ells=cls_lensed.ls)
		else:
			cov_order0_tot[index_block].data = gaussian_variance(cl11=cl_len_XX[1],cl22=cl_len_YY[1],
												cl12=cl_len_XY[1],ells=cls_lensed.ls)


	#############################################################
	## Contribution of the trispectrum to the covariance: O(clpp)
	#############################################################
	for block in blocks:
		index_block = blocks.index(block)
		flavor = 'cl%s%s'%(block[0].lower(),block[2].lower())

		if block in ['BBBB','TTBB','EEBB','TETE','TTEE','TTTE','EETE','TEBB']:
			## We do not consider the contribution of those terms (although you can)
			continue
		cov_order1 = covmat(0,lmax)

		if rank==0: print 'Order O(clpp): doing block %s (%s)\n'%(block,flavor)

		## Load spins
		spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner('clte')

		## Load weights (lensed spectra, and their noisy version)
		cl_unlen_TT, cl_unlen_EE, cl_unlen_TE = lib_spectra.load_weights(cls_unlensed, 'clte',
										noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		cl_unlen_vec = np.array([cl_unlen_TT[0], cl_unlen_EE[0], np.zeros_like(cl_unlen_TT[0]), cl_unlen_TE[0]])

		## Load weights (unlensed spectra)
		uup = block[0] + block[2]
		vvp = block[1] + block[3]
		uvp = block[0] + block[3]
		vup = block[1] + block[2]

		## Define range of ells, and distribute over procs.
		n_tot = comm.size
		l2range = range(lmin+rank,lmax+1,n_tot)
		l2dim = len(l2range)-1

		## Compute this term
		cov_order1.data = loop_lensing.covariance_cmbxcmb_order1_uvupvp(cl_unlen_vec,clpp,l2range,
									uup,vvp,uvp,vup,lmin,spinl2_x,spinl3_x,spinl2_y,spinl3_y,l2dim,lmax)
		comm.Barrier()

		## Reduce the results on the root
		comm.Reduce([cov_order1.data, MPI.DOUBLE],[cov_order1_tot[index_block].data, MPI.DOUBLE],op = MPI.SUM,root = 0)

		## Done for this block
		comm.Barrier()

	#############################################################
	## Compute dC^{CMB}/dC^{phiphi}
	## Here, you have two ways of computing the derivatives:
	## 		* Using series-expansion. Quick but less accurate.
	## 		* Using correlation functions. Less quick, but extra accurate.
	#############################################################
	file_manager_derivatives_CMB = util.file_manager('dCMB_over_dcpp_tot', exp, spec='v1', lmax=lmax,
													force_recomputation=False, folder=folder_cache,rank=rank)
	if file_manager_derivatives_CMB.FileExist is True:
		if rank==0:
			dCMB_over_dcpp_tot, names_dCMB_over_dcpp_tot = file_manager_derivatives_CMB.data
	else:
		for position_block,block in enumerate(names_dCMB_over_dcpp_tot):
			flavor = 'cl%s%s'%(block[0].lower(),block[1].lower())
			if rank==0: print 'Pre-compute derivatives for block %s (%s)\n'%(block,flavor)

			if not use_corrfunc:
				if rank == 0: print 'Use series-expansion to compute derivative (may not be exact)'
				if block == 'BBBB':
					## BB takes clee as unlensed weights (noiseless!)
					cl_unlen_XX, cl_unlen_YY, cl_unlen_XY = lib_spectra.load_weights(cls_unlensed, 'clee', 0.0,
																0.0, 2*lmax, extra='_long')
				else:
					## noiseless!
					cl_unlen_XX, cl_unlen_YY, cl_unlen_XY = lib_spectra.load_weights(cls_unlensed, flavor,
																0.0, 0.0, 2*lmax, extra='_long')

				## Define range of ells, and distribute over procs.
				n_tot = comm.size
				l2range = range(lmin+rank,lmax+1,n_tot)
				l2dim = len(l2range)-1

				## Load spins
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner(flavor)

				## Change order of spins. Why? do not know... but it works :D
				derivatives = loop_lensing.compute_derivatives_dcttdcpp_mpi(cl_unlen_XY[0],l2range,flavor,lmin,spinl3_x,
												spinl2_x,spinl3_y,spinl2_y,l2dim,lmax)

				## Reduce on the root
				comm.Reduce([derivatives, MPI.DOUBLE],[dCMB_over_dcpp_tot[position_block], MPI.DOUBLE],op = MPI.SUM,root = 0)
			else:
				if rank == 0: print 'Use correlation functions to compute derivative'
				if block == 'BBBB':
					## BB takes clee as unlensed weights
					## noiseless!
					cl_unlen_XX, cl_unlen_YY, cl_unlen_XY = lib_spectra.load_weights(cls_unlensed, 'clee',
																0.0, 0.0, 2*lmax, extra='_long')
				else:
					## noiseless!
					cl_unlen_XX, cl_unlen_YY, cl_unlen_XY = lib_spectra.load_weights(cls_unlensed, flavor,
																0.0, 0.0, 2*lmax, extra='_long')

				clpp_long = cls_unlensed.clpp_long
				## Define container used to reduce results on root
				derivatives_tot_tmp = np.zeros((lmax+1,lmax+1))

				## Define range of ells, and distribute over procs.
				n_tot = comm.size
				l2range = range(lmin+rank,lmax+1,n_tot)
				l2dim = len(l2range)-1

				## Compute.
				## Use cl_unlen_XY which is TT for TT, EE for EE, EE for BB, and TE for TE.
				dxim,dxip,dm,dp = correlation_functions.derivative_dclcmbdclpp_corr_func(ACCURACY_BOOST*lmax,flavor,
											accurate_lensing,cl_unlen_XY[0][0:lmax+EXTRA_MULTIPOLES+1],
											clpp_long[0:lmax+EXTRA_MULTIPOLES+1],l2range,lmin,lmax,l2dim,lmax+EXTRA_MULTIPOLES)

				if block in ['TTTT','TETE']:
					derivatives = 2*np.pi*np.dot(dxim.T,dm)
				if block=='EEEE':
					derivatives = 2*np.pi*(np.dot(dxim.T,dm) + np.dot(dxip.T,dp))/2.0
				if block=='BBBB':
					derivatives = 2*np.pi*(np.dot(dxip.T,dp) - np.dot(dxim.T,dm))/2.0

				## Reduce on root
				comm.Reduce([derivatives, MPI.DOUBLE],[dCMB_over_dcpp_tot[position_block], MPI.DOUBLE],op = MPI.SUM, root = 0)

				## Free a bit to save memory
				del dxip,dxim,dm,dp

			## Done for this block
			comm.Barrier()

		## /!\ /!\ /!\ The fortran routine returns the transposed! /!\ /!\ /!\
		array_to_save = [dCMB_over_dcpp_tot, names_dCMB_over_dcpp_tot]
		## /!\ /!\ /!\ The fortran routine returns the transposed! /!\ /!\ /!\

		if rank==0: file_manager_derivatives_CMB.save_data_on_disk(array_to_save)
	comm.Barrier()

	#############################################################
	## Compute dlensedCMB/dunlensedCMB
	## For the moment, we only consider the case for BB
	## All the others are set to 1
	#############################################################
	for block in ['BBBB']:
		flavor = 'cl%s%s'%(block[0].lower(),block[2].lower())
		## Define range of ells, and distribute over procs.
		n_tot = comm.size
		l2range = range(lmin+rank,lmax+1,n_tot)
		l2dim = len(l2range)-1

		## Compute
		derivatives = loop_lensing.compute_derivatives_dcttdcpp_mpi(clpp,l2range,flavor,lmin,2,-2,l2dim,lmax)

		## Reduce
		comm.Reduce([derivatives, MPI.DOUBLE],[dBB_over_dEE_tot, MPI.DOUBLE],op = MPI.SUM,root = 0)

		comm.Barrier()

	#############################################################
	## At this stage, we computed all necessary objects and the
	## Gaussian and O(clpp) terms have already been updated.
	## Loop over blocks to update terms with O(clpp^2) contributions.
	#############################################################
	if rank==0:
		gaussvar_phiphi_nonoise = np.diag(gaussian_variance(cl11=clpp,ells=cls_lensed.ls))
		for block in blocks:
			index_block = blocks.index(block)
			if rank == 0: print 'Updating block %s'%block

			if block == 'TTTT':
				index_derivative = names_dCMB_over_dcpp_tot.index('TTTT')
				cov_order2_tot[index_block].data = matvecmat(dCMB_over_dcpp_tot[index_derivative].T,gaussvar_phiphi_nonoise,
												dCMB_over_dcpp_tot[index_derivative],lstart=lmin)
			elif block == 'EEEE':
				index_derivative = names_dCMB_over_dcpp_tot.index('EEEE')
				cov_order2_tot[index_block].data = matvecmat(dCMB_over_dcpp_tot[index_derivative].T,gaussvar_phiphi_nonoise,
												dCMB_over_dcpp_tot[index_derivative],lstart=lmin)
			elif block == 'TTEE':
				index_derivative_TT = names_dCMB_over_dcpp_tot.index('TTTT')
				index_derivative_EE = names_dCMB_over_dcpp_tot.index('EEEE')
				cov_order2_tot[index_block].data = matvecmat(dCMB_over_dcpp_tot[index_derivative_TT].T,
												gaussvar_phiphi_nonoise,dCMB_over_dcpp_tot[index_derivative_EE],lstart=lmin)
			elif block == 'BBBB':
				index_derivative = names_dCMB_over_dcpp_tot.index('BBBB')
				cov_order2_tot[index_block].data += matvecmat(dCMB_over_dcpp_tot[index_derivative].T,
												gaussvar_phiphi_nonoise,dCMB_over_dcpp_tot[index_derivative],lstart=lmin)
				## Second term in O(clpp^2) uses unlensed EE
				cl_unlen_XX, cl_unlen_YY, cl_unlen_XY = lib_spectra.load_weights(cls_unlensed, 'clee',
														noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
				gaussvar_EE_unlensed = np.diag(gaussian_variance(cl11=cl_unlen_XX[0],ells=cls_lensed.ls))
				cov_order2_tot[index_block].data += matvecmat(dBB_over_dEE_tot.T,gaussvar_EE_unlensed,dBB_over_dEE_tot,lstart=lmin)
			elif block == 'EEBB':
				index_derivative_EE = names_dCMB_over_dcpp_tot.index('EEEE')
				index_derivative_BB = names_dCMB_over_dcpp_tot.index('BBBB')
				cov_order2_tot[index_block].data += matvecmat(dCMB_over_dcpp_tot[index_derivative_EE].T,
												gaussvar_phiphi_nonoise,dCMB_over_dcpp_tot[index_derivative_BB],lstart=lmin)
				## Second term in O(clpp^2) uses unlensed EE
				cl_unlen_XX, cl_unlen_YY, cl_unlen_XY = lib_spectra.load_weights(cls_unlensed, 'clee',
															noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
				gaussvar_EE_unlensed = np.diag(gaussian_variance(cl11=cl_unlen_XX[0],ells=cls_lensed.ls))
				cov_order2_tot[index_block].data += vecmat(gaussvar_EE_unlensed,dBB_over_dEE_tot,lstart=lmin)
			elif block == 'TTBB':
				index_derivative_TT = names_dCMB_over_dcpp_tot.index('TTTT')
				index_derivative_BB = names_dCMB_over_dcpp_tot.index('BBBB')
				cov_order2_tot[index_block].data += matvecmat(dCMB_over_dcpp_tot[index_derivative_TT].T,
											gaussvar_phiphi_nonoise,dCMB_over_dcpp_tot[index_derivative_BB],lstart=lmin)
				## Second term in O(clpp^2) uses unlensed TE
				cl_unlen_XX, cl_unlen_YY, cl_unlen_XY = lib_spectra.load_weights(cls_unlensed, 'clte',
														noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
				gaussvar_TE_unlensed = np.diag(cross_gaussian_variance(cl1=cl_unlen_XY[0],cl2=cl_unlen_XY[0],ells=cls_lensed.ls))
				cov_order2_tot[index_block].data += vecmat(gaussvar_TE_unlensed,dBB_over_dEE_tot,lstart=lmin)
			elif block == 'TETE':
				index_derivative = names_dCMB_over_dcpp_tot.index('TETE')
				cov_order2_tot[index_block].data = matvecmat(dCMB_over_dcpp_tot[index_derivative].T,
												gaussvar_phiphi_nonoise,dCMB_over_dcpp_tot[index_derivative],lstart=lmin)
			elif block == 'TTTE':
				index_derivative_TT = names_dCMB_over_dcpp_tot.index('TTTT')
				index_derivative_TE = names_dCMB_over_dcpp_tot.index('TETE')
				cov_order2_tot[index_block].data = matvecmat(dCMB_over_dcpp_tot[index_derivative_TT].T,
												gaussvar_phiphi_nonoise,dCMB_over_dcpp_tot[index_derivative_TE],lstart=lmin)
			elif block == 'EETE':
				index_derivative_EE = names_dCMB_over_dcpp_tot.index('EEEE')
				index_derivative_TE = names_dCMB_over_dcpp_tot.index('TETE')
				cov_order2_tot[index_block].data = matvecmat(dCMB_over_dcpp_tot[index_derivative_EE].T,
												gaussvar_phiphi_nonoise,dCMB_over_dcpp_tot[index_derivative_TE],lstart=lmin)
			elif block == 'TEBB':
				index_derivative_TE = names_dCMB_over_dcpp_tot.index('TETE')
				index_derivative_BB = names_dCMB_over_dcpp_tot.index('BBBB')
				cov_order2_tot[index_block].data += matvecmat(dCMB_over_dcpp_tot[index_derivative_TE].T,
												gaussvar_phiphi_nonoise,dCMB_over_dcpp_tot[index_derivative_BB],lstart=lmin)
				## Second term in O(clpp^2) uses unlensed EE and TE
				cl_unlen_XX, cl_unlen_YY, cl_unlen_XY = lib_spectra.load_weights(cls_unlensed, 'clte',
															noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
				gaussvar_EETE_unlensed = np.diag(cross_gaussian_variance(cl1=cl_unlen_YY[1],cl2=cl_unlen_XY[1],ells=cls_lensed.ls))
				cov_order2_tot[index_block].data += vecmat(gaussvar_EETE_unlensed,dBB_over_dEE_tot,lstart=lmin)

	return cov_order0_tot, cov_order1_tot, cov_order2_tot, blocks

def analytic_covariances_phixCMB(cls_unlensed,cls_lensed,lmin=2,noise_uK_arcmin=0.0,fwhm_arcmin=0.0,
				MPI=None,use_corrfunc=False,exp='CMB-S4',folder_cache='cache',path_to_N1matrix=''):
	'''
	This routine computes the covariance between lensed CMB spectrum and the lensing potential power spectrum.
	This is a perturbative derivation which focuses the 4 main parts of the covariance:
		- noise
		- signal
		- trispectrum contributions (Type A and B).
	This is an implementation of the Eqs. 35, 36, 39, and 40 in the paper respectively.
	Input
		* cls_unlensed: camb_clfile object, containing unlensed power spectra.
		* cls_lensed: camb_clfile object, containing lensed power spectra.
		* lmin: int, minimum multipole.
		* blocks: list of strings, corresponds to phi_XY phi_WZ x UV (e.g. [TTTT_TT, EEEE_EE, BBBB_BB, EBEB_EE, EBEB_BB])
		* noise_uK_arcmin: float, white noise level to be used (in uK.arcmin).
		* fwhm_arcmin: float, beam FWHM to be used.
		* MPI: module, module to pass if you want to parallelize the computation. If None, the computation is done serially.
		* use_corrfunc: boolean, if True use correlation function method to compute derivatives.
								Use series-expansion otherwise.
		* exp: string, name of the experiment you chose according to misc.py
		* folder_cache: string, where to save the data.
		* path_to_N1matrix: string, where to find the derivative of N1 wrt Cphiphi (external, from Biases_n1mat.f90)
	Output
		* cov_MV, cov_MV_signal, cov_MV_noise, cov_MV_trispA, cov_MV_trispB, combinations_CMB: 2D array,
			covariance matrices between lensed CMB spectrum and the lensing potential power spectrum:
				- total contribution
				- signal contribution
				- noise contribution
				- Type A trispectrum contribution
				- Type B trispectrum contribution
	'''
	try:
		comm = MPI.COMM_WORLD
		rank = comm.rank
		n_tot = comm.size
		Barrier = comm.Barrier
	except:
		## No MPI
		print 'Serial mode for the cross-covariances \n'
		print 'I hope you are not computing too big objects :D Otherwise, use the parallelized version.'
		comm = None
		rank = 0
		n_tot = 1
		Barrier = lambda: -1

	## Maximum multipole
	lmax = int(cls_unlensed.lmax)

	## Input lensing potential power spectrum
	clpp = cls_unlensed.clpp
	clpp_long = cls_unlensed.clpp_long

	## If you want to add more multipole for the
	## computation of dN^{(0)} over dC^CMB (not necessary)
	EXTRA_MULTIPOLES_lensing = 0

	## Choice of method for the computation of derivatives
	## By default, it uses the correlation function method (see App. B in the paper)
	## Note that if the derivatives are already computed, the code won't recompute them.
	if use_corrfunc:
		ACCURACY_BOOST 		 = 4 # 4 is conservative. oversampling of the theta grid
		accurate_lensing 	 = True # Gauss-Legendre (true) vs simple integration (False)
		EXTRA_MULTIPOLES_CMB = 0 # Add more multipole for the computation of dC^{CMB} over dC^{\phi \phi} (not necessary)
	else:
		EXTRA_MULTIPOLES = 0 # Add more multipole for the computation of dC^{CMB} over dC^{\phi \phi} (not necessary)

	## We have 4 derivatives to compute: TTTT, EEEE, BBBB, TETE
	blocks_for_CMB_derivatives = ['TTTT','EEEE','BBBB','TETE']
	dCMB_over_dcpp_tot = np.array([np.zeros((lmax+1,lmax+1)) for i in range(len(blocks_for_CMB_derivatives))])

	## Initialization of the N0 vectors
	blocks_for_N0 = ['TTTT_TT','EEEE_EE','BBBB_BB',
					 'EBEB_EE','EBEB_BB',
					 'TETE_TE','TETE_TT','TETE_EE',
					 'TTEE_TE',
					 'TTTE_TT','TTTE_TE',
					 'EETE_EE','EETE_TE']
	N0_tot = np.array([np.zeros(lmax+1) for i in range(len(blocks_for_N0))])

	## Load N0
	file_manager_N0 = util.file_manager('N0', exp, spec='v1', lmax=lmax, force_recomputation=False, folder=folder_cache,rank=rank)
	if file_manager_N0.FileExist is True:
		if rank==0: print 'N0 part loaded from %s.pkl'%file_manager_N0.fn
		N0, subblocks = file_manager_N0.data
		for position_block,block in enumerate(blocks_for_N0):
			index_N0 = subblocks.index(block[0:4])
			N0_tot[position_block] = N0[index_N0]
	else:
		if rank ==0: print 'You need to compute N0 first!'
		sys.exit()


	#############################################################
	# Computation of dN^{(0)}/dC^{CMB}
	#############################################################
	## Initialization of derivatives of N0 wrt lensed CMB spectra.
	blocks_for_phi_derivatives = ['TTTT_TT','EEEE_EE','BBBB_BB',
								  'EBEB_EE','EBEB_BB',
								  'TETE_TE','TETE_TT','TETE_EE',
								  'TTEE_TE',
								  'TTTE_TT','TTTE_TE',
								  'EETE_EE','EETE_TE']
	dN_over_dccmb_tot = np.array([np.zeros((lmax+1,lmax+1)) for i in range(len(blocks_for_phi_derivatives))])

	file_manager_derivatives = util.file_manager('dN_over_dccmb_tot', exp, spec='v1', lmax=lmax,
										force_recomputation=False, folder=folder_cache,rank=rank)
	if file_manager_derivatives.FileExist is True:
		if rank==0:
			dN_over_dccmb_tot, blocks_for_phi_derivatives = file_manager_derivatives.data
	else:
		### Precompute objects: derivatives of RDN0 wrt lensed CMB spectra and N0
		for block in blocks_for_phi_derivatives:
			flavor = 'cl%s%s'%(block[0].lower(),block[1].lower())
			index_derivative = blocks_for_phi_derivatives.index(block)

			if rank==0: print 'Pre-compute derivatives: doing block (%s)\n'%(block)

			## Lensed spectra used for the weights
			## Load spins
			if block == 'TETE_TE':
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner('clte')
				flavor = 'cltete_te' ## used for the derivative
			elif block == 'TETE_TT':
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner('clte')
				flavor = 'cltete_tt' ## used for the derivative
			elif block == 'TTTE_TT':
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner('clte')
				flavor = 'clttte_tt' ## used for the derivative
			elif block == 'TTTE_TE':
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner('clte')
				flavor = 'clttte_te' ## used for the derivative
			elif block == 'EETE_EE':
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner('clte')
				flavor = 'cleete_ee' ## used for the derivative
			elif block == 'EETE_TE':
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner('clte')
				flavor = 'cleete_te' ## used for the derivative
			elif block == 'TETE_EE':
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner('clte')
				flavor = 'cltete_ee' ## used for the derivative
			elif block == 'TTEE_TE':
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner('clte')
				flavor = 'clttee_te' ## used for the derivative
			elif block == 'EBEB_EE':
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, flavor, noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner(flavor)
				flavor = 'clebeb_ee'
			elif block == 'EBEB_BB':
				## Swap  the values
				cl_len_YY, cl_len_XX, cl_len_XY = lib_spectra.load_weights(cls_lensed, flavor, noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner(flavor)
				flavor = 'clebeb_ee'
			else:
				cl_len_XX, cl_len_YY, cl_len_XY = lib_spectra.load_weights(cls_lensed, flavor, noise_uK_arcmin,
														fwhm_arcmin, 2*lmax, extra='_long')
				spinl2_x, spinl3_x, spinl2_y, spinl3_y = lib_spectra.load_spin_values_wigner(flavor)

			weight_len_XX = np.zeros((2,2*(lmax+EXTRA_MULTIPOLES_lensing)+1),dtype=float)
			weight_len_XX[0] = cl_len_XX[0][0:2*(lmax+EXTRA_MULTIPOLES_lensing)+1]
			weight_len_XX[1] = cl_len_XX[1][0:2*(lmax+EXTRA_MULTIPOLES_lensing)+1]
			weight_len_YY = np.zeros((2,2*(lmax+EXTRA_MULTIPOLES_lensing)+1),dtype=float)
			weight_len_YY[0] = cl_len_YY[0][0:2*(lmax+EXTRA_MULTIPOLES_lensing)+1]
			weight_len_YY[1] = cl_len_YY[1][0:2*(lmax+EXTRA_MULTIPOLES_lensing)+1]
			weight_len_XY = np.zeros((2,2*(lmax+EXTRA_MULTIPOLES_lensing)+1),dtype=float)
			weight_len_XY[0] = cl_len_XY[0][0:2*(lmax+EXTRA_MULTIPOLES_lensing)+1]
			weight_len_XY[1] = cl_len_XY[1][0:2*(lmax+EXTRA_MULTIPOLES_lensing)+1]

			derivatives_tot_tmp = np.zeros((lmax+1+EXTRA_MULTIPOLES_lensing,lmax+1+EXTRA_MULTIPOLES_lensing))
			l2range = range(lmin+rank,lmax+1+EXTRA_MULTIPOLES_lensing,n_tot)
			l2dim = len(l2range)-1

			## Computation of the derivatives
			Bl2l3 = loop_lensing.xx_phiphi_disconnected_mpi(weight_len_XX,weight_len_YY,weight_len_XY,l2range,flavor,
										lmin,spinl2_x,spinl3_x,spinl2_y,spinl3_y,l2dim,lmax+EXTRA_MULTIPOLES_lensing)
			comm.Reduce([Bl2l3, MPI.DOUBLE],[derivatives_tot_tmp, MPI.DOUBLE],op = MPI.SUM,root = 0)
			if rank == 0:
				## Keep only elements up to lmax
				dN_over_dccmb_tot[index_derivative] = np.array([derivatives_tot_tmp[i][0:lmax+1] for i in range(0,lmax+1)])
				## Because fortran routine returns the transposed
				dN_over_dccmb_tot[index_derivative] = dN_over_dccmb_tot[index_derivative].T
			Barrier()

			if rank==0:
				if block[0:4] == 'EBEB':
					index1 = blocks_for_N0.index('EBEB_EE')
					index2 = blocks_for_N0.index('EBEB_EE')
					vecAmplitude = 1.0 * N0_tot[index1] * N0_tot[index2] / (2*cls_lensed.ls + 1.0)
				elif block in ['TETE_TE','TTTE_TT', 'TTTE_TE', 'EETE_EE', 'EETE_TE']:
					index1 = blocks_for_N0.index(block[0:2]+block[0:2]+'_'+block[0:2])
					index2 = blocks_for_N0.index(block[2:4]+block[2:4]+'_'+block[2:4])
					vecAmplitude = 2.0 * N0_tot[index1] * N0_tot[index2] / (2*cls_lensed.ls + 1.0)
				elif block in ['TETE_TT','TETE_EE']:
					index1 = blocks_for_N0.index(block[0:2]+block[0:2]+'_'+block[0:2])
					index2 = blocks_for_N0.index(block[2:4]+block[2:4]+'_'+block[2:4])
					vecAmplitude = 1.0 * N0_tot[index1] * N0_tot[index2] / (2*cls_lensed.ls + 1.0)
				else:
					index1 = blocks_for_N0.index(block[0:2]+block[0:2]+'_'+block[0:2])
					index2 = blocks_for_N0.index(block[2:4]+block[2:4]+'_'+block[2:4])
					vecAmplitude = 4.0 * N0_tot[index1] * N0_tot[index2] / (2*cls_lensed.ls + 1.0)

				dN_over_dccmb_tot[index_derivative] = vecmat(vecAmplitude,dN_over_dccmb_tot[index_derivative])

			Barrier()
		array_to_save = [dN_over_dccmb_tot, blocks_for_phi_derivatives]
		if rank==0: file_manager_derivatives.save_data_on_disk(array_to_save)
	## End of the loop. All products are computed.
	Barrier()


	#############################################################
	# Load cov(CMB,CMB)
	#############################################################
	file_manager_CMBxCMB = util.file_manager('covariances_CMBxCMB', exp, spec='v1', lmax=lmax,
								force_recomputation=False, folder=folder_cache,rank=rank)
	if file_manager_CMBxCMB.FileExist is True:
		if rank==0: print 'CMBxCMB part already computed. Loading from %s'%file_manager_CMBxCMB.fn
		cov_order0_CMBxCMB, cov_order1_CMBxCMB, cov_order2_CMBxCMB, blocks_CMBxCMB = file_manager_CMBxCMB.data
	else:
		if rank==0: print 'Computing CMBxCMB part'
		blocks_CMBxCMB = ['TTTT', 'EEEE', 'BBBB', 'EEBB', 'TTEE', 'TTBB', 'TETE', 'TTTE', 'EETE', 'TEBB']
		cov_order0_CMBxCMB, cov_order1_CMBxCMB, cov_order2_CMBxCMB, blocks_CMBxCMB = analytic_covariances_CMBxCMB(cls_unlensed,cls_lensed,
				lmin=lmin,blocks=blocks_CMBxCMB, noise_uK_arcmin=noise_uK_arcmin,
				fwhm_arcmin=fwhm_arcmin, MPI=MPI,exp=exp,folder_cache=folder_cache)


	#############################################################
	# Load dC^{CMB}/dC^{phiphi} (computation is done in cov(CMB,CMB))
	#############################################################
	file_manager_derivatives_CMB = util.file_manager('dCMB_over_dcpp_tot', exp, spec='v1', lmax=lmax,
									force_recomputation=False, folder=folder_cache,rank=rank)
	if file_manager_derivatives_CMB.FileExist is True:
		if rank==0:
			dCMB_over_dcpp_tot, blocks_for_CMB_derivatives = file_manager_derivatives_CMB.data
	else:
		if rank==0:
			print 'You need to have the derivatives of the lensed CMB spectra wrt Cphiphi!'
			print 'Have a look at analytic_covariances_CMBxCMB()'
			sys.exit()
	Barrier()

	#############################################################
	# Actual computation of cov(CMB,phiphi)
	#############################################################
	## Done! We have all objects!
	combinations_CMB = ['TT', 'EE', 'TE', 'BB']
	cov_MV = np.array([covmat(0,lmax) for i in range( len( combinations_CMB ) ) ] )
	cov_MV_noise = np.array([covmat(0,lmax) for i in range( len( combinations_CMB ) ) ] )
	cov_MV_trispA = np.array([covmat(0,lmax) for i in range( len( combinations_CMB ) ) ] )
	cov_MV_trispB = np.array([covmat(0,lmax) for i in range( len( combinations_CMB ) ) ] )
	cov_MV_signal = np.array([covmat(0,lmax) for i in range( len( combinations_CMB ) ) ] )
	## We need to loop over blocks and update covariance matrices.
	## Done on the root (can be optimized)
	if rank==0:
		gaussvar_phiphi_nonoise = np.diag(gaussian_variance(cl11=clpp,ells=cls_lensed.ls))

		## Names of spectra
		combinations = ['TT', 'EE', 'BB', 'TE', 'EB', 'TB']
		all_blocks = ['%s%s_%s'%(i,j,k) for i in combinations for j in combinations for k in combinations_CMB]

		## Form N0 matrix
		file_manager_N0 = util.file_manager('N0', exp, spec='v1', lmax=lmax,
							force_recomputation=False, folder=folder_cache,rank=rank)
		if file_manager_N0.FileExist is True:
			if rank==0: print 'Loading N0 from %s.pkl'%file_manager_N0.fn
			N0, names_for_N0 = file_manager_N0.data
			## Form minimum variance for N0 and compute weights
			minimum_variance_n0, weights_for_MV, names_for_weights = lib_spectra.compute_minimum_variance_weights(N0,names_for_N0)
		else:
			if rank==0: print 'You need N0!'
			sys.exit()

		## Load trispB data
		fns_trispB = glob.glob(os.path.join(folder_cache,'trispB_v1_%d_%s_*.pkl'%(lmax,exp)))
		data_trispB = []
		names_trispB = []
		for fn in fns_trispB:
			dataB, junk = np.load(fn)['data']
			data_trispB.append(dataB)
			name_trispB = fn.split('.')[0][-4:]
			names_trispB.append(name_trispB)

		for block in all_blocks:
			print 'Compute block %s (%d/%d)'%(block, all_blocks.index(block)+1,len(all_blocks))
			## Combination of terms
			names_noise, names_trispA, name_signal = combination_noise_signal_trispA(block)

			## CMB block for MV
			index_CMB = combinations_CMB.index(block[5:])

			## Weights
			wXY_index = names_for_weights.index(block[0:2])
			wZW_index = names_for_weights.index(block[2:4])
			vec1 = weights_for_MV[wXY_index] * weights_for_MV[wZW_index]

			## noise term
			mat_noise = covmat(0,lmax)
			for name_noise in names_noise:
				if name_noise[0] in blocks_for_phi_derivatives:
					index_N0_derivative = blocks_for_phi_derivatives.index(name_noise[0])
					derivative = dN_over_dccmb_tot[index_N0_derivative]
				elif name_noise[0][2:4]+name_noise[0][:2]+name_noise[0][4:] in blocks_for_phi_derivatives:
					index_N0_derivative = blocks_for_phi_derivatives.index(name_noise[0][2:4]+name_noise[0][:2]+name_noise[0][4:])
					derivative = dN_over_dccmb_tot[index_N0_derivative]
				else:
					continue

				## Take only the non-zero Gaussian-variance
				subblocks_CMBxCMB = ['TTTT', 'EEEE', 'BBBB', 'TTEE', 'TETE', 'TTTE', 'EETE']
				if name_noise[1] in subblocks_CMBxCMB:
					index_CMBxCMB = blocks_CMBxCMB.index(name_noise[1])
					cov = np.diag(cov_order0_CMBxCMB[index_CMBxCMB].data)
				elif name_noise[1][2:]+name_noise[1][:2] in subblocks_CMBxCMB:
					index_CMBxCMB = blocks_CMBxCMB.index(name_noise[1][2:]+name_noise[1][:2])
					cov = np.diag(cov_order0_CMBxCMB[index_CMBxCMB].data.T)
				else:
					continue

				## Update the cross-covariance for the MV
				mat_noise.data += matvec(derivative,cov,lstart=lmin)
			prod = vecmat(vec1,mat_noise.data,lstart=lmin)
			cov_MV[index_CMB].data += prod
			cov_MV_noise[index_CMB].data += prod

			## Trispectrum type-A
			mat_trispA = covmat(0,lmax)
			for name_trispA in names_trispA:
				if name_trispA[0] in blocks_for_phi_derivatives:
					index_N0_derivative = blocks_for_phi_derivatives.index(name_trispA[0])
					derivative = dN_over_dccmb_tot[index_N0_derivative]
				elif name_trispA[0][2:4]+name_trispA[0][:2]+name_trispA[0][4:] in blocks_for_phi_derivatives:
					index_N0_derivative = blocks_for_phi_derivatives.index(name_trispA[0][2:4]+name_trispA[0][:2]+name_trispA[0][4:])
					derivative = dN_over_dccmb_tot[index_N0_derivative]
				else:
					continue

				## Take all CMB blocks
				if name_trispA[1] in blocks_CMBxCMB:
					index_CMBxCMB = blocks_CMBxCMB.index(name_trispA[1])
					cov = cov_order1_CMBxCMB[index_CMBxCMB].data + cov_order2_CMBxCMB[index_CMBxCMB].data
				elif name_trispA[1][2:]+name_trispA[1][:2] in blocks_CMBxCMB:
					index_CMBxCMB = blocks_CMBxCMB.index(name_trispA[1][2:]+name_trispA[1][:2])
					cov = cov_order1_CMBxCMB[index_CMBxCMB].data.T + cov_order2_CMBxCMB[index_CMBxCMB].data.T
				else:
					continue

				## Update the cross-covariance for the MV
				mat_trispA.data += util.dot(derivative,cov)
			prod = vecmat(vec1,mat_trispA.data,lstart=lmin)
			cov_MV[index_CMB].data += prod
			cov_MV_trispA[index_CMB].data += prod

			## Trispectrum Type B
			mat_trispB = covmat(0,lmax)
			combinations_trispB = combination_trispectrumB(block)
			for combination_trispB in combinations_trispB:
				cl_name, n0_name, mat_name = combination_trispB

				## Load the pre-computed trispB matrix
				fn = os.path.join(folder_cache,'trispB_v1_%d_%s_%s.pkl'%(lmax,exp,mat_name))
				if not os.path.isfile(fn):
					print 'No file %s'%fn
					continue
				mat = data_trispB[names_trispB.index(mat_name)]

				## Load the cl
				junk, junk2, cl2_len = lib_spectra.load_weights(cls_lensed, 'cl'+cl_name.lower(),
												noise_uK_arcmin, fwhm_arcmin, lmax, extra='')

				## Load N0
				index_N0 = names_for_N0.index(n0_name)

				vec_left = clpp * N0[index_N0] / (2.0*cls_lensed.ls + 1.0)
				vec_right = cl2_len[1] / (2.0*cls_lensed.ls + 1.0)

				mat_trispB.data += vecmat(vec_left,matvec(mat,vec_right,lstart=lmin), lstart=lmin)

			## Update trispB part
			prod = vecmat(vec1,mat_trispB.data,lstart=lmin)
			cov_MV[index_CMB].data += prod
			cov_MV_trispB[index_CMB].data += prod

			## Signal
			index_CMB_derivative = blocks_for_CMB_derivatives.index(name_signal)

			## Load corresponding dN1_over_dCpp matrix
			tag = '_analytical'

			# N1 matrix from code and samplings
			try:
				Lout, Ls, dLs, dN1_over_dCpp_binned = readN1_matrix(ntag='All',atag=tag,
														matrixtag = block[:4], path_to_N1matrix='N1/%s'%exp)
			except:
				Lout, Ls, dLs, dN1_over_dCpp_binned = readN1_matrix(ntag='All',atag=tag,
														matrixtag = block[2:4]+block[:2], path_to_N1matrix='N1/%s'%exp)

			## N1 matrix operate on and return (L(L+1)]^2C^pp /2pi normalized quantities, so we need to undo this.
			dN1_over_dCpp_binned = np.array([[
						dN1_over_dCpp_binned[list(Lout).index(j)][list(Ls).index(i)] * (i*(i+1.))**2 / (j*(j+1.))**2 for i in Ls] \
						 for j in Lout])
			dN1_over_dCpp_binned = np.nan_to_num(dN1_over_dCpp_binned)

			for i,_ in enumerate(Lout):
				dN1_over_dCpp_binned[i,:] /= dLs

			## Interpolate to all \ells
			matsp = interpolate.RectBivariateSpline(Lout, Ls, dN1_over_dCpp_binned,kx=5,ky=5)
			dN1_over_dCpp = matsp(cls_unlensed.ls,cls_unlensed.ls, grid = True)

			## /!\ /!\ /!\ The fortran routine returns the transposed! /!\ /!\ /!\
			## In fact, it should be dCMB_over_dcpp_tot.T,
			## but as the fortran routine returned the transpose, we keep the matrix as it is.
			mat_simple = vecmat(gaussvar_phiphi_nonoise,dCMB_over_dcpp_tot[index_CMB_derivative],lstart=lmin)
			mat_N1 = util.dot(dN1_over_dCpp,np.diag(gaussvar_phiphi_nonoise),dCMB_over_dcpp_tot[index_CMB_derivative])

			## Update Signal part
			prod = vecmat(vec1,mat_simple+mat_N1,lstart=lmin)
			cov_MV[index_CMB].data += prod
			cov_MV_signal[index_CMB].data += prod

	## well done.
	return cov_MV, cov_MV_signal, cov_MV_noise, cov_MV_trispA, cov_MV_trispB, combinations_CMB

def analytic_covariances_phixphi(cls_unlensed,cls_lensed,lmin=2,noise_uK_arcmin=0.0,fwhm_arcmin=0.0,MPI=None,
						exp='CMB-S4',folder_cache='cache',fn_n1=None):
	'''
	This routine computes the auto covariance (correlation) of the lensing potential power spectrum
	Implementation of Eqs. 40, A19, A20 (plus the connected 8pt function)
	The routine is parallelized in a dump way (and within python), by splitting evenly the ells of the first sum
	within processors. For ell > 3000, better to use many procs (few minutes on 1 core).
	/!\ spin conventions are hardcoded in order to fit with fortran routines (see loop_lensing.f90).
	If you modify the fortran routines, the values of spin should be changed change. /!\
	Input
		* cls_unlensed: camb_clfile object, containing unlensed power spectra.
		* cls_lensed: camb_clfile object, containing lensed power spectra.
		* lmin: int, minimum multipole.
		* noise_uK_arcmin: float, white noise level to be used (in uK.arcmin).
		* fwhm_arcmin: float, beam FWHM to be used.
		* MPI: module, module to pass if you want to parallelize the computation. If None, the computation is done serially.
		* exp: string, name of the experiment you chose according to misc.py
		* folder_cache: string, where to save the data.
		* fn_n1: string, where to find N1 (external, from Biases_n1mat.f90)
	Output
		* cov_MV, cov_MV_RDN0: 2D arrays, auto-covariance matrix for
				the reconstructed lensing potential power spectrum (w/o and w/ RDN0 subtraction).
	'''
	try:
		comm = MPI.COMM_WORLD
		rank = comm.rank
		n_tot = comm.size
		Barrier = comm.Barrier
	except:
		## No MPI
		comm = None
		rank = 0
		n_tot = 1
		Barrier = lambda: -1

	## Maximum multipole
	lmax = int(cls_unlensed.lmax)

	## Input lensing potential power spectrum
	clpp = cls_unlensed.clpp

	## Range of multipoles for processors
	l2range = range(lmin+rank,lmax+1,n_tot)
	l2dim = len(l2range)-1

	#############################################################
	## Load N0 temporarily
	#############################################################
	blocks_for_N0 = ['TTTT_TT','EEEE_EE','BBBB_BB',
					 'EBEB_EE','EBEB_BB',
					 'TETE_TE','TETE_TT','TETE_EE',
					 'TTEE_TE',
					 'TTTE_TT','TTTE_TE',
					 'EETE_EE','EETE_TE']
	N0_tot = np.array([np.zeros(lmax+1) for i in range(len(blocks_for_N0))])
	file_manager_N0 = util.file_manager('N0', exp, spec='v1', lmax=lmax, force_recomputation=False, folder=folder_cache,rank=rank)
	if file_manager_N0.FileExist is True:
		if rank==0: print 'N0 part loaded from %s.pkl'%file_manager_N0.fn
		N0, subblocks = file_manager_N0.data
		for position_block,block in enumerate(blocks_for_N0):
			index_N0 = subblocks.index(block[0:4])
			N0_tot[position_block] = N0[index_N0]
	else:
		if rank ==0: print 'You need N0!'
		sys.exit()

	#############################################################
	## Load derivatives dN^{(0)}/dC^{CMB}
	#############################################################
	blocks_for_derivatives = ['TTTT_TT','EEEE_EE','BBBB_BB',
								  'EBEB_EE','EBEB_BB',
								  'TETE_TE','TETE_TT','TETE_EE',
								  'TTEE_TE',
								  'TTTE_TT','TTTE_TE',
								  'EETE_EE','EETE_TE']
	dN_over_dccmb_tot = np.array([np.zeros((lmax+1,lmax+1)) for i in range(len(blocks_for_derivatives))])

	file_manager_derivatives = util.file_manager('dN_over_dccmb_tot', exp, spec='v1', lmax=lmax, force_recomputation=False, folder=folder_cache,rank=rank)
	if file_manager_derivatives.FileExist is True:
		if rank==0:
			dN_over_dccmb_tot, blocks_for_derivatives = file_manager_derivatives.data
	else:
		if rank==0:
			print 'You need to precompute the derivative of N^{(0)} wrt lensed CMB!'
			sys.exit()
	Barrier()

	#############################################################
	## Load 6-pt functions (cross-covariance matrices) required for the
	## computation of the remaining 4-pt + 6-ptx2-pt function contributions.
	#############################################################
	file_manager_phixCMB = util.file_manager('covariances_phixCMB', exp, spec='v1', lmax=lmax, force_recomputation=False, folder=folder_cache,rank=rank)
	if file_manager_phixCMB.FileExist is True:
		if rank==0: print 'phixCMB part already computed. Loading from %s'%file_manager_phixCMB.fn
		cov_phiCMB_MV, cov_phiCMB_MV_signal, cov_phiCMB_MV_noise, cov_phiCMB_MV_trispA, cov_phiCMB_MV_trispB, combinations_CMB = file_manager_phixCMB.data
	else:
		if rank==0:
			print 'You need to precompute the cov(phi,CMB)!'
			sys.exit()

	#############################################################
	## Loop over pre-computed products, and form covariance matrices.
	#############################################################
	cov_MV = covmat(0,lmax)
	cov_MV_RDN0 = covmat(0,lmax)
	blocks_used = []
	if rank==0:
		## Names of spectra
		combinations = {'TT':0,'EE':1,'EB':2,'TE':3,'TB':4,'BB':5}
		step = len(combinations.keys())
		all_blocks = ['%s%s_%s%s'%(i,j,k,l) for i in combinations.keys() for j in combinations.keys() for k in combinations.keys() for l in combinations.keys()]

		## Form N0 matrix
		file_manager_N0 = util.file_manager('N0', exp, spec='v1', lmax=lmax, force_recomputation=False, folder=folder_cache,rank=rank)
		if file_manager_N0.FileExist is True:
			if rank==0: print 'Loading N0 from %s.pkl'%file_manager_N0.fn
			N0, names_for_N0 = file_manager_N0.data
		else:
			if rank==0: print 'You need N0!'
			sys.exit()

		## Load N1 if necessary
		if fn_n1 is None:
			if rank == 0: print 'No N1 used'
			n1_MV = np.zeros_like(cls_unlensed.ls)
		else:
			if rank == 0: print 'Loading N1 from %s'%fn_n1
			N1 = np.loadtxt(fn_n1).T
			n1_MV = lib_spectra.compute_minimum_variance_N1(cls_unlensed.ls,None,N1,N0,names_for_N0)

		## Form minimum variance for N0 and compute weights
		minimum_variance_n0, weights_for_MV, names_for_weights = lib_spectra.compute_minimum_variance_weights(N0,names_for_N0)

		gaussvar_MV = gaussian_variance(clpp+minimum_variance_n0+n1_MV,ells=cls_unlensed.ls)
		cov_MV.data += gaussvar_MV
		cov_MV_RDN0.data += gaussvar_MV

		#############################################################
		## Connected 8-pt (see App A2 in the paper)
		#############################################################
		cl_len_TT, cl_len_EE, cl_len_TE = lib_spectra.load_weights(cls_lensed, 'clte', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		cl_len_BB, junk1, junk2 = lib_spectra.load_weights(cls_lensed, 'clbb', noise_uK_arcmin, fwhm_arcmin, 2*lmax, extra='_long')
		dic_spectra = {'TT':cl_len_TT[1], 'EE':cl_len_EE[1], 'TE':cl_len_TE[1], 'BB':cl_len_BB[1]}
		for comb in ['TTTT_TTTT','EEEE_EEEE','BBBB_BBBB','TETE_TETE']:
			## Indices for blocks
			index_N0 = blocks_for_N0.index(comb[0:7])
			index_derivative = blocks_for_derivatives.index(comb[0:7])

			## SNR and MV vectors
			vec_amplitude = clpp / N0[index_N0]
			vec1 = weights_for_MV[names_for_weights.index(comb[0:2])] *  weights_for_MV[names_for_weights.index(comb[5:7])]

			## Build the cov
			gaussvarCMB = cross_gaussian_variance(cl1=dic_spectra[comb[0:2]],cl2=dic_spectra[comb[0:2]],ells=cls_lensed.ls)
			individual_covariance = util.dot(dN_over_dccmb_tot[index_derivative],gaussvarCMB,dN_over_dccmb_tot[index_derivative].T)
			prod = 8. * vecmat(vec1*vec_amplitude,matvec(individual_covariance,vec1*vec_amplitude,lstart=lmin),lstart=lmin)

			## Update the MV
			cov_MV.data += prod
			cov_MV_RDN0.data += prod

		for comb in ['EBEB_EBEB']:
			## Indices for blocks
			index_N0 = blocks_for_N0.index('EBEB_EE')
			index_derivative = blocks_for_derivatives.index('EBEB_EE')

			## SNR and MV vectors
			vec_amplitude = clpp / N0[index_N0]
			vec1 = weights_for_MV[names_for_weights.index(comb[0:2])] *  weights_for_MV[names_for_weights.index(comb[5:7])]

			## Build the cov
			gaussvarCMB = cross_gaussian_variance(cl1=dic_spectra['EE'],cl2=dic_spectra['EE'],ells=cls_lensed.ls)
			individual_covariance = util.dot(dN_over_dccmb_tot[index_derivative],gaussvarCMB,dN_over_dccmb_tot[index_derivative].T)
			prod = 8. * vecmat(vec1*vec_amplitude,matvec(individual_covariance,vec1*vec_amplitude,lstart=lmin),lstart=lmin)

			## Update
			cov_MV.data += prod
			cov_MV_RDN0.data += prod


		#############################################################
		## disconnected 8-pt (see Eq. 44 in the paper)
		#############################################################
		print 'Noise'
		for c1 in combinations.keys(): ## XY
			XY_index = combinations[c1]
			wXY_index = names_for_weights.index(c1)
			for c2 in combinations.keys(): ## ZW
				ZW_index = combinations[c2]
				wZW_index = names_for_weights.index(c2)

				## Weights (XY * ZW)
				vec1 = weights_for_MV[wXY_index] *  weights_for_MV[wZW_index]

				for CMB in ['TT','EE','TE','BB']:
					index_CMB = combinations_CMB.index(CMB)
					der_name = '%s%s_%s'%(c1,c2, CMB)
					if der_name in blocks_for_derivatives:
						index_derivative = blocks_for_derivatives.index(der_name)
					elif der_name[2:4]+der_name[0:2]+der_name[4:] in blocks_for_derivatives:
						index_derivative = blocks_for_derivatives.index(der_name[2:4]+der_name[0:2]+der_name[4:])
					else:
						continue

					print der_name, CMB

					## Just the noise
					cov_nonoise = cov_phiCMB_MV_noise[index_CMB].data

					mat1 = matvec(util.dot(cov_nonoise,dN_over_dccmb_tot[index_derivative].T),vec1,lstart=lmin)

					## Update the cross-covariance for the MV
					cov_MV.data += mat1

		#############################################################
		## Remaining 4-pt + 6-ptx2-pt function contribution (See eq. 42 in the paper)
		#############################################################
		print '6x2pt contribution'
		for c1 in combinations.keys(): ## XY
			XY_index = combinations[c1]
			wXY_index = names_for_weights.index(c1)
			for c2 in combinations.keys(): ## ZW
				ZW_index = combinations[c2]
				wZW_index = names_for_weights.index(c2)
				for CMB in [c1[0]+c2[0],c1[0]+c2[1],c1[1]+c2[0],c1[1]+c2[1]]:
					if CMB in ['EB','BE','TB','BT']: continue
					if CMB == 'ET': CMB = 'TE'

					index_CMB = combinations_CMB.index(CMB)
					der_name = '%s%s_%s'%(c1,c2, CMB)
					if der_name in blocks_for_derivatives:
						index_derivative = blocks_for_derivatives.index(der_name)
					elif der_name[2:4]+der_name[0:2]+der_name[4:] in blocks_for_derivatives:
						index_derivative = blocks_for_derivatives.index(der_name[2:4]+der_name[0:2]+der_name[4:])
					else:
						continue

					print der_name, CMB

					## Weights (XY * ZW)
					vec1 = weights_for_MV[wXY_index] *  weights_for_MV[wZW_index]

					## phixCMB without the noise
					cov_nonoise = cov_phiCMB_MV_signal[index_CMB].data + cov_phiCMB_MV_trispA[index_CMB].data + cov_phiCMB_MV_trispB[index_CMB].data

					mat1 = vecmat(vec1,util.dot(dN_over_dccmb_tot[index_derivative],cov_nonoise.T),lstart=lmin)
					mat2 = vecmat(vec1,util.dot(dN_over_dccmb_tot[index_derivative],cov_nonoise.T),lstart=lmin).T

					## Update the cross-covariance for the MV
					cov_MV.data += mat1 + mat2

	Barrier()
	return cov_MV, cov_MV_RDN0, np.unique(blocks_used)

########################################################################
# Manipulate manually matrices and vectors

def phi_to_kappa(mat,lstart,lstop,ellrange=None,spectrum=True):
	'''
	Perform change of variable for matrix elements from phi to kappa (factor l*(l+1)).
	Input:
		* mat: 2D-array, the matrix containing elements related to phi
		* lstart: int, minimum multipole
		* lstop: int, maximum multipole
		* spectrum: boolean, if True it applies [l*(l+1)]**2 instead
	Output:
		* data_norm: 2D-array, the inital matrix normalized
	'''
	if spectrum is True:
		norm_phi = lambda l: (l*(l+1))**2
	else:
		norm_phi = lambda l: l*(l+1)

	if ellrange is None:
		data_norm = np.array([mat[j-lstart] * norm_phi(j)*norm_phi(np.arange(lstart,lstop)) for j in range(lstart,lstop)])
	else:
		data_norm = np.array([mat[pos] * norm_phi(j)*norm_phi(ellrange) for pos,j in enumerate(ellrange)])
	return data_norm

def manual_correlation_matrix(mat, vec1, vec2, lstart=0,lstop=None):
	'''
	Compute manually correlation matrix from covariance matrix.
	See covmat class for automatic usage
	The code makes use of C to speed up the computation.
	Input:
		* lstart: int, minimum multipole
		* lstop: int, maximum multipole
		* mat: 2D array, matrix
		* vec1: 1D array
		* vec2: 1D array
	Output
		* corrmat: 2D array, correlation matrix.

	'''
	if lstop is None:
		lstop = len(vec1)

	corrmat = np.zeros_like(mat)

	code=r'''
	int l1,l2;
	for(l1 = lstart; l1 < lstop; l1++) {
		for(l2 = lstart; l2 < lstop; l2++) {
			corrmat(l1,l2) = mat(l1,l2) / sqrt(vec1(l1)*vec2(l2));
		}
	}
	'''

	inline(code,['mat','corrmat','vec1','vec2','lstart','lstop'], headers=["<math.h>"],type_converters = converters.blitz)

	return corrmat

def simple_matvec_c(vec,mat,fortran_conv=0):
	'''
	Simply compute the vector-matrix product (which is in fact a mat-mat product) vec_l2 * mat_l2l3 in C
	Input
		* vec: 1D array, the vector
		* mat: 2D array, the matrix
		* fortran_conv: boolean, use fortran convention (col major). TBI
	Output
		* matvec: 2D array, the product vec_l2 * mat_l2l3
	'''
	size = len(vec)
	matvec = np.zeros_like(mat)
	code = r'''
	int l2,l3;
	for(l2 = 2; l2 < size; l2++) {
		for(l3 = 2; l3 < size; l3++) {
			matvec(l2,l3) = vec(l2)*mat(l2,l3);
		}
	}
	'''
	inline(code,['matvec','vec','mat','size','fortran_conv'], headers=["<math.h>"],type_converters = converters.blitz)
	return matvec

def vecmat(vec,mat,lstart=False):
	'''
	Performs V_l * M_{ll'}
	Input:
		* vec: 1D array, input vector
		* mat: 2D array, input matrix
		* lstart: int, (optional) the starting ell.
	Output:
		* expanded_vec * mat: 2D array
	'''
	size = len(vec)
	expanded_vec = np.tile(vec,(size,1)).T
	if lstart:
		expanded_vec[0:lstart]=0.0
	return expanded_vec * mat

def matvec(mat,vec,lstart=False):
	'''
	Performs M_{ll'} * V_l'
	Input:
		* mat: 2D array, input matrix
		* vec: 1D array, input vector
		* lstart: int, (optional) the starting ell.
	Output:
		* expanded_vec * mat: 2D array
	'''
	size = len(vec)
	# expanded_vec = np.tile(vec,(size,1))
	expanded_vec = np.tile(vec,(size,1)).T
	if lstart:
		expanded_vec[0:lstart]=0.0
	return expanded_vec.T * mat

def matvecmat(mat1,vec,mat2,lstart=None):
	'''
	Computes M[i][j] = sum_l M1[i][l]*V[l]*M2[l][j]
	'''
	vecmat_prod = vecmat(vec,mat2,lstart=lstart)
	matvecmat_prod = np.dot(mat1,vecmat_prod)

	return matvecmat_prod

def compute_covmat_auto_sims(datavec):
	'''
	Auto-covariance for sims.
	Input:
		* datavec: 2D-array, contains the sims [N_sims , N_bins]
	Output:
		* covmat: 2D-array, covariance matrix
	'''
	nsims = datavec.shape[0]
	size  = datavec.shape[1]

	meanvec = np.mean(datavec,axis=0)

	X = datavec - meanvec
	covmat = np.cov(X.T)

	return covmat

def compute_covmat_cross_sims(lensingvec,CMBvec):
	'''
	Cross-covariance for sims.
	Input:
		* lensingvec: 2D-array, contains the sims of the lensing potential [N_sims , N_bins]
		* CMBvec: 2D-array, contains the sims of the lensed CMB [N_sims , N_bins]
	Output:
		* covmat: 2D-array, cross-covariance matrix
	'''
	nsims = lensingvec.shape[0]
	size_CMB  = CMBvec.shape[1]
	size_lensing  = lensingvec.shape[1]

	meanlensing = np.mean(lensingvec,axis=0)
	meanCMB = np.mean(CMBvec,axis=0)
	covmat_cross = np.zeros((size_lensing,size_CMB))

	code=r'''
	int s,i,j;
	for(s = 0; s < nsims; s++) {
		for(i = 0; i < size_lensing; i++) {
			for(j = 0; j < size_CMB; j++) {
				covmat_cross(i,j) = covmat_cross(i,j) + (lensingvec(s,i) - meanlensing(i))*(CMBvec(s,j) - meanCMB(j));
			}
		}
	}
	'''

	inline(code,['covmat_cross','lensingvec','meanlensing','CMBvec','meanCMB','nsims','size_lensing','size_CMB'], headers=["<math.h>"],type_converters = converters.blitz)

	return covmat_cross / (nsims-1)

def is_covmat(obj):
	""" check (by ducktyping) if obj is a covmat """
	if not hasattr(obj, 'lmax'):
		return False
	if not hasattr(obj, 'size'):
		return False
	return set(obj.__dict__.keys()).issubset( set( ['lmax', 'lmin', 'size', 'data' ] ) )

def is_covmatbinned(obj):
	""" check (by ducktyping) if obj is a covmat """
	if not hasattr(obj, 'lbins'):
		return False
	if not hasattr(obj, 'lboundaries'):
		return False
	return set(obj.__dict__.keys()).issubset( set( ['lbins', 'lboundaries', 'size', 'data' ] ) )

########################################################################
# Class for covariance matrices

class covmat(object):
	""" Class to handle covariance matrices """
	def __init__(self, lmin, lmax):
		"""
		Initialization of the covmat class
		Input
			 * lmax: int, maximum multipole
			 * lmin: int, minimum multipole
		"""

		self.lmax = lmax
		self.lmin = lmin
		self.size = lmax - lmin + 1
		self.data = np.zeros((self.size,self.size))

	def copy(self):
		"""
		Clone this object.
		"""
		return copy.deepcopy(self)

	def __add__(self, other):
		"""
		Sum two covmat objects. size must be the same for both.
		Input
			* other: 2D array, the matrix that you want to add
		"""
		if is_covmat(other):
			assert( self.size == other.size )
			ret = self.copy()
			zs  = np.zeros((self.size,self.size))
			for attr in ['data']:
				if (hasattr(self, attr) or hasattr(other, attr)):
					setattr(ret, attr, getattr(self, attr, zs) + getattr(other, attr, zs) )
			return ret
		else:
			assert(0)

	def plot(self, p=pl.imshow, **kwargs):
		"""
		Plot the matrix (use imshow by default)
		Input
			* p: function, plotting function to use p(x,y,**kwargs)
		"""
		p(self.data, **kwargs )

	def correlation_matrix(self, lmin=2, remove_diag=False):
		'''
		Compute correlation matrix from covariance matrix.
		The code makes use of C to speed up the computation.
		Input:
			* lmin: int, minimum multipole
			* remove_diag: boolean, return the correlation matrix with zero on the diagonal
		Output
			* corrmat: 2D array, correlation matrix. If remove_diag is True, the diagonal elements are subtracted.

		'''
		covmat = self.data
		size = self.size

		corrmat = np.zeros_like(covmat)

		code=r'''
		int l1,l2;
		for(l1 = lmin; l1 < size; l1++) {
			for(l2 = lmin; l2 < size; l2++) {
				corrmat(l1,l2) = covmat(l1,l2) / sqrt(covmat(l1,l1)*covmat(l2,l2));
			}
		}
		'''

		inline(code,['covmat','corrmat','size','lmin'], headers=["<math.h>"],type_converters = converters.blitz)

		if remove_diag:
			toremove = np.eye(size)
			toremove[0:lmin]=0.0
			return corrmat - toremove
		else:
			return corrmat

	def cross_correlation_matrix(self, vec1, vec2, lmin=2):
		'''
		Compute cross-correlation matrix from cross-covariance matrix between 1 and 2.
		The self.data vector should contains the covmat(1,2).
		The code makes use of C to speed up the computation.
		Input:
			* lmin: int, minimum multipole
			* vec1: 1D array, the gaussian variance of 1
			* vec2: 1D array, the gaussian variance of 2
		Output
			* corrmat: 2D array, cross-correlation matrix.

		'''
		covmat = self.data
		size = self.size

		corrmat = np.zeros_like(covmat)

		code=r'''
		int l1,l2;
		for(l1 = lmin; l1 < size; l1++) {
			for(l2 = lmin; l2 < size; l2++) {
				corrmat(l1,l2) = covmat(l1,l2) / sqrt(vec1(l1)*vec2(l2));
			}
		}
		'''

		inline(code,['covmat','corrmat','vec1','vec2','size','lmin'], headers=["<math.h>"],type_converters = converters.blitz)

		return corrmat

	def svd_decomposition(self,mat=None,transpose=False):
		'''
		Perform SVD decomposition of the covariance matrix (or correlation matrix)
		By default, it performs the SVD of the covariance matrix (self.data), but you
		can specified the covariance (or any other matrix) via the keyword mat.
		You can reconstruct mat = np.dot(U, np.dot(S, V)), where S = np.diag(s).
		Input:
			(optional) mat: 2-D array, a matrix
		Output:
			U: 2-D array, left singular vectors
			s: 1-D array, The singular values for every singular vector, sorted in descending order
			V: 2-D array, right singular vectors
		'''
		if mat is None:
			matrix = self.data
		else:
			matrix = mat

		if transpose:
			matrix = matrix.T

		U, s, V = np.linalg.svd(matrix)

		return U, s, V

	def normalize(self, transpose=True):
		'''
		Normalize with l(l+1)/2pi and (l(l+1))**2/2pi
		'''
		data_norm = np.zeros_like(self.data)
		norm_CMB 		= lambda l: l*(l+1)/(2.*np.pi)
		norm_phi 		= lambda l: (l*(l+1))**2/(2.*np.pi)

		if transpose:
			data_norm = np.array([self.data.T[j] * norm_CMB(j)*norm_phi(np.arange(self.size)) for j in range(self.size)])
		else:
			data_norm = np.array([self.data[j] * norm_CMB(j)*norm_phi(np.arange(self.size)) for j in range(self.size)])
		return data_norm

	def denormalize(self, transpose=True):
		'''
		Normalize with l(l+1)/2pi and (l(l+1))**2/2pi
		'''
		data_norm = np.zeros_like(self.data)
		norm_CMB 		= lambda l: l*(l+1)/(2.*np.pi)
		norm_phi 		= lambda l: (l*(l+1))**2/(2.*np.pi)

		if transpose:
			data_norm = np.array([self.data.T[j] / norm_CMB(j)/norm_phi(np.arange(self.size)) for j in range(self.size)])
		else:
			data_norm = np.array([self.data[j] / norm_CMB(j)/norm_phi(np.arange(self.size)) for j in range(self.size)])
		return data_norm

class covmatbinned(object):
	""" Class to handle binned covariance matrices (for simulations) """
	def __init__(self, lbins, lboundaries,lbins2=None,lboundaries2=None):
		"""
		Initialization of the covmatbinned class
		Input
			 * lbins: 1D array, center of the bins
			 * lboundaries: 1D array, boundaries of the bins [l0min, l0max=l1min, ...,lN-1max=lN-1min,lNmax]
		"""

		if lbins2 is not None:
			if lboundaries2 is None:
				print 'You have to specify boudaries for the 2nd dimension'
				sys.exit()
			self.lbins1 = lbins
			self.lbins2 = lbins2
			self.lboundaries1 = lboundaries
			self.lboundaries2 = lboundaries2
			self.size1 = len(lbins)
			self.size2 = len(lbins2)
			self.data = np.zeros((self.size1,self.size2))
		else:
			self.lbins1 = lbins
			self.lbins2 = lbins
			self.lboundaries1 = lboundaries
			self.lboundaries2 = lboundaries
			self.size1 = len(lbins)
			self.size2 = len(lbins)
			self.data = np.zeros((self.size1,self.size2))

	def copy(self):
		"""
		Clone this object.
		"""
		return copy.deepcopy(self)

	def __add__(self, other):
		"""
		Sum two covmatbinned objects. size must be the same for both.
		Input
			* other: 2D array, the matrix that you want to add
		"""
		if is_covmat(other):
			assert( self.size == other.size )
			ret = self.copy()
			zs  = np.zeros((self.size,self.size))
			for attr in ['data']:
				if (hasattr(self, attr) or hasattr(other, attr)):
					setattr(ret, attr, getattr(self, attr, zs) + getattr(other, attr, zs) )
			return ret
		else:
			assert(0)

	def plot(self, p=pl.imshow, **kwargs):
		"""
		Plot the matrix (use imshow by default)
		Input
			* p: function, plotting function to use p(x,y,**kwargs)
		"""
		p(self.data, **kwargs )

	def correlation_matrix(self, startbin=0, remove_diag=False):
		'''
		Compute correlation matrix from covariance matrix.
		The code makes use of C to speed up the computation.
		Input:
			* startbin: int, the first index for the calculation
			* remove_diag: boolean, return the correlation matrix with zero on the diagonal
		Output
			* corrmat: 2D array, correlation matrix. If remove_diag is True, the diagonal elements are subtracted.

		'''
		covmat = self.data
		size1 = self.size1
		size2 = self.size2
		if size1 != size2:
			print 'Matrix has to be square!'
			sys.exit()

		corrmat = np.zeros_like(covmat)

		code=r'''
		int l1,l2;
		for(l1 = startbin; l1 < size1; l1++) {
			for(l2 = startbin; l2 < size2; l2++) {
				corrmat(l1,l2) = covmat(l1,l2) / sqrt(covmat(l1,l1)*covmat(l2,l2));
			}
		}
		'''

		inline(code,['covmat','corrmat','size1','size2','startbin'], headers=["<math.h>"],type_converters = converters.blitz)

		if remove_diag:
			toremove = np.eye(size1)
			# toremove[0:lmin]=0.0
			return corrmat - toremove
		else:
			return corrmat

	def cross_correlation_matrix(self, vec1, vec2, startbin=0):
		'''
		Compute cross-correlation matrix from cross-covariance matrix between 1 and 2.
		The self.data vector should contains the covmat(1,2).
		The code makes use of C to speed up the computation.
		Input:
			* vec1: 1D array, the gaussian variance of 1
			* vec2: 1D array, the gaussian variance of 2
			* startbin: int, the first index for the calculation
		Output
			* corrmat: 2D array, cross-correlation matrix.

		'''
		covmat = self.data
		size1 = self.size1
		size2 = self.size2

		corrmat = np.zeros_like(covmat)

		code=r'''
		int l1,l2;
		for(l1 = startbin; l1 < size1; l1++) {
			for(l2 = startbin; l2 < size2; l2++) {
				corrmat(l1,l2) = covmat(l1,l2) / sqrt(vec1(l1)*vec2(l2));
			}
		}
		'''

		inline(code,['covmat','corrmat','vec1','vec2','size1','size2','startbin'], headers=["<math.h>"],type_converters = converters.blitz)

		return corrmat

	def svd_decomposition(self,mat=None,transpose=False):
		'''
		Perform SVD decomposition of the covariance matrix (or correlation matrix)
		By default, it performs the SVD of the covariance matrix (self.data), but you
		can specified the covariance (or any other matrix) via the keyword mat.
		You can reconstruct mat = np.dot(U, np.dot(S, V)), where S = np.diag(s).
		Input:
			(optional) mat: 2-D array, a matrix
		Output:
			U: 2-D array, left singular vectors
			s: 1-D array, The singular values for every singular vector, sorted in descending order
			V: 2-D array, right singular vectors
		'''
		if mat is None:
			matrix = self.data
		else:
			matrix = mat

		if transpose:
			matrix = matrix.T

		U, s, V = np.linalg.svd(matrix)

		return U, s, V

	def normalize(self, lbins_phi,lbins_CMB, transpose=True,allphi=False,denorm=False):
		'''
		Normalize with l(l+1)/2pi and (l(l+1))**2/2pi
		Input:
			* lbins: 1D array, center of the bins
		'''
		data_norm = np.zeros_like(self.data)
		if denorm:
			norm_CMB 		= lambda l: 1./(l*(l+1)/(2.*np.pi))
			norm_phi 		= lambda l: 1./((l*(l+1))**2/(2.*np.pi))
		else:
			norm_CMB 		= lambda l: l*(l+1)/(2.*np.pi)
			norm_phi 		= lambda l: (l*(l+1))**2/(2.*np.pi)

		if allphi:
			norm_CMB = norm_phi

		if transpose:
			data_norm = np.array([self.data.T[pos] * norm_CMB(j)*norm_phi(lbins_phi) for pos,j in enumerate(lbins_CMB)])
		else:
			data_norm = np.array([self.data[pos] * norm_CMB(j)*norm_phi(lbins_phi) for pos,j in enumerate(lbins_CMB)])
		return data_norm

	def denormalize(self, lbins_phi,lbins_CMB, quantity='CMB', transpose=True):
		'''
		Normalize with l(l+1)/2pi and (l(l+1))**2/2pi
		Input:
			* lbins: 1D array, center of the bins
		'''
		data_norm = np.zeros_like(self.data)
		if quantity == 'CMB':
			norm1 		= lambda l: l*(l+1)/(2.*np.pi)
			norm2 		= lambda l: l*(l+1)/(2.*np.pi)
		elif quantity == 'lensing':
			norm1 		= lambda l: (l*(l+1))**2/(2.*np.pi)
			norm2 		= lambda l: (l*(l+1))**2/(2.*np.pi)
		elif quantity == 'cross':
			norm1 		= lambda l: l*(l+1)/(2.*np.pi)
			norm2 		= lambda l: (l*(l+1))**2/(2.*np.pi)

		if transpose:
			data_norm = np.array([self.data.T[pos] / norm1(j)/norm2(lbins_phi) for pos,j in enumerate(lbins_CMB)])
		else:
			data_norm = np.array([self.data[pos] / norm1(j)/norm2(lbins_phi) for pos,j in enumerate(lbins_CMB)])
		return data_norm

########################################################################
# Misc

def readN1_matrix(ntag, atag='_S4_Noise1.5_Beam3', matrixtag = None,path_to_N1matrix=''):
	'''
	Read the derivative of N1 wrt Cpp from Antony Lewis's N1 derivative code.
	Input:
		* ntag: string, tag from Antony's code
		* atag: string, tag containing experiment name, noise and beam level
		* matrixtag: string, name of the combination (e.g. TTTT)
		* path_to_N1matrix: string, path to the data + 'N1_'
	Output:
		* Lout: 1D array, the outputed bins
		* Ls: 1D array, the input bins
		* Lin[:,1], 1D array, binsize
		* mat: 2D array, matrix of size (Lout, Ls)
	'''
	if matrixtag is None:
		matrixtag = ntag

	## Load the matrix
	mat = np.loadtxt(os.path.join(path_to_N1matrix,'N1_' + matrixtag + atag + '_matrix.dat'))[1:,1:]

	## Load the bins
	Lout = np.loadtxt(os.path.join(path_to_N1matrix,'N1_' + ntag + atag + '_Lout.dat')).astype(int)

	## Load the input bins and the binsize
	Lin = np.loadtxt(os.path.join(path_to_N1matrix,'N1_' + ntag + atag + '_Lin.dat'))
	Ls = Lin[:, 0].astype(int)

	return Lout, Ls, Lin[:, 1], mat
