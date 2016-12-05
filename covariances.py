# Copyright (C) 2016 Peloton
#########################
# Main script to compute CMB and lensing covariances
# author: julien@sussex
# See 1611.01446
#########################
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pylab as pl
pl.ioff()

import sys,os,glob
import argparse

import lib_covariances, lib_spectra
import misc, util

try:
	import LensingBiases as LB
except:
	print 'You need to install the lensingbiases package'
	print 'to compute N1. See '
	print 'https://github.com/JulienPeloton/lensingbiases'

font = {'family' : 'monospace',
			  'weight' : 'bold',
			  'size'   : 15}
pl.rc('font', **font)

## Define
DEBUG = False

def addargs(parser):
	''' Parse command line arguments '''
	parser.add_argument('-input_unlensed_spectra',dest='input_unlensed_spectra',help='Input unlensed spectra (CAMB)',required=True)
	parser.add_argument('-input_lensed_spectra',dest='input_lensed_spectra',help='Input lensed spectra (CAMB)',required=True)
	parser.add_argument('-exp',dest='exp',help='Instrumental configuration: Planck, Core++, CMB-S4, Ideal',required=True)
	parser.add_argument('-runmode',dest='runmode',help='What do you want to compute?',required=True)
	parser.add_argument('--mpi',dest='mpi',help='Run using mpi (mpi4py required)',action='store_true')

def grabargs(args_param=None):
	''' Parse command line arguments '''
	parser = argparse.ArgumentParser(description='Main script to compute auto- and cross-covariances.')
	addargs(parser)
	args = parser.parse_args(args_param)
	return args

def main():
	args_param = None
	args   = grabargs(args_param)
	covariances_main(args)

def covariances_main(args):
	'''
	Example script which uses functions from the package.
	'''
	if args.mpi:
		from mpi4py import MPI
		comm = MPI.COMM_WORLD
		rank = comm.rank
		barrier = comm.barrier
		if DEBUG: print "Hello! I'm rank %d from %d running in total..." % (comm.rank, comm.size)
	else:
		MPI = None
		rank = 0
		barrier = lambda: -1

	MODE = args.runmode
	exp = args.exp
	noise_uK_arcmin, fwhm_arcmin, lmin, lmax, folder_cache = misc.get_exp_configuration(exp)

	if rank == 0:
		if DEBUG: print args
		print '+-----------------------------+'
		print '+ Exp:			   ',exp
		print '+ Noise (uK.arcmin): ',noise_uK_arcmin
		print '+ FWHM (arcmin):	 ',fwhm_arcmin
		print '+ lmin:			  ',lmin
		print '+ lmax:			  ',lmax
		print '+ MODE:			  ',MODE
		print '+-----------------------------+'

	## Initialization of spectra
	cls_unlensed = lib_spectra.get_camb_cls(fname=args.input_unlensed_spectra,lmax=lmax)
	cls_lensed = lib_spectra.get_camb_cls(fname=args.input_lensed_spectra,lmax=lmax)

	########################################################################
	# List of different options

	if MODE == 'N0':
		'''
		Compute N0 bias.
		For more informations, see lib_spectra.compute_N0_XYXY()
		The blocks have to be of the form XXYY, where X,Y = TT, EE, BB, TE, TB, EB.
		If XX != YY, the code computes the general formula for N0 (see Eq. 17 in 1611.01446).
		Otherwise, it uses the fact that N0 = A^{XX} (for optimal weights).
		'''
		## The available blocks in the code.
		blocks = ['TTTT','EEEE','BBBB','TETE','TBTB','EBEB','TTEE','TTTE','EETE','EBTB']

		## Initialization of file manager
		file_manager = util.file_manager(MODE, exp, 'v1', lmax, force_recomputation=False, folder=folder_cache,rank=rank)

		if file_manager.FileExist is True:
			## Load data instead of recomputing everything.
			N0, blocks = file_manager.data
		else:
			N0, blocks = lib_spectra.compute_N0_XYWZ(cls_lensed, lmin=lmin, blocks=blocks,
				noise_uK_arcmin=noise_uK_arcmin, fwhm_arcmin=fwhm_arcmin, MPI=MPI)
			array_to_save = [N0, blocks]

	if MODE == 'N1':
		'''
		Compute N1 bias.
		It uses the lensingbiases package (see https://github.com/JulienPeloton/lensingbiases)
		The output is saved on the disk, and used later on.
		'''

		## Initialization of file manager
		file_manager = util.file_manager(MODE, exp, 'v1', lmax, force_recomputation=False, folder=folder_cache,rank=rank)

		LB.checkproc_py()

		## We need to compute N0 for internal purposes (flat-sky, so not used afterwards)
		bins, phiphi, n0_mat, indices = LB.compute_n0_py(from_args=None,phifile=args.input_unlensed_spectra,
							lensedcmbfile=args.input_lensed_spectra,
							FWHM=fwhm_arcmin,noise_level=noise_uK_arcmin,
							lmin=lmin,lmaxout=lmax,lmax=lmax,lmax_TT=lmax,
							tmp_output='N1/%s'%exp)

		## Compute N1, and derivatives of N1 wrt lensing potential power-spectrum
		LB.compute_n1_derivatives_py(from_args=None,phifile=args.input_unlensed_spectra,
							lensedcmbfile=args.input_lensed_spectra,
							FWHM=fwhm_arcmin,noise_level=noise_uK_arcmin,
							lmin=lmin,lmaxout=lmax,lmax=lmax,lmax_TT=lmax,
							tmp_output='N1/%s'%exp)
		sys.exit()

	elif MODE == 'covariances_CMBxCMB':
		'''
		Compute different parts of the CMB auto-covariance.
		'''
		## The available blocks in the code.
		blocks = ['TTTT', 'EEEE', 'BBBB', 'EEBB', 'TTEE', 'TTBB', 'TETE', 'TTTE', 'EETE', 'TEBB']
		file_manager = util.file_manager(MODE, exp, spec='v1', lmax=lmax,
									force_recomputation=False, folder=folder_cache,rank=rank)

		if file_manager.FileExist is True:
			cov_order0_tot, cov_order1_tot, cov_order2_tot, junk = file_manager.data
		else:
			cov_order0_tot, cov_order1_tot, cov_order2_tot, junk = lib_covariances.analytic_covariances_CMBxCMB(cls_unlensed,
					cls_lensed,lmin=lmin,blocks=blocks,
					noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin, MPI=MPI,
					use_corrfunc=True,exp=exp,folder_cache=folder_cache)
			array_to_save = [cov_order0_tot, cov_order1_tot, cov_order2_tot, blocks]

	elif MODE == 'covariances_phixCMB':
		'''
		Compute different parts of the cross-covariance
		between the lensing potential power-spectrum and the lensed CMB.
		'''
		## Initialization of file manager
		path_to_N1matrix = 'N1/%s'%(exp)
		file_manager = util.file_manager(MODE, exp, spec='v1', lmax=lmax,
									force_recomputation=False, folder=folder_cache,rank=rank)

		if file_manager.FileExist is True:
			cov_MV, cov_MV_signal, cov_MV_noise, cov_MV_trispA, cov_MV_trispB, combinations_CMB = file_manager.data
		else:
			cov_MV, cov_MV_signal, cov_MV_noise, cov_MV_trispA, cov_MV_trispB, combinations_CMB = \
							lib_covariances.analytic_covariances_phixCMB(cls_unlensed,cls_lensed,lmin=lmin,
							noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,MPI=MPI,use_corrfunc=True,exp=exp,
							folder_cache=folder_cache,path_to_N1matrix=path_to_N1matrix)
			array_to_save = [cov_MV, cov_MV_signal, cov_MV_noise, cov_MV_trispA, cov_MV_trispB, combinations_CMB]

	elif MODE == 'covariances_phixphi':
		'''
		Compute different parts of the auto-covariance
		for the reconstructed lensing potential power spectrum.
		'''
		## Initialization of file manager
		fn_n1 = 'N1/%s/N1_All_analytical.dat'%(exp)

		file_manager = util.file_manager(MODE, exp, spec='v1', lmax=lmax,
								force_recomputation=False, folder=folder_cache,rank=rank)

		if file_manager.FileExist is True:
			cov_MV, cov_RDN0_MV, blocks = file_manager.data
		else:
			cov_MV, cov_RDN0_MV, blocks = lib_covariances.analytic_covariances_phixphi(cls_unlensed,cls_lensed,lmin=lmin,
							noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,
							MPI=MPI,exp=exp,folder_cache=folder_cache,fn_n1=fn_n1)

			array_to_save = [cov_MV, cov_RDN0_MV, blocks]

	elif MODE == 'trispB':
		'''
		Section to compute trispectrum B in a different manner.
		We compute here the terms sum_l g_XY f_ZW. Those terms are then combined in lib_covariances
		as defined in Eq. 37 in 1611.01446.
		'''
		## Initialization of file manager
		file_manager = util.file_manager(MODE, exp, spec='v1', lmax=lmax,
								force_recomputation=False, folder=folder_cache,rank=rank)

		## This is just the upper triangle.
		## Transposed terms are computed as well.
		combinations = ['TTTT', 'TTEE', 'TTBB', 'TTTE', 'TTET',
						'EEEE', 'EEBB', 'EETE', 'EEET',
						'BBBB', 'BBTE', 'BBET',
						'TETE', 'TEET',
						'ETET',
						'EBEB', 'EBBE',
						'BEBE',
						'TBTB', 'TBBT',
						'BTBT',
						'EBTB','BETB','BEBT']
		for combination in combinations:
			trispB = lib_covariances.precompute_trispB(cls_unlensed,cls_lensed,combination=combination,lmin=lmin,
								noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,
								MPI=MPI,exp=exp,folder_cache=folder_cache)

			array_to_save = [trispB, combination]

			if rank==0:
				file_manager.save_data_on_disk(array_to_save,name=os.path.join(folder_cache,
									file_manager.basefn + '_%s.pkl'%combination))

			## Compute transposed terms as well
			if combination[0:2] != combination[2:4]:
				trispB = lib_covariances.precompute_trispB(cls_unlensed,cls_lensed,
									combination=combination[2:4]+combination[:2],lmin=lmin,
									noise_uK_arcmin=noise_uK_arcmin,fwhm_arcmin=fwhm_arcmin,
									MPI=MPI,exp=exp,folder_cache=folder_cache)

				array_to_save = [trispB, combination]

				if rank==0:
					file_manager.save_data_on_disk(array_to_save,name=os.path.join(folder_cache,
											file_manager.basefn + '_%s.pkl'%(combination[2:4]+combination[:2])))

			MPI.COMM_WORLD.Barrier()
		sys.exit()

	if file_manager.FileExist is False and rank==0:
		file_manager.save_data_on_disk(array_to_save)

if __name__ == "__main__":
		main()
