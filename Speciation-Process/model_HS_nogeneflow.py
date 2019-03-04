#!/usr/bin/python
# -*- coding: utf-8 -*- 

import dadi
import matplotlib.pyplot as plt
import numpy as np
import itertools as itt
import pandas as pd
import sys
import os
import time

#################################
### Creating input data files ###

dd = dadi.Misc.make_data_dict('opt1_output_2.txt')

fs_arcgarmac = dadi.Spectrum.from_data_dict(dd,pop_ids=['arc','gar','mac'],projections=[44,65,50],polarized=True)
fs_arcgarmac = fs_arcgarmac.project([10,10,10])
fs_arcgarmac.to_file('Coenonympha_arcgarmac_10.fs')

fs_arcgardar = dadi.Spectrum.from_data_dict(dd,pop_ids=['arc','gar','dar'],projections=[44,65,20],polarized=True)
fs_arcgardar = fs_arcgardar.project([10,10,10])
fs_arcgardar.to_file('Coenonympha_arcgardar_10.fs')


#################################################################################
## scenario hybrid speciation with gene flow between hybrid species and parentals

def admix(params, ns, pts):
    ## parse params
    N1,N2,N3,T2,T1,f = params
    ## create a search grid
    xx = dadi.Numerics.default_grid(pts)
    ## make ancestral pop that splits into two
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
    ## allow drift to occur along each of these branches
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=N1, nu2=N2, m12=0., m21=0.)
    ## create pop 3 from a mixture of 1 and 2
    phi = dadi.PhiManip.phi_2D_to_3D_admix(phi, f, xx, xx, xx)
    ## allow drift and migration to occur along these branches
    ## cuba population shrinks in size after divergence
    phi = dadi.Integration.three_pops(phi, xx, T1,
                                      nu1=N1, nu2=N2, nu3=N3,
                                      m12=0.0, m13=0.0,
                                      m21=0.0, m23=0.0,
                                      m31=0.0, m32=0.0)
    ## simulate the fs
    fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx,xx))
    return fs

###############################################################################

fs_list = ['Coenonympha_arcgarmac_10.fs','Coenonympha_arcgardar_10.fs']
pop_ids = [['arc','gar','mac'],['arc','gar','dar']]

i = 0

while i < 2:
	print fs_list[i]

	data = dadi.Spectrum.from_file(fs_list[i])	

	## sample sizes
	ns = data.sample_sizes

	## Popultation IDs
	## pop_ids = Pop_ids[i]

	## points used for exrapolation
	pts_l = [40,50,60]

	## starting values for params
	N1 = N2 = N3 = 1.0
	T2 = 2.0
	T1 = 2.0
	m13 = m31 = m23 = m32 = 0.1
	f  = 0.4999

	## create starting parameter sets
	params_admix = np.array([N1,N2,N3,T2,T1,f])

	## search limits
	upper_admix = [100, 100, 100,  10,  2, 0.99999]
	lower_admix = [1e-2, 1e-2, 1e-2, 1e-3, 1e-3, 0.00001]

	###############################################################################

	func = admix
	maxiters = 2

	## Make the extrapolating version of our demographic model function
	func_ex = dadi.Numerics.make_extrap_log_func(func)	

	## Do the optimization
	p0 = dadi.Misc.perturb_params(params_admix, fold=2.,
                              upper_bound=upper_admix,
                              lower_bound=lower_admix),

	popt = dadi.Inference.optimize_log_lbfgsb(p0, data, func_ex, pts_l, lower_admix, upper_admix, 20, maxiters)

	## Calculate the best-fit model AFS.
	mod = func_ex(popt, ns, pts_l)
	ll_opt = dadi.Inference.ll_multinom(mod, data)

	## Add from Rougeux et al. 2017
	theta = dadi.Inference.optimal_sfs_scaling(mod, data)
	AIC = 2*len(params_admix)-2*ll_opt

	## Print results
	print 'Optimized parameters', repr(popt)
	print 'Optimized log-likelihood:', ll_opt
	print 'theta:', theta
	
	## Write results
	outputname = fs_list[i] + "_" + repr(time.localtime()[0]) + "_" + repr(time.localtime()[1]) + "_" + repr(time.localtime()[2]) + "_" + repr(time.localtime()[3])	
	os.mkdir("./RESULTS_HS_nogeneflow/" + outputname)
	output_file = open(("./RESULTS_HS_nogeneflow/" + outputname + "/" + outputname + ".txt"), "w")
	line = ("\n" + "HS_nogeneflow" + "\n" + "\n" "Optimization : " + repr("BFGS") + "\n"  "Optimized parameters: " + repr(popt) + "\n" + "Optimized log-likelihood: " + repr(ll_opt) + "\n" + "theta: " + repr(theta) + "\n" + "AIC: " + repr(AIC) + "\n")
	output_file.write(line)

	## Plot a comparison of the resulting fs with the data.
	#import pylab
	#fig = pylab.figure()
	#dadi.Plotting.plot_3d_comp_multinom(mod, data, vmin=0.01, resid_range=3, pop_ids = pop_ids[i])
	#pylab.savefig(("./RESULTS_HS_nogeneflow/" + outputname + "/" + outputname + "_" + "HS_nogeneflow" + ".png"), dpi=300)
	#pylab.close(fig)

	i = i + 1

	output_file.close()

