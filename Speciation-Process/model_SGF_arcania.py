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

def IM_split1(params, ns, pts):
    ## parse params
    N1,N2,N3,T2,T1,m13,m31,m23,m32 = params
    ## create a search grid
    xx = dadi.Numerics.default_grid(pts)
    ## create ancestral pop
    phi = dadi.PhiManip.phi_1D(xx)
    ## split ancestral pop into two species
    phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
    ## allow drift to occur along each of these branches
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=N1, nu2=N2, m12=0., m21=0.)
    ## split pop1 into pops 1 and 3
    phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)
    ## allow drift and migration to occur along these branches
    phi = dadi.Integration.three_pops(phi, xx, T1,
                                      nu1=N1, nu2=N2, nu3=N3,
                                      m12=0.0, m13=m13,
                                      m21=0.0, m23=m23,
                                      m31=m31, m32=m32)
    ## simulate the fs
    fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx,xx))
    return fs

###############################################################################

fs_list = ['Coenonympha_arcgarmac_10.fs', 'Coenonympha_arcgardar_10.fs']
pop_ids = [['arc','gar','mac'],['arc','gar','dar']]

i = 1

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

	## create starting parameter sets
	params_SGF = np.array([N1,N2,N3,T2,T1,m13,m31,m23,m32])

	## search limits
	upper_SGF = [100, 100, 100,  10,  2, 40, 40, 40, 40]
	lower_SGF = [1e-2, 1e-2, 1e-2, 1e-3, 1e-3, 0, 0, 0, 0]

	###############################################################################

	func = IM_split1
	maxiters = 2

	## Make the extrapolating version of our demographic model function
	func_ex = dadi.Numerics.make_extrap_log_func(func)	

	## Do the optimization
	p0 = dadi.Misc.perturb_params(params_SGF, fold=2.,
                              upper_bound=upper_SGF,
                              lower_bound=lower_SGF),

	popt = dadi.Inference.optimize_log_lbfgsb(p0, data, func_ex, pts_l, lower_SGF, upper_SGF, 20, maxiters)

	## Calculate the best-fit model AFS.
	mod = func_ex(popt, ns, pts_l)
	ll_opt = dadi.Inference.ll_multinom(mod, data)

	## Add from Rougeux et al. 2017
	theta = dadi.Inference.optimal_sfs_scaling(mod, data)
	AIC = 2*len(params_SGF)-2*ll_opt

	## Print results
	print 'Optimized parameters', repr(popt)
	print 'Optimized log-likelihood:', ll_opt
	print 'theta:', theta
	
	## Write results
	outputname = fs_list[i] + "_" + repr(time.localtime()[0]) + "_" + repr(time.localtime()[1]) + "_" + repr(time.localtime()[2]) + "_" + repr(time.localtime()[3]) + repr(time.localtime()[4])	
	os.mkdir("./RESULTS_SGF_arcania/" + outputname)
	output_file = open(("./RESULTS_SGF_arcania/" + outputname + "/" + outputname + ".txt"), "w")
	line = ("\n" + "SGF_arcania" + "\n" + "\n" "Optimization : " + repr("BFGS") + "\n"  "Optimized parameters: " + repr(popt) + "\n" + "Optimized log-likelihood: " + repr(ll_opt) + "\n" + "theta: " + repr(theta) + "\n" + "AIC: " + repr(AIC) + "\n")
	output_file.write(line)

	## Plot a comparison of the resulting fs with the data.
	#import pylab
	#fig = pylab.figure()
	#dadi.Plotting.plot_3d_comp_multinom(mod, data, vmin=0.01, resid_range=3, pop_ids = pop_ids[i])
	#pylab.savefig(("./RESULTS_HS_geneflow_arcania/" + outputname + "/" + outputname + "_" + "HS_geneflow_arcania" + ".png"), dpi=300)
	#pylab.close(fig)

	i = i + 1

	output_file.close()

