##### Title: Script to run DaDi inference
##### Input: MST 
##### Output: 

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 14:27:13 2021

@author: ksliao
"""

import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
import pandas as pd
import numpy
from multiprocessing import Process
import os
import subprocess
import sys

#Import genomewide annotated MST vcf 
all_chr = pd.read_csv("/net/wonderland/home/ksliao/zoellner_research/scripts/output/genomewide_fullanno.vcf", delimiter="\t", header=None)
all_chr.columns = ["CHR", "REF","ALT", "AA", "AC_correct", "ANNO", "MT", "KMER", "MST"]

MST = sys.argv[1]
print(MST)

#Define demographic history model
#At time TF + TB in past, equilibirium pop goes through bottleneck of depth nuB recovering to nuF
def bottleneck(params, ns, pts):
    nuB, nuF, TB, TF = params
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TB, nuB)
    phi = Integration.one_pop(phi, xx, TF, nuF)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

#Subset vcf to just the particular MST
all_chr_MST = all_chr[all_chr['MST'] == MST]


#Create site frequency spectrum using frequency counts
item_list = range(7113);

all_chr_MST_afs = pd.DataFrame(all_chr_MST['AC_correct'].value_counts().reindex(item_list, fill_value = 0))
all_chr_MST_afs.columns = ['ac']

all_chr_MST_afs.sort_index(inplace=True)

#Create an array from dataframe ac column values
all_chr_MST_ac_array = all_chr_MST_afs['ac'].values

#Create data objects and run dadi
data = dadi.Spectrum(all_chr_MST_ac_array)
ns = data.sample_sizes

print(ns)

theta_w = data.Watterson_theta()
theta_pi = data.pi()

pts_l = [7200, 7210, 7220]

#Provide upper and lower bounds of optimization
upper_bound = [1, 500, 10, 10]
lower_bound = [1e-3, 1e-3, 0, 0]

#Starting values
p0 = [0.5, 20, 1, .2]

func_ex = dadi.Numerics.make_extrap_log_func(bottleneck)

p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)


print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, verbose=len(p0), maxiter=30)

    # The verbose argument controls how often progress of the optimizer should be
    # printed. It's useful to keep track of optimization process.
print('Finshed optimization **************************************************')

print("here is popt")
print(popt)

model = func_ex(popt, ns, pts_l)

ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))

theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))


#Import jeds mutation rates
#mut_rates = pd.read_csv("/net/wonderland/home/ksliao/zoellner_research/demographic_inf/data/Jeds mut rates simplified.csv", delimiter=",", header=0)
#mut_rates.drop(['Unnamed: 0'], axis=1)
#MST_mutRate = mut_rates[mut_rates['MST'] == MST].iloc[0]['ERV_rel_rate']

#Create dataframe to store optimized theta, MST counts, jeds mutation rates, and normalized theta
theta_df = pd.DataFrame([[theta]], columns=['theta'])
theta_df['max_ll'] = ll_model
theta_df['nuB'] = popt[0]
theta_df['nuF'] = popt[1]
theta_df['TB'] = popt[2]
theta_df['TF'] = popt[3]
theta_df['MST'] = MST
                               
 
print(theta_df.head())

theta_df.to_csv(r"/net/wonderland/home/ksliao/zoellner_research/demographic_inf/output/bottleneck_7_1_21/" + MST + ".txt", sep=' ', index=None, header=None)
