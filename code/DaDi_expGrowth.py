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

#Define exponential growth model
#Nu is relative size of ancestral population to the reference population. Often reference population is ancestry so nu defaults to 1
#T is the integration time
def growth(params , ns, pts):
    nu,T = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    nu_func = lambda t: numpy.exp(numpy.log(nu) * t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)
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
upper_bound = [750, 10]
lower_bound = [1e-3, 1e-3]

p0 = [150, 3]

func_ex = dadi.Numerics.make_extrap_log_func(growth)

p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, verbose=len(p0), maxiter=15)

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

#Create dataframe to store optimized theta, MST counts, jeds mutation rates, and normalized theta
theta_df = pd.DataFrame([[theta]], columns=['theta'])
theta_df['max_ll'] = ll_model
theta_df['nu'] = popt[0]
theta_df['T'] = popt[1]
theta_df['MST'] = MST
theta_df['theta_w'] = theta_w
theta_df['theta_pi'] = theta_pi

print(theta_df.head())

theta_df.to_csv(r"/net/wonderland/home/ksliao/zoellner_research/demographic_inf/output/exp_growth_10_6_19/" + MST + ".txt", sep=' ', index=None, header=None)

