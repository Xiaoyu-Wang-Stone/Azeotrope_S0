import os, sys
import numpy as np
import scipy
from scipy import signal
import matplotlib.pyplot as plt
from math import isinf
from math import pi
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.optimize import curve_fit

# desrciption of input data 
name_dict={'AA':'Met-Met', 'AB':'Met-Water','BB':'Water-Water'}

mol_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


def fit_OrnsteinZernike(x, s0, xi):
    # S(k) = S(0)/(1+xi*k^2)
    return s0/(1.+xi*x**2)

for ii in ['II','IW','WW']:
    Sk[ii] = {}
    for mol in mol_list:
        Sk[ii][mol] = np.loadtxt(f"{mol}_met_water/Sk-{ii}-real-avg.list", \
            skiprows=1)[:,:]

#
#left, bottom, width, height = 0, 0, 1, 0.8
#fig = plt.figure(figsize=(6,3),dpi=500)
#axs0 = fig.add_axes([left, bottom, width, height])# main axes
#axs0.margins(-0.0, x=None, y=None, tight=True)


# calculate S0
S0_dict = {}
k_sq_cut = 0.01
for mol in mol_list:
    S0_dict[mol] = {}
    for ii in ['II','IW','WW']:
        Sk_cut = np.asarray([sk for sk in Sk[ii][mol][1:] if sk[0]**2.+sk[1]**2.+sk[2]**2. < k_sq_cut 
                            and not np.isnan(sk[5])])
        ksqr_now = [ sk[0]**2.+sk[1]**2.+sk[2]**2. for sk in Sk_cut]
        sk_now = Sk_cut[:,3]
        sk_error_now = Sk_cut[:,5]
   
        popt, pcov = curve_fit(fit_OrnsteinZernike, ksqr_now, sk_now)
        perr = np.sqrt(np.diag(pcov)) # error 
            
        if isinf(perr[0]): 
            popt = [float('nan'), float('nan')]  # why
            
        S0_dict[mol][ii]={'S0': popt[0], 'S0_error': perr[0], 'value': popt, 'error': perr}
        


    print(mol, S0_dict[mol]["II"]["S0"],S0_dict[mol]["II"]["S0_error"], \
        S0_dict[mol]["IW"]["S0"],S0_dict[mol]["IW"]["S0_error"], \
        S0_dict[mol]["WW"]["S0"],S0_dict[mol]["WW"]["S0_error"])


