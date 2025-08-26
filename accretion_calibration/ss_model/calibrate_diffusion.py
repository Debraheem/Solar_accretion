# coding: utf-8
import numpy as np
import pandas as pd
import scipy as sp
import subprocess
import os
from scipy import optimize
import pickle
import sys
from functools import lru_cache, wraps

def np_cache(function):
    @lru_cache()
    def cached_wrapper(hashable_array):
        array = np.array(hashable_array)
        return function(array)

    @wraps(function)
    def wrapper(array):
        return cached_wrapper(tuple(array))

    # copy lru_cache attributes over too
    wrapper.cache_info = cached_wrapper.cache_info
    wrapper.cache_clear = cached_wrapper.cache_clear

    return wrapper

np.set_printoptions(precision=10)
np.random.seed(42) # for reproducibility 

Z_X_solar  =  0.02293 #0.02307 

# Directories and save files 
save_dir = "LOGS"
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

age = 4.568
if len(sys.argv) > 2:
    age = sys.argv[2]

X_names = ['Y', 'Z', 'a']
X = [0.269, 0.0187, 2.22]
X_var   = [0.005,        0.005,           0.01]
bounds  = [(0.25, 0.29), (0.012, 0.022), (1.5, 2.5)]
"""
X_names = [ 'Y',   'a']
X       = [0.26, 1.665]
X_var   = [0.005, 0.01]
bounds  = [(0.25, 0.29), (1.5, 2.5)]
print(X_names)
"""

P = lambda x: x#(x-X)/X_var
R = lambda x: x#x*X_var+X
lower = P(np.array([bound[0] for bound in bounds]))
upper = P(np.array([bound[1] for bound in bounds]))

def get_flags(names, args):
    return ' -'.join([''] + list(map(lambda x,y: x+' '+str(y), names, args)))

@np_cache
def calibrate(theta, fast=True, single=False):
    _theta = R(np.copy(theta))
    print("parameters:", _theta)
    
    Y, Z, alpha = _theta
    #Z = (1-Y) * Z_X_solar / (1+Z_X_solar) # [Fe/H] = 0 = log10(Z/X/0.02293); Z=1-Y-X
    
    if np.any(theta < lower) or np.any(theta > upper):
        print('out of bounds')
        if single: 
            return 1
        return np.ones(len(theta))
    
    bash_cmd = get_flags(X_names + ['Z'], list(_theta) + [Z])
    print(bash_cmd)
    
    full = "./mesa.sh"+bash_cmd
    if fast:
        full = full + ' -f'
    print(full)
    with open('temp.txt', 'w') as output:
        process = subprocess.Popen(full.split(), 
            shell=False, stdout=output)
        process.wait(timeout=60000)
    
    # process the results 
    hist_file = os.path.join(save_dir, 'history.data')
    pro_file  = os.path.join(save_dir, 'solar.data')
    if not os.path.isfile(hist_file) or not os.path.isfile(pro_file): 
        print('No output')
        if single: 
            return 1
        return np.ones(len(theta))
    DF = pd.read_table(hist_file, sep='\s+', skiprows=5)
    prof = pd.read_table(pro_file, sep='\s+', skiprows=5)
    
    logR = DF['log_R'].values[-1]
    logL = DF['log_L'].values[-1]
    Fe_H = np.log10(prof.z.values[0] / prof.x.values[0] / Z_X_solar)
    print('log R:', logR)
    print('log L:', logL)
    print('[Fe/H]:', Fe_H)
    
    if np.abs(logR) < 1e-7 and np.abs(logL) < 1e-7 and np.abs(Fe_H) < 1e-4:
        print('good enough!')
        sys.exit()
    
    if single:
        return np.log10(np.abs(logR)) + \
               np.log10(np.abs(logL)) + \
               np.log10(np.abs(Fe_H)) 
    return np.array([logR, logL, Fe_H])

eps = np.finfo(float).eps
result = optimize.root(calibrate, x0=P(np.array(X)))

print(result)
