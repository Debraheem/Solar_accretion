#!/usr/bin/env python3
# inlist_calibrate_diffusion.py

# Step 1: Imports & setup
import numpy as np
import pandas as pd
import scipy as sp
import subprocess
import os
from scipy import optimize
import pickle
import sys
from functools import lru_cache, wraps

# Step 2: Caching decorator for speed
def np_cache(function):
    @lru_cache()
    def cached_wrapper(hashable_array):
        array = np.array(hashable_array)
        return function(array)
    @wraps(function)
    def wrapper(array):
        return cached_wrapper(tuple(array))
    wrapper.cache_info = cached_wrapper.cache_info
    wrapper.cache_clear = cached_wrapper.cache_clear
    return wrapper

# Step 3: Global constants & output dir
np.set_printoptions(precision=10)
np.random.seed(42)               # reproducibility
Z_X_solar = 0.0181# agss09 !0.02293 for gs98              # solar Z/X
save_dir  = "LOGS"
os.makedirs(save_dir, exist_ok=True)

# Step 4: Parse optional age argument
age = 4.572
if len(sys.argv) > 1:
    age = float(sys.argv[1])

# Step 5: Calibration variables, names, and bounds
X_names = ['Y', 'Z', 'a']        # helium, metals, αMLT
X       = [0.269, 0.0187, 2.22]
X_var   = [0.005,        0.005,           0.01]
bounds  = [(0.25, 0.29), (0.012, 0.022), (1.5, 2.5)]

lower = np.array([b[0] for b in bounds])
upper = np.array([b[1] for b in bounds])

# Step 6: Flag‐builder for mesa.sh
def get_flags(names, args):
    """
    Builds a string like:  -Y 0.27 -Z 0.0187 -a 2.22
    """
    pairs = ["{} {}".format(n, v) for n, v in zip(names, args)]
    return ' ' + ' '.join(['-' + p for p in pairs])

# Step 7: Core calibration function
@np_cache
def calibrate(theta, fast=True, single=False):
    # 7.1 Rescale back to physical parameters
    Y, Z, alpha = theta
    
    # 7.2 Check bounds
    if np.any(theta < lower) or np.any(theta > upper):
        print("→ Out of bounds:", theta)
        return 1.0 if single else np.ones(len(theta))
    
    # 7.3 Build & run mesa.sh
    flags = get_flags(X_names, [Y, Z, alpha])
    cmd   = "./mesa.sh" + flags
    if fast:
        cmd += " -f"
    print("→ Running:", cmd)
    with open('temp.txt', 'w') as out:
        proc = subprocess.Popen(cmd.split(), stdout=out)
        proc.wait(timeout=60000)
    
    # 7.4 Read results
    hist = os.path.join(save_dir, 'history.data')
    prof = os.path.join(save_dir, 'solar.data')
    if not (os.path.isfile(hist) and os.path.isfile(prof)):
        print("→ No output files found")
        return 1.0 if single else np.ones(len(theta))
    
    DF   = pd.read_table(hist, sep='\s+', skiprows=5)
    PR   = pd.read_table(prof, sep='\s+', skiprows=5)
    
    # 7.5 Compute residuals
    logR = DF['log_R'].iloc[-1]
    logL = DF['log_L'].iloc[-1]
    Fe_H = np.log10(PR.z.values[0] / PR.x.values[0] / Z_X_solar)
    print(f"  logR={logR:.3e}, logL={logL:.3e}, [Fe/H]={Fe_H:.3e}")
    
    # 7.6 Convergence check
    if abs(logR) < 1e-7 and abs(logL) < 1e-7 and abs(Fe_H) < 1e-4:
        print("→ Converged!")
        sys.exit(0)
    
    # 7.7 Return either scalar or vector of residuals
    if single:
        return np.log10(abs(logR)) + np.log10(abs(logL)) + np.log10(abs(Fe_H))
    else:
        return np.array([logR, logL, Fe_H])

# Step 8: Run the optimizer
if __name__ == "__main__":
    result = optimize.root(calibrate, x0=np.array(X))
    print("Optimization result:", result)
