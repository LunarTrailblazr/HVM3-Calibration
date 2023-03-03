#!/home/drt/src/anaconda3/bin/python
# C
# David R Thompson
# david.r.thompson@jpl.nasa.gov
# 
import sys, os
import pylab as plt
from scipy.ndimage import median_filter
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from spectral.io import envi
from glob import glob
import numpy as np
import os.path

base = os.path.dirname(os.path.realpath(__file__)) 
config = base + 'l1/config/hvm3.json'


with open('input_files.txt','r') as fin:
  for line in fin.readlines():

     dark, raw, mode = line.split()

     # Average the dark data into a single dark frame
     dark_average = dark.replace('.raw','') + '_avg'
     if not os.path.exists(dark_average):
         cmd = 'python %s/l1/utils/hvm3_make_dark.py %s %s' % (base, dark, dark_average)
         print(cmd)
         os.system(cmd)

     # Dark subtraction
     darksub = raw.replace('.raw','') + '_darksub'
     if not os.path.exists(darksub):
         cmd = 'python %s/l1/utils/hvm3_subtract_dark.py %s %s %s' % (base, raw, dark_average, darksub)
         print(cmd)
         os.system(cmd)

     # Pedestal shift
     pedestal = raw.replace('.raw','') + '_pedestal'
     if not os.path.exists(pedestal):
         cmd = 'python %s/l1/utils/hvm3_fix_pedestal.py %s %s %s' % (base, darksub, config, pedestal)
         print(cmd)
         os.system(cmd)

     # Flat field 
     flatted = raw.replace('.raw','') + '_flat'
     if not os.path.exists(flatted):
         cmd = 'python %s/l1/utils/hvm3_apply_flat.py %s %s %s' % (base, pedestal, config, flatted)
         print(cmd)
         os.system(cmd)

     # Find and fix bad elements
     badfix = raw.replace('.raw','') + '_badfix'
     badmap = raw.replace('.raw','') + '_badmap'
     if not os.path.exists(badfix):
         cmd = 'srun -n 1 -N 1 -c 40 --mem 180GB python %s/l1/utils/hvm3_fix_bad.py %s --npca 3 --output_bad %s %s' % (base, flatted, badmap, badfix)
         #cmd = 'python %s/l1/utils/hvm3_fix_bad.py %s --npca 3 --num_cpus 1 --output_bad %s %s' % (base, flatted, badmap, badfix)
         print(cmd)
         os.system(cmd)

     # Apply RCCs
     radiance = raw.replace('.raw','') + '_radiance'
     if not os.path.exists(radiance):
         cmd = 'python %s/l1/utils/hvm3_apply_rccs.py %s --mode %s %s %s' % (base, badfix, mode, config, radiance)
         print(cmd)
         os.system(cmd)
