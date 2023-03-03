#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# HVM3 Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys, os.path
import scipy as sp
import numpy as np
from spectral.io import envi
from datetime import datetime, timezone
from scipy import linalg, polyfit, polyval
import json
import logging
import argparse
import multiprocessing
import ray
import pylab as plt

# Import some HVM3-specific functions
my_directory, my_executable = os.path.split(os.path.abspath(__file__))
sys.path.append(my_directory + '/utils/')

from hvm3_config import HVM3Config


header_template = """ENVI
description = {{HVM3 L1B calibrated spectral radiance (units: uW nm-1 cm-2 sr-1)}}
samples = {columns}
lines = {lines}
bands = {rows}
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bil
byte order = 0
wavelength units = Nanometers
wavelength = {{{wavelength_string}}}
fwhm = {{{fwhm_string}}}
band names = {{{band_names_string}}}
"""


def bin_frame(frame,factor):
  out = np.zeros((frame.shape[0], int(frame.shape[1]/factor)))
  for i in range(out.shape[1]):
     out[:,i] = np.mean(frame[:,(i*factor):((i+1)*factor)], axis=1)
  return out


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')

     

def main():

    description = "Spectroradiometric Calibration"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--mode', default = 'directsolar_wide')
    parser.add_argument('input_file', default='')
    parser.add_argument('config_file', default='')
    parser.add_argument('output_file', default='')
    args = parser.parse_args()
    config = HVM3Config(args.config_file, mode=args.mode)

    ray.init()

    # Load radiometry
    _, rccs, uncert = sp.loadtxt(config.radiometric_coefficient_file).T

    # Load wavelengths
    _, wl_full, fwhm_full = sp.loadtxt(config.spectral_calibration_file).T * 1000

    infile = envi.open(find_header(args.input_file))

    if int(infile.metadata['data type']) == 2:
        dtype = np.int16
    elif int(infile.metadata['data type']) == 12:
        dtype = np.uint16
    elif int(infile.metadata['data type']) == 4:
        dtype = np.float32
    else:
        raise ValueError('Unsupported data type')
    if infile.metadata['interleave'] != 'bil':
        raise ValueError('Unsupported interleave')

    rows = int(infile.metadata['bands'])
    columns = int(infile.metadata['samples'])
    if rows != config.frame_size[0] or columns != config.frame_size[1]:
        raise ValueError('Image size does not match configuration')
    lines = int(infile.metadata['lines'])
    nframe = rows * columns
    lines_analyzed = 0
    noises = []

    if hasattr(config, 'distributed_columns'):
         distributed_columns = np.arange(config.distributed_columns[0],
                                         config.distributed_columns[1], 
                                         dtype=int)
    else:
         distributed_columns = np.arange(lines, dtype=int)
    if hasattr(config, 'distributed_rows'):
         distributed_rows = np.arange(config.distributed_rows[0],
                                         config.distributed_rows[1], 
                                         dtype=int)
    else:
         distributed_rows = np.arange(rows, dtype=int)
    
    with open(args.input_file,'rb') as fin:
       with open(args.output_file,'wb') as fout:

            raw = sp.fromfile(fin, count=nframe, dtype=dtype)
            jobs = []
               
            while len(raw)>0:

                # Read a frame of data 
                raw = np.array(raw, dtype=sp.float32)
                frame = raw.reshape((rows,columns))

                # Absolute radiometry
                frame = (frame.T * rccs).T
          
                # Clip the FPA 
                frame = frame[:,distributed_columns]
                frame = frame[distributed_rows,:]

                # Flip short-to-long
                frame = np.flip(frame, axis=0)

                # Crosstrack binning
                frame = bin_frame(frame, int(config.crosstrack_binning))
    
                # Write to file
                np.asarray(frame, dtype=sp.float32).tofile(fout)
            
                # Read next chunk
                raw = sp.fromfile(fin, count=nframe, dtype=dtype)

    # Form output metadata strings
    wl = np.flip(wl_full[distributed_rows])
    fwhm = np.flip(fwhm_full[distributed_rows])

    # Prepare wavelength metadata
    band_names_string = ','.join(['channel_'+str(i) \
       for i in range(len(wl))])
    fwhm_string =  ','.join([str(w) for w in fwhm])
    wavelength_string = ','.join([str(w) for w in wl])

    # Adjust rows and columns to account for the clipping
    rows = len(distributed_rows)
    columns = int(len(distributed_columns)/config.crosstrack_binning)
    
    # Write the header
    with open(args.output_file+'.hdr','w') as fout:
        fout.write(header_template.format(**locals()))


if __name__ == '__main__':

    main()
