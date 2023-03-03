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
description = {{HVM3 L1B Flat Field}}
samples = {columns}
lines = {lines}
bands = {rows}
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bil
byte order = 0
wavelength units = Nanometers
"""


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

    # Load flat field
    flat_field = sp.fromfile(config.flat_field_file, dtype = sp.float32)
    flat_field = flat_field.reshape((config.frame_size))
    flat_field[np.logical_not(np.isfinite(flat_field))] = 0

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
    
    with open(args.input_file,'rb') as fin:
        with open(args.output_file,'wb') as fout:

           raw = sp.fromfile(fin, count=nframe, dtype=dtype)
               
           while len(raw)>0:

               # Read a frame of data 
               raw = np.array(raw, dtype=sp.float32)
               frame = raw.reshape((rows,columns))

               # Apply flat field
               frame = frame * flat_field
    
               # Write to file
               np.asarray(frame, dtype=sp.float32).tofile(fout)
            
               # Read next chunk
               raw = sp.fromfile(fin, count=nframe, dtype=dtype)

    # Place all calibration parameters in header metadata
    params = {'lines': lines}

    # Write the header
    params.update(**locals())
    with open(args.output_file+'.hdr','w') as fout:
        fout.write(header_template.format(**params))


if __name__ == '__main__':

    main()
