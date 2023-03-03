#! /usr/bin/env python
#
#  Copyright 2022 California Institute of Technology
#
# HVM3 Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys
import numpy as np
from spectral.io import envi
import pylab as plt
import json
import logging
import argparse
from numba import jit
from math import pow
from scipy.signal import medfilt
from scipy.interpolate import interp1d
from hvm3_config import HVM3Config
from sklearn.linear_model import RANSACRegressor
import subprocess
import ray


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


@ray.remote
def fix_emission(frame, factors):
  est_emission = np.zeros_like(frame)
  model = RANSACRegressor()
  for i in range(frame.shape[1]):
    y = frame[170:,i:(i+1)]
    x = factors[170:,np.newaxis]
    if all(y==y[0]):
      continue
    model.fit(x,y)
    est_emission[:,i] = model.predict(factors[:,np.newaxis])[:,0]
  est_emission_smooth = medfilt(est_emission, 5)
  new = frame - est_emission_smooth
  return new
  


def main():

    description = "Remove spectrometer self-emission"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('config')
    parser.add_argument('dark')
    parser.add_argument('--mode',default='directsolar_wide')
    parser.add_argument('--max_jobs',default=40)
    parser.add_argument('output')
    args = parser.parse_args()
    ray.init()

    config = HVM3Config(args.config, mode=args.mode)
    dark = np.squeeze(envi.open(args.dark+'.hdr').load()[:,:,0])
    factors = np.median(dark[:,30:-30],axis=1)
    
    factors = factors - np.median(dark[:,config.masked_cols], axis=1)

    infile = envi.open(find_header(args.input))

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
    lines = int(infile.metadata['lines'])
    nframe = rows * columns
    metadata = infile.metadata.copy()
    metadata['data type'] = 4

    envi.write_envi_header(args.output+'.hdr',metadata)

    with open(args.input,'rb') as fin:
      with open(args.output,'wb') as fout:

        frame = np.fromfile(fin, count=nframe, dtype=dtype)
        jobs = []
        lines_analyzed =  0
               
        while len(frame)>0:

            # Read a frame of data 
            frame = np.array(frame.reshape((rows, columns)),dtype=dtype)

            if lines_analyzed%10==0:
                logging.info('Calibrating line '+str(lines_analyzed))
                print(lines_analyzed)

            #jobs.append(fix_emission(frame, factors))
            jobs.append(fix_emission.remote(frame, factors))
            lines_analyzed = lines_analyzed + 1

            if len(jobs) == args.max_jobs:
                
                # Write to file
                result = ray.get(jobs)
                for frame in result:
                    np.asarray(frame, dtype=np.float32).tofile(fout)
                jobs = []
        
            # Read next chunk
            frame = np.fromfile(fin, count=nframe, dtype=dtype)

        # Do any final jobs
        result = ray.get(jobs)
        for frame in result:
            np.asarray(frame, dtype=np.float32).tofile(fout)

    print('done') 

if __name__ == '__main__':

    main()

  
