#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# HVM3 Radiometric Calibration code
# Author: David R Thompson, david.r.thompson@jpl.nasa.gov

import scipy.linalg
import os, sys
import numpy as np
from spectral.io import envi
import json
import logging
import argparse


header_string = """ENVI
description = {HVM3 Dark Frame}
samples = %i
lines = %i
bands = 2
header offset = 0
file type = ENVI Standard
data type = 4
interleave = bsq
byte order = 0"""


def find_header(infile):
  if os.path.exists(infile+'.hdr'):
    return infile+'.hdr'
  elif os.path.exists('.'.join(infile.split('.')[:-1])+'.hdr'):
    return '.'.join(infile.split('.')[:-1])+'.hdr'
  else:
    raise FileNotFoundError('Did not find header file')


def subtract_dark(frame, dark):
    if frame.shape[0] != dark.shape[0] or \
       frame.shape[1] != dark.shape[1]: 
         logging.error('Mismatched dark and frame sizes')
    return frame-dark


def emission_scale(x, profile, spectrum, startchan):
    err = abs(x[0]+x[1]*profile - spectrum)
    return sum(err[startchan:])


def main():

    description = "Subtract dark frame from a data cube"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input')
    parser.add_argument('--max_column',type=int,default=9999)
    parser.add_argument('--min_column',type=int,default=-9999)
    parser.add_argument('--masked_columns_right',type=int,default=5)
    parser.add_argument('--masked_columns_left',type=int,default=20)
    parser.add_argument('--exclude_emission',type=str,default='')
    parser.add_argument('dark')
    parser.add_argument('output')
    args = parser.parse_args()

    infile = envi.open(find_header(args.input))
    darkfile = envi.open(find_header(args.dark))
    dark = np.squeeze(darkfile.load()[:,:,0])

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

    if len(args.exclude_emission)>1:
        emission_profile = envi.open(args.exclude_emission+hdr).load()
        emission_profile = np.squeeze(emisison_profile)
        for i in range(dark.shape[1]):
            xbest = minimize(emission_scale, [0,1], args=(profile,dark[i,:],150))
            dark[i,:] = dark[i,:] - xbest.x[1] * emission_profile
    

    rows = int(infile.metadata['bands'])
    columns = int(infile.metadata['samples'])
    lines = int(infile.metadata['lines'])
    nframe = rows * columns

    metadata = infile.metadata.copy()
    metadata['data type'] = 4
    envi.write_envi_header(args.output+'.hdr', metadata)

    with open(args.input,'rb') as fin:
      with open(args.output,'wb') as fout:

        for line in range(lines):

            # Read a frame of data
            if line%10==0:
                logging.info('Line '+str(line))
            frame = np.fromfile(fin, count=nframe, dtype=dtype)
            frame = np.array(frame.reshape((rows, columns)),dtype=np.float32)
            frame = subtract_dark(frame, dark)
            frame[:,args.max_column:(-args.masked_columns_right)] = 0
            frame[:,args.masked_columns_left:args.min_column] = 0
            np.array(frame, dtype=np.float32).tofile(fout)

    print('done') 

if __name__ == '__main__':

    main()
