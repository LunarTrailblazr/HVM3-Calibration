# David R Thompson
import pylab as plt
import numpy as np
from glob import glob
from datetime import datetime, timedelta 
from collections import OrderedDict
from spectral.io import envi
import os, sys

def celcius2kelvin(c):
    return c+273.15

timestamp_dir = '../data/temperature_logs/'
files = glob(timestamp_dir+'2*.txt')
temperatures = []
files.sort()
for fi in files:
  with open(fi,'r') as fin:
    for line in fin.readlines()[1:]:
      tokens = line.split()
      time_s = float(tokens[0])
      spectrometer_T = celcius2kelvin(float(tokens[2]))
      bus_T = celcius2kelvin(float(tokens[3]))
      fpa_T = celcius2kelvin(float(tokens[1]))
      epoch = datetime(year=1904,month=1,day=1,hour=0,minute=0,second=0)
      offset = timedelta(seconds=int(time_s))
      time = epoch+offset
      print(fi,time,spectrometer_T)
      temperatures.append((time,spectrometer_T))
temperatures.sort()
    
tol = timedelta(seconds=20)

def utc2temp(temperatures, query, lo=None, hi=None, margin=5):
  if lo is None:
     lo = 0
     hi = int(len(temperatures)-1)
  test = int((hi-lo)/2+lo)
  current = temperatures[test][0]
  if abs(current-query) < tol: 
     return temperatures[test][1]
  elif test<margin:
     return None
  elif test >= len(temperatures)-margin:
     return None
  elif current > query:
     hi = test
  else: 
     lo = test
  return utc2temp(temperatures, query, lo=lo, hi=hi)

test = datetime(year=1900,month=1,day=1)
base = '/Users/drt/data/22HVM3/darks/'

temps = {}
with open(timestamp_dir+'timestamps_modified.txt') as fin:
  for line in fin.readlines():
     timestr, filename = line.split()    
     if not filename.endswith('.hdr'): 
       continue
     filename = base+filename.replace('.hdr','_avg.hdr')
     outname = filename.replace('.hdr','_emission.hdr')
     if not os.path.exists(filename):
       print('did not find',filename)
       continue
     print(filename)
     dt = test.strptime(timestr,'%Y%m%d%H%M%S') 
     dt = dt - timedelta(minutes=4)
     temperature = utc2temp(temperatures, query=dt)
      
     I = envi.open(filename).load()
     dk = np.median(np.squeeze(I[:,:9,0]),axis=1)
     bg = np.median(np.squeeze(I[:,11:19,0]),axis=1)
     outname = filename.replace('.hdr','_emission_%4.1f.hdr'%temperature)
     envi.save_image(outname,np.reshape(bg-dk,(1,1,len(bg))),force=True,ext='')

     dk = np.median(np.squeeze(I[148:153,0:6,0]))
     bg = np.median(np.squeeze(I[148:153,:,0]),axis=0)
     profilename = filename.replace('.hdr','_zone3_%4.1f.hdr'%temperature)
     envi.save_image(profilename,np.reshape(bg-dk,(1,1,len(bg))),force=True,ext='')

     dk = np.median(np.squeeze(I[200:210,0:6,0]))
     bg = np.median(np.squeeze(I[200:210,:,0]),axis=0)
     profilename = filename.replace('.hdr','_zone2_%4.1f.hdr'%temperature)
     envi.save_image(profilename,np.reshape(bg-dk,(1,1,len(bg))),force=True,ext='')

     dk = np.median(np.squeeze(I[275:295,0:6,0]))
     bg = np.median(np.squeeze(I[275:295,:,0]),axis=0)
     profilename = filename.replace('.hdr','_zone1_%4.1f.hdr'%temperature)
     envi.save_image(profilename,np.reshape(bg-dk,(1,1,len(bg))),force=True,ext='')
     #plt.plot(bg-dk)

plt.show()
   
