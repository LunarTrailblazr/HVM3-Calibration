#!/home/drt/src/anaconda3/bin/python
# David R Thompson
import sys, os
import pylab as plt
from scipy.ndimage import median_filter
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from spectral.io import envi
from glob import glob
import numpy as np
import os.path

header_template = '''ENVI
description = {LIVEVIEW raw export file, 1 frame mean per grab}
samples = 640
lines   = %i 
bands   = 1
header offset = 0
file type = ENVI Standard
data type = 12
interleave = bil
sensor type = Unknown
byte order = 0
wavelength units = Unknown
'''

profile = np.loadtxt('profile.txt')
profilex = np.arange(len(profile))
config = '~/src/hvm3/l1/config/hvm3.json'
config_right = '~/src/hvm3/l1/config/hvm3_rightside.json'


def fit_shift(I, shift0=0, subsamp=4):
  subsamp = int(I.shape[0]/1000)
  print(subsamp)
  filt = [1,3]
  N = median_filter(I[0:I.shape[0]:subsamp,:],filt)
  def shift_err(x):
    I2 = np.zeros_like(N)
    I2[0,:] = N[0,:]
    for i in range(1,N.shape[0]):
      I2[i,:]=interp1d(np.arange(640)-x,I2[i-1,:],fill_value=0,bounds_error=False)(np.arange(640))
    err=np.nanmedian(abs(median_filter(I2,filt)-N))
    print(x,err)
    return err
  x = minimize(shift_err,shift0*subsamp,method='Nelder-Mead',tol=0.01)
  return x.x/subsamp


def forward(x,i):
   ctr = x[0]+x[1]*i
   my_profile = interp1d(profilex+ctr, profile, bounds_error=False, fill_value=0)(np.arange(640))
   my_profile = my_profile * x[2] + x[3]
   return my_profile, ctr


def fiterr(x,I):
  errs = []
  for i in range(I.shape[0]):
     my_profile, ctr = forward(x,i)
     use = my_profile[i]>1e-6
     errs.append(sum(abs(I[i,use]-my_profile[use])))
  return np.median(errs)


def extract(dark, dark_extract):
    lines = 0
    with open(dark_extract,'wb') as fout:
        with open(dark,'rb') as fin:
            x = np.fromfile(fin,count=640*321,dtype=np.uint16)
            while len(x)>0:
               x = x.reshape(321,640)
               lines=lines+1
               fout.write(x[19,:])
               x = np.fromfile(fin,count=640*321,dtype=np.uint16)
    with open(dark_extract+'.hdr','w') as hout:
       hout.write(header_template%lines)


def make_mask(dark_extract, dark_mask, threshold):
    I = np.squeeze(envi.open(dark_extract+'.hdr').load())
    shift = fit_shift(I, shift0*shiftdir)
    max_row = np.argmax(np.sum(I>threshold,axis=1))
    segment = np.where(I[max_row,:]>threshold)[0]
    ctr = (max(segment)-min(segment))/2+min(segment)
    offset = ctr+(max_row*shift)
    print('max_row',max_row,'ctr',ctr,'offset',offset,'shift',shift)
    mask = np.zeros_like(I)
    margin = 25
    for i in range(I.shape[0]):
       ctr = int(offset-i*shift)
       lo = max([int(ctr-margin),0])
       hi = min([int(ctr+margin),640])
       if hi>lo:
          mask[i,lo:hi]=1
    envi.save_image(dark_mask+'.hdr',mask.reshape((mask.shape[0],mask.shape[1],1)),ext='',force=True)


for temp in [325,250,275,282,285,290,300,350,375,400,425]:
   for shiftdir,direction in [(-1,'LR'),(1,'RL')]:
   #for shiftdir,direction in [(1,'RL')]:
       #for integration,threshold,shift0 in [('227ms',5000,0.057)]:
       for integration,threshold,shift0 in [('227ms',5000,0.057),('44ms',1500,0.0115)]:

          pattern = '*UTC*_BB%i_%s_Dark_%s/*.raw'%(temp,integration,direction)
          print(pattern)
          hits = glob(pattern)
          if len(hits)<1:
              continue
          dark = hits[0]
          light = glob('*UTC*_BB%i_%s_%s/*.raw'%(temp,integration,direction))[0]

          dark_extract = dark.replace('.raw','_b20')
          if not os.path.exists(dark_extract):
              print('extracting')
              extract(dark, dark_extract)

          dark_mask = dark_extract.replace('_b20','_mask')
          if not os.path.exists(dark_mask):
              make_mask(dark_extract, dark_mask, threshold)

          avg = dark.replace('.raw','_avg')
          if not os.path.exists(avg):
              cmd = 'python /home/drt/src/hvm3/l1/utils/makedark.py %s --mask %s %s' % (dark,dark_mask,avg)
              print(cmd)
              os.system(cmd)

          emiss = dark.replace('.raw','_emission')
          if not os.path.exists(emiss):
              cmd = 'python /home/drt/src/hvm3/l1/utils/dark2emission.py %s %s %s' % (avg,config,emiss)
              print(cmd)
              os.system(cmd)

          if temp==425:
              threshold=8000
          light_extract = light.replace('.raw','_b20')
          if not os.path.exists(light_extract):
              extract(light, light_extract)

          light_mask = light_extract.replace('_b20','_mask')
          if not os.path.exists(light_mask):
              make_mask(light_extract, light_mask, threshold)

          dksub = light.replace('.raw','_darksub')
          if not os.path.exists(dksub):
              if 'LR' in dksub:
                  cmd = 'python /home/drt/src/hvm3/l1/utils/darksubtract.py %s --max_column 370 %s %s' % (light, avg, dksub)
              elif 'RL' in dksub:
                  cmd = 'python /home/drt/src/hvm3/l1/utils/darksubtract.py %s --min_column 284 %s %s' % (light, avg, dksub)
              print(cmd)
              os.system(cmd)

          pedestal = dksub.replace('_darksub','_pedestal')
          if not os.path.exists(pedestal):
              cmd = 'python /home/drt/src/hvm3/l1/utils/hvm3_pedestal.py %s %s %s' % (dksub, config, pedestal)
              print(cmd)
              os.system(cmd)

          noemiss = dksub.replace('_darksub','_noemis')
          if not os.path.exists(noemiss):
              cmd = 'python /home/drt/src/hvm3/l1/utils/fix_emission.py %s %s %s %s' % (pedestal, config, avg, noemiss)
              print(cmd)
              os.system(cmd)

          badfix = dksub.replace('_darksub','_badfix')
          badmap = dksub.replace('_darksub','_badmap')

          if not os.path.exists(badfix):
              cmd = 'srun -n 1 -N 1 -c 40 --mem 180GB python /home/drt/src/hvm3/l1/utils/findbad.py %s --npca 3 --output_bad %s %s' % (noemiss, badmap, badfix)
              print(cmd)
              os.system(cmd)

          flat = noemiss.replace('_noemis','_flat')
          if not os.path.exists(flat):
              cmd = 'python /home/drt/src/emit-sds-l1b/utils/makeflat.py %s --config %s --mask_image %s %s' % (noemiss, config, light_mask, flat)
              print(cmd)
              os.system(cmd)
