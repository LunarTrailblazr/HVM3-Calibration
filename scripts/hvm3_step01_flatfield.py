# David R Thompson
import os
import numpy as np
from spectral.io import envi



bgfile = '/beegfs/scratch/drt/20220419_HVM3_ColdAlignment/20220510_232328_FlatField_Background/20220510_232328_UTC_FlatField_Fields614--30.raw'
bgdarkfile = '/beegfs/scratch/drt/20220419_HVM3_ColdAlignment/20220510_232328_FlatField_Background/20220510_235023_UTC_FlatField_dark.raw'
bgavgfile = '/beegfs/scratch/drt/20220419_HVM3_ColdAlignment/20220510_232328_FlatField_Background/20220510_235023_UTC_FlatField_dark_avg'
bgdksubfile = bgfile.replace('.raw','_bgsub')
bgflatfile = bgdksubfile + '_backgroundflat'
bgfile = bgflatfile.replace('backgroundflat','background')

fgfile = '/beegfs/scratch/drt/20220419_HVM3_ColdAlignment/20220510_225238_FlatField/20220510_225238_UTC_FlatField_Fields614--30.raw'
fgdarkfile = '/beegfs/scratch/drt/20220419_HVM3_ColdAlignment/20220510_225238_FlatField/20220510_231934_UTC_FlatField_dark.raw'
fgavgfile = '/beegfs/scratch/drt/20220419_HVM3_ColdAlignment/20220510_225238_FlatField/20220510_231934_UTC_FlatField_dark_avg'
fgdksubfile = fgfile.replace('.raw','_bgsub')
flatfile = fgfile.replace('.raw','_flat')

emitdir = '/home/drt/src/emit-sds-l1b/'
configfile = '/home/drt/src/hvm3/l1/config/hvm3.json'
outfile = bgfile + '_flat'
darkfile = bgfile + '_bg'

if False:
  cmd = 'python %s/utils/emit2dark.py %s %s' %(emitdir,bgdarkfile,bgavgfile)
  print(cmd)
  os.system(cmd)
  cmd = 'python %s/utils/darksubtract.py %s %s %s' %(emitdir,bgfile,bgavgfile,bgdksubfile)
  print(cmd)
  os.system(cmd)
  cmd = 'python %s/utils/emit2dark.py %s %s' %(emitdir,fgdarkfile,fgavgfile)
  print(cmd)
  os.system(cmd)
  cmd = 'python %s/utils/darksubtract.py %s %s %s' %(emitdir,fgfile,fgavgfile,fgdksubfile)
  print(cmd)
  os.system(cmd)

if True:
  cmd = 'python %s/utils/makeflat.py --config %s %s %s' %(emitdir, configfile, fgdksubfile, flatfile)
  print(cmd)
  os.system(cmd)

if False:
  I = envi.open(bgflatfile+'.hdr')
  X = np.squeeze(I.load())
  print(I.metadata)
  dn = np.array([float(f) for f in I.metadata['average_dns']])
  bg = (X.T*dn).T
  print(bg.shape)
  bg = np.asarray(bg.reshape((321,640)),dtype=np.float32)
  with open(bgfile,'wb') as fout:
      bg.tofile(fout)

if False:
  cmd = 'python %s/utils/makeflat.py --config %s %s %s' %(emitdir, configfile, fgdksubfile, flatfile)
  print(cmd)
  os.system(cmd)

fgfile = '/beegfs/scratch/drt/20220419_HVM3_ColdAlignment/20220510_225238_FlatField/20220510_225238_UTC_FlatField_Fields614--30.raw'
dksubfile = fgfile + '_darksub'
flatfile = dksubfile + '_flat'

if False:
  cmd = 'python %s/utils/darksubtract.py %s %s %s' %(emitdir,fgfile,darkfile,dksubfile)
  print(cmd)
  os.system(cmd)




testfile = dksubfile+'_uniform'
if False:
  ff = np.squeeze(envi.open(flatfile+'.hdr').load())
  with open(dksubfile,'rb') as fin:
    with open(testfile,'wb') as fout:
     for line in range(2000):
       print(line)
       X = np.fromfile(fin,count=321*640,dtype=np.float32)
       X = X.reshape((321,640))
       X = X * ff
       np.array(X,dtype=np.float32).tofile(fout)


     
dkfile = '/beegfs/scratch/drt/20220419_HVM3_ColdAlignment/20220510_225238_FlatField/20220510_225238_UTC_FlatField_dark_avg'
dksub2 = fgfile +'_traditionaldarksub'
if False:
  cmd = 'python %s/utils/darksubtract.py %s %s %s' %(emitdir,fgfile,dkfile,dksub2)
  print(cmd)
  os.system(cmd)


if False:
  cmd = 'python %s/utils/makeflat.py --config %s %s %s' %(emitdir, configfile, dksub2, flatfile)
  print(cmd)
  os.system(cmd)

testfile = dksubfile+'_uniform'
if False:
  ff = np.squeeze(envi.open(flatfile+'.hdr').load())
  with open(dksub2,'rb') as fin:
    with open(testfile,'wb') as fout:
     for line in range(2000):
       print(line)
       X = np.fromfile(fin,count=321*640,dtype=np.float32)
       X = X.reshape((321,640))
       X = X * ff
       np.array(X,dtype=np.float32).tofile(fout)
 
