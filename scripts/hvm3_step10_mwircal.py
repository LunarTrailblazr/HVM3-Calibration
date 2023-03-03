# David R Thompson
import numpy as np
from glob import glob
from datetime import datetime, timedelta 
from collections import OrderedDict
import pylab as plt
from spectral.io import envi
from scipy.interpolate import interp1d

q,wl,fwhm = np.loadtxt('../data/HVM3_Wavelengths_20221016.txt').T * 1000.0

# Timestamps and temperature records
def celcius2kelvin(c):
    return c+273.15

timestamp_dir = '../data/temperature_logs/'
files = glob(timestamp_dir+'202211*.txt')
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
      bb_ctr_left_T = celcius2kelvin(float(tokens[37]))
      bb_ctr_right_T = celcius2kelvin(float(tokens[37]))
      bb_T = np.mean([bb_ctr_left_T,bb_ctr_right_T])
      epoch = datetime(year=1904,month=1,day=1,hour=0,minute=0,second=0)
      offset = timedelta(seconds=int(time_s))
      time = epoch+offset
      temperatures.append((time,bb_T))
temperatures.sort()
    
tol = timedelta(seconds=20)

# translate UTC to blackbody temperature
def utc2temp(temperatures, query, lo=None, hi=None, margin=5):
  if query > temperatures[-1][0] or query < temperatures[0][0]:
     raise ValueError('out of bounds time!')
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

# Black body emitted radiance in uW sr-1 cm-2 nm-1
def planck(emissivity, T, wl):
    # Radiance of a surface due to emission"""
    c_1 = 1.88365e32/np.pi
    c_2 = 14387690
    J_per_eV = 1.60218e-19
    wl_um = wl / 1000.0
    ph_per_sec_cm2_sr_nm = c_1/(wl**4)/(np.exp(c_2/wl/T)-1.0) * emissivity
    # photon energy in eV
    eV_per_sec_cm2_sr_nm = 1.2398 * ph_per_sec_cm2_sr_nm/wl_um
    W_per_cm2_sr_nm = J_per_eV * eV_per_sec_cm2_sr_nm
    uW_per_cm2_sr_nm = W_per_cm2_sr_nm*1e6
    return uW_per_cm2_sr_nm

def resample(wl_old, spectrum, method='linear'):
  p = interp1d(wl_old, spectrum, kind=method, fill_value='extrapolate', bounds_error=False)
  return p(wl)

bwl, emiss = np.loadtxt('../data/MWIRCal/bb_emissivity.txt').T
emiss = 0.0*resample(bwl, emiss) +0.998

test = datetime(year=1900,month=1,day=1)
base = '/Users/drt/data/22HVM3/mwircal/'
electrons_per_DN = 42

for integration_time in ['227ms','44ms']: 
  legends,ratio_food,names = [],[],[]
  for flatpath in sorted(glob(base+'/*'+integration_time+'*flat')):

    direction = 'RL'
    if 'LR' in flatpath:
        direction = 'LR' 

    img = envi.open(flatpath+'.hdr')
    meta = img.metadata.copy()
    dns = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    flat = img.load().copy()
        
    T = float(flatpath.split('/')[-1].split('_')[-4].replace('BB',''))
    
    timestr = ''.join(flatpath.split('/')[-1].split('_')[:2])
    dt = test.strptime(timestr,'%Y%m%d%H%M%S') 

   #if dt > test.strptime('20221110','%Y%m%d'):
   #    continue

    temperature = utc2temp(temperatures, query=dt)
    if not temperature:
        temperature = T
    print(T,temperature)

    rdn = planck(emiss, temperature, wl)
    factors = rdn / dns

    np.savetxt('mwir_radiance_%iK.txt'%T,np.c_[wl,rdn],fmt='%8.6f')

    pattern = base+'/*BB%i*%s_Dark_%s_emission'%(T,integration_time,direction)
    emission = glob(pattern)
    print(pattern)
    if len(emission)>0:
        emitted_dns = envi.open(emission[0]+'.hdr').metadata['average_dns']
        emitted_dns = np.array([float(q) for q in emitted_dns])
        outfile = '../data/MWIRCal/HVM3_DN2RDN_%s_BB%i_%s_MWIR.txt' %(integration_time,T,direction)
        np.savetxt(outfile,np.c_[wl,emitted_dns,dns,rdn], fmt='%10.8f',header='Wavelength (nm), Self-emission signal (DN), Scene signal (DN), Radiance (uW nm-1 sr-1 cm-2)')

   #forfitting = np.logical_and(wl>2500,np.logical_and(wl<3400,np.logical_or(wl<2600,wl>3100)))
   #ctm = np.polyfit(wl[forfitting],factors[forfitting],3)
   #plt.plot(wl,factors)
   #plt.plot(wl,np.polyval(ctm,wl))
   #plt.show()
    
   #plt.plot(wl,np.polyval(ctm,wl)/factors)
   #plt.ylim([0,2])
   #plt.show()
   
    outfile = '../data/MWIRCal/HVM3_RadiometricCoeffs_%s_BB%i_%s_MWIR.txt' %(integration_time,T,direction)
    channels = np.arange(len(wl))
    factors_uncert = np.ones(len(wl))
    np.savetxt(outfile,np.c_[channels,factors,factors_uncert], fmt='%10.8f')
    legends.append(outfile.split('/')[-1])
    if integration_time == '227ms':
        plt.title('227 ms integration') 
    else:
        plt.title('44 ms integration')
    factors[factors<0]=np.nan
    plt.semilogy(wl,factors)
    if T<350:
       ratio_food.append(factors)
       names.append(flatpath)

  plt.legend(legends)
  plt.show()
  legends = []
  for i in range(100):
     a=np.random.randint(0,high=len(ratio_food))
     b = a
     while a == b:
         b=np.random.randint(0,high=len(ratio_food))
     legends.append('%i,%i'%(a,b))
     plt.plot(wl,ratio_food[a]/ratio_food[b])
  plt.legend(legends)
 #for i,name in enumerate(names):
 #   print(i,name)  
  plt.ylim([0.5,1.5])
  plt.show()
