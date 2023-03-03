# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi
from glob import glob
from scipy.optimize import minimize

c,wl,fwhm = np.loadtxt('/Users/drt/src/hvm3/l1/data/HVM3_Wavelengths_20221016.txt').T * 1000.0 
wl = np.arange(600,3810,10)

k          = 1.380648e-23    # Boltzmann constant, J K-1
q          = 1.60217733e-19  # elementary charge, in Coulombs
c_1        = 1.88365e32/np.pi # first rad. constant, photon radiance
c_2        = 14387690.0        # second radiometric constant
c_1watt    = 3.7417749e16/np.pi # first radiometric constant for W/cm2/sr/nm
hc         = 1.986e-16 # J nm

# blackbody emission in oh / (sec cm2 sr nm)
def blackbody(wvl, T):
  T = float(T)
  em = np.array([c_1 / pow(float(w),4) / (np.exp(c_2 / float(w) / float(T))-1.0) for w in wvl])
  return em
  #return c_1 / pow(wvl,4) / (s.exp(c_2 / wvl / T)-1.0)

cutoff = 3650

def spectrometer_emission(spectrometer_temp, gain, t_intp):
       spectrometer_emissivity = 1.0 
       detector_pitch = 0.003 # in centimeters
       del_wl = abs(wl[1]-wl[0])
       emission = blackbody(wl, spectrometer_temp) * \
                           spectrometer_emissivity * \
                           del_wl * t_intp 
       hemisphere     = np.pi
       detector_area  = detector_pitch**2
       qe             = 1.0
       alpha_spec_int = emission * (hemisphere * detector_area)    
       alpha_spec_int[wl>(cutoff)] = 0
       sig = sum(alpha_spec_int * qe)
       return sig * gain

def err(x,T,X,t_intp):
    gain = x[0]
    error = 0
    for temp,emis in zip(T,X):
       mdl = spectrometer_emission(temp,gain,t_intp) 
       error = error + pow(emis-mdl,2)
    return error

framerate = '4p4'
framerate = '22'

files = glob('/Users/drt/data/22HVM3/darks/*emission_*hdr')
X,T = [],[]
for fi in files:
  fname = fi.split('/')[-1]
  if not framerate in fname:
    continue
  temperature = fi.split('/')[-1].split('_')[-1].replace('.hdr','')
  print(fname,temperature)
  emis = np.squeeze(envi.open(fi).load())
  X.append(emis)
  T.append(temperature)
X = np.array(X)
T = np.array(T)


if framerate == '22':
  t_intp = 1/22
  x0 = [0.001]
  with open('HVM3_Spectrometer_Emission_44ms.txt','w') as fout:
    for i in range(X.shape[1]):
      xbest = minimize(err,x0,args=(T,X[:,i],t_intp),method='Nelder-Mead').x
      print(xbest,err(xbest,T,X[:,i],t_intp))
     #for j in range(len(T)):
     #    plt.plot(T[j], X[j,i],'ro')
     #    plt.plot(T[j], spectrometer_emission(T[j],xbest[0],t_intp),'k.')
     #plt.ylim([0,5000])
     #plt.show() 
      fout.write('%8.6f\n'%(xbest[0]))

elif framerate == '4p4':
  t_intp = 1/4.4
  x0 = [0.001]
  with open('HVM3_Spectrometer_Emission_227ms.txt','w') as fout:
    for i in range(X.shape[1]):
      xbest = minimize(err,x0,args=(T,X[:,i],t_intp),method='Nelder-Mead').x
      print(xbest,err(xbest,T,X[:,i],t_intp))
     #for j in range(len(T)):
     #    plt.plot(T[j], X[j,i],'ro')
     #    plt.plot(T[j], spectrometer_emission(T[j],xbest[0],t_intp),'k.')
     #plt.show() 
      fout.write('%8.6f\n'%(xbest[0]))


