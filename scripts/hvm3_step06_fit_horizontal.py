# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi
from scipy.interpolate import splrep, splev
from glob import glob


zone1 = envi.open('/Users/drt/data/22HVM3/darks/RFF_4p4Dark6_avg_zone1_189.8.hdr').load()
zone1 = np.squeeze(zone1)
x = np.arange(len(zone1))
sp = splrep(x, zone1, s=2e6)
plt.plot(zone1,'k.')
zone1 = splev(x, sp)
plt.plot(zone1,'r')

zone2 = envi.open('/Users/drt/data/22HVM3/darks/RFF_4p4Dark6_avg_zone2_189.8.hdr').load()
zone2 = np.squeeze(zone2)
sp = splrep(x, zone2, s=2e6)
plt.plot(zone2,'k.')
zone2 = splev(x, sp)
plt.plot(zone2,'r')


zone3 = envi.open('/Users/drt/data/22HVM3/darks/RFF_4p4Dark6_avg_zone3_189.8.hdr').load()
zone3 = np.squeeze(zone3)
zone3[178:183] = (zone3[177]+zone3[183])/2
sp = splrep(x, zone3, s=2e6)
plt.plot(zone3,'k.')
zone3 = splev(x, sp)
plt.plot(zone3,'r')
plt.show()

print(zone1.shape)
coeffs = np.zeros((321,640))
# 156, 256
coeffs[256:,:] = zone1/np.median(zone1[11:19])
coeffs[156:256,:] = zone2/np.median(zone2[11:19])
coeffs[:156,:] = zone3/np.median(zone3[11:19])

spectr = np.loadtxt('HVM3_Spectrometer_Emission_4p4.txt')
spectr[0]=0
x = np.arange(len(spectr))
sp = splrep(x, spectr, s=0.000005)
plt.plot(spectr,'k.')
spectr = splev(x, sp)
plt.plot(spectr,'r')
plt.show()

total = (coeffs.T * spectr).T
outfile = '/Users/drt/src/hvm3/l1/data/HVM3_EmissionFactors_20221026.hdr'
envi.save_image(outfile,np.array(total,dtype=np.float32),ext='',force=True)
