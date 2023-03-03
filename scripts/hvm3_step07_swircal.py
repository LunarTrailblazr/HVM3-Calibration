# David R Thompson
import numpy as np
import pylab as plt
from spectral.io import envi
from scipy.interpolate import interp1d, splev, splrep


q,wl,fwhm = np.loadtxt('../data/HVM3_Wavelengths_20221016.txt').T * 1000.0

def resample(wl_old, spectrum, method='linear'):
  p = interp1d(wl_old, spectrum, kind=method, fill_value='extrapolate', bounds_error=False)
  return p(wl)



fivepercent = False

if fivepercent:
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_22Cal10_05_flat.hdr')
    meta = img.metadata.copy()
    dns_22a = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_22Cal11_05_flat.hdr')
    meta = img.metadata.copy()
    dns_22b = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_4p4Cal10_05_flat.hdr')
    meta = img.metadata.copy()
    dns_4p4a = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_4p4Cal11_05_flat.hdr')
    meta = img.metadata.copy()
    dns_4p4b = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)

    # Spectralon reflectance
    wl_spec, spec_rfl =\
         np.loadtxt('../data/VSWIRCal/Five_pct_Panel/Reflectance_SRT-05-100_51314-1-1_1nm.txt',skiprows=1).T  
    b2500 = np.argmin(abs(wl_spec-2500))
    spectralon_rfl = splev(wl,splrep(wl_spec, spec_rfl,s=0.1))
    plt.plot(wl,spectralon_rfl)
    plt.plot(wl_spec,spec_rfl)
    plt.show()

    dns_4p4 = np.median(np.array([dns_4p4a, dns_4p4b]),axis=0)   
    dns_22 = np.median(np.array([dns_22a, dns_22b]),axis=0)   

else:
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_22Cal4_10_flat.hdr')
    meta = img.metadata.copy()
    dns_22a = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_22Cal5_10_flat.hdr')
    meta = img.metadata.copy()
    dns_22b = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)

    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_22Cal6_10_flat.hdr')
    meta = img.metadata.copy()
    dns_22c = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_22Cal7_10_flat.hdr')
    meta = img.metadata.copy()
    dns_22d = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_4p4Cal4_10_flat.hdr')
    meta = img.metadata.copy()
    dns_4p4a = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_4p4Cal5_10_flat.hdr')
    meta = img.metadata.copy()
    dns_4p4b = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_4p4Cal6_10_flat.hdr')
    meta = img.metadata.copy()
    dns_4p4c = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)
    
    img = envi.open('/Users/drt/data/22HVM3/swircal/RFF_4p4Cal7_10_flat.hdr')
    meta = img.metadata.copy()
    dns_4p4d = np.array([float(q) for q in meta['average_dns']], dtype=np.float32)

    # Spectralon reflectance
    wl_spec, spec_rfl =\
         np.loadtxt('../data/VSWIRCal/Ten_pct_Panel/Reflectance_SRT-10-100_10AA01-0622-2521_1nm.txt',skiprows=1).T  

    spectralon_rfl = splev(wl,splrep(wl_spec, spec_rfl,s=0.1))
    #spectralon_rfl = resample(wl_spec, spec_rfl)

    plt.plot(wl_spec,spec_rfl)

    dns_22 = np.median(np.array([dns_22a, dns_22b, dns_22c, dns_22d]),axis=0)   
    dns_4p4 = np.median(np.array([dns_4p4a, dns_4p4b, dns_4p4c, dns_4p4d]),axis=0)   

plot = True

# Load irradiance, translate uncertainty from percenmt to one sigma irradiance
wl_irr, irr = np.loadtxt('../data/VSWIRCal/NIST_Lamp_1357/S1357_22.std', skiprows=1).T 
irradiance = resample(wl_irr, irr)

# Convert to microwatts
irradiance = 1e6 * irradiance

if plot:
    plt.plot(wl,spectralon_rfl)
    plt.show()

wwl, wxmit = np.loadtxt('../data/VSWIRCal/sapphire.txt').T
window_transmittance = resample(wwl,wxmit) * 0.01

plt.plot(wl,window_transmittance)
plt.plot(wl,spectralon_rfl)
plt.show()

rdn = irradiance / np.pi *  spectralon_rfl * window_transmittance

all_factors = []
#for DN, name in [(dns_22,'44ms'),(dns_4p4,'227ms')]:
#for DN, name in [(dns_4p4,'44')]:
for DN, name in [(dns_22,'227ms')]:
    factors = rdn / DN

    plt.semilogy(wl,factors)
    degree = 2

    # Interpolate over water vapor absorption at 1.88 microns
    a = np.argmin(abs(wl-1770))
    b = np.argmin(abs(wl-1790))
    c = np.argmin(abs(wl-1940))
    d = np.argmin(abs(wl-1950))
    edges = np.concatenate((np.arange(d,c), np.arange(b,a)), axis=0)
    interior = np.arange(c,b)
    channels = np.arange(len(factors))
    model = np.polyfit(channels[edges],factors[edges],degree)
    factors[interior] = np.polyval(model, channels[interior])

    # Interpolate over water vapor absorption at 1.38 microns
    a = np.argmin(abs(wl-1310))
    b = np.argmin(abs(wl-1330))
    c = np.argmin(abs(wl-1420))
    d = np.argmin(abs(wl-1430))
    edges = np.concatenate((np.arange(d,c), np.arange(b,a)), axis=0)
    interior = np.arange(c,b)
    channels = np.arange(len(factors))
    model = np.polyfit(channels[edges],factors[edges],degree)
    factors[interior] = np.polyval(model, channels[interior])

    # Show the interpolated result
    if plot:
        plt.semilogy(wl,factors)
        plt.ylim([-0.01,0.1])
        plt.show()
    all_factors.append(factors)

    factors_uncert = np.ones(len(factors))#rdn_uncert / DN

    channels = np.arange(len(rdn))
    outfile = '../data/HVM3_RadiometricCoeffs'+name+'SWIR_20221028.txt'
    np.savetxt(outfile,np.c_[channels,factors,factors_uncert], fmt='%10.8f')
    print('writing',outfile)

    rdn = DN*factors

    if plot:
        plt.plot(wl,DN*factors,'k')
        plt.box(False)
        plt.grid(True)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Measured radiance ($\mu$W sr$^{-1}$ nm$^{-1}$ cm$^{-2}$)')
        np.savetxt('swir_radiance.txt',np.c_[wl,rdn],fmt='%8.6f')
        
        plt.show()
        plt.plot(wl,rdn/interp1d([1780,1980],[rdn[np.argmin(abs(wl-1780))],rdn[np.argmin(abs(wl-1980))]],fill_value='extrapolate',bounds_error=False)(wl))
        plt.show()    

if plot and len(all_factors)>1:
    plt.plot(all_factors[0]/all_factors[1])
    plt.show()

