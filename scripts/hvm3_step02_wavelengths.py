# David R Thompson
import pylab as plt
import numpy as np
from scipy.stats import norm
from spectral.io import envi
from scipy.interpolate import interp1d
import astropy.modeling as modeling
from sklearn.linear_model import RANSACRegressor
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter
from scipy.optimize import minimize
from scipy.signal import medfilt
from numpy.random import randn
import json

# Fit a gaussian to a single peak in an ordered series of data
def find_peak(x, plot=False):
    fitter = modeling.fitting.LevMarLSQFitter()
    model = modeling.models.Gaussian1D(amplitude=np.max(x),
                                       mean=np.argmax(x),
                                       stddev=1.0/2.35)   # depending on the data you need to give some initial values
    fitted_model = fitter(model, np.arange(len(x)), x)
    if plot:
         print(fitted_model.mean[0], fitted_model.amplitude[0], fitted_model.stddev[0])
         plt.plot(x)
         plt.plot(fitted_model(np.arange(len(x))))
         plt.show()
    return fitted_model.mean[0], fitted_model.amplitude[0], fitted_model.stddev[0]


# We will analyze the laser sphere data from TVAC. 
basedir = '/Users/drt/data/22HVM3/'
I = envi.open(basedir+'20221006_215546_UTC_LaS_Fields470-173_darksub_avg.hdr')
I = I.load()

# Laser info 
# Wavelengths from Holly.  We start with an initial first guess of the channels
# based on her mesurements.
wavelengths = np.array([3391,1899,1308,631])
channel_guess = np.array([33.04,182.24,241.34,308.96])
nlasers = len(wavelengths)

# Change to refractive wavelength of vacuum
index_of_refraction = np.array([1,1,1,1])  # already done
wavelengths = wavelengths * index_of_refraction
 
margin = 10
ncols = 640
nrows = 321
mid = int(ncols/2)
use = np.arange(170,475)

# our list of laser fits, one sublist per laser
fits = [[] for c in channel_guess]

# Find a the spatial location of the lasers, and fit each laser peak
# put them as column, row, wavelength triplets into the "observed" list
frame = np.squeeze(I[:,:,0])
robust_model = RANSACRegressor()
for i, chn in zip(range(nlasers), channel_guess):
  for col in range(0,ncols):

    # Extract the segment with the laser line
    idx = np.arange(int(chn-margin),int(chn+margin+1), dtype=int)
    n = len(idx)
    meas = frame[idx,int(round(col))]

    # Subtract the continuum
    x = np.arange(n)[:,np.newaxis]
    robust_model.fit(x,meas)
    ctm = robust_model.predict(x)
    meas = meas - ctm

    # peak fitting
    row,_,_ = find_peak(meas)#, plot=(col==300))
    row = row+idx[0] 
    fits[i].append([col,row]) 

  fits[i] = np.array(fits[i])

# Now plot the result, and save the laser fits to a file for later plotting
X = np.zeros((nlasers,ncols))
Y = np.zeros((nlasers,ncols))

for i in range(nlasers):

    # columns and rows for a single laser
    cols, rows = fits[i].T

    # fit the line, centering the horizontal axis on zero
    # so that the intercept is the shift at column 320
    robust_model.fit(cols[use,np.newaxis]-mid,rows[use])
    x = np.arange(ncols)[:,np.newaxis] - mid
    rows_est = robust_model.predict(x)
    intercept = robust_model.estimator_.intercept_
    slope = robust_model.estimator_.coef_[0]
    print('slope',slope,'intercept',intercept)

    plt.plot(cols,rows-intercept,'.')
    plt.plot(np.arange(ncols),rows_est-intercept,'k')
    np.savetxt('../data/plots/HVM3_Laser_%i_ColRow.txt'%wavelengths[i],
            np.c_[cols,rows],fmt='%8.6f')

    X[i,:] = np.squeeze(rows_est)
    Y[i,:] = wavelengths[i]
     
plt.show()


# Perform the fit for each column
ctrs = np.zeros((nrows, ncols, 2))
for i in range(ncols):
  p = np.polyfit(X[:,i],Y[:,i],1)
  wl = np.polyval(p, np.arange(nrows))
  err = max(Y[:,i]-np.polyval(p,X[:,i]))
  print('error for column',i,':',err)
  ctrs[:,i,0] = wl

# Now simulate wavelength error due to uncertain wavelengths
errs = []
ref = mid
uncert = 0.5
for trial in range(1000):
    a = X[:,ref]
    y = wavelengths
    p = np.polyfit(a,y,1)
    y2 = y + randn(len(y)) * uncert
    p2 = np.polyfit(a,y2,1)
    err = np.polyval(p,np.arange(nrows)) - np.polyval(p2,np.arange(nrows))
    errs.append(err)
errs = abs(np.array(errs)).mean(axis=0)
plt.plot(errs)
plt.show()
for c in range(ncols):
    ctrs[:,c,1] = errs

# save out wavelength center matrix
envi.save_image('../data/HVM3_WavelengthCenters_20221016.hdr',np.array(ctrs,dtype=np.float32),ext='',force=True)

# save out ascii text file, no FWHMs, in microns
wvl = np.median(np.squeeze(ctrs[:,:,0]),axis=1)
fwhm = np.zeros(nrows)
chn = np.arange(nrows)
np.savetxt('../data/HVM3_Wavelengths_20221016.txt',np.c_[chn,wvl/1000.0,fwhm], fmt='%10.8f')
