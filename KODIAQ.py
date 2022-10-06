#!/usr/bin/python3

from astropy.io import fits
from astropy import units as u
import numpy as np
import pandas
from matplotlib import pyplot as plt
import glob
from astropy.visualization import quantity_support
import statistics
quantity_support()  # for getting units on the axes below 
    
class Spectra():
    
        def __init__(self, catalog):
            'input the directory where you have the KODIAQ data'
            self.catalog = catalog 
            
        def values(self, qsoname):
            self.qsoname =qsoname
            qsocat = pandas.read_csv(self.catalog)
            index = qsocat[qsocat['QSO'] == self.qsoname] 
            number = index.index[0]
            #finds the redshift value for that qso and converts it to a string
            redshift = index['z'][number]
            redshiftname = str(redshift)
            self.redshiftname = redshiftname
    
        #finds the directory name in the column
            dirname = index['File Paths']
            
    
        #this is the for loop to find the fits data out of the fits file
            for fits_file_names in dirname:
    #opening the fits files 
                hdu1 = fits.open(fits_file_names)
    #finding the header
                hdr = hdu1[0].header
    #find spec flux 
                flux = hdu1[0].data
                flux_w_units = flux*(10**16)
        #this is used to get rid of a problimatic area of the flux spectrum that is seen in all spectra
                fluxarray = np.array(flux_w_units[730:])
                self.flux = fluxarray
    #finding error spec
                errorspec = hdu1[1].data
                errorspec_w_units = errorspec*(10**16)
        #this is used to get rid of a problimatic area of the flux spectrum that is seen in all spectra
                errorspecarray = np.array(errorspec_w_units[730:])
                self.error = errorspecarray
    #find cont flux 
                cont = hdu1[3].data
                cont_w_units = cont*(10**16)
        #this is used to get rid of a problimatic area of the flux spectrum that is seen in all spectra
                contarray = np.array(cont_w_units[730:])
                self.cont = contarray
    #find wave lenght range
            CVAL= hdr['CRVAL1']
            CDEL = hdr['CDELT1']
        
    #find wave lennght with correct units 
            wave = 10**(CVAL + np.arange(hdr['NAXIS1'])*CDEL)
        #doing the same thing to match the wavelenghts to area we removed in spectra
            wave_array =wave[730:]
            self.obswave = wave_array
        
        #finding the value of the rest frame wavelength
            restwave =wave_array/(redshift+1)
            self.restwave = restwave
        #finds the first index with 0 so we can find the last flux value
            zeroindex = int(np.argwhere(fluxarray == 0)[0:1])
        #min and max of graph in rest fram
            minrest = min(wave_array)/(redshift+1)
            self.minrest = minrest
            maxrest = max(wave_array[:zeroindex])/(redshift+1)
            self.maxrest = maxrest
            minobs = min(wave_array)
            self.minobs = minobs
            maxobs = max(wave_array[:zeroindex])
            self.maxobs = maxobs
        #this is the max flux so we can plot the graph with good bounds
            fluxmax = max(fluxarray)
            self.fluxmax = fluxmax
        #finding the flux valuse in this range, which is from 1100 to 1150 Angstroms
            fluxrange = fluxarray[2620:3960]
            avgflux = statistics.mean(fluxrange)
        #normalizing
            normflux = fluxarray/avgflux
            self.normflux= normflux
            normcont = contarray/avgflux
            self.normcont = normcont
            normerror = errorspecarray/avgflux
            self.normerror = normerror 
        #max bound for normplot
            normmax = max(normflux)
            self.normmax= normmax
        #min and max lya range
            minlya = 1020
            maxlya =1240
            self.minlya = minlya 
            self.maxlya = maxlya
            return(self.flux, self.cont, self.error, self.restwave, self.obswave,self.minrest,
                   self.maxrest ,self.minobs,self.maxobs,self.fluxmax, self.redshiftname,self.normflux,
                  self.normcont,self.normerror, self.normmax,self.minlya,self.maxlya)
        
        def restplot(self):
            return(plt.figure(figsize = (15,6)),
               plt.ylim(0,self.fluxmax),
               plt.xlim(self.minrest, self.maxrest),
               plt.xlabel('Rest Wavelength ($\AA$)'),
               plt.ylabel('Normalized Flux'),
               plt.title(self.qsoname),
               plt.title('z = ' + self.redshiftname, loc ='right'),
               plt.plot(self.restwave, self.flux, lw= .2, color = 'blue', label = 'Spectral Flux') ,
               plt.plot(self.restwave, self.cont, lw= 1, color = 'r', label = 'Continuum'),
               plt.plot(self.restwave, self.error,  lw = .4, color = 'y', label = 'Measurment Error'),
               plt.legend())
        def obsplot(self):
            return(plt.figure(figsize = (15,6)),
               plt.ylim(0,self.fluxmax),
               plt.xlim(self.minobs, self.maxobs),
               plt.xlabel('Observed Wavelength ($\AA$)'),
               plt.ylabel('Normalized Flux'),
               plt.title(self.qsoname),
               plt.title('z = ' + self.redshiftname, loc ='right'),
               plt.plot(self.obswave, self.flux, lw= .2, color = 'blue', label = 'Spectral Flux') ,
               plt.plot(self.obswave, self.cont, lw= 1, color = 'r', label = 'Continuum'),
               plt.plot(self.obswave, self.error,  lw = .4, color = 'y', label = 'Measurment Error'),
               plt.legend())
        def lyarestnorm(self):
            return(plt.figure(figsize = (15,6)),
               plt.ylim(0,self.normmax),
               plt.xlim(self.minlya, self.maxlya),
               plt.xlabel('Rest Wavelength ($\AA$)'),
               plt.ylabel('Normalized Flux'),
               plt.title(self.qsoname),
               plt.title('z = ' + self.redshiftname, loc ='right'),
               plt.plot(self.restwave, self.normflux, lw= .2, color = 'blue', label = 'Spectral Flux') ,
               plt.plot(self.restwave, self.normcont, lw= 1, color = 'r', label = 'Continuum'),
               plt.plot(self.restwave, self.normerror,  lw = .4, color = 'y', label = 'Measurment Error'),
               plt.legend())
            