# First pass redshift finder of redmonster.  Able to search entire
# parameter-space and redshift-space of models in increments
# of d(loglam)/d(pixel)
#
# Tim Hutchinson, University of Utah @ IAC, May 2014
# t.hutchinson@utah.edu
from __future__ import division
from os import environ, makedirs
from os.path import join, exists
from time import gmtime, strftime

import numpy as n
from astropy.io import fits
from scipy import linalg
from scipy.optimize import nnls
import matplotlib as m
from matplotlib import pyplot as p

import multiprocessing
import sys
import time

from redmonster.datamgr.ssp_prep import SSPPrep
from redmonster.physics.misc import poly_array, two_pad
from redmonster.datamgr.io2 import read_ndArch, write_chi2arr

# Assumes all templates live in $REDMONSTER_DIR/templates/

def _zchi2(arg) :
    return zchi2_single_template(**arg)

def _zchi2_no_poly(arg) :
    return zchi2_single_template_no_poly(**arg)


def zchi2_single_template(j,poly_fft, t_fft, t2_fft, data_fft, ivar_fft, pmat_pol, bvec_pol, chi2_0, chi2_null, num_z, npixstep, zminpix, flag_val_neg_model) :
    npoly=poly_fft.shape[0] # degree+1 of polynomial
    fftnaxis1=t_fft.shape[0] # number of redshift bins
    pmat = pmat_pol.copy() # precomputed block for polynomial terms , n.zeros( (npoly+1, npoly+1, fftnaxis1), dtype=float)
    bvec = bvec_pol.copy() # precomputed vector block for polynomial terms , n.zeros( (npoly+1, fftnaxis1), dtype=float)

    # fill matrix
    pmat[0,0] = n.fft.ifft(t2_fft * ivar_fft.conj()).real
    bvec[0]   = n.fft.ifft(t_fft * data_fft.conj()).real
    for ipos in range(npoly):
        pmat[ipos+1,0] = pmat[0,ipos+1] = n.fft.ifft(t_fft*poly_fft[ipos].conj()).real

    # solve for each z

    zchi2arr=n.zeros((num_z))
    zwarning=n.zeros((num_z))
    for l in n.arange(num_z)*npixstep:
        try : # try to solve for this redshift
            f = linalg.solve(pmat[:,:,l+zminpix],bvec[:,l+zminpix])

            zchi2arr[(l//npixstep)] = chi2_0 - n.dot(n.dot(f,pmat[:,:,l+zminpix]),f) # is this true ?????
            if f[0]<0 :
                zwarning[(l//npixstep)] = int(zwarning[(l//npixstep)]) | flag_val_neg_model
                zchi2arr[(l//npixstep)] = chi2_null
                try:
                    n.dot(n.dot(f,pmat[:,:,l+zminpix]),f)
                except Exception as e:
                    print("Except: %r" % e)
                    zchi2arr[(l//npixstep)] = chi2_null
        except : # failure
            #print "failure"
            #print sys.exc_info()
            zchi2arr[(l//npixstep)] = chi2_null
    return j,zchi2arr,zwarning

def zchi2_single_template_no_poly(j,t_fft, t2_fft, data_fft, ivar_fft, chi2_0, num_z, npixstep, zminpix, flag_val_neg_model) :

    a = n.fft.ifft(t2_fft * ivar_fft.conj()).real
    b = n.fft.ifft(t_fft * data_fft.conj()).real

    ii = zminpix+n.arange(num_z)*npixstep
    f = (a[ii]!=0)*b[ii]/(a[ii]+(a[ii]==0))
    zchi2arr = chi2_0 - a[ii]*f**2 # is this true ?????

    zwarning=n.zeros((num_z))
    zchi2arr[f<0] = chi2_0
    zwarning[f<0] = flag_val_neg_model
    return j,zchi2arr,zwarning



class ZFinder:
    def __init__(self, fname=None, group=[0], npoly=None, zmin=None, zmax=None, nproc=1):
        self.fname = fname
        if type(group) == list:
            self.group = group
        elif type(group) == int:
            self.group = [group,]
        else:
            self.group = eval(group)
        self.npoly = npoly if npoly is not None else 4
        self.zmin = float(zmin)
        self.zmax = float(zmax)
        self.nproc = nproc
        self.pixoffset = None
        self.zchi2arr = None

        try:
            self.templatesdir = environ['REDMONSTER_TEMPLATES_DIR']
        except KeyError as e:
            self.templatesdir = None
            print("Enviromental variable 'REDMONSTER_TEMPLATES_DIR' not set: \
                    %r" % e)
        self.read_template()
        self.npars = len(self.templates.shape) - 1
        self.templates_flat = n.reshape(self.templates, (-1,self.fftnaxis1))
        self.t_fft = n.fft.fft(self.templates_flat)
        self.t2_fft = n.fft.fft(self.templates_flat**2)
        self.f_nulls = []
        self.chi2_null = []
        self.sn2_data = []


    def read_template(self):
        self.templates, self.baselines, self.infodict = \
                read_ndArch(join(self.templatesdir,self.fname))
        self.type = self.infodict['class']
        self.origshape = self.templates.shape
        self.ntemps = self.templates[...,0].size
        self.fftnaxis1 = two_pad(self.origshape[-1])
        self.tempwave = 10**(self.infodict['coeff0'] + \
                             n.arange(self.infodict['nwave']) * \
                             self.infodict['coeff1'])
        templates_pad = n.zeros( self.origshape[:-1]+(self.fftnaxis1,) )
        templates_pad[...,:self.origshape[-1]] = self.templates
        self.templates = templates_pad


    def create_z_baseline(self, loglam0):
        self.zbase = ((10**loglam0)/self.tempwave) - 1


    def conv_zbounds(self):
        zmaxpix = n.where( abs((self.zbase-self.zmin)) == \
                          n.min(abs(self.zbase-self.zmin)) )[0][0]
        zminpix = n.where( abs((self.zbase-self.zmax)) == \
                          n.min(abs(self.zbase-self.zmax)) )[0][0]
        return zminpix, zmaxpix


    def zchi2(self, specs, specloglam, ivar, npixstep=1, chi2file=False,
              plate=None, mjd=None, fiberid=None):
        self.chi2file = chi2file
        self.npixstep = npixstep
        self.zwarning = n.zeros(specs.shape[0])
        flag_val_unplugged = int('0b10000000',2)
        flag_val_neg_model = int('0b1000',2)
        self.create_z_baseline(specloglam[0])
        if (self.zmin != None) and (self.zmax != None) and \
                (self.zmax > self.zmin):

            zminpix, zmaxpix = self.conv_zbounds()
            self.pixoffset = zminpix
            num_z = int(n.floor( (zmaxpix - zminpix) / npixstep ))
            zinds = zminpix + n.arange(num_z)*npixstep
            self.zbase = self.zbase[zinds]
        else:
            zminpix = 0
            # Number of pixels to be fitted in redshift
            num_z = int(n.floor( (zself.origshape[-1] - specs.shape[-1]) /
                                npixstep ))

        # Create arrays for use in routine
        # Create chi2 array of shape (# of fibers, template_parameter_1,
        # ..., template_parameter_N, # of redshifts)
        zchi2arr = n.zeros((specs.shape[0], self.templates_flat.shape[0],
                            num_z))
        temp_zwarning = n.zeros(zchi2arr.shape)

        # Pad data and SSPs to a power of 2 for faster FFTs
        data_pad = n.zeros(specs.shape[:-1] + (self.fftnaxis1,), dtype=float)
        data_pad[...,:specs.shape[-1]] = specs
        ivar_pad = n.zeros(ivar.shape[:-1] + (self.fftnaxis1,), dtype=float)
        ivar_pad[...,:specs.shape[-1]] = ivar

        # Pre-compute FFTs for use in convolutions
        data_fft = n.fft.fft(data_pad * ivar_pad)
        ivar_fft = n.fft.fft(ivar_pad)

        if self.npoly>0 :
            # Compute poly terms, noting that they will stay fixed with
            # the data - assumes data is passed in as shape (nfibers, npix)
            polyarr = poly_array(self.npoly, specs.shape[1])
            pmat = n.zeros( (self.npoly+1, self.npoly+1, self.fftnaxis1),
                            dtype=float)
            bvec = n.zeros( (self.npoly+1, self.fftnaxis1), dtype=float)

            # Pad to a power of 2 for faster FFTs
            poly_pad = n.zeros((self.npoly, self.fftnaxis1), dtype=float)
            poly_pad[...,:polyarr.shape[-1]] = polyarr

            # Pre-compute FFTs for use in convolutions
            poly_fft = n.zeros((ivar_pad.shape[0], self.npoly, self.fftnaxis1),dtype=complex)
            for i in range(self.npoly):
                poly_fft[:,i,:] = n.fft.fft(poly_pad[i] * ivar_pad)





        # Compute z for all fibers

        for i in range(specs.shape[0]): # Loop over fibers

            start=time.time()

            #print 'INFO Fitting fiber %s of %s for template %s' % \
            #        (i+1, specs.shape[0], self.fname)

            # If flux is all zeros, flag as unplugged according to BOSS
            # zwarning flags and don't bother with doing fit
            if len(n.where(specs[i] != 0.)[0]) == 0:
                self.zwarning[i] = int(self.zwarning[i]) | flag_val_unplugged
            else: # Otherwise, go ahead and do fit

                self.sn2_data.append (n.sum( (specs[i]**2)*ivar[i] ) )

                if self.npoly>0 :
                    for ipos in range(self.npoly):
                        bvec[ipos+1] = n.sum( poly_pad[ipos] * data_pad[i] *
                                              ivar_pad[i])
                    for ipos in range(self.npoly):
                        for jpos in range(self.npoly):
                            pmat[ipos+1,jpos+1] = n.sum( poly_pad[ipos] *
                                                         poly_pad[jpos] *ivar_pad[i]) # CAN GO FASTER HERE (BUT NOT LIMITING = 0.001475

                    f_null = linalg.solve(pmat[1:,1:,0],bvec[1:,0])
                    self.f_nulls.append( f_null )
                    self.chi2_null.append( self.sn2_data[i] -
                                           n.dot(n.dot(f_null,pmat[1:,1:,0]),f_null))
                else :
                    self.chi2_null.append( self.sn2_data[i] )
                # print 'INFO Chi^2_null value is %s' % self.chi2_null[i]
                # Loop over templates
                # multiprocessing

                func_args = []
                for j in range(self.templates_flat.shape[0]):
                    if self.npoly>0 :
                        arguments = {"j":j,"poly_fft":poly_fft[i], "t_fft":self.t_fft[j], "t2_fft":self.t2_fft[j], "data_fft":data_fft[i], "ivar_fft":ivar_fft[i], "pmat_pol":pmat, "bvec_pol":bvec, "chi2_0":self.sn2_data[i], "chi2_null":self.chi2_null[i],"num_z":num_z, "npixstep":self.npixstep, "zminpix":zminpix,"flag_val_neg_model":flag_val_neg_model}
                    else :
                        arguments = {"j":j,"t_fft":self.t_fft[j], "t2_fft":self.t2_fft[j], "data_fft":data_fft[i], "ivar_fft":ivar_fft[i], "chi2_0":self.sn2_data[i], "num_z":num_z, "npixstep":self.npixstep, "zminpix":zminpix, "flag_val_neg_model":flag_val_neg_model}
                    func_args.append(arguments)

                results = None
                if self.nproc > 1:
                    pool = multiprocessing.Pool(self.nproc)
                    if self.npoly>0 :
                        results = pool.map(_zchi2, func_args)
                    else :
                        results = pool.map(_zchi2_no_poly, func_args)
                    pool.close()
                    pool.join()
                else:
                    if self.npoly>0 :
                        results = [ _zchi2(x) for x in func_args ]
                    else :
                        results = [ _zchi2_no_poly(x) for x in func_args ]

                for result in results :
                    j                  = result[0]
                    zchi2arr[i,j]      = result[1]
                    temp_zwarning[i,j] = result[2]

                stop=time.time()

                print("INFO fitted fiber %d/%d, chi2_null=%f, %d templates in %s, npoly=%d, using %d procs in %f sec"%(i+1, specs.shape[0],self.chi2_null[i],self.templates_flat.shape[0],self.fname,self.npoly,self.nproc,stop-start))


        # Use only neg_model flag from best fit model/redshift and add
        # it to self.zwarning
        for i in range(self.zwarning.shape[0]):
            minpos = ( n.where(zchi2arr[i] == n.min(zchi2arr[i]))[0][0],
                      n.where(zchi2arr[i] == n.min(zchi2arr[i]))[1][0] )
            self.zwarning[i] = int(self.zwarning[i]) | \
                    int(temp_zwarning[i,minpos[0],minpos[1]])
        zchi2arr = n.reshape(zchi2arr, (specs.shape[0],) + self.origshape[:-1] +
                             (num_z,) )
        bestl = n.where(zchi2arr == n.min(zchi2arr))[-1][0]
        thisz = ((10**(specloglam[0]))/self.tempwave[bestl+zminpix])-1

        #return zchi2arr
        self.zchi2arr = zchi2arr
        #self.store_models(specs, ivar)
        if self.chi2file is True:
            if (plate is not None) & (mjd is not None) & (fiberid is not None):
                write_chi2arr(plate, mjd, fiberid, self.zchi2arr)
            else:
                print('WARNING Plate/mjd/fiberid not given - unable to write chi2 file!')
        else:
            print('INFO Not writing chi2')

    def store_models(self, specs, ivar):
        self.models = n.zeros( (specs.shape) )
        for i in range(self.models.shape[0]):
            minloc = n.unravel_index( self.zchi2arr[i].argmin(),
                                     self.zchi2arr[i].shape )
            pmat = n.zeros( (specs.shape[-1],self.npoly+1) )
            this_temp = self.templates[minloc[:-1]]
            pmat[:,0] = this_temp[(minloc[-1]*self.npixstep) +
                                  self.pixoffset:(minloc[-1]*self.npixstep) +
                                  self.pixoffset+specs.shape[-1]]
            polyarr = poly_array(self.npoly, specs.shape[-1])
            pmat[:,1:] = n.transpose(polyarr)
            ninv = n.diag(ivar[i])
            try: # Some eBOSS spectra have ivar[i] = 0 for all i
                '''
                f = linalg.solve( n.dot(n.dot(n.transpose(pmat),ninv),pmat),
                                n.dot( n.dot(n.transpose(pmat),ninv),specs[i]) )
                '''
                f = nnls( n.dot(n.dot(n.transpose(pmat),ninv),pmat),
                         n.dot( n.dot(n.transpose(pmat),ninv),specs[i]) );\
                                f = n.array(f)[0]
                self.models[i] = n.dot(pmat, f)
            except Exception as e:
                self.models[i] = n.zeros(specs.shape[-1])
                print("Exception: %r" % r)


    def write_chi2arr(self, plate, mjd, fiberid):
        prihdu = fits.PrimaryHDU(self.zchi2arr)
        thdulist = fits.HDUList([prihdu])
        try:
            rsr = environ['REDMONSTER_SPECTRO_REDUX']
            run2d = environ['RUN2D']
            run1d = environ['RUN1D']
            if (rsr is not None) & (run2d is not None) & (run1d is not None):
                testpath = join(rsr, run2d, '%s' % plate, run1d)
                if exists(testpath):
                    dest = testpath
                else:
                    try:
                        makedirs(testpath)
                        dest = testpath
                    except Exception as e:
                        print("Exception: %r" % e)
                        dest = None
        except Exception as e:
            print("Exception: %r" % e)
            dest = None
        if dest is not None:
            try:
                thdulist.writeto(join(dest, '%s' % 'chi2arr-%s-%s-%s-%03d.fits'
                                      % (self.type, plate, mjd, fiberid)),
                                 clobber=True)
                print('Writing chi2 file to %s' % \
                        join(dest, '%s' % 'chi2arr-%s-%s-%s-%03d.fits' %
                             (self.type, plate, mjd, fiberid)))
            except Exception as e:
                print("Exception: %r" % e)
                print('Environment variables not set or path does not exist - \
                        not writing chi2 file!')
        else:
            print('Environment variables not set or path does not exist - \
                    not writing chi2 file!')
