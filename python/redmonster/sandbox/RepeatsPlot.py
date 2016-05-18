import pyfits
import glob
import numpy as N
import pylab as P
import os
import shutil
import sys

class Repeatability:

    directory = '/uufs/astro.utah.edu/common/home/u0814744/compute/scratch/repeatability'

    @staticmethod
    def makeSpAllRepeats():

        hdu = pyfits.open('/uufs/chpc.utah.edu/common/home/sdss02/ebosswork/eboss/spectro/redux/v5_8_0/spAll-v5_8_0.fits')

        print 'Selecting LRGs'
        a = hdu[1].data
        w = N.where( a.EBOSS_TARGET1 & 2 > 0)[0]
        a = a[w]


        wr = list()

        for i in range(len(a)):
            sys.stderr.write('\r '+str(i)+' '+str(len(a)))
            t = a.THING_ID[i]
            w = N.where( a.THING_ID == t )[0]
            if len(w) > 1:
                wr.extend(w)
        print len(wr)

        wr = N.unique(N.array(wr))
        newa = a[wr]
        print len(N.unique( newa.THING_ID )), 'repeats'
        newhdu = pyfits.BinTableHDU(data=newa)
        newhdu.writeto(Repeatability.directory+'/spAll-v5_8_0-repeats_lrg.fits')

    @staticmethod
    def getDeltaVexistent():

        a = pyfits.open(Repeatability.directory+'/spAll-v5_8_0-repeats_lrg.fits')[1].data

        data = N.empty( 0, dtype= [ ('thing_id', int), ('dv', float), ('z1', float), \
                                    ('z2', float), ('dc1', float), ('dc2', float), \
                                    ('dc_min', float)])
        for i in range(len(a)):
            t = a.THING_ID[i]
            if t in data['thing_id']:
                continue
            w = N.where( a.THING_ID == t)[0]
            if a.RCHI2DIFF_NOQSO[w[0]] < a.RCHI2DIFF_NOQSO[w[1]]:
                j1 = w[1]
                j2 = w[0]
            else:
                j1 = w[0]
                j2 = w[1]

            z1 = a.Z_NOQSO[j1]
            z2 = a.Z_NOQSO[j2]
            dc1 = a.RCHI2DIFF_NOQSO[j1]
            dc2 = a.RCHI2DIFF_NOQSO[j2]

            c_kms = 299792.458
            dv = abs(z1-z2)*c_kms/(1+min([z1, z2]))
            dc_min = min([dc1, dc2])

            data = N.append(data, N.array( (t, dv, z1, z2, dc1, dc2, dc_min), dtype=data.dtype))

        return data

    @staticmethod
    def getDeltaVsplits(guy=0):

        dir = '/uufs/chpc.utah.edu/common/home/sdss02/ebosswork/eboss/spectro/redux/test/bautista'
        if guy:
            a = pyfits.open(dir+'/v5_8_guy_split1/spAll-v5_8_guy_split1.fits')[1].data
            b = pyfits.open(dir+'/v5_8_guy_split2/spAll-v5_8_guy_split2.fits')[1].data
        else:
            a = pyfits.open(dir+'/split1/spAll-split1.fits')[1].data
            b = pyfits.open(dir+'/split2/spAll-split2.fits')[1].data

        w = N.where( a.EBOSS_TARGET1 & 2 > 0)[0]
        ngals = len(w)
        a = a[w]
        b = b[w]

        data =   N.empty( 0, dtype= [ ('thing_id', int), ('dv', float), ('z1', float), \
                                    ('z2', float), ('dc1', float), ('dc2', float), \
                                    ('dc_min', float)])
        for i in range(ngals):
            t = a.THING_ID[i]
            if b.THING_ID[i] != t:
                print 'Mismatch in thing_id!', t, b.THING_ID[i]

            if a.RCHI2DIFF_NOQSO[i] < b.RCHI2DIFF_NOQSO[i]:
                z1 = b.Z_NOQSO[i]
                z2 = a.Z_NOQSO[i]
                dc1 = b.RCHI2DIFF_NOQSO[i]
                dc2 = a.RCHI2DIFF_NOQSO[i]
            else:
                z1 = a.Z_NOQSO[i]
                z2 = b.Z_NOQSO[i]
                dc1 = a.RCHI2DIFF_NOQSO[i]
                dc2 = b.RCHI2DIFF_NOQSO[i]

            c_kms = 299792.458
            dv = abs(z1-z2)*c_kms/(1+ min([z1, z2]))
            dc_min = min([dc1, dc2])

            data = N.append(data, N.array( (t, dv, z1, z2, dc1, dc2, dc_min), dtype=data.dtype))


        return data

    @staticmethod
    def getRepeatsData():
        data1 = Repeatability.getDeltaVexistent()
        data2 = Repeatability.getDeltaVsplits()
        data = N.append(data1, data2)
        return data

    @staticmethod
    def makeComp(save=0):

        data = Repeatability.getRepeatsData()

        dc = data['dc_min']
        dv = data['dv']

        P.figure()
        P.plot( dc, dv, 'k.', alpha=0.4)
        P.ylim(0.1, 1e6)
        P.xlim(1e-6, 1)
        P.xscale('log')
        P.yscale('log')
        ylim = P.ylim()
        P.plot( [1e-2, 1e-2], ylim, 'b--', lw=2)
        P.plot( [5e-3, 5e-3], ylim, 'r--', lw=2)
        P.plot( [1e-6, 1], [1000, 1000], 'm--', lw=2)
        P.xlabel(r'$\Delta \chi^2/dof$')
        P.ylabel(r'$\Delta v$ (km/s)')
        P.tight_layout()
        if save:
            P.savefig(Repeatability.directory+'/plots/Repeat_all_rchi2_vel.pdf',\
                         bbox_inches='tight')

        ngals = len(dv)*1.
        print 'Total galaxies in plot', ngals
        print '   dchi2 < 0.01', N.sum( dc< 0.01),  N.sum( dc< 0.01)/ngals
        print '   dchi2 < 0.005', N.sum( dc< 0.005),  N.sum( dc< 0.005)/ngals
        print '   dv > 1000 km/s', N.sum( dv > 1000.), N.sum( dv > 1000.)/ngals

        return data



