from redmonster.application import zfind

plates = [3686]
mjds = [55268]
fiberid = [100,101] #[i for i in xrange(1000)]
templates = ['ndArch-all-CAP-grids.fits', 'ndArch-ssp_em_galaxy-v000.fits',
             'ndArch-QSO-V003.fits']
zmin = [-.005,-.01,.4]
zmax = [.005, 1.2, 3.5]
npoly = [4,4,4]
npixstep = [1,2,4]
dest = '/Users/timhutchinson'


zf = zfind.Zfind(clobber=True)
zf.set_templates(templates=templates, zmin=zmin, zmax=zmax, npoly=npoly,
                 npixstep=npixstep)
#inifile= '/uufs/astro.utah.edu/common/home/u0814744/software/redmonster/\
        master/conf/zfind.ini'
#zf = zfind.Zfind(inifile=inifile, clobber=False)
for i in xrange(len(plates)):
    zf.reduce_plate_mjd(plate=plates[i], mjd=mjds[i], fiberid=fiberid)