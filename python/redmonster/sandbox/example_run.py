from redmonster.application import zfind

plates = [8123]
mjds = [56931]
fiberid = [i for i in xrange(1000)]
templates = ['ndArch-all-CAP-grids.fits', 'ndArch-ssp_em_galaxy-v000.fits', 'ndArch-QSO-V003.fits']
zmin = [-.005,-.01,.4]
zmax = [.005, 1.2, 3.5]
npoly = [4,4,4]
npixstep = [1,2,4]
dest = '/Users/boltonlab3/scratch'


for i in xrange(len(plates)):
    zfind.Zfind(plate=plates[i], mjd=mjds[i], templates=templates, fiberid=fiberid, zmin=zmin, zmax=zmax, npoly=npoly, npixstep=npixstep, clobber=False)