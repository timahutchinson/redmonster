from redmonster.sandbox import main as M

plate = 3686
mjd = 55268
templates = ['ndArch-ssp_em_galaxy-v000.fits', 'ndArch-all-CAP-grids.fits', 'ndArch-QSO-V003.fits']
fiberid = [i for i in xrange(1000)]
zmin = [.01, -.005, .4]
zmax = [1.2, .005, 3.5]
npixstep = [1,1,4]
dest = '/Users/boltonlab3/scratch'

M.main(plate=plate, mjd=mjd, templates=templates, fiberid=fiberid, zmin=zmin, zmax=zmax, npixstep=npixstep)