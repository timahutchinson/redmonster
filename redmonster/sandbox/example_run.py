from redmonster.sandbox import main

plates = [3686]
mjds = [55268]
fiberid = [100,101,102]
templates = ['ndArch-all-CAP-grids.fits', 'ndArch-ssp_em_galaxy-v000.fits', 'ndArch-QSO-V003.fits']
zmin = [-.005,-.01,.4]
zmax = [.005, 1.2, 3.5]
npoly = [4,4,4]
npixstep = [1,4,4]
dest = '/Users/timhutchinson/compute/scratch'

for i in xrange(len(plates)):
    main.main(plate=plates[i], mjd=mjds[i], templates=templates, fiberid=fiberid, zmin=zmin, zmax=zmax, npoly=npoly, npixstep=npixstep, dest=dest)