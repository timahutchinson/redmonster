import numpy as n
from redmonster.sandbox import yanny as y

# Read yanny file
x = y.yanny(filename='/uufs/astro.utah.edu/common/home/u0814744/boss/spInspect_alltest_bolton.par.txt', np=True)

# Get fibers, zpipe, zperson for each plate
args = n.where(x['BOSSOBJECT']['plate'] == 3686)[0]
fibers3686 = []
zpipe3686 = []
zperson3686 = []
for i in args:
    fibers3686.append( x['BOSSOBJECT'][i][2])
    zpipe3686.append( x['BOSSOBJECT'][i][5])
    zperson3686.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3687)[0]
fibers3687 = []
zpipe3687 = []
zperson3687 = []
for i in args:
    fibers3687.append( x['BOSSOBJECT'][i][2])
    zpipe3687.append( x['BOSSOBJECT'][i][5])
    zperson3687.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3804)[0]
fibers3804 = []
zpipe3804 = []
zperson3804 = []
for i in args:
    fibers3804.append( x['BOSSOBJECT'][i][2])
    zpipe3804.append( x['BOSSOBJECT'][i][5])
    zperson3804.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3805)[0]
fibers3805 = []
zpipe3805 = []
zperson3805 = []
for i in args:
    fibers3805.append( x['BOSSOBJECT'][i][2])
    zpipe3805.append( x['BOSSOBJECT'][i][5])
    zperson3805.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3853)[0]
fibers3853 = []
zpipe3853 = []
zperson3853 = []
for i in args:
    fibers3853.append( x['BOSSOBJECT'][i][2])
    zpipe3853.append( x['BOSSOBJECT'][i][5])
    zperson3853.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3855)[0]
fibers3855 = []
zpipe3855 = []
zperson3855 = []
for i in args:
    fibers3855.append( x['BOSSOBJECT'][i][2])
    zpipe3855.append( x['BOSSOBJECT'][i][5])
    zperson3855.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3856)[0]
fibers3856 = []
zpipe3856 = []
zperson3856 = []
for i in args:
    fibers3856.append( x['BOSSOBJECT'][i][2])
    zpipe3856.append( x['BOSSOBJECT'][i][5])
    zperson3856.append( x['BOSSOBJECT'][i][6])

args = n.where(x['BOSSOBJECT']['plate'] == 3860)[0]
fibers3860 = []
zpipe3860 = []
zperson3860 = []
for i in args:
    fibers3860.append( x['BOSSOBJECT'][i][2])
    zpipe3860.append( x['BOSSOBJECT'][i][5])
    zperson3860.append( x['BOSSOBJECT'][i][6])