#Modifying code for future laziness purposes, Brittany Dygert 3 February 2026
#On JASMIN, 
#1. "module load ants" in command line
#2. in this file, change opfile, rsel & totalems for region and amount Tg/yr
#3. python3 bdd_ancil.py
#4. "ancil_2anc.py --output <ancillary_file_name.anc> --grid-staggering 6 <pp_file_name.pp>" in command line
#5. From local terminal, scp [jasmin server]:'[path to file]' [monsoon server/etc]:'[path to folder]/.'

#bdd: region select (scroll down to R1-R15 for defined regions)
#Put in as list in case you want multiple regions
#rsel=[16]
rsel=list(range(1, 17)) #for activating all regions
#bdd: Total desired emissions in Tg/yr
totalems=50;
numreg=len(rsel) #number of regions selected
# Define output filename:
# ======================
opfile = 'test2.pp'

import cf_units
import iris
import iris.analysis.cartography as iac
import numpy as np
import warnings

# Set values for output directory ("dir", ending in a slash) and 
# output filename ("opfile", ending ".pp"). Then set the emission
# amounts for each region (Tg/yr). Once you're happy with things,
# un-comment the final line and it will save the file (as a 
# pp-file) to your specified location. It can then be converted
# to the required ancillary-file format as follows (on VDI):
# 
# => module load ants
# => ancil_2anc.py --output <ancillary_file_name.anc> --grid-staggering 6 <pp_file_name.pp>
# 
# As well as giving you the ancillary file with the .anc suffix it 
# will also give you a netCDF file (.anc.nc) which you can use to
# plot the emissions field if you want to look at it.
#  
# The ancillary file can then be transferred to the HPC using the
# scp command.


warnings.filterwarnings("ignore", category=UserWarning, message="Collapsing a non-contiguous coordinate.")
warnings.filterwarnings("ignore", category=UserWarning, message="Unable to create instance of HybridHeightFactory.")
warnings.filterwarnings("ignore", category=UserWarning, message="has_year_zero kwarg ignored for idealized calendars")

dir = '~/MCB_ancils/'



# Define region limits [W, S, E, N]
# =================================
# All longitudes are specified in degrees East.

# R1 (North-East Pacfic): 
R1 = [210, 0, 250, 30]

# R2 (South-East Pacific): 
R2 = [250, -30, 290, 0]

# R3 (North Pacific): 
R3 = [170, 30, 240, 50]

# R4 (South Pacific): 
R4 = [190, -50, 270, -30]

# R5-6 (South-East Atlantic - straddles the Greenwich meridian):
R5 = [335, -30, 359.5, 0]
R6 = [0.5, -30, 15, 0]

# R7 (Western North Pacific):
R7 = [140, 30, 210, 50]

# R8 (North-West Pacific)
R8 = [120, 0, 160, 30]

# R9 (South-West Pacific)
R9 = [150, -30, 190, 0]

# R10-11 (South Atlantic - straddles the Greenwich meridian):
R10 = [305, -50, 359.5, -30]
R11 = [0.5, -50, 15, -30]

# R12 (North-West Atlantic)
R12 = [290, 0, 335, 30]

# R13 (North Atlantic)
R13 = [290, 30, 360, 50]

# R14 (South-East Indian Ocean)
R14 = [70, -30, 110, 0]

# R15 (Northern Indian Ocean)
R15 = [45, 0, 100, 30]

# R16 (Northern Oceans - extended):
R16 = [0, 50, 359.5, 80]

#Redo so you can just pick a region out of an array
#bdd
rs=[R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R16]
regsel=[rs[i-1] for i in rsel]
rates=np.zeros((numreg,1))
areas=np.zeros((numreg,1))
rlats=[]
rlons=[]

#-----------------------------------------------------------------------------------------
# Read in a field of surface temperature to use as a pattern
# to construct the emissions data:
#infile = '/data/users/alex.mason/MCB_50Tg/cq295/temp_at_1.5m_mm/cq295a.p52049sep.pp'
#ss_emiss = iris.load_cube(infile, iris.AttributeConstraint(STASH='m01s52i024'))
#infile = '/data/users/hadna/Controller/MCB/mod_seasalt_test/cp109a.pm_2040_jan_to_dec_00024.pp'
infile = 'cp109a.pm_2040_jan_to_dec_00024.pp'
ss_emiss = iris.load_cube(infile, iris.AttributeConstraint(STASH='m01s00i024'))


# Zero the data:
ss_emiss.data[:,:,:] = 0.0          # Dimensions are [time, lat, lon]


# Set STASH code to that for 2D user ancil:
ss_emiss.attributes['STASH'] = iris.fileformats.pp.STASH(1, 00, 301)


# Change its name:
ss_emiss.rename('Sea-salt emissions')


# Set the units:
ss_emiss.units = cf_units.Unit('kg m-2 s-1')


# Get lats and lons (as "DimCoord" thingies):
lats = ss_emiss.coord('latitude')
lons = ss_emiss.coord('longitude')


# Get gridbox areas (as 2-D numpy array):
cube = ss_emiss[0].copy()
if not cube.coord('latitude').has_bounds():
   cube.coord('latitude').guess_bounds()
if not cube.coord('longitude').has_bounds():
   cube.coord('longitude').guess_bounds()
grid_areas = iac.area_weights(cube)
#bdd, make array to check if a grid area has already been used
#alternate way to avoid overlap
#I would just zero out the grid_area but it would mess up validation later
overlap = np.ones(np.shape(grid_areas))


# Get land-fraction field:
landfile = 'aw310a.land_fraction.pp'
land_frac = iris.load_cube(landfile)


# Get indices of points in regions and associated areas:
#bdd
for k in range(0,numreg):
    rlats.append(np.where((lats.points >= regsel[k][1]) & (lats.points <= regsel[k][3]))[0])
    rlons.append(np.where((lons.points >= regsel[k][0]) & (lons.points <= regsel[0][2]))[0])
    for i in rlats[k]:
        for j in rlons[k]:
            if land_frac.data[i, j] == 0.0:
                areas[k] = areas[k] + ((grid_areas[i, j])*(overlap[i, j]))
                #Set to zero so that the next the index is hit, it doesn't contribute to area
                overlap[i, j] = 0.0
#bdd
#now that we have area, calculate emission rate
#I'm assuming we want equal emissions in all regions?
#if not just put in a diff array with specified emissions
#bdd
areatot= np.sum(areas[:])
rates=rates+(totalems*1.0e9/(60.*60.*24.*360.*areatot))
print(rates)
print(areatot)

#bdd loop back around to put in injection amounts
# Insert the required injection amount in the appropriate areas:
# NOTE that region R7 (Western North Pacific) partially overlaps
#      region R3 (North Pacific)
#bdd, ^ this is what I added overlap catch for.  My way is a little inefficient\# in that it still loops through every grid point and will just overwrite 
#a grid point it's already hit.
#this is ok because I set the code to make every region to have the same rate 
#(rate * area = desired total emission)
#probably a problem if you want different rates for different regions: 
#if problem: just use Andy Jones' original code create_ancil.py
for k in range(0,numreg):
    for i in rlats[k]:
        for j in rlons[k]:
            if land_frac.data[i, j] == 0.0:
                ss_emiss.data[:, i, j] = rates[k]
#end bdd

# Calculate global total using data from a single month (they're all the same):
ss_emiss_total = ss_emiss[0].collapsed(['longitude','latitude'],iris.analysis.SUM, weights=grid_areas).data


# Convert from kg/s to Tg/yr for printout purposes:
ss_emiss_Tg_yr = ss_emiss_total *60.*60.*24.*360.*1e-9


# Print out total area used for injection:
print('\nTotal area used for injection (million km2) = '+str(areatot * 1e-12))


# Print out individual emission rates (per unit area):
print('  ')
for k in range(0,numreg):
    print('Emissions for R'+str(rsel[k]),' (NEP), kg m-2 sec-1: '+str((rates[k] * 1.0e9 / (60.*60.*24.*360.))))
# Print out total:
print('  ')
print('Total emissions (Tg[sea-salt]/yr) = '+str(ss_emiss_Tg_yr))
print('  ')

# Save emissions to the specified output file (if desired):
iris.save(ss_emiss, opfile)
