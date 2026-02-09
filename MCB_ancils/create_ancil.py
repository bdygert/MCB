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

# R5 (South-East Atlantic - straddles the Greenwich meridian):
R5_W = [335, -30, 359.5, 0]
R5_E = [0.5, -30, 15, 0]

# Northern Ocean ("R6") has been superseded by NOx ("R15").

# R7 (Western North Pacific):
R7 = [140, 30, 210, 50]

# R8 (North-West Pacific)
R8 = [120, 0, 160, 30]

# R9 (South-West Pacific)
R9 = [150, -30, 190, 0]

# R10 (South Atlantic - straddles the Greenwich meridian):
R10_W = [305, -50, 359.5, -30]
R10_E = [0.5, -50, 15, -30]

# R11 (North-West Atlantic)
R11 = [290, 0, 335, 30]

# R12 (North Atlantic)
R12 = [290, 30, 360, 50]

# R13 (South-East Indian Ocean)
R13 = [70, -30, 110, 0]

# R14 (Northern Indian Ocean)
R14 = [45, 0, 100, 30]

# R15 (Northern Oceans - extended):
R15 = [0, 50, 359.5, 80]



# Define output filename:
# ======================
opfile = dir + 'test_AM04_seasalt_ems.pp'


# Define sea-salt injection rates for each region
# ===============================================
R1_ss_rate =  0.000             # Tg year-1    
R2_ss_rate =  0.000 #11.000
R3_ss_rate =  0.000
R4_ss_rate =  0.000
R5_ss_rate =  0.000
R7_ss_rate =  0.000
R8_ss_rate =  0.000 
R9_ss_rate =  0.000
R10_ss_rate =  0.000
R11_ss_rate =  0.000
R12_ss_rate =  0.000
R13_ss_rate =  0.000
R14_ss_rate =  0.000
R15_ss_rate =  11.000 #0.000


#-----------------------------------------------------------------------------------------
# Read in a field of surface temperature to use as a pattern
# to construct the emissions data:
#infile = '/data/users/alex.mason/MCB_50Tg/cq295/temp_at_1.5m_mm/cq295a.p52049sep.pp'
#ss_emiss = iris.load_cube(infile, iris.AttributeConstraint(STASH='m01s52i024'))
#infile = '/data/users/hadna/Controller/MCB/mod_seasalt_test/cp109a.pm_2040_jan_to_dec_00024.pp'
infile = dir + 'cp109a.pm_2040_jan_to_dec_00024.pp'
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


# Get land-fraction field:
landfile = dir +'aw310a.land_fraction.pp'
land_frac = iris.load_cube(landfile)


# Get indices of points in regions:
R1_lat_indices = np.where((lats.points >= R1[1]) & (lats.points <= R1[3]))[0]
R1_lon_indices = np.where((lons.points >= R1[0]) & (lons.points <= R1[2]))[0]

R2_lat_indices = np.where((lats.points >= R2[1]) & (lats.points <= R2[3]))[0]
R2_lon_indices = np.where((lons.points >= R2[0]) & (lons.points <= R2[2]))[0]

R3_lat_indices = np.where((lats.points >= R3[1]) & (lats.points <= R3[3]))[0]
R3_lon_indices = np.where((lons.points >= R3[0]) & (lons.points <= R3[2]))[0]

R4_lat_indices = np.where((lats.points >= R4[1]) & (lats.points <= R4[3]))[0]
R4_lon_indices = np.where((lons.points >= R4[0]) & (lons.points <= R4[2]))[0]

R5_lat_indices = np.where( ((lats.points >= R5_W[1]) & (lats.points <= R5_W[3])) \
                          |((lats.points >= R5_E[1]) & (lats.points <= R5_E[3])) )[0]
R5_lon_indices = np.where( ((lons.points >= R5_W[0]) & (lons.points <= R5_W[2])) \
                          |((lons.points >= R5_E[0]) & (lons.points <= R5_E[2])) )[0]

R7_lat_indices = np.where((lats.points >= R7[1]) & (lats.points <= R7[3]))[0]
R7_lon_indices = np.where((lons.points >= R7[0]) & (lons.points <= R7[2]))[0]

R8_lat_indices = np.where((lats.points >= R8[1]) & (lats.points <= R8[3]))[0]
R8_lon_indices = np.where((lons.points >= R8[0]) & (lons.points <= R8[2]))[0]

R9_lat_indices = np.where((lats.points >= R9[1]) & (lats.points <= R9[3]))[0]
R9_lon_indices = np.where((lons.points >= R9[0]) & (lons.points <= R9[2]))[0]

R10_lat_indices = np.where( ((lats.points >= R10_W[1]) & (lats.points <= R10_W[3])) \
                           |((lats.points >= R10_E[1]) & (lats.points <= R10_E[3])) )[0]
R10_lon_indices = np.where( ((lons.points >= R10_W[0]) & (lons.points <= R10_W[2])) \
                           |((lons.points >= R10_E[0]) & (lons.points <= R10_E[2])) )[0]

R11_lat_indices = np.where((lats.points >= R11[1]) & (lats.points <= R11[3]))[0]
R11_lon_indices = np.where((lons.points >= R11[0]) & (lons.points <= R11[2]))[0]

R12_lat_indices = np.where((lats.points >= R12[1]) & (lats.points <= R12[3]))[0]
R12_lon_indices = np.where((lons.points >= R12[0]) & (lons.points <= R12[2]))[0]

R13_lat_indices = np.where((lats.points >= R13[1]) & (lats.points <= R13[3]))[0]
R13_lon_indices = np.where((lons.points >= R13[0]) & (lons.points <= R13[2]))[0]

R14_lat_indices = np.where((lats.points >= R14[1]) & (lats.points <= R14[3]))[0]
R14_lon_indices = np.where((lons.points >= R14[0]) & (lons.points <= R14[2]))[0]

R15_lat_indices = np.where((lats.points >= R15[1]) & (lats.points <= R15[3]))[0]
R15_lon_indices = np.where((lons.points >= R15[0]) & (lons.points <= R15[2]))[0]


# Get area of open ocean within specified regions:
R1_ocean_area = 0.0
for i in R1_lat_indices:
   for j in R1_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R1_ocean_area = R1_ocean_area + grid_areas[i, j]

R2_ocean_area = 0.0
for i in R2_lat_indices:
   for j in R2_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R2_ocean_area = R2_ocean_area + grid_areas[i, j]

R3_ocean_area = 0.0
for i in R3_lat_indices:
   for j in R3_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R3_ocean_area = R3_ocean_area + grid_areas[i, j]

R4_ocean_area = 0.0
for i in R4_lat_indices:
   for j in R4_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R4_ocean_area = R4_ocean_area + grid_areas[i, j]

R5_ocean_area = 0.0
for i in R5_lat_indices:
   for j in R5_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R5_ocean_area = R5_ocean_area + grid_areas[i, j]

R7_ocean_area = 0.0
for i in R7_lat_indices:
   for j in R7_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R7_ocean_area = R7_ocean_area + grid_areas[i, j]

R8_ocean_area = 0.0
for i in R8_lat_indices:
   for j in R8_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R8_ocean_area = R8_ocean_area + grid_areas[i, j]

R9_ocean_area = 0.0
for i in R9_lat_indices:
   for j in R9_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R9_ocean_area = R9_ocean_area + grid_areas[i, j]

R10_ocean_area = 0.0
for i in R10_lat_indices:
   for j in R10_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R10_ocean_area = R10_ocean_area + grid_areas[i, j]

R11_ocean_area = 0.0
for i in R11_lat_indices:
   for j in R11_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R11_ocean_area = R11_ocean_area + grid_areas[i, j]

R12_ocean_area = 0.0
for i in R12_lat_indices:
   for j in R12_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R12_ocean_area = R12_ocean_area + grid_areas[i, j]

R13_ocean_area = 0.0
for i in R13_lat_indices:
   for j in R13_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R13_ocean_area = R13_ocean_area + grid_areas[i, j]

R14_ocean_area = 0.0
for i in R14_lat_indices:
   for j in R14_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R14_ocean_area = R14_ocean_area + grid_areas[i, j]

R15_ocean_area = 0.0
for i in R15_lat_indices:
   for j in R15_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R15_ocean_area = R15_ocean_area + grid_areas[i, j]

# The following calculates the combined area of the R3 & R7 regions 
# allowing for their overlap; it's not needed for calculating emissions,
# it's just for adding up the total injection area. The longitude indices
# use R3 for W and R7 for E; either can be used for N/S.

R3_R7_lat_indices = np.where((lats.points >= R3[1]) & (lats.points <= R3[3]))[0]
R3_R7_lon_indices = np.where((lons.points >= R3[0]) & (lons.points <= R7[2]))[0]

R3_R7_ocean_area = 0.0
for i in R3_R7_lat_indices:
   for j in R3_R7_lon_indices:
      if land_frac.data[i, j] == 0.0:
         R3_R7_ocean_area = R3_R7_ocean_area + grid_areas[i, j]




# Insert the required injection amount in the appropriate areas:
# NOTE that region R7 (Western North Pacific) partially overlaps
#      region R3 (North Pacific) so, as R7 comes after R3, there's
#      a slightly different approach for R7:

for i in R1_lat_indices:
   for j in R1_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R1_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R1_ocean_area

for i in R2_lat_indices:
   for j in R2_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R2_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R2_ocean_area

for i in R3_lat_indices:
   for j in R3_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R3_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R3_ocean_area

for i in R4_lat_indices:
   for j in R4_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R4_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R4_ocean_area

for i in R5_lat_indices:
   for j in R5_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R5_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R5_ocean_area

for i in R7_lat_indices:
   for j in R7_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = ss_emiss.data[:, i, j] + ((R7_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R7_ocean_area)

for i in R8_lat_indices:
   for j in R8_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R8_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R8_ocean_area

for i in R9_lat_indices:
   for j in R9_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R9_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R9_ocean_area

for i in R10_lat_indices:
   for j in R10_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R10_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R10_ocean_area

for i in R11_lat_indices:
   for j in R11_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R11_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R11_ocean_area

for i in R12_lat_indices:
   for j in R12_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R12_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R12_ocean_area

for i in R13_lat_indices:
   for j in R13_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R13_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R13_ocean_area

for i in R14_lat_indices:
   for j in R14_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R14_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R14_ocean_area

for i in R15_lat_indices:
   for j in R15_lon_indices:
      if land_frac.data[i, j] == 0.0:
         ss_emiss.data[:, i, j] = (R15_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R15_ocean_area


# Calculate global total using data from a single month (they're all the same):
ss_emiss_total = ss_emiss[0].collapsed(['longitude','latitude'],iris.analysis.SUM, weights=grid_areas).data


# Convert from kg/s to Tg/yr for printout purposes:
ss_emiss_Tg_yr = ss_emiss_total *60.*60.*24.*360.*1e-9


# Print out total area used for injection:
tot_area = R1_ocean_area + R2_ocean_area + R3_R7_ocean_area + R4_ocean_area + R5_ocean_area \
         + R8_ocean_area + R9_ocean_area + R10_ocean_area + R11_ocean_area + R12_ocean_area \
         + R13_ocean_area + R14_ocean_area + R15_ocean_area
print('\nTotal area used for injection (million km2) = '+str(tot_area * 1e-12))


# Print out individual emission rates (per unit area):
print('  ')
print('Emissions for R1 (NEP), kg m-2 sec-1: '+str((R1_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R1_ocean_area))
print('Emissions for R2 (SEP), kg m-2 sec-1: '+str((R2_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R2_ocean_area))
print('Emissions for R3 (NP), kg m-2 sec-1 : '+str((R3_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R3_ocean_area))
print('Emissions for R4 (SP), kg m-2 sec-1 : '+str((R4_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R4_ocean_area))
print('Emissions for R5 (SEA), kg m-2 sec-1 : '+str((R5_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R5_ocean_area))
print('Emissions for R7 (WNP), kg m-2 sec-1 : '+str((R7_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R7_ocean_area))
print('Emissions for R8 (NWP), kg m-2 sec-1 : '+str((R8_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R8_ocean_area))
print('Emissions for R9 (SWP), kg m-2 sec-1 : '+str((R9_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R9_ocean_area))
print('Emissions for R10 (SA), kg m-2 sec-1 : '+str((R10_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R10_ocean_area))
print('Emissions for R11 (NWA), kg m-2 sec-1 : '+str((R11_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R11_ocean_area))
print('Emissions for R12 (NA), kg m-2 sec-1 : '+str((R12_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R12_ocean_area))
print('Emissions for R13 (SEI), kg m-2 sec-1 : '+str((R13_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R13_ocean_area))
print('Emissions for R14 (NI), kg m-2 sec-1 : '+str((R14_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R14_ocean_area))
print('Emissions for R15 (NOx), kg m-2 sec-1 : '+str((R15_ss_rate * 1.0e9 / (60.*60.*24.*360.)) / R15_ocean_area))
# Print out total:
print('  ')
print('Total emissions (Tg[sea-salt]/yr) = '+str(ss_emiss_Tg_yr))
print('  ')

# Save emissions to the specified output file (if desired):
#iris.save(ss_emiss, opfile)
