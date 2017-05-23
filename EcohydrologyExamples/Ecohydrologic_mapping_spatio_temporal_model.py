
# coding: utf-8

# <a href="http://landlab.github.io"><img style="float: left" src="https://raw.githubusercontent.com/landlab/tutorials/master/landlab_header.png"></a>

# # Modeling ecohydrological response to a storm event

# This Landlab driver illustrates the use of Landlab ecohydrology components to model semi-arid ecohydrological dynamics driven by a storm pulse and solar radiation. Components (names given in parenthesis) we will use are:
# * Solar radiation (Radiation)
# * Potential Evapotranspiration (PotentialEvapotranspiration)
# * Soil Moisture (SoilMoisture)
# * Vegetation (Vegetation)
# A digital elevation model (DEM) of a headwater region in central New Mexico (latitude 34N) will be used as input. 
# Components will be introduced step by step. First, we will start with mapping solar radiation and potential evapotranspiration (PET). Note that some of the commands used are only to provide information about the in/outputs of components and can be deleted or not run. We will then run soil moisture and vegetation modules and show how to write outputs in a file.
# 
# Let’s being importing python libraries, Landlab components, and python plotting tools 

# In[ ]:

import numpy as np
from landlab.components import (PrecipitationDistribution,
                                Radiation, PotentialEvapotranspiration,
                                SoilMoisture, Vegetation)
from landlab.io import read_esri_ascii
from landlab.plot import imshow_grid
get_ipython().magic(u'matplotlib inline')
from landlab.plot.imshow import imshow_grid
import matplotlib.pyplot as plt
plt.show()


# Read an existing esri grid as watershed and map the elevation field

(watershed, z)=read_esri_ascii('DEM_10m.asc', name='topographic__elevation')
imshow_grid(watershed, 'topographic__elevation')
plt.show()

# Inputs for spatio-temporal Ecohydrology Model
n = 150   # Number of storms
current_time = 0            # Initial time in years
doy__start_of_monsoon = 182  # Record the start day of monsoon (Julian day)
doy__end_of_monsoon = 273    # Record the ending day of monsoon (Julian day)
vegetation_type=0                 # grass=0, shrub=1, tree=2, bare=3
vegetation_cover=0.5;             # initial vegetation ground cover fraction (0-1) keep this fixed, updated when the vegetation model is used
initial_soil_saturation=0.3;      # initial soil saturation in the root zone (0-1)
initial_live_biomass=30.        # Initial live biomass in g DM m^-2
initial_dead_biomass=30.        # Initial dead biomass in g DM m^-2

# Initialize fields for SoilMoisture and Vegetation components
watershed.at_cell['vegetation__plant_functional_type']= (
    vegetation_type * np.ones(watershed.number_of_cells, dtype=int))
watershed.at_cell['soil_moisture__initial_saturation_fraction']=(
    initial_soil_saturation*np.ones(watershed.number_of_cells))
watershed.at_cell['vegetation__cover_fraction']=(
    vegetation_cover*np.ones(watershed.number_of_cells))

# Create empty arrays to hold spatio-temporal data
precip = np.empty([n])
storm_dt = np.empty([n])
inter_storm_dt = np.empty([n])
sm_saturation_fraction = np.empty([n, watershed.number_of_cells])
veg_live_leaf_area_index = np.empty([n, watershed.number_of_cells])
veg_live_biomass = np.empty([n, watershed.number_of_cells])

# Instantiate classes
# We shall create storms with different characteristics
# for dry season and monsoon seasoon
# Instantiating PrecipitationDistribution component for dry season
precip_dry = PrecipitationDistribution(mean_storm_duration=2.016,
                                       mean_interstorm_duration=159.36,
                                       mean_storm_depth=3.07)
# Instantiating PrecipitationDistribution component for monsoon season
precip_wet = PrecipitationDistribution(mean_storm_duration=1.896,
                                       mean_interstorm_duration=84.24,
                                       mean_storm_depth=4.79)
# Instantiate radiation component
rad = Radiation(watershed, method='Grid', latitude=34.)
# Instantiate potential evapotranspiration component and
# use 'Cosine' method
pet_grass = PotentialEvapotranspiration(watershed, method='Cosine',
                                        MeanTmaxF=4.96, delta_d=7.)
soil_moisture = SoilMoisture(watershed, runon=1, f_bare=0.7)
vegetation = Vegetation(watershed, Blive_init=initial_live_biomass,
                        Bdead_init=initial_dead_biomass)

# Run time loop - storm loop
for i in range(0, n):
    # Calculate the Julian day to identify the season
    julian = np.int(np.floor((current_time - np.floor(current_time)) * 365.))
    # Generate seasonal storms
    # for Dry season
    if doy__start_of_monsoon <= julian <= doy__end_of_monsoon:
        precip_wet.update()
        precip[i] = precip_wet.get_storm_depth()
        storm_dt[i] = precip_wet.get_precipitation_event_duration()
        inter_storm_dt[i] = precip_wet.get_interstorm_event_duration()
    # Wet Season—Jul to Sep—NA Monsoon
    else:
        precip_dry.update()
        precip[i] = precip_dry.get_storm_depth()
        storm_dt[i] = precip_dry.get_precipitation_event_duration()
        inter_storm_dt[i] = precip_dry.get_interstorm_event_duration()
    watershed.at_cell['rainfall__daily_depth'] = (precip[i] *
        np.ones(watershed.number_of_cells))
    rad.update(current_time)
    pet_grass.update(current_time)
    watershed.at_cell['surface__potential_evapotranspiration_30day_mean'] = (
        watershed['cell']['surface__potential_evapotranspiration_rate'])
    current_time = soil_moisture.update(current_time, Tb=inter_storm_dt[i],
                                        Tr=storm_dt[i])
    sm_saturation_fraction[i] = (
        watershed.at_cell['soil_moisture__saturation_fraction'])
    vegetation.update(Tb=inter_storm_dt[i], Tr=storm_dt[i])
    veg_live_leaf_area_index[i] = (
        watershed.at_cell['vegetation__live_leaf_area_index'])
    veg_live_biomass[i] = (
        watershed.at_cell['vegetation__live_biomass'])


# # Plot figures

# Lets find two random cells, one on north facing slope and one on south facing cells to plot model outputs for.

# Calculate slopes and aspects - and categorize
slope, aspect = watershed.calculate_slope_aspect_at_nodes_burrough(vals='topographic__elevation')
slope = np.degrees(slope)
aspect = np.degrees(aspect)
aspect = np.remainder(aspect, 360.)
# Aspect categorization
north_slopes = np.where(np.logical_or(aspect>=315., aspect<45.) == True)[0]
south_slopes = np.where(np.logical_and(aspect>=135, aspect<225) == True)[0]

# Check how many cells are north-facing
north_slopes.shape

# Select a random north-facing cell
random_north_facing_cell = north_slopes[500]


# Check how many cells are south-facing
south_slopes.shape

# Select a random south-facing cell
random_south_facing_cell = south_slopes[500]

# Plot Soil Moisture Saturation Fraction (Theta)
plt.figure()
plt.plot(sm_saturation_fraction[:, random_north_facing_cell], 'b',
         label='North Facing Cell')
plt.plot(sm_saturation_fraction[:, random_south_facing_cell], 'r',
         label='South Facing Cell')
plt.xlabel('Storm Number')
plt.ylabel('Soil Moisture Saturation Fraction (Theta)')
plt.xlim(xmin=0, xmax=n)
plt.ylim(ymin=0, ymax=1)
plt.legend(loc=1)
plt.show()


# Now lets plot the figures!

# Plot live leaf area index
plt.figure()
plt.plot(veg_live_leaf_area_index[:, random_north_facing_cell], 'b',
         label='North Facing Cell')
plt.plot(veg_live_leaf_area_index[:, random_south_facing_cell], 'r',
         label='South Facing Cell')
plt.xlabel('Storm Number')
plt.ylabel('Live Leaf Area Index')
plt.xlim(xmin=0, xmax=n)
plt.ylim(ymin=0, ymax=1.5)
plt.legend(loc=1)
plt.show()

# Plot live biomass (g m^-2 d^-1)
plt.figure()
plt.plot(veg_live_biomass[:, random_north_facing_cell], 'b',
         label='North Facing Cell')
plt.plot(veg_live_biomass[:, random_south_facing_cell], 'r',
         label='South Facing Cell')
plt.xlabel('Storm Number')
plt.ylabel('Live Biomass')
plt.xlim(xmin=0, xmax=n)
plt.ylim(ymin=0)
plt.legend(loc=1)
plt.show()

