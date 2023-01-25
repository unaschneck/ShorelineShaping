# ShorelineShaping
 Code accompaniment for Schneck et al. - The Shoreline Shaping Capability of Waves and Tides at Titan's Lakes

## Figures
Figure 1 - SAR Mosaics by [S. Birch and A. Hayes](https://geomorph-sbirch.com/data-products/)

Figure 2 - model_tidal_shields.m

Figure 3 - km_oscillatory_shields.m

Figure 4 - entrainment_depth_bathtub.m

Figure 5 - entrainment_depth_bathtub.m

Figure 6 - entrainment_depth_SAR.m

Figure 7 - entrainment_depth_bathtub.m

Figure 8 - entrainment_depth_SAR.m

Figure 9 - fw_oscillatory_shields.m

Figure S1-S4 - test_wave_models.m

Figure S6 - model_tidal_shields.m

Figure S7 – flow_w_depth.m

Figure S8 - make_bathtub_lakes.m

Figure S9 - [TADPOLE](https://github.com/MITGeomorph/Tadpole) with Ligeia Mare bathymetry

## Directory ReadMe
```
ShorelineShaping
|    |     |______________tides
|    |                    |____flow_w_depth.m : unidirectional flow with depth in nearshore and strait
|    |                    |____model_tidal_shields.m : Shield diagram for grain entrainment and suspension due to tidal flows at strait and nearshore
|    |                           |______make_shields_diagram.m : (function) makes Shields diagram
|    |                           |______tidal_shields.m : (function) finds particle Reynolds number and Shields parameter for unidirectional flow
|    |
|    |____________________waves
|                         |_____entrainment_depth_bathtub.m: finds entrainment depth for waves and makes contour plots
|                         |______fw_oscillatory_shields.m : finds Shields parameter for waves using wave friction factor
|                         |______km_oscillatory_shields.m : finds Shields parameter for waves using Komar & Miller threshold
|                         |______user_defined_functions
|                                |_______make_wave_COE.m : (function) grows waves towards maturity 
|                                |_______shoal_wave_COE : (function) shoals Airy waves towards breaking
|                                |_______test_wave_model.m : tests make_wave_COE and shoal_wave_COE functions
|
|_________________________resources
                          |____helper_Functions
                          |    |_______figure_settings.m : makes figures pretty
                          |
                          |_____lake_bathymetries
                          |     |______bathtub_bathy : lake bathymetry data
                          |     |      |______make_bathtub_lake.m : (function) makes bathtub model of lakes
                          |     |______SAR_bathy_cleaned : clean version of SAR lake bathymetry 
                          |     |      |______clean_lake_bathy.m : (function) cleak the SAR lake data and get perimeter
                          |     |______SAR_bathy_original : original data pulled from cub files from A. Hayes
                          |     |______SAR_bathy_smoothed : SAR lake bathymetry that has been cleaned and smoothed using TADPOLE
                          |
                          |
                          |_____max_fetch : maximum fetch for each point on the shoreline
                          |     |_____get_max_fetch_lakes : (function) gets maximum fetch for all points on shoreline using VisiLibity
                          |
                          |_____outside_data : terrestrial gravel-bedded river data from S. Birch
                          |     |___________ fetch_VisiLibity.m : (function) finds visible points at a point on shoreline (adapted from R. Palermo)
                          |     |___________ Pint.m : (function) helper function for fetch_VisiLibity (adapted from R. Palermo)
                          |     |___________ visilibity_polygon.m : (function) helper function for fetch_Visilibity (adapted from R. Palermo)
                          |
                          |_____shoreline_coord : coordinates for the perimeters of the lakes from ArcGIS
                                
 ```                               

