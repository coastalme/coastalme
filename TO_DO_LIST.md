# TO DO LIST

068 is max (also using github permalinks to code lines)

## Table of contents
- [User input](#user-input)
- [Error handling](#error-handling)
- [Theory/Efficiency](#theory-efficiency)
- [Output](#output)
- [Completed](#completed)



## User input
- TODO 000 [Should user input be split in two main files]
- TODO 011 Should this constant be a user input? 
		[1](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/cme.h#L600-L602)
		[2](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/simulation.cpp#L211)
		[3](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/simulation.cpp#L371)
- TODO 036 [Read in changed deep water wave values](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/init_grid.cpp#L177) 
- TODO 030 [Do we also need to be able to input landform sub-categories?] [code](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/gis_raster.cpp#L750)
- TODO 027 [Sort out GDAL problem with raster reference units](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/gis_raster.cpp#L537)
- TODO 022 [Get intervention update working](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_intervention.cpp#L32)
- TODO 041 [Read in SWL per-timestep](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/read_input.cpp#L1977)
- TODO 042 [Should we have a smallest valid input for KLS in the CERC equation?] [code](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/read_input.cpp#L2378-L2380)
- TODO 045 [Method of getting depth of closure value needs to be a user input] [code](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/simulation.h#L771)
- TODO 047 [Where is the GDAL description for the deep water wave stations vector file?](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/simulation.h#L1198)
- TODO 048 [Where is the GDAL description for the flood input locations point or vector file?](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/simulation.h#L1220)
- TODO 049 [Handle other command line parameters e.g. path to .ini file, path to datafile](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/utils.cpp#L125)
- TODO 035 [Also handle other EPSG for vector spatial reference systems](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/gis_vector.cpp#L492)
- TODO 054 [Choose more files to omit from "usual" raster output]

##   Error handling
-  TODO 004 [Improve error handling of situation where we have a valid shadow zone but cannot find a neighbouring cell which is 'under' the coastline](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/calc_shadow_zones.cpp#L490)
-  TODO 006 [Check GDALGridCreate() with only start-of-coast or an end-of-coast profiles](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/calc_waves.cpp#L214)
-  TODO 009 [Decide what to do when we have eroded down to basement](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/calc_waves.cpp#L214)
-  TODO 015 [Improve situation where profile hits another profile which belongs to a different coast object (will certainly need this for estuaries)](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/create_profiles.cpp#L1478) [2](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/create_profiles.cpp#L1537)
-  TODO 017 [Extra safety check needed, make sure that each point is within valid grid](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_beach_within_polygon.cpp#L146)
-  TODO 018 [Improve situation where new landwards point on parallel profile is not within the raster grid](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_beach_within_polygon.cpp#L222) [2](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_beach_within_polygon.cpp#L394)
-  TODO 019 [Improve situation where Dean profile has a near-zero elevation difference](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_beach_within_polygon.cpp#L286) [2](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_beach_within_polygon.cpp#L394)
-  TODO 020 [Check calculation of elevation of coast point of Dean parallel profile](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_beach_within_polygon.cpp#L949) [2](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_beach_within_polygon.cpp#L1435)
-  TODO 021 [Improve situation where all layers have zero thickness](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_cliff_collapse.cpp#L819)
-  TODO 025 [Improve situation where this point has only zero thickness layers](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_shore_platform_erosion.cpp#L195) [2](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_shore_platform_erosion.cpp#L508)
-  TODO 026 [Check situation where cell in parallel profile is not in a polygon](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_shore_platform_erosion.cpp#L495)
-  TODO 028 [Give a warning if raster input layer has several bands](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/gis_raster.cpp#L605)
-  TODO 038 [Do better error handling if insufficient memory](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/gis_raster.cpp#L71)
-  TODO 053 [Improve handling of situation where landward elevation of profile is -ve](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/calc_waves.cpp#L1581)
-  TODO 055 [Maybe add a safety check?]

##   Theory/Efficiency
-  TODO 002 [Do we really need D50 for drift landform class? What do we need for drift?](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/assign_landforms.cpp#L346)
-  TODO 005 [Maybe give every coast point a value for end-of-profile wave height and direction instead of for deep water wave height and direction](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/calc_waves.cpp#L58)
-  TODO 010 [Do we also need to update the active zone cells?](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/calc_waves.cpp#L1828)
-  TODO 012 [Change finding of adjacent polygons, and calculation of the length of shared normals, when we make polygon seaward length determined by depth of closure](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/create_polygons.cpp#L504)
-  TODO 013 [Change calculation (need user input?) of coastline smoothing convexity threshold](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/create_profiles.cpp#L119)
-  TODO 014 [Profile spacing, vould try gradually increasing the profile spacing with increasing concavity, and decreasing the profile spacing with increasing convexity](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/create_profiles.cpp#L396)
-  TODO 016 [Check mass balance for recirculating unconsolidated sediment option](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_beach_sediment_movement.cpp#L539)
-  TODO 023 [Only calculate shore platform erosion if cell is in a polygon](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_shore_platform_erosion.cpp#L51)
-  TODO 024 [Should we calculate platform erosion on a profile that has hit dry land?](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/do_shore_platform_erosion.cpp#L124)
-  TODO 037 [Need more info on nFindIndex](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/interpolate.cpp#L103-L121)
-  TODO 044 [Estuaries :-)](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/simulation.cpp#L933)
-  TODO 050 Update for recent versions of Windows
-  TODO 051 [Implement other ways of calculating depth of closure, see TODO 045](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/utils.cpp#L2482-L2493)
-  TODO 056 [Check this please Andres]
-  TODO 057 [Check this please Manuel]
-  TODO 059 [Implement dune landform class]
-  TODO 060 [Remove 'magic numbers' from code here]
-  TODO 061 [Is this safety check to depth of breaking a reasonable thing to do?]
-  TODO 066 [Should this be for all layers? Check]
-  TODO 067 [Is this ever non-zero? Check]

##  Output
-  TODO 065 [Get GPKG output working] 
-  TODO 063 Add NetCDF support, see https://trac.osgeo.org/gdal/wiki/NetCDF
-  TODO 064 Add support for grids that are not oriented N-S and W-E, but are still rectangular (will need to add a transformation in the reading and writing process, the first to bring it to the local base and the second to save it in global coordinates)
-  TODO 031 [Get raster slice output working with multiple slices](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/gis_raster.cpp#L1122-L1125)
-  TODO 032 [Improve output scaling for DBL_NODATA situation](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/gis_raster.cpp#L1617)
-  TODO 033 [Also test and configure (e.g. by passing open() options) other vector output file formats](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/gis_utils.cpp#L893)
-  TODO 034 [Also test and configure (e.g. by passing open() options) other raster output file formats](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/gis_utils.cpp#L1642)
-  TODO 043 [When outputting profiles, how do we deal with randomness of profile spacing (since profile location is determined by curvature)?]
-  TODO 052 [Improve saving of profiles and parallel profiles]
-  TODO 062 [Show end-of-iteration number of cells with sediment somewhere]
-  TODO 068 Only show output in log file that is relevant to processes being simulated
   
## Completed
-  TODO 003 Make coastline curvature moving window size a user input *** FIXED in 1.1.22
-  TODO 046 Why is cliff collapse eroded during deposition (three size classes) no longer calculated? *** FIXED in 1.1.22
-  TODO 058 Dave to check this *** FIXED in 1.1.22
-  TODO 039 [Rewrite reading of multiple random number seeds](https://github.com/coastalme/coastalme/blob/730a0be2274de02125681ab3e14424dbbfadbaeb/src/read_input.cpp#L535) FIXED in 1.2.1 (8 Nov 2024)

   