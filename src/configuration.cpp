/*!

   \file configuration.cpp
   \brief Implementation of unified configuration class for CoastalME
   \details Provides default values and initialization for simulation parameters
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License

*/

/* ==============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include "cme.h"
#include "configuration.h"
// #include "simulation.h"
#include <algorithm>

//===============================================================================================================================
//! Constructor
//===============================================================================================================================
CConfiguration::CConfiguration()
{
   InitializeDefaults();
}

//===============================================================================================================================
//! Destructor
//===============================================================================================================================
CConfiguration::~CConfiguration()
{
}

//===============================================================================================================================
//! Initialize all parameters with default values
//===============================================================================================================================
void CConfiguration::InitializeDefaults()
{
   // Run Information
   m_strRunName = "";
   m_nLogFileDetail = 1;
   m_bCSVPerTimestepResults = true;

   // Simulation timing
   m_strStartDateTime = "00-00-00 01/01/2000";
   m_strDuration = "365 days";
   m_strTimestep = "1 day";
   m_vecSaveTimes.clear();
   m_nRandomSeed = 0;
   m_bUseSystemTimeForSeed = true;

   // GIS Output
   m_nMaxSaveDigits = 3;
   m_strSaveDigitsMode = "sequential";
   m_vecRasterFiles.clear();
   m_vecRasterFiles.push_back("usual");
   m_strRasterFormat = "";
   m_bWorldFile = false;
   m_bScaleValues = false;
   m_vecSliceElevations.clear();
   m_vecVectorFiles.clear();
   m_vecVectorFiles.push_back("all");
   m_strVectorFormat = "ESRI Shapefile";
   m_vecTimeSeriesFiles.clear();
   m_vecTimeSeriesFiles.push_back("all");

   // Grid and Coastline
   m_nCoastlineSmoothing = 0;
   m_nCoastlineSmoothingWindow = 7;
   m_nPolynomialOrder = 4;
   m_strOmitGridEdges = "";
   m_nProfileSmoothingWindow = 0;
   m_dMaxLocalSlope = 1.0;
   m_dMaxBeachElevation = 10.0;

   // Layers and Files
   m_nNumLayers = 1;
   m_strBasementDEMFile = "";
   m_vecUnconsFineFiles.clear();
   m_vecUnconsSandFiles.clear();
   m_vecUnconsCoarseFiles.clear();
   m_vecConsFineFiles.clear();
   m_vecConsSandFiles.clear();
   m_vecConsCoarseFiles.clear();
   m_strSuspendedSedFile = "";
   m_strLandformFile = "";
   m_strInterventionClassFile = "";
   m_strInterventionHeightFile = "";

   // Hydrology
   m_nWavePropagationModel = 1;      // CShore
   m_dSeawaterDensity = 1029.0;
   m_dInitialWaterLevel = 0.0;
   m_dFinalWaterLevel = 0.0;
   m_bHasFinalWaterLevel = false;
   m_dDeepWaterWaveHeight = 1.0;
   m_strWaveHeightTimeSeries = "";
   m_dDeepWaterWaveOrientation = 270.0;
   m_dWavePeriod = 10.0;
   m_strTideDataFile = "";
   m_dBreakingWaveRatio = 0.8;

   // Sediment and Erosion
   m_bCoastPlatformErosion = true;
   m_dPlatformErosionResistance = 2e6;
   m_bBeachSedimentTransport = true;
   m_nBeachTransportAtEdges = 1;      // open
   m_nBeachErosionEquation = 0;       // CERC
   m_dFineMedianSize = 0.0;
   m_dSandMedianSize = 0.0;
   m_dCoarseMedianSize = 0.0;
   m_dSedimentDensity = 2650.0;
   m_dBeachSedimentPorosity = 0.4;
   m_dFineErosivity = 1.0;
   m_dSandErosivity = 0.7;
   m_dCoarseErosivity = 0.3;
   m_dTransportKLS = 0.4;
   m_dKamphuis = 5.0;
   m_dBermHeight = 0.25;

   // Cliff parameters
   m_bCliffCollapse = true;
   m_dCliffErosionResistance = 2.5e8;
   m_dNotchOverhang = 0.5;
   m_dNotchBase = 0.3;
   m_dCliffDepositionA = 0.0;
   m_dTalusWidth = 15.0;
   m_dMinTalusLength = 10.0;
   m_dMinTalusHeight = 0.5;

   // Flood parameters
   m_bFloodInput = false;
   m_strFloodCoastline = "";
   m_strRunupEquation = "";
   m_strFloodLocations = "";
   m_strFloodInputLocation = "";

   // Sediment input parameters
   m_bSedimentInput = false;
   m_strSedimentInputLocation = "";
   m_strSedimentInputType = "";
   m_strSedimentInputDetails = "";

   // Physics and Geometry
   m_dGravitationalAcceleration = 9.81;
   m_dNormalSpacing = 0.0;
   m_dRandomFactor = 0.25;
   m_dNormalLength = 130.0;
   m_dStartDepthRatio = 30.0;

   // Profile and Output Options
   m_bSaveProfileData = false;
   m_vecProfileNumbers.clear();
   m_vecProfileTimesteps.clear();
   m_bSaveParallelProfiles = false;
   m_bOutputErosionPotential = false;
   m_nCurvatureWindow = 11;

   // Cliff Edge Processing
   m_nCliffEdgeSmoothing = 1;
   m_nCliffEdgeSmoothingWindow = 33;
   m_nCliffEdgePolynomialOrder = 4;
   m_dCliffSlopeLimit = 0.3;
}

//===============================================================================================================================
//! Get raster files with keyword expansion support
//===============================================================================================================================
vector<string> CConfiguration::GetRasterFiles() const
{
   // Case 11: Raster GIS files to output - expand "all" and "usual" keywords
   vector<string> expandedFiles;

   for (string const &fileSpec : m_vecRasterFiles)
   {
      string fileSpecLower = fileSpec;
      std::transform(fileSpecLower.begin(), fileSpecLower.end(),
                     fileSpecLower.begin(), ::tolower);

      if (fileSpecLower == "all")
      {
         // Add all possible raster outputs (Case 11 "all" mode)
         expandedFiles.insert(expandedFiles.end(),
                              {"suspended_sediment",
                               "avg_suspended_sediment",
                               "fine_uncons",
                               "fine_cons",
                               "sand_uncons",
                               "sand_cons",
                               "coarse_uncons",
                               "coarse_cons",
                               "sediment_top_elevation",
                               "top_elevation",
                               "sea_depth",
                               "wave_height",
                               "wave_orientation",
                               "wave_period",
                               "potential_platform_erosion",
                               "actual_platform_erosion",
                               "total_potential_platform_erosion",
                               "total_actual_platform_erosion",
                               "potential_beach_erosion",
                               "actual_beach_erosion",
                               "total_potential_beach_erosion",
                               "total_actual_beach_erosion",
                               "beach_deposition",
                               "total_beach_deposition",
                               "landform",
                               "local_cons_sediment_slope",
                               "slope",
                               "cliff",
                               "avg_sea_depth",
                               "avg_wave_height",
                               "avg_wave_orientation",
                               "beach_protection",
                               "basement_elevation",
                               "coastline",
                               "coast_normal",
                               "active_zone",
                               "cliff_collapse",
                               "total_cliff_collapse",
                               "cliff_collapse_deposition",
                               "total_cliff_collapse_deposition",
                               "polygon",
                               "potential_platform_erosion_mask",
                               "sea_mask",
                               "beach_mask",
                               "shadow_zone_codes",
                               "deep_water_wave_angle",
                               "deep_water_wave_height",
                               "deep_water_wave_period",
                               "polygon_uncons_sediment_up_or_down_drift",
                               "polygon_uncons_sediment_gain_or_loss"});
      }
      else if (fileSpecLower == "usual")
      {
         // Add usual/standard raster outputs (Case 11 "usual" mode)
         expandedFiles.insert(expandedFiles.end(),
                              {"suspended_sediment",
                               "avg_suspended_sediment",
                               "fine_uncons",
                               "fine_cons",
                               "sand_uncons",
                               "sand_cons",
                               "coarse_uncons",
                               "coarse_cons",
                               "sediment_top_elevation",
                               "top_elevation",
                               "sea_depth",
                               "wave_height",
                               "wave_orientation",
                               "potential_platform_erosion",
                               "actual_platform_erosion",
                               "total_potential_platform_erosion",
                               "total_actual_platform_erosion",
                               "potential_beach_erosion",
                               "actual_beach_erosion",
                               "total_potential_beach_erosion",
                               "total_actual_beach_erosion",
                               "beach_deposition",
                               "total_beach_deposition",
                               "landform",
                               "local_cons_sediment_slope",
                               "slope",
                               "avg_wave_height",
                               "avg_wave_orientation",
                               "beach_protection",
                               "basement_elevation",
                               "active_zone",
                               "cliff_collapse",
                               "total_cliff_collapse",
                               "cliff_collapse_deposition",
                               "total_cliff_collapse_deposition",
                               "polygon",
                               "shadow_zone_codes",
                               "polygon_uncons_sediment_up_or_down_drift",
                               "polygon_uncons_sediment_gain_or_loss",
                               "coast_normal",
                               "coastline",
                               "fine_uncons",
                               "fine_cons",
                               "sand_uncons",
                               "sand_cons",
                               "coarse_uncons",
                               "coarse_cons"});
      }
      else if (fileSpecLower == "Cmetools")
      {
         // Add usual/standard raster outputs (Case 11 "usual" mode)
         expandedFiles.insert(expandedFiles.end(),
                              {"fine_uncons",
                               "fine_cons",
                               "sand_uncons",
                               "sand_cons",
                               "coarse_uncons",
                               "coarse_cons",
                               "top_elevation",
                               "sea_depth",
                               "wave_height",
                               "total_actual_platform_erosion",
                               "total_actual_beach_erosion",
                               "total_beach_deposition",
                               "landform",
                               "basement_elevation",
                               "active_zone",
                               "cliff_collapse",
                               "total_cliff_collapse",
                               "cliff_collapse_deposition",
                               "total_cliff_collapse_deposition",
                               "polygon",
                               "coast_normal",
                               "coastline"});
      }
      else
      {
         // Regular file specification - add as-is
         expandedFiles.push_back(fileSpec);
      }
   }

   return expandedFiles;
}

//===============================================================================================================================
//! Get vector files with keyword expansion support
//===============================================================================================================================
vector<string> CConfiguration::GetVectorFiles() const
{
   // Case 16: Vector GIS files to output - expand "all" and "usual" keywords
   vector<string> expandedFiles;

   for (string const &fileSpec : m_vecVectorFiles)
   {
      string fileSpecLower = fileSpec;
      std::transform(fileSpecLower.begin(), fileSpecLower.end(),
                     fileSpecLower.begin(), ::tolower);

      if (fileSpecLower == "all")
      {
         // Add all possible vector outputs (Case 16 "all" mode)
         expandedFiles.insert(expandedFiles.end(), {"coast", "cliff_edge", "wave_angle", "normals", "invalid_normals", "avg_wave_angle", "wave_energy", "mean_wave_energy", "breaking_wave_height", "coast_curvature", "polygon_node", "polygon", "cliff_notch", "shadow_boundary", "downdrift_boundary", "deep_water_wave_angle", "wave_setup", "storm_surge", "run_up", "flood_line"});
      }
      else if (fileSpecLower == "usual")
      {
         // Add usual/standard vector outputs (Case 16 "usual" mode)
         expandedFiles.insert(
            expandedFiles.end(),
            {"coast", "cliff_edge", "wave_angle", "normals", "invalid_normals",
             "avg_wave_angle", "wave_energy", "mean_wave_energy",
             "breaking_wave_height", "polygon", "cliff_notch",
             "shadow_boundary", "downdrift_boundary", "deep_water_wave_angle"});
      }
      else
      {
         // Regular file specification - add as-is
         expandedFiles.push_back(fileSpec);
      }
   }

   return expandedFiles;
}

string CConfiguration::GetOmitGridEdges() const
{
   // This needs to be lower case
   std::string my_text{m_strOmitGridEdges};
   std::transform(my_text.begin(), my_text.end(), my_text.begin(), ::tolower);
   return my_text;
}
