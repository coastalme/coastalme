/*!
 *
 * \file write_output.cpp
 * \brief Writes non-GIS output files
 * \details TODO 001 A more detailed description of this routine.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2024
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include <assert.h>

#include <ctime>
using std::localtime;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;
using std::noshowpos;
using std::showpos;

#include <iomanip>
using std::put_time;
using std::resetiosflags;
using std::setfill;
using std::setiosflags;
using std::setprecision;
using std::setw;

#include <sstream>
using std::stringstream;

#include <string>
using std::to_string;

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "interpolate.h"

//===============================================================================================================================
//! Writes beginning-of-run information to Out and Log files
//===============================================================================================================================
void CSimulation::WriteStartRunDetails(void)
{
   // Set the Out file output format to fixed point
   OutStream << std::fixed;

   // Start outputting stuff
   OutStream << PROGRAM_NAME << " for " << PLATFORM << " " << strGetBuild() << " on " << strGetComputerName() << endl
             << endl;

   LogStream << PROGRAM_NAME << " for " << PLATFORM << " " << strGetBuild() << " on " << strGetComputerName() << endl
             << endl;

   // ----------------------------------------------- Run Information ----------------------------------------------------------
   OutStream << "RUN DETAILS" << endl;
   OutStream << " Name                                                      \t: " << m_strRunName << endl;
   OutStream << " Run started                                               \t: " << put_time(localtime(&m_tSysStartTime), "%T %A %d %B %Y") << endl;

   // Same info. for Log file
   LogStream << m_strRunName << " run started at " << put_time(localtime(&m_tSysStartTime), "%T on %A %d %B %Y") << endl
             << endl;

   // Contine with Out file
   OutStream << " Initialization file                                       \t: "
#ifdef _WIN32
             << pstrChangeToForwardSlash(&m_strCMEIni) << endl;
#else
             << m_strCMEIni << endl;
#endif

   OutStream << " Input data read from                                      \t: "
#ifdef _WIN32
             << pstrChangeToForwardSlash(&m_strDataPathName) << endl;
#else
             << m_strDataPathName << endl;
#endif

   OutStream << " Main output file (this file)                              \t: "
#ifdef _WIN32
             << pstrChangeToForwardSlash(&m_strOutFile) << OUTEXT << endl;
#else
             << m_strOutFile << OUTEXT << endl;
#endif

   LogStream << "Main output file                                          \t: "
#ifdef _WIN32
             << pstrChangeToForwardSlash(&m_strOutFile) << OUTEXT << endl;
#else
             << m_strOutFile << OUTEXT << endl;
#endif

   OutStream << " Log file                                                  \t: "
#ifdef _WIN32
             << pstrChangeToForwardSlash(&m_strOutFile) << LOGEXT << endl;
#else
             << m_strOutFile << LOGEXT << endl;
#endif

   LogStream << "Log file (this file)                                      \t: "
#ifdef _WIN32
             << pstrChangeToForwardSlash(&m_strOutFile) << LOGEXT << endl;
#else
             << m_strOutFile << LOGEXT << endl;
#endif

   OutStream << " Level of Log detail                                       \t: ";
   if (m_nLogFileDetail == NO_LOG_FILE)
      OutStream << "0 (least detail)";
   else if (m_nLogFileDetail == LOG_FILE_LOW_DETAIL)
      OutStream << "1 (least detail)";
   else if (m_nLogFileDetail == LOG_FILE_MIDDLE_DETAIL)
      OutStream << "2 (medium detail)";
   else if (m_nLogFileDetail == LOG_FILE_HIGH_DETAIL)
      OutStream << "3 (high detail)";
   else if (m_nLogFileDetail == LOG_FILE_ALL)
      OutStream << "4 (everything)";
   OutStream << endl;

   LogStream << "Level of Log detail                                       \t: ";
   if (m_nLogFileDetail == LOG_FILE_LOW_DETAIL)
      LogStream << "1 (least detail)";
   else if (m_nLogFileDetail == LOG_FILE_MIDDLE_DETAIL)
      LogStream << "2 (medium detail)";
   else if (m_nLogFileDetail == LOG_FILE_HIGH_DETAIL)
      OutStream << "3 (high detail)";
   else if (m_nLogFileDetail == LOG_FILE_ALL)
      OutStream << "4 (everything)";
   LogStream << endl << endl;

   OutStream << " Simulation start date/time                                \t: ";
   // hh:mm:ss dd/mm/yyyy
   char cPrev = OutStream.fill('0');
   OutStream << setw(2) << m_nSimStartHour << COLON << setw(2) << m_nSimStartMin << COLON << setw(2) << m_nSimStartSec << SPACE << setw(2) << m_nSimStartDay << SLASH << setw(2) << m_nSimStartMonth << SLASH << setw(2) << m_nSimStartYear << endl;
   OutStream.fill(cPrev);

   OutStream << " Duration of simulation                                    \t: ";
   OutStream << strDispSimTime(m_dSimDuration) << endl;
   if (m_bSaveRegular)
   {
      // Saves at regular intervals
      OutStream << " Time between saves                                        \t: ";
      OutStream << strDispSimTime(m_dRegularSaveInterval) << endl;
   }
   else
   {
      // Saves at user-defined intervals
      OutStream << " Saves at                                                  \t: ";
      string strTmp;
      for (int i = 0; i < m_nUSave; i++)
      {
         strTmp.append(strDispSimTime(m_dUSaveTime[i]));
         strTmp.append(", ");
      }

      // Also at end of run
      strTmp.append(strDispSimTime(m_dSimDuration));
      OutStream << strTmp << endl;
   }
   OutStream << " Raster GIS output format                                  \t: " << m_strGDALRasterOutputDriverLongname << endl;
   OutStream << " Maximum number of GIS Save Number digits                  \t: " << m_nGISMaxSaveDigits << endl;
   OutStream << " GIS Save Numbers sequential (S) or iteration number (I)   \t: " << (m_bGISSaveDigitsSequential ? "S" : "I") << endl;
   OutStream << " Random number seeds                                       \t: ";
   {
      for (int i = 0; i < NRNG; i++)
         OutStream << m_ulRandSeed[i] << '\t';
   }
   OutStream << endl;

   OutStream << "*First random numbers generated                            \t: " << ulGetRand0() << '\t' << ulGetRand1() << endl;
   OutStream << " Raster GIS output format                                  \t: " << m_strGDALRasterOutputDriverLongname << endl;
   OutStream << " Raster output values scaled (if needed)                   \t: " << (m_bScaleRasterOutput ? "Y" : "N") << endl;
   OutStream << " Raster world files created (if needed)                    \t: " << (m_bWorldFile ? "Y" : "N") << endl;
   OutStream << " Raster GIS files saved                                    \t: " << strListRasterFiles() << endl;
   if (m_bSliceSave)
   {
      OutStream << std::fixed << setprecision(3);
      OutStream << " Elevations for 'slice' raster output files                \t: ";
      for (int i = 0; i < static_cast<int>(m_VdSliceElev.size()); i++)
         OutStream << m_VdSliceElev[i] << " ";
      OutStream << endl;
   }

   OutStream << " Vector GIS output format                                  \t: " << m_strVectorGISOutFormat << endl;
   OutStream << " Vector GIS files saved                                    \t: " << strListVectorFiles() << endl;
   OutStream << " Output file (this file)                                   \t: "
#ifdef _WIN32
             << pstrChangeToForwardSlash(&m_strOutFile) << endl;
#else
             << m_strOutFile << endl;
#endif
   OutStream << " Log file                                                  \t: "
#ifdef _WIN32
             << pstrChangeToForwardSlash(&m_strLogFile) << endl;
#else
             << m_strLogFile << endl;
#endif

   OutStream << " Optional time series files saved                          \t: " << strListTSFiles() << endl;

   OutStream << " Coastline vector smoothing algorithm                      \t: ";
   switch (m_nCoastSmooth)
   {
   case SMOOTH_NONE:
   {
      OutStream << "none";
      break;
   }

   case SMOOTH_RUNNING_MEAN:
   {
      OutStream << "running mean";
      break;
   }

   case SMOOTH_SAVITZKY_GOLAY:
   {
      OutStream << "Savitzky-Golay";
      break;
   }
   }
   OutStream << endl;

   OutStream << " Grid edge(s) to omit when searching for coastline         \t: " << (m_bOmitSearchNorthEdge ? "N" : "") << (m_bOmitSearchSouthEdge ? "S" : "") << (m_bOmitSearchWestEdge ? "W" : "") << (m_bOmitSearchEastEdge ? "E" : "") << endl;

   if (m_nCoastSmooth != SMOOTH_NONE)
   {
      OutStream << " Size of coastline vector smoothing window                 \t: " << m_nCoastSmoothWindow << endl;

      if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
         OutStream << " Savitzky-Golay coastline smoothing polynomial order       \t: " << m_nSavGolCoastPoly << endl;
   }
   OutStream << " Size of profile slope smoothing window                    \t: " << m_nProfileSmoothWindow << endl;
   OutStream << " Max local slope on profile (m/m)                          \t: " << m_dProfileMaxSlope << endl;
   OutStream << " Vertical tolerance for beach to be included in smoothing  \t: " << m_dMaxBeachElevAboveSWL << " m" << endl;
   OutStream << endl;

   // --------------------------------------------------- Raster GIS stuff -------------------------------------------------------
   OutStream << "Raster GIS Input Files" << endl;
   OutStream << " Basement DEM file                                         \t: "
#ifdef _WIN32
             << pstrChangeToForwardSlash(&m_strInitialBasementDEMFile) << endl;
#else
             << m_strInitialBasementDEMFile << endl;
#endif
   OutStream << " Basement DEM driver code                                  \t: " << m_strGDALBasementDEMDriverCode << endl;
   OutStream << " GDAL basement DEM driver description                      \t: " << m_strGDALBasementDEMDriverDesc << endl;
   OutStream << " GDAL basement DEM projection                              \t: " << m_strGDALBasementDEMProjection << endl;
   OutStream << " GDAL basement DEM data type                               \t: " << m_strGDALBasementDEMDataType << endl;
   OutStream << " Grid size (X by Y)                                        \t: " << m_nXGridMax << " by " << m_nYGridMax << endl;
   OutStream << resetiosflags(ios::floatfield);
   OutStream << std::fixed << setprecision(1);
   OutStream << "*Coordinates of NW corner of grid (external CRS)           \t: " << m_dNorthWestXExtCRS << ", " << m_dNorthWestYExtCRS << endl;
   OutStream << "*Coordinates of SE corner of grid (external CRS)           \t: " << m_dSouthEastXExtCRS << ", " << m_dSouthEastYExtCRS << endl;
   OutStream << "*Cell size                                                 \t: " << m_dCellSide << " m" << endl;
   OutStream << "*Grid area                                                 \t: " << m_dExtCRSGridArea << " m^2" << endl;
   OutStream << std::fixed << setprecision(2);
   OutStream << "*Grid area                                                 \t: " << m_dExtCRSGridArea * 1e-6 << " km^2" << endl;

   if (! m_strInitialLandformFile.empty())
   {
      OutStream << " Initial Landform Class file                               \t: " << m_strInitialLandformFile << endl;
      OutStream << " GDAL Initial Landform Class file driver code              \t: " << m_strGDALLDriverCode << endl;
      OutStream << " GDAL Initial Landform Class file driver description       \t: " << m_strGDALLDriverDesc << endl;
      OutStream << " GDAL Initial Landform Class file projection               \t: " << m_strGDALLProjection << endl;
      OutStream << " GDAL Initial Landform Class file data type                \t: " << m_strGDALLDataType << endl;
      OutStream << endl;
   }

   if (! m_strInterventionClassFile.empty())
   {
      OutStream << " Intervention Class file                                   \t: " << m_strInterventionClassFile << endl;
      OutStream << " GDAL Intervention Class file driver code                  \t: " << m_strGDALICDriverCode << endl;
      OutStream << " GDAL Intervention Class file driver description           \t: " << m_strGDALICDriverDesc << endl;
      OutStream << " GDAL Intervention Class file projection                   \t: " << m_strGDALICProjection << endl;
      OutStream << " GDAL Intervention Class file data type                    \t: " << m_strGDALICDataType << endl;
      OutStream << endl;
   }

   if (! m_strInterventionHeightFile.empty())
   {
      OutStream << " Intervention Height file                                  \t: " << m_strInterventionHeightFile << endl;
      OutStream << " GDAL Intervention Height file driver code                 \t: " << m_strGDALIHDriverCode << endl;
      OutStream << " GDAL Intervention Height file driver description          \t: " << m_strGDALIHDriverDesc << endl;
      OutStream << " GDAL Intervention Height file projection                  \t: " << m_strGDALIHProjection << endl;
      OutStream << " GDAL Intervention Height file data type                   \t: " << m_strGDALIHDataType << endl;
      OutStream << endl;
   }

   if (! m_strInitialSuspSedimentFile.empty())
   {
      OutStream << " Initial Susp Sediment file                                \t: " << m_strInitialSuspSedimentFile << endl;
      OutStream << " GDAL Initial Susp Sediment file driver code               \t: " << m_strGDALISSDriverCode << endl;
      OutStream << " GDAL Initial Susp Sediment file driver description        \t: " << m_strGDALISSDriverDesc << endl;
      OutStream << " GDAL Initial Susp Sediment file projection                \t: " << m_strGDALISSProjection << endl;
      OutStream << " GDAL Initial Susp Sediment file data type                 \t: " << m_strGDALISSDataType << endl;
      OutStream << endl;
   }

   for (int i = 0; i < m_nLayers; i++)
   {
      if (m_nLayers == 1)
         OutStream << " Only one layer" << endl;
      else
         OutStream << " Layer " << i << (i == 0 ? "(Top)" : "") << (i == m_nLayers - 1 ? "(Bottom)" : "") << endl;

      if (! m_VstrInitialFineUnconsSedimentFile[i].empty())
      {
         OutStream << "    Initial Fine Uncons Sediment file                      \t: " << m_VstrInitialFineUnconsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Fine Uncons Sediment file driver code     \t: " << m_VstrGDALIUFDriverCode[i] << endl;
         OutStream << "    GDAL Initial Fine Uncons Sediment file driver desc     \t: " << m_VstrGDALIUFDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Fine Uncons Sediment file projection      \t: " << m_VstrGDALIUFProjection[i] << endl;
         OutStream << "    GDAL Initial Fine Uncons Sediment file data type       \t: " << m_VstrGDALIUFDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialSandUnconsSedimentFile[i].empty())
      {
         OutStream << "    Initial Sand Uncons Sediment file                      \t: " << m_VstrInitialSandUnconsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Sand Uncons Sediment file driver code     \t: " << m_VstrGDALIUSDriverCode[i] << endl;
         OutStream << "    GDAL Initial Sand Uncons Sediment file driver desc     \t: " << m_VstrGDALIUSDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Sand Uncons Sediment file projection      \t: " << m_VstrGDALIUSProjection[i] << endl;
         OutStream << "    GDAL Initial Sand Uncons Sediment file data type       \t: " << m_VstrGDALIUSDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialCoarseUnconsSedimentFile[i].empty())
      {
         OutStream << "    Initial Coarse Uncons Sediment file                    \t: " << m_VstrInitialCoarseUnconsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Coarse Uncons Sediment file driver code   \t: " << m_VstrGDALIUCDriverCode[i] << endl;
         OutStream << "    GDAL Initial Coarse Uncons Sediment file driver desc   \t: " << m_VstrGDALIUCDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Coarse Uncons Sediment file projection    \t: " << m_VstrGDALIUCProjection[i] << endl;
         OutStream << "    GDAL Initial Coarse Uncons Sediment file data type     \t: " << m_VstrGDALIUCDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialFineConsSedimentFile[i].empty())
      {
         OutStream << "    Initial Fine Cons Sediment file                        \t: " << m_VstrInitialFineConsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Fine Cons Sediment file driver code       \t: " << m_VstrGDALICFDriverCode[i] << endl;
         OutStream << "    GDAL Initial Fine Cons Sediment file driver desc       \t: " << m_VstrGDALICFDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Fine Cons Sediment file projection        \t: " << m_VstrGDALICFProjection[i] << endl;
         OutStream << "    GDAL Initial Fine Cons Sediment file data type         \t: " << m_VstrGDALICFDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialSandConsSedimentFile[i].empty())
      {
         OutStream << "    Initial Sand Cons Sediment file                        \t: " << m_VstrInitialSandConsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Sand Cons Sediment file driver code       \t: " << m_VstrGDALICSDriverCode[i] << endl;
         OutStream << "    GDAL Initial Sand Cons Sediment file driver desc       \t: " << m_VstrGDALICSDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Sand Cons Sediment file projection        \t: " << m_VstrGDALICSProjection[i] << endl;
         OutStream << "    GDAL Initial Sand Cons Sediment file data type         \t: " << m_VstrGDALICSDataType[i] << endl;
         OutStream << endl;
      }

      if (! m_VstrInitialCoarseConsSedimentFile[i].empty())
      {
         OutStream << "    Initial Coarse Cons Sediment file                      \t: " << m_VstrInitialCoarseConsSedimentFile[i] << endl;
         OutStream << "    GDAL Initial Coarse Cons Sediment file driver code     \t: " << m_VstrGDALICCDriverCode[i] << endl;
         OutStream << "    GDAL Initial Coarse Cons Sediment file driver desc     \t: " << m_VstrGDALICCDriverDesc[i] << endl;
         OutStream << "    GDAL Initial Coarse Cons Sediment file projection      \t: " << m_VstrGDALICCProjection[i] << endl;
         OutStream << "    GDAL Initial Coarse Cons Sediment file data type       \t: " << m_VstrGDALICCDataType[i] << endl;
         OutStream << endl;
      }
   }
   //   OutStream << endl;

   // ---------------------------------------------------- Vector GIS stuff ------------------------------------------------------
   OutStream << "Vector GIS Input Files" << endl;

   if (m_bSingleDeepWaterWaveValues)
      OutStream << " None" << endl;
   else
   {

      OutStream << " Deep water wave stations shapefile                        \t: " << m_strDeepWaterWaveStationsShapefile << endl;
      OutStream << " GDAL/OGR deep water wave stations shapefile driver code   \t: " << m_strOGRDWWVDriverCode << endl;
      OutStream << " GDAL/OGR deep water wave stations shapefile data type     \t: " << m_strOGRDWWVDataType << endl;
      OutStream << " GDAL/OGR deep water wave stations shapefile geometry      \t: " << m_strOGRDWWVGeometry << endl;
      OutStream << " Deep water wave values file                               \t: " << m_strDeepWaterWavesTimeSeriesFile << endl;

      if (m_dWaveDataWrapHours > 0)
         OutStream << " Deep water wave values will wrap every " << m_dWaveDataWrapHours << " hours" << endl;
   }
   OutStream << endl;

   // -------------------------------------------------------- Other data --------------------------------------------------------
   OutStream << "Other Input Data" << endl;

   OutStream << " Wave propagation model                                    \t: ";
   if (m_nWavePropagationModel == WAVE_MODEL_COVE)
      OutStream << "COVE";
   else if (m_nWavePropagationModel == WAVE_MODEL_CSHORE)
      OutStream << "CShore";
   OutStream << endl;
   OutStream << " Density of sea water                                     \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(0) << m_dSeaWaterDensity << " kg/m^3" << endl;
   OutStream << " Initial still water level                                 \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(1) << m_dOrigSWL << " m" << endl;
   OutStream << " Final still water level                                   \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(1) << m_dFinalSWL << " m" << endl;
   if (m_bSingleDeepWaterWaveValues)
   {
      OutStream << " Deep water wave height                                    \t: " << m_dAllCellsDeepWaterWaveHeight << " m" << endl;
      OutStream << " Deep water wave orientation                               \t: " << m_dAllCellsDeepWaterWaveAngle << " degrees" << endl;
      OutStream << " Wave period                                               \t: " << m_dAllCellsDeepWaterWavePeriod << " s" << endl;
   }
   else
   {
      OutStream << " Maximum User input Deep water wave height                 \t: " << m_dMaxUserInputWaveHeight << " m" << endl;
      OutStream << " Maximum User input Deep waterWave period                  \t: " << m_dMaxUserInputWavePeriod << " s" << endl;
   }
   OutStream << " Start depth for wave calcs (*deep water wave height)      \t: " << m_dWaveDepthRatioForWaveCalcs << endl;
   OutStream << "*Depth of closure                                          \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(3) << m_dDepthOfClosure << " m" << endl;
   OutStream << " Tide data file                                            \t: " << m_strTideDataFile << endl;
   OutStream << " Do coast platform erosion?                                \t: " << (m_bDoShorePlatformErosion ? "Y" : "N") << endl;
   OutStream << " Coast platform resistance to erosion                      \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(3) << m_dR << endl;
   OutStream << " Do beach sediment transport?                              \t: " << (m_bDoBeachSedimentTransport ? "Y" : "N") << endl;
   OutStream << " Handling of beach sediment at grid edges                  \t: ";
   if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
      OutStream << "closed";
   else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
      OutStream << "open";
   else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
      OutStream << "recirculate";
   OutStream << endl;
   OutStream << " Beach potential erosion/deposition equation               \t: ";
   if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_CERC)
      OutStream << "CERC";
   else if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_KAMPHUIS)
      OutStream << "Kamphuis";
   OutStream << endl;
   OutStream << " Median particle size of fine sediment                     \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dD50Fine << " mm" << endl;
   OutStream << " Median particle size of sand sediment                     \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dD50Sand << " mm" << endl;
   OutStream << " Median particle size of coarse sediment                   \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dD50Coarse << " mm" << endl;
   OutStream << " Beach sediment density                                    \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dBeachSedimentDensity << " kg/m^3" << endl;
   OutStream << " Beach sediment porosity                                   \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dBeachSedimentPorosity << endl;
   OutStream << " Fine-sized sediment relative erodibility                  \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(1) << m_dFineErodibility << endl;
   OutStream << " Sand-sized sediment relative erodibility                  \t: " << resetiosflags(ios::floatfield) << m_dSandErodibility << endl;
   OutStream << " Coarse-sized sediment relative erodibility                \t: " << m_dCoarseErodibility << endl;
   if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_CERC)
      OutStream << " Transport parameter KLS for CERC equation                 \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(3) << m_dKLS << endl;
   if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_KAMPHUIS)
      OutStream << " Transport parameter for Kamphuis equation                 \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(3) << m_dKamphuis << endl;
   OutStream << " Height of Dean profile start above SWL                    \t: " << resetiosflags(ios::floatfield) << std::fixed << setprecision(1) << m_dDeanProfileStartAboveSWL << " m" << endl;
   OutStream << " Sediment input at a point                                 \t: " << (m_bSedimentInput ? "Y" : "N") << endl;
   if (m_bSedimentInput)
   {
      OutStream << " Sediment input shapefile                                  \t: " << m_strSedimentInputEventShapefile << endl;
      OutStream << " Sediment input type                                       \t: ";
      if (m_bSedimentInputAtPoint)
         OutStream << "at point";
      else if (m_bSedimentInputAtCoast)
         OutStream << "in block on coast";
      else if (m_bSedimentInputAlongLine)
         OutStream << "where line interests with coast";
      OutStream << endl;
      OutStream << " Sediment input time series file                           \t: " << m_strSedimentInputEventTimeSeriesFile << endl;
   }
   OutStream << " Do cliff collapse?                                        \t: " << (m_bDoCliffCollapse ? "Y" : "N") << endl;
   OutStream << " Cliff resistance to erosion                               \t: " << m_dCliffErosionResistance << endl;
   OutStream << " Notch overhang to initiate collapse                       \t: " << m_dNotchDepthAtCollapse << " m" << endl;
   OutStream << " Notch base below SWL                                      \t: " << m_dNotchBaseBelowSWL << " m" << endl;
   OutStream << " Scale parameter A for cliff deposition                    \t: ";
   if (bFPIsEqual(m_dCliffDepositionA, 0.0, TOLERANCE))
      OutStream << "auto";
   else
      OutStream << m_dCliffDepositionA << "  m^(1/3)";
   OutStream << endl;
   OutStream << " Planview width of cliff deposition talus                  \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dCliffDepositionPlanviewWidth << " m" << endl;
   OutStream << " Planview length of cliff deposition talus                 \t: " << m_dCliffTalusMinDepositionLength << " m" << endl;
   OutStream << " Min height of land-end talus (fraction of cliff elevation)\t: " << m_dMinCliffTalusHeightFrac << endl;
   OutStream << " Do riverine flooding?                                     \t: " << (m_bDoRiverineFlooding ? "Y" : "N") << endl;
   if (m_bDoRiverineFlooding)
   {
      // BUG 002 Need more info on this
      OutStream << " FloodSWLSetupLine                                         \t: " << (m_bFloodSWLSetupLine ? "Y" : "N") << endl;
      OutStream << " FloodSWLSetupSurgeLine                                    \t: " << (m_bFloodSWLSetupSurgeLine ? "Y" : "N") << endl;
      OutStream << " m_bFloodSWLSetupSurgeRunupLine                            \t: " << (m_bFloodSWLSetupSurgeRunupLine ? "Y" : "N") << endl;
   }
   OutStream << " Gravitational acceleration                                \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dG << " m^2/s" << endl;
   OutStream << " Minimum spacing of coastline normals                      \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dCoastNormalAvgSpacing << " m" << endl;
   OutStream << " Random factor for spacing of normals                      \t: " << resetiosflags(ios::floatfield) << std::fixed << m_dCoastNormalRandSpacingFactor << endl;
   OutStream << " Length of coastline normals                               \t: " << m_dCoastNormalLength << " m" << endl;
   OutStream << " Maximum number of 'cape' normals                          \t: " << m_nNaturalCapeNormals << endl;
   OutStream << endl;
   /*
      OutStream << std::fixed << setprecision(8);
      OutStream << " Erosion potential shape function:" << endl;
      OutStream << "\tDepth over DB\tErosion potential\tFirst derivative of erosion potential" << endl;
      for (int i = 0; i < m_VdDepthOverDB.size(); i++)
         OutStream << "\t" << m_VdDepthOverDB[i] << "\t\t" << m_VdErosionPotential[i] << "\t\t" << m_VdErosionPotentialFirstDeriv[i] << endl;
      OutStream << endl;
   */
   // ------------------------------------------------------ Testing only --------------------------------------------------------
   OutStream << "Testing only" << endl;

   OutStream << " Output profile data?                                      \t: " << (m_bOutputProfileData ? "Y" : "N") << endl;
   OutStream << " Profile numbers to be saved                               \t: ";
   for (unsigned int i = 0; i < m_VnProfileToSave.size(); i++)
      OutStream << m_VnProfileToSave[i] << SPACE;
   OutStream << endl;
   OutStream << " Timesteps when profiles are saved                         \t: ";
   for (unsigned int i = 0; i < m_VulProfileTimestep.size(); i++)
      OutStream << m_VulProfileTimestep[i] << SPACE;
   OutStream << endl;
   OutStream << " Output parallel profile data?                             \t: " << (m_bOutputParallelProfileData ? "Y" : "N") << endl;
   OutStream << " Output erosion potential look-up data?                    \t: " << (m_bOutputErosionPotentialData ? "Y" : "N");
   if (m_bOutputErosionPotentialData)
      OutStream << " (see " << m_strOutPath << EROSION_POTENTIAL_LOOKUP_FILE << ")";
   OutStream << endl;
   OutStream << " Erode coast in alternate directions?                      \t: " << (m_bErodeShorePlatformAlternateDirection ? "Y" : "N") << endl;
   OutStream << " Size of moving window for calculating coastline curvature \t: " << m_nCoastCurvatureMovingWindowSize << endl;

   OutStream << endl
             << endl;

   // -------------------------------------------------- Per-iteration output ----------------------------------------------------
   OutStream << std::fixed << setprecision(3);

   // Write per-timestep headers to .out file
   OutStream << PER_ITER_HEAD << endl;
   OutStream << "Sea depth in metres. All erosion and deposition values in millimetres" << endl;
   OutStream << "GISn = GIS files saved as <filename>n." << endl;
   OutStream << endl;

   OutStream << PER_ITER_HEAD1 << endl;
   OutStream << PER_ITER_HEAD2 << endl;
   OutStream << PER_ITER_HEAD3 << endl;
   OutStream << PER_ITER_HEAD4 << endl;
   OutStream << PER_ITER_HEAD5 << endl;
}

//===============================================================================================================================
//! Write the results for this timestep to the .out file
//===============================================================================================================================
bool CSimulation::bWritePerTimestepResults(void)
{
   OutStream << resetiosflags(ios::floatfield);
   OutStream << std::fixed << setprecision(0);

   // Output timestep and simulated time info ===================================================================================
   OutStream << setw(4) << m_ulIter;
   OutStream << setw(7) << m_dSimElapsed; // In hours
   OutStream << resetiosflags(ios::floatfield);
   OutStream << std::fixed << setprecision(0);
   OutStream << setw(7) << m_dSimElapsed / (24 * 365.25); // In years

   // Output average sea depth (m) per sea cell =================================================================================
   OutStream << resetiosflags(ios::floatfield);
   OutStream << std::fixed << setprecision(2);
   double dAvgSeaDepth = m_dThisIterTotSeaDepth / static_cast<double>(m_ulThisIterNumSeaCells);
   OutStream << setw(6) << dAvgSeaDepth;
   OutStream << " ";

   // Output the this-timestep % of sea cells with potential shore platform erosion =============================================
   OutStream << std::fixed << setprecision(0);
   OutStream << setw(6) << 100 * static_cast<double>(m_ulThisIterNumPotentialPlatformErosionCells) / static_cast<double>(m_ulThisIterNumSeaCells);

   // Output per-timestep potential shore platform erosion in m (average for all sea cells)
   OutStream << std::fixed << setprecision(1);
   OutStream << setw(6) << 1000 * m_dThisIterPotentialPlatformErosion / static_cast<double>(m_ulThisIterNumSeaCells);

   // Output per-timestep potential shore platform erosion in m (average for all cells with potential shore platform erosion)
   OutStream << std::fixed << setprecision(1);
   if (m_ulThisIterNumPotentialPlatformErosionCells > 0)
      OutStream << setw(6) << 1000 * m_dThisIterPotentialPlatformErosion / static_cast<double>(m_ulThisIterNumPotentialPlatformErosionCells);
   else
      OutStream << setw(6) << SPACE;

   // Output the this-timestep % of sea cells with actual shore platform erosion ================================================
   OutStream << std::fixed << setprecision(0);
   OutStream << setw(6) << 100 * static_cast<double>(m_ulThisIterNumActualPlatformErosionCells) / static_cast<double>(m_ulThisIterNumSeaCells);

   // Output per-timestep actual shore platform erosion in m (average for all sea cells)
   OutStream << std::fixed << setprecision(1);
   double dThisIterActualPlatformErosion = m_dThisIterActualPlatformErosionFineCons + m_dThisIterActualPlatformErosionSandCons + m_dThisIterActualPlatformErosionCoarseCons;
   OutStream << setw(6) << 1000 * dThisIterActualPlatformErosion / static_cast<double>(m_ulThisIterNumSeaCells);

   // Output per-timestep actual shore platform erosion in m (average for all cells with actual shore platform erosion)
   OutStream << std::fixed << setprecision(1);
   if (m_ulThisIterNumActualPlatformErosionCells > 0)
      OutStream << setw(5) << 1000 * dThisIterActualPlatformErosion / static_cast<double>(m_ulThisIterNumActualPlatformErosionCells);
   else
      OutStream << setw(5) << SPACE;

   // Output per-timestep actual shore platform erosion in m (average for all sea cells)
   OutStream << std::fixed << setprecision(1);

   if (m_dThisIterActualPlatformErosionFineCons > 0)
      OutStream << setw(4) << 1000 * m_dThisIterActualPlatformErosionFineCons / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   if (m_dThisIterActualPlatformErosionSandCons > 0)
      OutStream << setw(4) << 1000 * m_dThisIterActualPlatformErosionSandCons / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   if (m_dThisIterActualPlatformErosionCoarseCons > 0)
      OutStream << setw(4) << 1000 * m_dThisIterActualPlatformErosionCoarseCons / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   // Output the this-timestep % of sea cells with potential beach erosion ======================================================
   OutStream << setw(7) << 100 * static_cast<double>(m_ulThisIterNumPotentialBeachErosionCells) / static_cast<double>(m_ulThisIterNumSeaCells);

   // Output per-timestep potential beach erosion in m (average for all sea cells)
   OutStream << std::fixed << setprecision(0);
   // assert(m_ulThisIterNumSeaCells > 0);
   double dTmp = 1000 * m_dThisIterPotentialBeachErosion / static_cast<double>(m_ulThisIterNumSeaCells);
   if (dTmp > 99999)
   {
      OutStream << setw(6) << std::scientific << setprecision(0) << dTmp;
      OutStream << std::fixed;
   }
   else
      OutStream << setw(6) << dTmp;

   // Output per-timestep potential beach erosion in m (average for all cells with potential beach erosion)
   OutStream << std::fixed << setprecision(1);
   if (m_ulThisIterNumPotentialBeachErosionCells > 0)
   {
      dTmp = 1000 * m_dThisIterPotentialBeachErosion / static_cast<double>(m_ulThisIterNumPotentialBeachErosionCells);
      if (dTmp > 99999)
      {
         OutStream << setw(6) << std::scientific << setprecision(0) << dTmp;
         OutStream << std::fixed;
      }
      else
         OutStream << setw(6) << dTmp;
   }
   else
      OutStream << setw(6) << SPACE;

   // This-timestep % of sea cells with actual beach erosion ====================================================================
   OutStream << std::fixed << setprecision(0);
   OutStream << setw(7) << 100 * static_cast<double>(m_ulThisIterNumActualBeachErosionCells) / static_cast<double>(m_ulThisIterNumSeaCells);

   // Output per-timestep actual beach erosion in m (average for all sea cells)
   double dThisIterActualBeachErosion = m_dThisIterBeachErosionFine + m_dThisIterBeachErosionSand + m_dThisIterBeachErosionCoarse;
   OutStream << setw(6) << 1000 * dThisIterActualBeachErosion / static_cast<double>(m_ulThisIterNumSeaCells);

   // Per-iteration actual beach erosion in m (average for all cells with actual beach erosion)
   OutStream << std::fixed << setprecision(1);
   if (m_ulThisIterNumActualBeachErosionCells > 0)
      OutStream << setw(7) << 1000 * dThisIterActualBeachErosion / static_cast<double>(m_ulThisIterNumActualBeachErosionCells);
   else
      OutStream << setw(7) << SPACE;

   // Per-iteration actual beach erosion in m (average for all sea cells)
   OutStream << std::fixed << setprecision(1);

   if (m_dThisIterBeachErosionFine > 0)
      OutStream << setw(4) << 1000 * m_dThisIterBeachErosionFine / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   if (m_dThisIterBeachErosionSand > 0)
      OutStream << setw(4) << 1000 * m_dThisIterBeachErosionSand / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   if (m_dThisIterBeachErosionCoarse > 0)
      OutStream << setw(4) << 1000 * m_dThisIterBeachErosionCoarse / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   // Output the this-timestep % of sea cells with beach deposition =============================================================
   OutStream << std::fixed << setprecision(0);
   OutStream << setw(7) << 100 * static_cast<double>(m_ulThisIterNumBeachDepositionCells) / static_cast<double>(m_ulThisIterNumSeaCells);

   // Per-iteration beach deposition in m (average for all sea cells)
   double dThisIterBeachDeposition = m_dThisIterBeachDepositionSand + m_dThisIterBeachDepositionCoarse;
   OutStream << setw(6) << 1000 * dThisIterBeachDeposition / static_cast<double>(m_ulThisIterNumSeaCells);

   // Per-iteration beach deposition in m (average for all cells with beach deposition)
   OutStream << std::fixed << setprecision(1);
   if (m_ulThisIterNumBeachDepositionCells > 0)
      OutStream << setw(9) << 1000 * dThisIterBeachDeposition / static_cast<double>(m_ulThisIterNumBeachDepositionCells);
   else
      OutStream << setw(9) << SPACE;

   // Per-iteration beach deposition in m (average for all sea cells)
   OutStream << std::fixed << setprecision(1);

   if (m_dThisIterBeachDepositionSand > 0)
      OutStream << setw(4) << 1000 * m_dThisIterBeachDepositionSand / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   if (m_dThisIterBeachDepositionCoarse > 0)
      OutStream << setw(4) << 1000 * m_dThisIterBeachDepositionCoarse / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   // Output the this-timestep sediment input in m ==============================================================================
   OutStream << std::fixed << setprecision(1);

   if (m_dThisiterUnconsFineInput > 0)
      OutStream << setw(4) << m_dThisiterUnconsFineInput;
   else
      OutStream << setw(4) << SPACE;

   if (m_dThisiterUnconsSandInput > 0)
      OutStream << setw(4) << m_dThisiterUnconsSandInput;
   else
      OutStream << setw(4) << SPACE;

   if (m_dThisiterUnconsCoarseInput > 0)
      OutStream << setw(4) << m_dThisiterUnconsCoarseInput;
   else
      OutStream << setw(4) << SPACE;

   // Per-iteration cliff collapse erosion (both cons and uncons) in m (average for all coast cells) ============================
   OutStream << std::fixed << setprecision(1);

   if ((m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons) > 0)
      OutStream << setw(4) << 1000 * (m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons) / static_cast<double>(m_ulThisIterNumCoastCells);
   else
      OutStream << setw(4) << SPACE;

   if ((m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons) > 0)
      OutStream << setw(4) << 1000 * (m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons) / static_cast<double>(m_ulThisIterNumCoastCells);
   else
      OutStream << setw(4) << SPACE;

   if ((m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) > 0)
      OutStream << setw(4) << 1000 * (m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) / static_cast<double>(m_ulThisIterNumCoastCells);
   else
      OutStream << setw(4) << SPACE;

   // Per-iteration cliff collapse deposition in m (average for all sea cells) ==================================================
   if (m_dThisIterUnconsSandCliffDeposition > 0)
      OutStream << setw(4) << 1000 * m_dThisIterUnconsSandCliffDeposition / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   if (m_dThisIterUnconsCoarseCliffDeposition > 0)
      OutStream << setw(4) << 1000 * m_dThisIterUnconsCoarseCliffDeposition / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(4) << SPACE;

   // Output per-timestep fine sediment going to suspension, in m (average for all sea cells) ==================================
   if (m_dThisIterFineSedimentToSuspension > 0)
      OutStream << setw(6) << 1000 * m_dThisIterFineSedimentToSuspension / static_cast<double>(m_ulThisIterNumSeaCells);
   else
      OutStream << setw(6) << SPACE;

   OutStream << " ";

   // Finally, set 'markers' for events that have occurred this timestep
   if (m_bSaveGISThisIter)
      OutStream << " GIS" << m_nGISSave;

   OutStream << endl;

   // Did a text file write error occur?
   if (OutStream.fail())
      return false;

   return true;
}

//===============================================================================================================================
//! Write the results for this timestep to the time series CSV files
//===============================================================================================================================
bool CSimulation::bWriteTSFiles(void)
{
   // Sea area
   if (m_bSeaAreaTSSave)
   {
      // Output in external CRS units
      SeaAreaTSStream << m_dSimElapsed << "\t,\t" << m_dExtCRSGridArea * static_cast<double>(m_ulThisIterNumSeaCells) / static_cast<double>(m_ulNumCells) << endl;

      // Did a time series file write error occur?
      if (SeaAreaTSStream.fail())
         return false;
   }

   // Still water level
   if (m_bStillWaterLevelTSSave)
   {
      // Output as is (m)
      StillWaterLevelTSStream << m_dSimElapsed << "\t,\t" << m_dThisIterSWL << endl;

      // Did a time series file write error occur?
      if (StillWaterLevelTSStream.fail())
         return false;
   }

   // Actual platform erosion (fine, sand, and coarse)
   if (m_bActualPlatformErosionTSSave)
   {
      // Output as is (m depth equivalent)
      PlatformErosionTSStream << m_dSimElapsed << "\t,\t" << m_dThisIterActualPlatformErosionFineCons << ",\t" << m_dThisIterActualPlatformErosionSandCons << ",\t" << m_dThisIterActualPlatformErosionCoarseCons << endl;

      // Did a time series file write error occur?
      if (PlatformErosionTSStream.fail())
         return false;
   }

   // Cliff collapse erosion (fine, sand, and coarse)
   if (m_bCliffCollapseErosionTSSave)
   {
      // Output as is (m depth equivalent)
      CliffCollapseErosionTSStream << m_dSimElapsed << "\t,\t" << m_dThisIterCliffCollapseErosionFineUncons << ",\t" << m_dThisIterCliffCollapseErosionSandUncons << ",\t" << m_dThisIterCliffCollapseErosionCoarseUncons << endl;

      // Did a time series file write error occur?
      if (CliffCollapseErosionTSStream.fail())
         return false;
   }

   // Cliff collapse deposition (sand and coarse)
   if (m_bCliffCollapseDepositionTSSave)
   {
      // Output as is (m depth equivalent)
      CliffCollapseDepositionTSStream << m_dSimElapsed << "\t,\t" << m_dThisIterUnconsSandCliffDeposition << ",\t" << m_dThisIterUnconsCoarseCliffDeposition << endl;

      // Did a time series file write error occur?
      if (CliffCollapseDepositionTSStream.fail())
         return false;
   }

   // Cliff collapse net
   if (m_bCliffCollapseNetTSSave)
   {
      // Output as is (m depth equivalent)
      CliffCollapseNetChangeTSStream << noshowpos << m_dSimElapsed << "\t,\t" << showpos << -m_dThisIterCliffCollapseFineErodedDuringDeposition + (m_dThisIterUnconsSandCliffDeposition - m_dThisIterCliffCollapseSandErodedDuringDeposition) + (m_dThisIterUnconsCoarseCliffDeposition - m_dThisIterCliffCollapseCoarseErodedDuringDeposition) << endl;

      // Did a time series file write error occur?
      if (CliffCollapseNetChangeTSStream.fail())
         return false;
   }

   // Beach erosion (fine, sand, and coarse)
   if (m_bBeachErosionTSSave)
   {
      // Output as is (m depth equivalent)
      BeachErosionTSStream << m_dSimElapsed << "\t,\t" << m_dThisIterBeachErosionFine << ",\t" << m_dThisIterBeachErosionSand << ",\t" << m_dThisIterBeachErosionCoarse << endl;

      // Did a time series file write error occur?
      if (BeachErosionTSStream.fail())
         return false;
   }

   // Beach deposition (sand and coarse)
   if (m_bBeachDepositionTSSave)
   {
      // Output as is (m depth equivalent)
      BeachDepositionTSStream << m_dSimElapsed << "\t,\t" << m_dThisIterBeachDepositionSand << ",\t" << m_dThisIterBeachDepositionCoarse << endl;

      // Did a time series file write error occur?
      if (BeachDepositionTSStream.fail())
         return false;
   }

   // Net change in beach sediment
   if (m_bBeachSedimentChangeNetTSSave)
   {
      // Output as is (m depth equivalent)
      BeachSedimentNetChangeTSStream << noshowpos << m_dSimElapsed << "\t,\t" << showpos << -m_dThisIterBeachErosionFine + (m_dThisIterBeachDepositionSand - m_dThisIterBeachErosionSand) + (m_dThisIterBeachDepositionCoarse - m_dThisIterBeachErosionCoarse) << endl;

      // Did a time series file write error occur?
      if (BeachSedimentNetChangeTSStream.fail())
         return false;
   }

   if (m_bSuspSedTSSave)
   {
      // Output as is (m depth equivalent)
      FineSedSuspensionTSStream << m_dSimElapsed << "\t,\t" << m_dThisIterFineSedimentToSuspension << endl;

      // Did a time series file write error occur?
      if (FineSedSuspensionTSStream.fail())
         return false;
   }

   if (m_bFloodSetupSurgeTSSave)
   {
      // Output as is (m depth equivalent)
      FloodSetupSurgeTSStream << m_dSimElapsed << "\t,\t" << m_dThisIterDiffWaveSetupSurgeWaterLevel << endl;

      // Did a time series file write error occur?
      if (FloodSetupSurgeTSStream.fail())
         return false;
   }

   if (m_bFloodSetupSurgeRunupTSSave)
   {
      // Output as is (m depth equivalent)
      FloodSetupSurgeRunupTSStream << m_dSimElapsed << "\t,\t" << m_dThisIterDiffWaveSetupSurgeRunupWaterLevel << endl;

      // Did a time series file write error occur?
      if (FloodSetupSurgeRunupTSStream.fail())
         return false;
   }

   return true;
}

//===============================================================================================================================
//! Output the erosion potential look-up values, for checking purposes
//===============================================================================================================================
void CSimulation::WriteLookUpData(void) const
{
   // Open the output file
   string strLookUpFile = m_strOutPath;
   strLookUpFile.append(EROSION_POTENTIAL_LOOKUP_FILE);
   ofstream LookUpOutStream;
   LookUpOutStream.open(strLookUpFile.c_str(), ios::out | ios::trunc);

   if (LookUpOutStream)
   {
      // File opened OK, so output the values
      LookUpOutStream << "DepthOverDB, \tErosionPotential" << endl;
      double dDepthOverDB = 0.0;
      while (dDepthOverDB <= m_dDepthOverDBMax)
      {
         double dErosionPotential = dGetInterpolatedValue(&m_VdDepthOverDB, &m_VdErosionPotential, dDepthOverDB, false);
         LookUpOutStream << dDepthOverDB << ",\t" << dErosionPotential << endl;
         dDepthOverDB += DEPTH_OVER_DB_INCREMENT;
      }
      LookUpOutStream << endl;

      // And close the file
      LookUpOutStream.close();
   }
}

//===============================================================================================================================
//! Save a coastline-normal profile
//===============================================================================================================================
int CSimulation::nSaveProfile(int const nProfile, int const nCoast, int const nProfSize, vector<double> const* pdVDistXY, vector<double> const* pdVZ, vector<double> const* pdVDepthOverDB, vector<double> const* pdVErosionPotentialFunc, vector<double> const* pdVSlope, vector<double> const* pdVRecessionXY, vector<double> const* pdVChangeElevZ, vector<CGeom2DIPoint> *const pPtVGridProfile, vector<double> const* pdVScapeXY) const
{
   // TODO 052 Make this more efficient, also give warnings if no profiles will be output
   for (unsigned int i = 0; i < m_VulProfileTimestep.size(); i++)
   {
      for (unsigned int j = 0; j < m_VnProfileToSave.size(); j++)
      {
         if ((m_ulIter == m_VulProfileTimestep[i]) && (nProfile == m_VnProfileToSave[j]))
         {
            if (! bWriteProfileData(nCoast, nProfile, nProfSize, pdVDistXY, pdVZ, pdVDepthOverDB, pdVErosionPotentialFunc, pdVSlope, pdVRecessionXY, pdVChangeElevZ, pPtVGridProfile, pdVScapeXY))
               return RTN_ERR_PROFILEWRITE;
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Writes values for a single profile, for checking purposes
//===============================================================================================================================
bool CSimulation::bWriteProfileData(int const nCoast, int const nProfile, int const nProfSize, vector<double> const* pdVDistXY, vector<double> const* pdVZ, vector<double> const* pdVDepthOverDB, vector<double> const* pdVErosionPotentialFunc, vector<double> const* pdVSlope, vector<double> const* pdVRecessionXY, vector<double> const* pdVChangeElevZ, vector<CGeom2DIPoint> *const pPtVGridProfile, vector<double> const* pdVScapeXY) const
{
   string strFName = m_strOutPath;
   stringstream ststrTmp;

   strFName.append("profile_");
   ststrTmp << FillToWidth('0', 3) << nProfile;
   strFName.append(ststrTmp.str());

   strFName.append("_timestep_");
   ststrTmp.clear();
   ststrTmp.str(string());
   ststrTmp << FillToWidth('0', 4) << m_ulIter;
   strFName.append(ststrTmp.str());

   strFName.append(".csv");

   ofstream OutProfStream;
   OutProfStream.open(strFName.c_str(), ios::out | ios::trunc);
   if (!OutProfStream)
   {
      // Error, cannot open file
      cerr << ERR << "cannot open " << strFName << " for output" << endl;
      return false;
   }

   OutProfStream << "\"Dist\", \"X\", \"Y\", \"Z (before erosion)\", \"Depth/DB\", \"Erosion Potential\", \"Slope\", \"Recession XY\", \"Change Elev Z\", \"Grid X\",  \"Grid Y\",  \"Weight\",  \"For profile " << nProfile << " from coastline " << nCoast << " at timestep " << m_ulIter << "\"" << endl;
   for (int i = 0; i < nProfSize; i++)
   {
      double dX = dGridCentroidXToExtCRSX(pPtVGridProfile->at(i).nGetX());
      double dY = dGridCentroidYToExtCRSY(pPtVGridProfile->at(i).nGetY());

      OutProfStream << pdVDistXY->at(i) << ",\t" << dX << ",\t" << dY << ",\t" << pdVZ->at(i) << ",\t" << pdVDepthOverDB->at(i) << ",\t" << pdVErosionPotentialFunc->at(i) << ",\t" << pdVSlope->at(i) << ",\t" << pdVRecessionXY->at(i) << ",\t" << pdVChangeElevZ->at(i) << ",\t" << pPtVGridProfile->at(i).nGetX() << ",\t" << pPtVGridProfile->at(i).nGetY() << ", \t" << pdVScapeXY->at(i) << endl;
   }

   OutProfStream.close();

   return true;
}

//===============================================================================================================================
//! Save a coastline-normal parallel profile
//===============================================================================================================================
int CSimulation::nSaveParProfile(int const nProfile, int const nCoast, int const nParProfSize, int const nDirection, int const nDistFromProfile, vector<double> const* pdVDistXY, vector<double> const* pdVZ, vector<double> const* pdVDepthOverDB, vector<double> const* pdVErosionPotentialFunc, vector<double> const* pdVSlope, vector<double> const* pdVRecessionXY, vector<double> const* pdVChangeElevZ, vector<CGeom2DIPoint> *const pPtVGridProfile, vector<double> const* pdVScapeXY) const
{
   // TODO 052 Make this more efficient, also give warnings if no profiles will be output
   for (unsigned int i = 0; i < m_VulProfileTimestep.size(); i++)
   {
      for (unsigned int j = 0; j < m_VnProfileToSave.size(); j++)
      {
         if ((m_ulIter == m_VulProfileTimestep[i]) && (nProfile == m_VnProfileToSave[j]))
         {
            if (! bWriteParProfileData(nCoast, nProfile, nParProfSize, nDirection, nDistFromProfile, pdVDistXY, pdVZ, pdVDepthOverDB, pdVErosionPotentialFunc, pdVSlope, pdVRecessionXY, pdVChangeElevZ, pPtVGridProfile, pdVScapeXY))
               return RTN_ERR_PROFILEWRITE;
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Writes values for a single parallel profile, for checking purposes
//===============================================================================================================================
bool CSimulation::bWriteParProfileData(int const nCoast, int const nProfile, int const nProfSize, int const nDirection, int const nDistFromProfile, vector<double> const* pdVDistXY, vector<double> const* pdVZ, vector<double> const* pdVDepthOverDB, vector<double> const* pdVErosionPotentialFunc, vector<double> const* pdVSlope, vector<double> const* pdVRecessionXY, vector<double> const* pdVChangeElevZ, vector<CGeom2DIPoint> *const pPtVGridProfile, vector<double> const* pdVScapeXY) const
{
   string strFName = m_strOutPath;
   stringstream ststrTmp;

   strFName.append("profile_");
   ststrTmp << FillToWidth('0', 3) << nProfile;
   strFName.append(ststrTmp.str());

   strFName.append("_parallel_");
   ststrTmp.clear();
   ststrTmp.str(string());
   ststrTmp << FillToWidth('0', 3) << nDistFromProfile;
   strFName.append(ststrTmp.str());

   strFName.append((nDirection == 0 ? "_F" : "_B"));

   strFName.append("_timestep_");
   ststrTmp.clear();
   ststrTmp.str(string());
   ststrTmp << FillToWidth('0', 4) << m_ulIter;
   strFName.append(ststrTmp.str());

   strFName.append(".csv");

   ofstream OutProfStream;
   OutProfStream.open(strFName.c_str(), ios::out | ios::trunc);
   if (!OutProfStream)
   {
      // Error, cannot open file
      cerr << ERR << "cannot open " << strFName << " for output" << endl;
      return false;
   }

   OutProfStream << "\"Dist\", \"X\", \"Y\", \"Z (before erosion)\", \"Depth/DB\", \"Erosion Potential\", \"Slope\", \"Recession XY\", \"Change Elev Z\", \"Grid X\",  \"Grid Y\",  \"Weight\",  \"For profile " << nProfile << " from coastline " << nCoast << " at timestep " << m_ulIter << "\"" << endl;
   for (int i = 0; i < nProfSize; i++)
   {
      double dX = dGridCentroidXToExtCRSX(pPtVGridProfile->at(i).nGetX());
      double dY = dGridCentroidYToExtCRSY(pPtVGridProfile->at(i).nGetY());

      OutProfStream << pdVDistXY->at(i) << ",\t" << dX << ",\t" << dY << ",\t" << pdVZ->at(i) << ",\t" << pdVDepthOverDB->at(i) << ",\t" << pdVErosionPotentialFunc->at(i) << ",\t" << pdVSlope->at(i) << ",\t" << pdVRecessionXY->at(i) << ",\t" << pdVChangeElevZ->at(i) << ",\t" << pPtVGridProfile->at(i).nGetX() << ",\t" << pPtVGridProfile->at(i).nGetY() << ", \t" << pdVScapeXY->at(i) << endl;
   }

   OutProfStream.close();

   return true;
}

//===============================================================================================================================
//! Writes end-of-run information to Out, Log and time-series files
//===============================================================================================================================
int CSimulation::nWriteEndRunDetails(void)
{
   // Final write to time series CSV files
   if (! bWriteTSFiles())
      return (RTN_ERR_TIMESERIES_FILE_WRITE);

   // Save the values from the RasterGrid array into raster GIS files
   if (! bSaveAllRasterGISFiles())
      return (RTN_ERR_RASTER_FILE_WRITE);

   // Save the vector GIS files
   if (! bSaveAllVectorGISFiles())
      return (RTN_ERR_VECTOR_FILE_WRITE);

   OutStream << " GIS" << m_nGISSave << endl;

   // Print out run totals etc.
   OutStream << PER_ITER_HEAD1 << endl;
   OutStream << PER_ITER_HEAD2 << endl;
   OutStream << PER_ITER_HEAD3 << endl;
   OutStream << PER_ITER_HEAD4 << endl;
   OutStream << PER_ITER_HEAD5 << endl;

   OutStream << std::fixed << setprecision(3);
   OutStream << endl
             << endl;

   // Write out hydrology grand totals etc.
   OutStream << ENDHYDROLOGYHEAD << endl;
   OutStream << "Minimum still water level = " << m_dMinSWL << endl;
   OutStream << "Maximum still water level = " << m_dMaxSWL << endl;
   OutStream << endl;

   // Now write out sediment movement grand totals etc.
   OutStream << ENDSEDIMENTHEAD << endl
             << endl;

   OutStream << "TOTAL PLATFORM EROSION" << endl;
   OutStream << "Potential platform erosion, all size classes           = " << m_ldGTotPotentialPlatformErosion * m_dCellArea << " m^3" << endl
             << endl;
   OutStream << "Actual platform erosion, fine                          = " << m_ldGTotFineActualPlatformErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual platform erosion, sand                          = " << m_ldGTotSandActualPlatformErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual platform erosion, coarse                        = " << m_ldGTotCoarseActualPlatformErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual platform erosion, all size classes              = " << (m_ldGTotFineActualPlatformErosion + m_ldGTotSandActualPlatformErosion + m_ldGTotCoarseActualPlatformErosion) * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "TOTAL CLIFF COLLAPSE EROSION" << endl;
   OutStream << "Cliff collapse, fine                                   = " << m_ldGTotCliffCollapseFine * m_dCellArea << " m^3" << endl;
   OutStream << "Cliff collapse, sand                                   = " << m_ldGTotCliffCollapseSand * m_dCellArea << " m^3" << endl;
   OutStream << "Cliff collapse, coarse                                 = " << m_ldGTotCliffCollapseCoarse * m_dCellArea << " m^3" << endl;
   OutStream << "Cliff collapse, all size classes                       = " << (m_ldGTotCliffCollapseFine + m_ldGTotCliffCollapseSand + m_ldGTotCliffCollapseCoarse + m_ldGTotCliffCollapseFineErodedDuringDeposition + m_ldGTotCliffCollapseSandErodedDuringDeposition + m_ldGTotCliffCollapseCoarseErodedDuringDeposition) * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "TOTAL DEPOSITION AND SUSPENSION OF CLIFF COLLAPSE TALUS" << endl;
   OutStream << "Cliff collapse to suspension, fine                     = " << m_ldGTotCliffTalusFineToSuspension * m_dCellArea << " m^3" << endl;
   OutStream << "Cliff collapse deposition, sand                        = " << m_ldGTotCliffTalusSandDeposition * m_dCellArea << " m^3" << endl;
   OutStream << "Cliff collapse deposition, coarse                      = " << m_ldGTotCliffTalusCoarseDeposition * m_dCellArea << " m^3" << endl;
   OutStream << "Cliff collapse deposition, sand and coarse             = " << (m_ldGTotCliffTalusSandDeposition + m_ldGTotCliffTalusCoarseDeposition) * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "TOTAL BEACH EROSION" << endl;
   OutStream << "Potential beach erosion, all size classes              = " << m_ldGTotPotentialBeachErosion * m_dCellArea << " m^3" << endl
             << endl;
   OutStream << "Actual fine beach erosion, fine                        = " << m_ldGTotActualFineBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual sand beach erosion, sand                        = " << m_ldGTotActualSandBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual coarse beach erosion, coarse                    = " << m_ldGTotActualCoarseBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual beach erosion, all size classes                 = " << (m_ldGTotActualFineBeachErosion + m_ldGTotActualSandBeachErosion + m_ldGTotActualCoarseBeachErosion) * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "TOTAL BEACH DEPOSITION" << endl;
   OutStream << "Beach deposition, sand                                 = " << m_ldGTotSandBeachDeposition * m_dCellArea << " m^3" << endl;
   OutStream << "Beach deposition, coarse                               = " << m_ldGTotCoarseBeachDeposition * m_dCellArea << " m^3" << endl;
   OutStream << "Beach deposition, sand and coarse                      = " << (m_ldGTotSandBeachDeposition + m_ldGTotCoarseBeachDeposition) * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "TOTAL SEDIMENT INPUT EVENTS" << endl;
   OutStream << "Sediment from sediment input events, fine              = " << m_ldGTotFineSedimentInput * m_dCellArea << " m^3" << endl;
   OutStream << "Sediment from sediment input events, sand              = " << m_ldGTotSandSedimentInput * m_dCellArea << " m^3" << endl;
   OutStream << "Sediment from sediment input events, coarse            = " << m_ldGTotCoarseSedimentInput * m_dCellArea << " m^3" << endl;
   OutStream << "Sediment from sediment input events, all size classes  = " << (m_ldGTotFineSedimentInput + m_ldGTotSandSedimentInput + m_ldGTotCoarseSedimentInput) * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "TOTAL SUSPENDED SEDIMENT" << endl;
   OutStream << "Suspended fine sediment                                = " << m_ldGTotSuspendedSediment * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "TOTAL LOST FROM GRID BY BEACH MOVEMENT" << endl;
   OutStream << "Potential sediment lost, all size classes              = " << m_ldGTotPotentialSedLostBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual sediment lost, fine                             = " << m_ldGTotActualFineLostBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual sediment lost, sand                             = " << m_ldGTotActualSandLostBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual sediment lost, coarse                           = " << m_ldGTotActualCoarseLostBeachErosion * m_dCellArea << " m^3" << endl;
   OutStream << "Actual sediment lost, all size classes                 = " << (m_ldGTotActualFineLostBeachErosion + m_ldGTotActualSandLostBeachErosion + m_ldGTotActualCoarseLostBeachErosion) * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "TOTAL LOST FROM GRID BY CLIFF COLLAPSE" << endl;
   OutStream << "Sediment lost, sand                                    = " << m_ldGTotSandSedLostCliffCollapse * m_dCellArea << " m^3" << endl;
   OutStream << "Sediment lost, coarse                                  = " << m_ldGTotCoarseSedLostCliffCollapse * m_dCellArea << " m^3" << endl;
   OutStream << endl;

   OutStream << "ALL-PROCESS TOTALS (all size classes)" << endl;
   long double ldFineEroded = m_ldGTotFineActualPlatformErosion + m_ldGTotCliffCollapseFine + m_ldGTotActualFineBeachErosion;
   OutStream << "Fine sediment eroded                                   = " << ldFineEroded * m_dCellArea << " m^3" << endl;
   OutStream << "Fine sediment to suspension                            = " << m_ldGTotSuspendedSediment * m_dCellArea << " m^3" << endl;
   if (! bFPIsEqual(ldFineEroded, m_ldGTotSuspendedSediment, 1.0L))
      OutStream << MASS_BALANCE_ERROR << endl;
   
   long double ldSandEroded = m_ldGTotSandActualPlatformErosion + m_ldGTotCliffCollapseSand + m_ldGTotActualSandBeachErosion;
   OutStream << "Sand sediment eroded                                   = " << ldSandEroded * m_dCellArea << " m^3" << endl;
   long double ldSandDeposited = m_ldGTotCliffTalusSandDeposition + m_ldGTotSandBeachDeposition;
   OutStream << "Sand sediment deposited                                = " << ldSandDeposited * m_dCellArea << " m^3" << endl;
   long double ldSandLost = m_ldGTotActualSandLostBeachErosion + m_ldGTotSandSedLostCliffCollapse;
   OutStream << "Sand sediment lost from grid                           = " << ldSandLost * m_dCellArea << " m^3" << endl;
   if (! bFPIsEqual(ldSandEroded, (ldSandDeposited + ldSandLost), 1.0L))
      OutStream << MASS_BALANCE_ERROR << endl;
   
   long double ldCoarseEroded = m_ldGTotCoarseActualPlatformErosion + m_ldGTotCliffCollapseCoarse + m_ldGTotActualCoarseBeachErosion;
   OutStream << "Coarse sediment eroded                                 = " << ldCoarseEroded * m_dCellArea << " m^3" << endl;
   long double ldCoarseDeposited = m_ldGTotCliffTalusCoarseDeposition + m_ldGTotCoarseBeachDeposition;
   OutStream << "Coarse sediment deposited                              = " << ldCoarseDeposited * m_dCellArea << " m^3" << endl;
   long double ldCoarseLost = m_ldGTotActualCoarseLostBeachErosion + m_ldGTotCoarseSedLostCliffCollapse;
   OutStream << "Coarse sediment lost from grid                         = " << ldCoarseLost * m_dCellArea << " m^3" << endl;
   if (! bFPIsEqual(ldCoarseEroded, (ldCoarseDeposited + ldCoarseLost), 1.0L))
      OutStream << MASS_BALANCE_ERROR << endl;
   
   OutStream << endl;
   
   long double ldActualTotalEroded = m_ldGTotFineActualPlatformErosion + m_ldGTotSandActualPlatformErosion + m_ldGTotCoarseActualPlatformErosion + m_ldGTotCliffCollapseFine + m_ldGTotCliffCollapseSand + m_ldGTotCliffCollapseCoarse + m_ldGTotCliffCollapseFineErodedDuringDeposition + m_ldGTotCliffCollapseSandErodedDuringDeposition + m_ldGTotCliffCollapseCoarseErodedDuringDeposition + m_ldGTotActualFineBeachErosion + m_ldGTotActualSandBeachErosion + m_ldGTotActualCoarseBeachErosion;
   OutStream << "Total sediment eroded (all processes)                  = " << ldActualTotalEroded * m_dCellArea << " m^3" << endl;

   long double ldTotalDepositedAndSuspension = m_ldGTotCliffTalusSandDeposition + m_ldGTotCliffTalusCoarseDeposition + m_ldGTotSandBeachDeposition + m_ldGTotCoarseBeachDeposition + m_ldGTotSuspendedSediment;
   OutStream << "Total sediment deposited/to suspension (all processes) = " << ldTotalDepositedAndSuspension * m_dCellArea << " m^3" << endl;

   long double ldTotalLost = m_ldGTotActualFineLostBeachErosion + m_ldGTotActualSandLostBeachErosion + m_ldGTotActualCoarseLostBeachErosion + m_ldGTotSandSedLostCliffCollapse + m_ldGTotCoarseSedLostCliffCollapse;
   OutStream << "Total sediment lost from grid (all processes)          = " << ldTotalLost * m_dCellArea << " m^3" << endl;
   OutStream << "                                                       = " << 24 * 365.25 * ldTotalLost * m_dCellArea / m_dSimDuration << " m^3/year" << endl;
   OutStream << "                                                       = " << 24 * ldTotalLost * m_dCellArea / m_dSimDuration << " m^3/day" << endl;
   OutStream << "                                                       = " << ldTotalLost * m_dCellArea / m_dSimDuration << " m^3/hour" << endl;
   OutStream << std::fixed << setprecision(6);
   OutStream << "                                                       = " << ldTotalLost * m_dCellArea / (m_dSimDuration * 3600) << " m^3/sec" << endl
             << endl;
   OutStream << std::fixed << setprecision(3);

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   {
      OutStream << "Grid edge option is ";
      if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
         OutStream << "CLOSED.";
      else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
         OutStream << "OPEN.";
      else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
         OutStream << "RE-CIRCULATING.";
      OutStream << endl << endl;
   }

   // Finally calculate performance details
   OutStream << PERFORMHEAD << endl;

   // Get the time that the run ended
   m_tSysEndTime = time(nullptr);

   OutStream << "Run ended at " << put_time(localtime(&m_tSysEndTime), "%T on %A %d %B %Y") << endl;
   OutStream << "Time simulated: " << strDispSimTime(m_dSimDuration) << endl
             << endl;

   // Write to log file
   LogStream << "END OF RUN ================================================================================================" << endl
             << endl;

   LogStream << "ALL-PROCESS TOTALS (all size classes)" << endl;
   LogStream << "Sediment added                                           = " << (m_ldGTotFineSedimentInput + m_ldGTotSandSedimentInput + m_ldGTotCoarseSedimentInput) * m_dCellArea << " m^3" << endl;
   LogStream << "Sediment eroded (all processes)                          = " << ldActualTotalEroded * m_dCellArea << " m^3" << endl;

   LogStream << "Sediment deposited and in suspension (all processes)     = " << ldTotalDepositedAndSuspension * m_dCellArea << " m^3" << endl;

   LogStream << "Sediment lost from grid (all processes)                  = " << ldTotalLost * m_dCellArea << " m^3" << endl;
   LogStream << "                                                         = " << 24 * 365.25 * ldTotalLost * m_dCellArea / m_dSimDuration << " m^3/year" << endl;
   LogStream << "                                                         = " << 24 * ldTotalLost * m_dCellArea / m_dSimDuration << " m^3/day" << endl;
   LogStream << "                                                         = " << ldTotalLost * m_dCellArea / m_dSimDuration << " m^3/hour" << endl;
   LogStream << "                                                         = " << setprecision(6) << ldTotalLost * m_dCellArea / (m_dSimDuration * 3600) << " m^3/sec" << endl;
   LogStream << endl;
   LogStream << std::fixed << setprecision(3);

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   {
      LogStream << "Grid edge option is ";
      if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_CLOSED)
         LogStream << "CLOSED.";
      else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_OPEN)
         LogStream << "OPEN.";
      else if (m_nUnconsSedimentHandlingAtGridEdges == GRID_EDGE_RECIRCULATE)
         LogStream << "RE-CIRCULATING.";
      LogStream << endl << endl;

      // Output averages for on-profile and between-profile potential shore platform erosion, ideally these are roughly equal
      LogStream << std::fixed << setprecision(6);
      LogStream << "On-profile average potential shore platform erosion      = " << (m_ulTotPotentialPlatformErosionOnProfiles > 0 ? m_dTotPotentialPlatformErosionOnProfiles / static_cast<double>(m_ulTotPotentialPlatformErosionOnProfiles) : 0) << " mm (n = " << m_ulTotPotentialPlatformErosionOnProfiles << ")" << endl;
      LogStream << "Between-profile average potential shore platform erosion = " << (m_ulTotPotentialPlatformErosionBetweenProfiles > 0 ? m_dTotPotentialPlatformErosionBetweenProfiles / static_cast<double>(m_ulTotPotentialPlatformErosionBetweenProfiles) : 0) << " mm (n = " << m_ulTotPotentialPlatformErosionBetweenProfiles << ")" << endl;
      LogStream << endl;
   }
   
#if !defined RANDCHECK
   // Calculate length of run, write in file (note that m_dSimDuration is in hours)
   CalcTime(m_dSimDuration * 3600);
#endif

   // Calculate statistics re. memory usage etc.
   CalcProcessStats();
   OutStream << endl
             << "END OF RUN" << endl;
   LogStream << endl
             << "END OF RUN" << endl;

   // Need to flush these here (if we don't, the buffer may not get written)
   LogStream.flush();
   OutStream.flush();

   return RTN_OK;
}

//===============================================================================================================================
//! Writes to the log file a table showing polygon to polygon shares of unconsolidated sediment transport, etc.
//===============================================================================================================================
void CSimulation::WritePolygonShareTable(int const nCoast)
{
   LogStream << endl << m_ulIter << ": Polygon-to-adjacent polygon shares (non-dimensional), per-polygon seawater volume (m^3), per-polygon D50 values (mm: a blank D50 value means that there is no unconsolidated sediment on that polygon)." << endl;
   
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------------------------------------" << endl;
   LogStream << strCentre("Polygon", 11) << "|" << strCentre("Coast", 11) << "|" << strCentre("Polygon", 11) << "|" << strCentre("Seawater", 14) << "|" << strCentre("Uncons d50", 14) << "| " << strCentre("(Dir'n Adj Share)...", 14) << endl;
   LogStream << strCentre("Global ID", 11) << "|" << strCentre("", 11) << "|" << strCentre("Coast ID", 11) << "|" << strCentre("Volume", 14) << "|" << strCentre("", 14) << "| " << strCentre("", 14) << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------------------------------------" << endl;

   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      LogStream << strIntRight(m_pVCoastPolygon[n]->nGetGlobalID(), 11) << "|" << strIntRight(nCoast, 11) << "|" << strIntRight(m_pVCoastPolygon[n]->nGetCoastID(), 11) << "|" << strDblRight(m_pVCoastPolygon[n]->dGetSeawaterVolume(), 0, 14) << "|" << strDblRight(m_pVCoastPolygon[n]->dGetAvgUnconsD50(), 0, 14) << "| ";

      for (int m = 0; m < m_pVCoastPolygon[n]->nGetNumUpCoastAdjacentPolygons(); m++)
      {
         if (! m_pVCoastPolygon[n]->bDownCoastThisIter())
            LogStream << "(UP  \t" << m_pVCoastPolygon[n]->nGetUpCoastAdjacentPolygon(m) << "\t" << m_pVCoastPolygon[n]->dGetUpCoastAdjacentPolygonBoundaryShare(m) << ")\t";
      }

      for (int m = 0; m < m_pVCoastPolygon[n]->nGetNumDownCoastAdjacentPolygons(); m++)
      {
         if (m_pVCoastPolygon[n]->bDownCoastThisIter())
            LogStream << "(DOWN\t" << m_pVCoastPolygon[n]->nGetDownCoastAdjacentPolygon(m) << "\t" << m_pVCoastPolygon[n]->dGetDownCoastAdjacentPolygonBoundaryShare(m) << ")\t";
      }
      LogStream << endl;
   }
   
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------------------------------------" << endl << endl;
}

//===============================================================================================================================
//! Writes to the log file a table showing per-polygon pre-existing unconsolidated sediment
//===============================================================================================================================
void CSimulation::WritePolygonPreExistingSediment(int const nCoast)
{
   double dTmpTot = 0;
   double dTmpFineTot = 0;
   double dTmpSandTot = 0;
   double dTmpCoarseTot = 0;

   m_dTotalFineConsInPolygons =
   m_dTotalSandConsInPolygons =
   m_dTotalCoarseConsInPolygons =
   m_dTotalFineUnconsInPolygons =
   m_dTotalSandUnconsInPolygons =
   m_dTotalCoarseUnconsInPolygons = 0;

   // ALSO PRINT m_dStartIterUnconsFineAllCells etc.

   LogStream << m_ulIter << ": Per-polygon pre-existing unconsolidated sediment. Note that this does not include pre-existing unconsolidated sediment outside the polygons.";
   if (m_ulIter > 1)
      LogStream << " Note that the all-polygon total will be slightly different from the all-polygon total at the end of the last timestep, since the coastline and polygons have been re-drawn.";
   LogStream << endl;
   
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << strCentre("Polygon", 11) << "|" << strCentre("Coast", 11) << "|" << strCentre("Polygon", 11) << "|" << strCentre("All", 14) << "|" << strCentre("Fine", 14) << "|" << strCentre("Sand", 14) << "|" << strCentre("Coarse", 14) << "|" << endl;
   LogStream << strCentre("Global ID", 11) << "|" << strCentre("", 11) << "|" << strCentre("Coast ID", 11) << "|" << strCentre("Sediment", 14) << "|" << strCentre("Sediment", 14) << "|" << strCentre("Sediment", 14) << "|" << strCentre("Sediment", 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
   
   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      // Add to this-iteration all-polygon totals of consolidated sediment within polygons
      m_dTotalFineConsInPolygons += m_pVCoastPolygon[n]->dGetPreExistingConsFine();
      m_dTotalSandConsInPolygons += m_pVCoastPolygon[n]->dGetPreExistingConsSand();
      m_dTotalCoarseConsInPolygons += m_pVCoastPolygon[n]->dGetPreExistingConsCoarse();

      // Now consider unconsolidated sediment
      double dThisFine = m_pVCoastPolygon[n]->dGetPreExistingUnconsFine();
      double dThisSand = m_pVCoastPolygon[n]->dGetPreExistingUnconsSand();
      double dThisCoarse = m_pVCoastPolygon[n]->dGetPreExistingUnconsCoarse();

      LogStream << strIntRight(m_pVCoastPolygon[n]->nGetGlobalID(), 11) << "|" << strIntRight(nCoast, 11) << "|" << strIntRight(m_pVCoastPolygon[n]->nGetCoastID(), 11) << "|" << strDblRight((dThisFine + dThisSand + dThisCoarse) * m_dCellArea, 0, 14) << "|" << strDblRight(dThisFine * m_dCellArea, 0, 14) << "|" <<  strDblRight(dThisSand * m_dCellArea, 0, 14) << "|" << strDblRight(dThisCoarse * m_dCellArea, 0, 14) << "|" << endl;

      dTmpFineTot += (dThisFine * m_dCellArea);
      dTmpSandTot += (dThisSand * m_dCellArea);
      dTmpCoarseTot += (dThisCoarse * m_dCellArea);
      dTmpTot += (dThisFine + dThisSand + dThisCoarse) * m_dCellArea;

      // Add to this-iteration all-polygon totals of unconsolidated sediment within polygons
      m_dTotalFineUnconsInPolygons += dThisFine;
      m_dTotalSandUnconsInPolygons += dThisSand;
      m_dTotalCoarseUnconsInPolygons += dThisCoarse;
   }

   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << "TOTAL                              |" << strDblRight(dTmpTot, 0, 14) << "|" << strDblRight(dTmpFineTot, 0, 14) << "|" << strDblRight(dTmpSandTot, 0, 14) << "|" << strDblRight(dTmpCoarseTot, 0, 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl << endl;
}

//===============================================================================================================================
//! Writes to the log file a table showing per-polygon unconsolidated sand/coarse sediment derived from erosion of the consolidated shore platform
//===============================================================================================================================
void CSimulation::WritePolygonShorePlatformErosion(int const nCoast)
{
   double dTmpTot = 0;
   double dTmpFineTot = 0;
   double dTmpSandTot = 0;
   double dTmpCoarseTot = 0;

   LogStream << endl << m_ulIter << ": Per-polygon unconsolidated sand/coarse sediment derived from erosion of the consolidated shore platform (all m^3). All fine sediment eroded from the shore platform goes to suspension." << endl;
   
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << strCentre("Polygon", 11) << "|" << strCentre("Coast", 11) << "|" << strCentre("Polygon", 11) << "|" << strCentre("All", 14) << "|" << strCentre("Fine", 14) << "|" << strCentre("Sand", 14) << "|" << strCentre("Coarse", 14) << "|" << endl;
   LogStream << strCentre("Global ID", 11) << "|" << strCentre("", 11) << "|" << strCentre("Coast ID", 11) << "|" << strCentre("Sediment", 14) << "|" << strCentre("Sediment", 14) << "|" << strCentre("Sediment", 14) << "|" << strCentre("Sediment", 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
   
   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      LogStream << strIntRight(m_pVCoastPolygon[n]->nGetGlobalID(), 11) << "|" << strIntRight(nCoast, 11) << "|" << strIntRight(m_pVCoastPolygon[n]->nGetCoastID(), 11) << "|" << strDblRight((m_pVCoastPolygon[n]->dGetPlatformErosionUnconsSand() + m_pVCoastPolygon[n]->dGetPlatformErosionUnconsCoarse()) * m_dCellArea, 0, 14) << "|" << strDblRight(0, 0, 14) << "|" <<  strDblRight(m_pVCoastPolygon[n]->dGetPlatformErosionUnconsSand() * m_dCellArea, 0, 14) << "|" << strDblRight(m_pVCoastPolygon[n]->dGetPlatformErosionUnconsCoarse() * m_dCellArea, 0, 14) << "|" << endl;
      
      dTmpTot += (m_pVCoastPolygon[n]->dGetPlatformErosionUnconsSand() + m_pVCoastPolygon[n]->dGetPlatformErosionUnconsCoarse()) * m_dCellArea;
      dTmpSandTot += (m_pVCoastPolygon[n]->dGetPlatformErosionUnconsSand() * m_dCellArea);
      dTmpCoarseTot += (m_pVCoastPolygon[n]->dGetPlatformErosionUnconsCoarse() * m_dCellArea);
   }
   
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << "TOTAL from shore platform          |" << strDblRight(dTmpTot, 0, 14) << "|" << strDblRight(dTmpFineTot, 0, 14) << "|" << strDblRight(dTmpSandTot, 0, 14) << "|" << strDblRight(dTmpCoarseTot, 0, 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl << endl;
}

//===============================================================================================================================
//! Writes to the log file a table showing per-polygon per-polygon cliff collapse
//===============================================================================================================================
void CSimulation::WritePolygonCliffCollapseErosion(int const nCoast)
{
   LogStream << endl << m_ulIter << ": Per-polygon cliff collapse (all m^3). Fine sediment derived from cliff collapse goes to suspension, sand/coarse sediment derived from cliff collapse becomes unconsolidated talus (DDPD = During Dean Profile Deposition)." << endl;
   
   LogStream << "-----------|-----------|-----------|--------------------------------------------|-----------------------------|--------------------------------------------|--------------------------------------------|" << endl;
   LogStream << strCentre("Polygon", 11) << "|" << strCentre("Coast", 11) << "|" << strCentre("Polygon", 11) << "|" << strCentre("All sediment", 44) << "|" << strCentre("Fine sediment", 29) << "|" << strCentre("Sand sediment", 44) << "|" << strCentre("Coarse sediment", 44) << "|" << endl;
   LogStream << strCentre("Global ID", 11) << "|" << strCentre("", 11) << "|" << strCentre("Coast ID", 11) << "|" << strCentre("Eroded Cliff", 15) << strCentre("Eroded DDPD", 15) << strCentre("Deposited", 14) << "|" << strCentre("Eroded", 14) << "|" << strCentre("Suspension", 14) << "|" << strCentre("Eroded Cliff", 15) << strCentre("Eroded DDPD", 15) << strCentre("Deposited", 14) << "|" << strCentre("Eroded Cliff", 15) << strCentre("Eroded DDPD", 15) << strCentre("Deposited", 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
   
   double dTmpErosionTot = 0;
   double dTmpErosionDDPDTot = 0;
   double dTmpDepositTot = 0;

   double dTmpErosionFineTot = 0;
   double dTmpSuspensionFineTot = 0;

   double dTmpErosionSandTot = 0;
   double dTmpErosionSandDDPDTot = 0;
   double dTmpDepositSandTot = 0;

   double dTmpErosionCoarseTot = 0;
   double dTmpErosionCoarseDDPDTot = 0;
   double dTmpDepositCoarseTot = 0;
   
   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      LogStream << strIntRight(m_pVCoastPolygon[n]->nGetGlobalID(), 11) << "|" << strIntRight(nCoast, 11) << "|" << strIntRight(m_pVCoastPolygon[n]->nGetCoastID(), 11)
         // All
         << "|" << strDblRight((m_pVCoastPolygon[n]->dGetCliffCollapseErosionFine() + m_pVCoastPolygon[n]->dGetCliffCollapseErosionSand() + m_pVCoastPolygon[n]->dGetCliffCollapseErosionCoarse()) * m_dCellArea, 0, 14) << " " << strDblRight((m_pVCoastPolygon[n]->dGetCliffCollapseSandErodedDeanProfile() + m_pVCoastPolygon[n]->dGetCliffCollapseCoarseErodedDeanProfile()) * m_dCellArea, 0, 14) << "|" << strDblRight((m_pVCoastPolygon[n]->dGetCliffCollapseToSuspensionFine() + m_pVCoastPolygon[n]->dGetCliffCollapseUnconsSandDeposition() + m_pVCoastPolygon[n]->dGetCliffCollapseUnconsCoarseDeposition()) * m_dCellArea, 0, 14) << "|"
         // Fine
         << strDblRight(m_pVCoastPolygon[n]->dGetCliffCollapseErosionFine(), 0, 14) << "|" << strDblRight(m_pVCoastPolygon[n]->dGetCliffCollapseToSuspensionFine(), 0, 14) << "|"
         // Sand
         << strDblRight(m_pVCoastPolygon[n]->dGetCliffCollapseErosionSand(), 0, 14) << " " << strDblRight(m_pVCoastPolygon[n]->dGetCliffCollapseSandErodedDeanProfile() * m_dCellArea, 0, 14) << "|" <<  strDblRight(m_pVCoastPolygon[n]->dGetCliffCollapseUnconsSandDeposition() * m_dCellArea, 0, 14) << "|"
         // Coarse
         << strDblRight(m_pVCoastPolygon[n]->dGetCliffCollapseErosionCoarse(), 0, 14) << " " << strDblRight(m_pVCoastPolygon[n]->dGetCliffCollapseCoarseErodedDeanProfile() * m_dCellArea, 0, 14)<< "|" << strDblRight(m_pVCoastPolygon[n]->dGetCliffCollapseUnconsCoarseDeposition() * m_dCellArea, 0, 14) << "|" << endl;
      
      dTmpErosionTot += ((m_pVCoastPolygon[n]->dGetCliffCollapseErosionFine() + m_pVCoastPolygon[n]->dGetCliffCollapseErosionSand() + m_pVCoastPolygon[n]->dGetCliffCollapseErosionCoarse()) * m_dCellArea);

      dTmpErosionDDPDTot += ((m_pVCoastPolygon[n]->dGetCliffCollapseSandErodedDeanProfile() + m_pVCoastPolygon[n]->dGetCliffCollapseCoarseErodedDeanProfile()) * m_dCellArea);

      dTmpDepositTot += ((m_pVCoastPolygon[n]->dGetCliffCollapseToSuspensionFine() + m_pVCoastPolygon[n]->dGetCliffCollapseUnconsSandDeposition() + m_pVCoastPolygon[n]->dGetCliffCollapseUnconsCoarseDeposition()) * m_dCellArea);

      dTmpErosionFineTot += (m_pVCoastPolygon[n]->dGetCliffCollapseErosionFine() * m_dCellArea);            
      dTmpSuspensionFineTot += (m_pVCoastPolygon[n]->dGetCliffCollapseToSuspensionFine() * m_dCellArea);

      dTmpErosionSandTot += ((m_pVCoastPolygon[n]->dGetCliffCollapseErosionSand()) * m_dCellArea);
      dTmpErosionSandDDPDTot += ((m_pVCoastPolygon[n]->dGetCliffCollapseSandErodedDeanProfile()) * m_dCellArea);
      dTmpDepositSandTot += (m_pVCoastPolygon[n]->dGetCliffCollapseUnconsSandDeposition() * m_dCellArea);

      dTmpErosionCoarseTot += ((m_pVCoastPolygon[n]->dGetCliffCollapseErosionCoarse()) * m_dCellArea);
      dTmpErosionCoarseDDPDTot += ((m_pVCoastPolygon[n]->dGetCliffCollapseCoarseErodedDeanProfile()) * m_dCellArea);
      dTmpDepositCoarseTot += m_pVCoastPolygon[n]->dGetCliffCollapseUnconsCoarseDeposition() * m_dCellArea;
   }
   
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << "TOTAL from cliff collapse          |" << strDblRight(dTmpErosionTot, 0, 14) << " " << strDblRight(dTmpErosionDDPDTot, 0, 14) << "|" << strDblRight(dTmpDepositTot, 0, 14) << "|"
   << strDblRight(dTmpErosionFineTot, 0, 14) << "|" << strDblRight(dTmpSuspensionFineTot, 0, 14) << "|"
   << strDblRight(dTmpErosionSandTot, 0, 14) << " " << strDblRight(dTmpErosionSandDDPDTot, 0, 14) << "|" << strDblRight(dTmpDepositSandTot, 0, 14) << "|"
   << strDblRight(dTmpErosionCoarseTot, 0, 14) << " " << strDblRight(dTmpErosionCoarseDDPDTot, 0, 14) << "|" << strDblRight(dTmpDepositCoarseTot, 0, 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|" << endl << endl;
}

//===============================================================================================================================
//! Writes to the log file a table showing per-polygon totals of stored unconsolidated beach sediment prior to polygon-to-polygon movement
//===============================================================================================================================
void CSimulation::WritePolygonSedimentBeforeMovement(int const nCoast)
{
   LogStream << endl << m_ulIter << ": Per-polygon totals of stored unconsolidated beach sediment prior to polygon-to-polygon movement (all m^3). Note that this does not include unconsolidated sediment stored outside the polygons." << endl;
   
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << strCentre("Polygon", 11) << "|" << strCentre("Coast", 11) << "|" << strCentre("Polygon", 11) << "|" << strCentre("All", 14) << "|" << strCentre("Fine", 14) << "|" << strCentre("Sand", 14) << "|" << strCentre("Coarse", 14) << "|" << endl;
   LogStream << strCentre("Global ID", 11) << "|" << strCentre("", 11) << "|" << strCentre("Coast ID", 11) << "|" << strCentre("Sediment", 14) << "|" << strCentre("Sediment", 14) << "|" << strCentre("Sediment", 14) << "|" << strCentre("Sediment", 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
   
   double dTmpTot = 0;
   double dTmpFineTot = 0;
   double dTmpSandTot = 0;
   double dTmpCoarseTot = 0;
   
   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      double dFine = m_pVCoastPolygon[n]->dGetPreExistingUnconsFine();
      double dSand = m_pVCoastPolygon[n]->dGetPreExistingUnconsSand();
      double dCoarse = m_pVCoastPolygon[n]->dGetPreExistingUnconsCoarse();
      
      LogStream << strIntRight(m_pVCoastPolygon[n]->nGetGlobalID(), 11) << "|" << strIntRight(nCoast, 11) << "|" << strIntRight(m_pVCoastPolygon[n]->nGetCoastID(), 11) << "|" << strDblRight((dFine + dSand + dCoarse) * m_dCellArea, 0, 14) << "|" << strDblRight(dFine * m_dCellArea, 0, 14) << "|" <<  strDblRight(dSand * m_dCellArea, 0, 14) << "|" << strDblRight(dCoarse * m_dCellArea, 0, 14) << "|" << endl;
      
      dTmpTot += (dFine + dSand + dCoarse) * m_dCellArea;
      dTmpFineTot += (dFine * m_dCellArea);
      dTmpSandTot += (dSand * m_dCellArea);
      dTmpCoarseTot += (dCoarse * m_dCellArea);         
   }
   
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << "TOTAL                              |" << strDblRight(dTmpTot, 0, 14) << "|" << strDblRight(dTmpFineTot, 0, 14) << "|" << strDblRight(dTmpSandTot, 0, 14) << "|" << strDblRight(dTmpCoarseTot, 0, 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl << endl;
}

//===============================================================================================================================
//! Writes to the log file a table showing per-polygon potential erosion of all size classes of unconsolidated beach sediment
//===============================================================================================================================
void CSimulation::WritePolygonPotentialErosion(int const nCoast)
{
   LogStream << m_ulIter << ": Per-polygon potential (i.e. not considering sediment availability) erosion of all size classes of unconsolidated beach sediment (-ve, all m^3), calculated with the ";
   if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_CERC)
      LogStream << "CERC";
   else if (m_nBeachErosionDepositionEquation == UNCONS_SEDIMENT_EQUATION_KAMPHUIS)
      LogStream << "Kamphuis";
   LogStream << " equation." << endl;
   
   LogStream << "-----------|-----------|-----------|--------------|" << endl;
   LogStream << strCentre("Polygon", 11) << "|" << strCentre("Coast", 11) << "|" << strCentre("Polygon", 11) << "|" << strCentre("Potential", 14) << "|" << endl;
   LogStream << strCentre("Global ID", 11) << "|" << strCentre("", 11) << "|" << strCentre("Coast ID", 11) << "|" << strCentre("Erosion", 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|" << endl;
   
   double dTmpTot = 0;
   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      LogStream << strIntRight(m_pVCoastPolygon[n]->nGetGlobalID(), 11) << "|" << strIntRight(nCoast, 11) << "|" << strIntRight(m_pVCoastPolygon[n]->nGetCoastID(), 11) << "|" << strDblRight(m_pVCoastPolygon[n]->dGetPotentialErosion() * m_dCellArea, 0, 14) << "|" << endl;
      
      dTmpTot += (m_pVCoastPolygon[n]->dGetPotentialErosion() * m_dCellArea);
   }      
   LogStream << "-----------|-----------|-----------|--------------|" << endl;
   LogStream << "TOTAL                              |" << strDblRight(dTmpTot, 0, 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|" << endl << endl;
}

// //===============================================================================================================================
// //! Writes to the log file a table showing per-polygon supply-limited erosion of unconsolidated beach sediment
// //===============================================================================================================================
// void CSimulation::WritePolygonUnconsErosion(int const nCoast)
// {
//    LogStream << m_ulIter << ": per-polygon supply-limited erosion of unconsolidated beach sediment (-ve, all m^3). All fine sediment eroded goes to suspension." << endl;
//
//    LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
//    LogStream << strCentre("Polygon", 11) << "|" << strCentre("Coast", 11) << "|" << strCentre("Polygon", 11) << "|" << strCentre("All", 14) << "|" << strCentre("Fine", 14) <<"|" << strCentre("Sand", 14) << "|" << strCentre("Coarse", 14) << "|" << endl;
//    LogStream << strCentre("Global ID", 11) << "|" << strCentre("", 11) << "|" << strCentre("Coast ID", 11) << "|" << strCentre("Estimated", 14) << "|" << strCentre("Estimated", 14) <<"|" << strCentre("Estimated", 14) << "|" << strCentre("Estimated", 14) << "|" << endl;
//    LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
//
//    double
//       dTmpTot = 0,
//       dTmpFineTot = 0,
//       dTmpSandTot = 0,
//       dTmpCoarseTot = 0;
//
//    for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
//    {
//       CGeomCoastPolygon const* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
//
//       LogStream << strIntRight(m_pVCoastPolygon[nPoly]->nGetGlobalID(), 11) << "|" << strIntRight(nCoast, 11) << "|" << strIntRight(m_pVCoastPolygon[nPoly]->nGetCoastID(), 11) << "|" << strDblRight((pPolygon->dGetBeachErosionUnconsFine() + pPolygon->dGetBeachErosionUnconsSand() + pPolygon->dGetBeachErosionUnconsCoarse()) * m_dCellArea, 0, 14) << "|" << strDblRight(pPolygon->dGetBeachErosionUnconsFine() * m_dCellArea, 0, 14) << "|" << strDblRight(pPolygon->dGetBeachErosionUnconsSand() * m_dCellArea, 0, 14) << "|" << strDblRight(pPolygon->dGetBeachErosionUnconsCoarse() * m_dCellArea, 0, 14) << "|" << endl;
//
//       dTmpTot += (pPolygon->dGetBeachErosionUnconsFine() + pPolygon->dGetBeachErosionUnconsSand() + pPolygon->dGetBeachErosionUnconsCoarse()) * m_dCellArea;
//       dTmpFineTot += (pPolygon->dGetBeachErosionUnconsFine() * m_dCellArea);
//       dTmpSandTot += (pPolygon->dGetBeachErosionUnconsSand() * m_dCellArea);
//       dTmpCoarseTot += (pPolygon->dGetBeachErosionUnconsCoarse() * m_dCellArea);
//    }
//
//    LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl;
//    LogStream << "TOTAL estimated erosion            |" << strDblRight(dTmpTot, 0, 14) << "|" << strDblRight(dTmpFineTot, 0, 14) << "|" << strDblRight(dTmpSandTot, 0, 14) << "|" << strDblRight(dTmpCoarseTot, 0, 14) << "|" << endl;
//    LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|" << endl << endl;
// }

// //===============================================================================================================================
// //! Writes to the log file a table showing the unsorted sequence of polygon processing
// //===============================================================================================================================
// void CSimulation::WritePolygonUnsortedSequence(int const nCoast, vector<vector<int> >& pnVVPolyAndAdjacent)
// {
//    LogStream << m_ulIter << ": unsorted sequence of polygon processing" << endl;
//    LogStream << strCentre("Coast", 14) << "|" << strCentre("From Poly", 14) << "|" << strCentre("Direction", 14) << "|" << strCentre("To Poly", 14) << "|" << endl;
//    LogStream << "--------------|--------------|--------------|--------------|" << endl;
//
//    for (int n = 0; n < static_cast<int>(pnVVPolyAndAdjacent.size()); n++)
//    {
//       LogStream << strIntRight(nCoast, 14) << "|";
//
//       for (int m = 0; m < static_cast<int>(pnVVPolyAndAdjacent[n].size()); m++)
//       {
//
//          if (m == 1)
//          {
//             if (pnVVPolyAndAdjacent[n][m] == true)
//                LogStream << strCentre("DOWN ", 14) << "|";
//             else
//                LogStream << strCentre("UP   ", 14) << "|";
//          }
//          else
//             LogStream << strIntRight(pnVVPolyAndAdjacent[n][m], 14) << "|";
//       }
//       LogStream << endl;
//    }
//    LogStream << "--------------|--------------|--------------|--------------|" << endl << endl;
// }

//===============================================================================================================================
//! Writes to the log file a table showing the sorted sequence of polygon processing, and any circularities
//===============================================================================================================================
void CSimulation::WritePolygonSortedSequence(int const nCoast, vector<vector<int> >& pnVVPolyAndAdjacent)
{
   // Show sorted order of polygon processing, and any circularities
   LogStream << m_ulIter << ": Sorted sequence of polygon processing (" << INT_NODATA << " = leaves grid), and any X -> Y -> X circularities" << endl;

   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|" << endl;
   LogStream << strCentre("From", 11) << "|" << strCentre("Coast", 11) << "|" << strCentre("From", 11) << "|" << strCentre("Direction", 14) << "|" << strCentre("To", 14) << "|" << strCentre("Circular", 14) << "|" << endl;
   LogStream << strCentre("Polygon", 11) << "|" << strCentre("", 11) << "|" << strCentre("Polygon", 11) << "|" << strCentre("", 14) << "|" << strCentre("Polygon", 14) << "|" << strCentre("-ity?", 14) << "|" << endl;
   LogStream << strCentre("Global ID", 11) << "|" << strCentre("", 11) << "|" << strCentre("Coast ID", 11) << "|" << strCentre("", 14) << "|" << strCentre("Coast ID", 14) << "|" << strCentre("", 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|" << endl;
   
   for (int nPoly = 0; nPoly < static_cast<int>(pnVVPolyAndAdjacent.size()); nPoly++)
   {
      const CGeomCoastPolygon* pPoly = m_VCoast[nCoast].pGetPolygon(pnVVPolyAndAdjacent[nPoly][0]);
      vector<int> VCirc = pPoly->VnGetCircularities();

      LogStream << strIntRight(m_pVCoastPolygon[pnVVPolyAndAdjacent[nPoly][0]]->nGetGlobalID(), 11) << "|" << strIntRight(nCoast, 11) << "|" << strIntRight(pnVVPolyAndAdjacent[nPoly][0], 11) << "|";
      
      string strTmp = "";
      for (int m = 0; m < static_cast<int>(pnVVPolyAndAdjacent[nPoly].size()); m++)
      {
         if (m == 1)
         {
            if (pnVVPolyAndAdjacent[nPoly][m] == true)
               LogStream << strCentre("DOWN ", 14) << "|";
            else
               LogStream << strCentre("UP   ", 14) << "|";
         }
         else if (m > 1)
         {
            // These are the "To" polygons
            strTmp += to_string(pnVVPolyAndAdjacent[nPoly][m]);
                  
            if (m < (static_cast<int>(pnVVPolyAndAdjacent[nPoly].size()) - 1))
               strTmp += ", ";
         }
      }
      LogStream << strRight(strTmp, 14) << "|";
      
      strTmp = "";
      
      // Now check for circularities                  
      if (! VCirc.empty())
      {                     
         // There is at least one circularity                     
         for (unsigned int i = 0; i < VCirc.size(); i++)
         {
            strTmp += to_string(VCirc[i]);
            
            if (i < (VCirc.size()-1))
               strTmp += ", ";
         }
      }
      LogStream << strCentre(strTmp, 14) << "|" << endl;
   }
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|" << endl;
}

//===============================================================================================================================
//! Writes to the log file a table showing per-polygon actual movement of unconsolidated beach sediment
//===============================================================================================================================
void CSimulation::WritePolygonActualMovement(int const nCoast, vector<vector<int> > const& pnVVPolyAndAdjacent)
{
   // Show estimated polygon-to-polygon movement
   LogStream << m_ulIter << ": per-polygon erosion and deposition of unconsolidated beach sediment, all m^3. Fine sediment is moved to suspension, not deposited (DDPD = During Dean Profile Deposition)." << endl;
   
   LogStream << "-----------|-----------|-----------|--------------------------------------------|-----------------------------|--------------------------------------------|--------------------------------------------|" << endl;
   LogStream << strCentre("Polygon", 11) << "|" << strCentre("Coast", 11) << "|" << strCentre("Polygon", 11) << "|" << strCentre("All", 44) << "|" << strCentre("Fine", 29) << "|" << strCentre("Sand", 44) << "|" << strCentre("Coarse", 44) << "|" << endl;
   LogStream << strCentre("Global ID", 11) << "|" << strCentre("", 11) << "|" << strCentre("Coast ID", 11) << "|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << strCentre("", 11) << "|" << strCentre("", 11) << "|" << strCentre("", 11) << "|" << strCentre("Erosion", 14) << " " << strCentre("Erosion DPDD", 14) << "|" << strCentre("Dep + Susp", 14) << "|" << strCentre("Erosion", 14) << "|" << strCentre("Suspension", 14) << "|" << strCentre("Erosion", 14) << " " << strCentre("Erosion DPDD", 14) << "|" << strCentre("Deposition", 14) << "|" << strCentre("Erosion", 14) << " " << strCentre("Erosion DPDD", 14) << "|" << strCentre("Deposition", 14) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;

   double dTmpTotAllErosion = 0;
   double dTmpTotAllErosionDDPD = 0;
   double dTmpTotAllDeposition = 0;

   double dTmpFineErosion = 0;
   double dTmpFineDeposition = 0;

   double dTmpSandErosion = 0;
   double dTmpSandErosionDDPD = 0;
   double dTmpSandDeposition = 0;

   double dTmpCoarseErosion = 0;
   double dTmpCoarseErosionDDPD = 0;
   double dTmpCoarseDeposition = 0;
      
   for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
   {
      int nPoly = pnVVPolyAndAdjacent[n][0];

      double dAllErosionNotDDPD = -m_pVCoastPolygon[nPoly]->dGeBeachErosionAllUncons() - (m_pVCoastPolygon[nPoly]->dGetBeachSandErodedDeanProfile() + m_pVCoastPolygon[nPoly]->dGetBeachCoarseErodedDeanProfile());
      double dSandErosionNotDDPD = -m_pVCoastPolygon[nPoly]->dGetBeachErosionUnconsSand() - m_pVCoastPolygon[nPoly]->dGetBeachSandErodedDeanProfile();
      double dCoarseErosionNotDDPD = -m_pVCoastPolygon[nPoly]->dGetBeachErosionUnconsCoarse() - m_pVCoastPolygon[nPoly]->dGetBeachCoarseErodedDeanProfile();

      LogStream << strIntRight(m_pVCoastPolygon[nPoly]->nGetGlobalID(), 11) << "|" << strIntRight(nCoast, 11) << "|" << strIntRight(m_pVCoastPolygon[nPoly]->nGetCoastID(), 11) << "|"
         // All
         << strDblRight(dAllErosionNotDDPD * m_dCellArea, 0, 14) << " " << strDblRight((m_pVCoastPolygon[nPoly]->dGetBeachSandErodedDeanProfile() + m_pVCoastPolygon[nPoly]->dGetBeachCoarseErodedDeanProfile()) * m_dCellArea, 0, 14) << "|" << strDblRight(m_pVCoastPolygon[nPoly]->dGetBeachDepositionAndSuspensionAllUncons() * m_dCellArea, 0, 14) << "|"
         // Fine
         << strDblRight(-m_pVCoastPolygon[nPoly]->dGetBeachErosionUnconsFine() * m_dCellArea, 0, 14) << "|" << strDblRight(-m_pVCoastPolygon[nPoly]->dGetBeachErosionUnconsFine() * m_dCellArea, 0, 14) << "|"
         // Sand
         << strDblRight(dSandErosionNotDDPD * m_dCellArea, 0, 14) << " " << strDblRight(m_pVCoastPolygon[nPoly]->dGetBeachSandErodedDeanProfile() * m_dCellArea, 0, 14) <<  "|" << strDblRight(m_pVCoastPolygon[nPoly]->dGetBeachDepositionUnconsSand() * m_dCellArea, 0, 14) << "|"
         // Coarse
         << strDblRight(dCoarseErosionNotDDPD * m_dCellArea, 0, 14) << " " << strDblRight(m_pVCoastPolygon[nPoly]->dGetBeachCoarseErodedDeanProfile() * m_dCellArea, 0, 14) << "|" << strDblRight(m_pVCoastPolygon[nPoly]->dGetToDoBeachDepositionUnconsCoarse() * m_dCellArea, 0, 14) << "|" << endl;
      
      dTmpTotAllErosion += (dAllErosionNotDDPD * m_dCellArea);
      dTmpTotAllErosionDDPD += ((m_pVCoastPolygon[nPoly]->dGetBeachSandErodedDeanProfile() + m_pVCoastPolygon[nPoly]->dGetBeachCoarseErodedDeanProfile()) * m_dCellArea);
      dTmpTotAllDeposition += (m_pVCoastPolygon[nPoly]->dGetBeachDepositionAndSuspensionAllUncons() * m_dCellArea);
      
      dTmpFineErosion += (m_pVCoastPolygon[nPoly]->dGetBeachErosionUnconsFine() * m_dCellArea);
      // dTmpFineDeposition +=
      
      dTmpSandErosion += (dSandErosionNotDDPD * m_dCellArea);
      dTmpSandErosionDDPD += (m_pVCoastPolygon[nPoly]->dGetBeachSandErodedDeanProfile() * m_dCellArea);
      dTmpSandDeposition += (m_pVCoastPolygon[nPoly]->dGetBeachDepositionUnconsSand() * m_dCellArea);

      dTmpCoarseErosion += (dCoarseErosionNotDDPD * m_dCellArea);
      dTmpCoarseErosionDDPD += (m_pVCoastPolygon[nPoly]->dGetBeachCoarseErodedDeanProfile() * m_dCellArea);
      dTmpCoarseDeposition += (m_pVCoastPolygon[nPoly]->dGetToDoBeachDepositionUnconsCoarse() * m_dCellArea);
   }         
   
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << "Lost from grid                     |" << strLeft("", 14) << "|" << strLeft("", 14) << "|" << strDblRight((m_dThisIterLeftGridUnconsSand + m_dThisIterLeftGridUnconsCoarse) * m_dCellArea, 0, 14) << "|" << strLeft("", 14) << "|" << strLeft("", 14) << "|" << strLeft("", 14) << "|" << strLeft("", 14) << "|" << strDblRight(m_dThisIterLeftGridUnconsSand * m_dCellArea, 0, 14) << "|" << strLeft("", 14) <<  "|" << strLeft("", 14) <<  "|" << strDblRight(m_dThisIterLeftGridUnconsCoarse * m_dCellArea, 0, 14) << "|" << endl;

   dTmpTotAllDeposition += ((m_dThisIterLeftGridUnconsSand + m_dThisIterLeftGridUnconsCoarse) * m_dCellArea);
   dTmpSandDeposition += (m_dThisIterLeftGridUnconsSand * m_dCellArea);
   dTmpCoarseDeposition += (m_dThisIterLeftGridUnconsCoarse * m_dCellArea);
   
   bool bShowZeroFine = false;
   bool bShowZeroSand = false;
   bool bShowZeroCoarse = false;
   
   if (! bFPIsEqual(dTmpFineErosion, 0.0, MASS_BALANCE_TOLERANCE))
      bShowZeroFine = true;

   if (! bFPIsEqual(dTmpSandErosion, 0.0, MASS_BALANCE_TOLERANCE))
      bShowZeroSand = true;

   if (! bFPIsEqual(dTmpCoarseErosion, 0.0, MASS_BALANCE_TOLERANCE))
      bShowZeroCoarse = true;

   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
   LogStream << "TOTAL                              |"
      // All
      << strDblRight(-dTmpTotAllErosion, 0, 14) << "|" << strDblRight(dTmpTotAllErosionDDPD, 0, 14) << "|" << strDblRight(dTmpTotAllDeposition, 0, 14) << "|"
      // Fine
      << strDblRight(-dTmpFineErosion, 0, 14, bShowZeroFine) << "|" << strDblRight(dTmpFineErosion, 0, 14, bShowZeroFine) << "|"
      // Sand
      << strDblRight(-dTmpSandErosion, 0, 14, bShowZeroSand) << "|" << strDblRight(dTmpSandErosionDDPD, 0, 14, bShowZeroSand) << "|" << strDblRight(dTmpSandDeposition, 0, 14, bShowZeroSand) << "|"
      // Coarse
      << strDblRight(-dTmpCoarseErosion, 0, 14, bShowZeroCoarse) << "|" << strDblRight(dTmpCoarseErosionDDPD, 0, 14, bShowZeroCoarse) << "|" << strDblRight(dTmpCoarseDeposition, 0, 14, bShowZeroCoarse) << "|" << endl;
   LogStream << "-----------|-----------|-----------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|" << endl;
}

//===============================================================================================================================
//! Update and print totals at the end of each timestep
//===============================================================================================================================
void CSimulation::DoTimestepTotals(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
   {
      // TODO 062 show these results somewhere
      // At end of timestep, get the numbers of cells with consolidated and unconsolidated sediment stored or in suspension, for all cells (both inside and outside polygons)
      int nSuspFineCellsAllCells = 0;
      int nUnconsFineCellsAllCells = 0;
      int nUnconsSandCellsAllCells = 0;
      int nUnconsCoarseCellsAllCells = 0;

      // TODO 062 show these results somewhere
      // Also get the numbers of cells with consolidated and unconsolidated sediment stored or in suspension, only for cells inside polygons
      int nSuspFineCellsInPolygons = 0;
      int nUnconsFineCellsInPolygons = 0;
      int nUnconsSandCellsInPolygons = 0;
      int nUnconsCoarseCellsInPolygons = 0;

      // Get the depths of consolidated and unconsolidated sediment stored or in suspension, for all cells (both inside and outside polygons)
      double dEndIterSuspFineAllCells = 0;
      double dEndIterUnconsFineAllCells = 0;
      double dEndIterUnconsSandAllCells = 0;
      double dEndIterUnconsCoarseAllCells = 0;
      double dEndIterConsFineAllCells = 0;
      double dEndIterConsSandAllCells = 0;
      double dEndIterConsCoarseAllCells = 0;

      // Also get the depths of consolidated and unconsolidated sediment stored or in suspension, only for cells inside polygons
      double dEndIterSuspFineInPolygons = 0;
      double dEndIterUnconsFineInPolygons = 0;
      double dEndIterUnconsSandInPolygons = 0;
      double dEndIterUnconsCoarseInPolygons = 0;
      double dEndIterConsFineInPolygons = 0;
      double dEndIterConsSandInPolygons = 0;
      double dEndIterConsCoarseInPolygons = 0;

      // Get depth of consolidated and unconsolidated (and suspended) sediment from each cell (both within and outside polygons)
      for (int nX = 0; nX < m_nXGridMax; nX++)
      {
         for (int nY = 0; nY < m_nYGridMax; nY++)
         {
            dEndIterConsFineAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetTotConsFineThickConsiderNotch();
            dEndIterConsSandAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetTotConsSandThickConsiderNotch();
            dEndIterConsCoarseAllCells += m_pRasterGrid->m_Cell[nX][nY].dGetTotConsCoarseThickConsiderNotch();

            double dSuspFine = m_pRasterGrid->m_Cell[nX][nY].dGetSuspendedSediment();
            if (dSuspFine > 0)
            {
               dEndIterSuspFineAllCells += dSuspFine;
               nSuspFineCellsAllCells++;
            }

            double dUnconsFine = m_pRasterGrid->m_Cell[nX][nY].dGetTotUnconsFine();
            if (dUnconsFine > 0)
            {
               dEndIterUnconsFineAllCells += dUnconsFine;
               nUnconsFineCellsAllCells++;
            }

            double dUnconsSand = m_pRasterGrid->m_Cell[nX][nY].dGetTotUnconsSand();
            if (dUnconsSand > 0)
            {
               dEndIterUnconsSandAllCells += dUnconsSand;
               nUnconsSandCellsAllCells++;
            }

            double dUnconsCoarse = m_pRasterGrid->m_Cell[nX][nY].dGetTotUnconsCoarse();
            if (dUnconsCoarse > 0)
            {
               dEndIterUnconsCoarseAllCells += dUnconsCoarse;
               nUnconsCoarseCellsAllCells++;
            }

            // Is this cell within a polygon?
            if (m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID() != INT_NODATA)
            {
               // It is within a polygon
               dEndIterConsFineInPolygons += m_pRasterGrid->m_Cell[nX][nY].dGetTotConsFineThickConsiderNotch();
               dEndIterConsSandInPolygons += m_pRasterGrid->m_Cell[nX][nY].dGetTotConsSandThickConsiderNotch();
               dEndIterConsCoarseInPolygons += m_pRasterGrid->m_Cell[nX][nY].dGetTotConsCoarseThickConsiderNotch();

               dSuspFine = m_pRasterGrid->m_Cell[nX][nY].dGetSuspendedSediment();
               if (dSuspFine > 0)
               {
                  dEndIterSuspFineInPolygons += dSuspFine;
                  nSuspFineCellsInPolygons++;
               }

               dUnconsFine = m_pRasterGrid->m_Cell[nX][nY].dGetTotUnconsFine();
               if (dUnconsFine > 0)
               {
                  dEndIterUnconsFineInPolygons += dUnconsFine;
                  nUnconsFineCellsInPolygons++;
               }

               dUnconsSand = m_pRasterGrid->m_Cell[nX][nY].dGetTotUnconsSand();
               if (dUnconsSand > 0)
               {
                  dEndIterUnconsSandInPolygons += dUnconsSand;
                  nUnconsSandCellsInPolygons++;
               }

               dUnconsCoarse = m_pRasterGrid->m_Cell[nX][nY].dGetTotUnconsCoarse();
               if (dUnconsCoarse > 0)
               {
                  dEndIterUnconsCoarseInPolygons += dUnconsCoarse;
                  nUnconsCoarseCellsInPolygons++;
               }
            }
         }
      }

      double dFineTmp = 0;
      double dSandTmp = 0;
      double dCoarseTmp = 0;

      LogStream << endl << m_ulIter << ": Consolidated sediment budget, all m^3. This includes sediment within and outside the polygons." << endl;

      // Stored consolidated sediment at start of timestep: sediment within and outside the polygons
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("", 8) << "|" << strCentre("In Polygons", 14) << "|" << strCentre("Outside", 14) << "|" << strCentre("ALL", 14) << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion (since polygons are not constructed if no beach erosion)
      double dTmp1 = (m_bDoBeachSedimentTransport ? m_dTotalFineConsInPolygons : 0);
      double dTmp2 = (m_bDoBeachSedimentTransport ? (m_dStartIterConsFineAllCells - m_dTotalFineConsInPolygons) : 0);

      LogStream << strLeft("At start of timestep " + to_string(m_ulIter) + ", stored consolidated sediment", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(dTmp1 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp2 * m_dCellArea, 0, 14) << "|" << strDblRight(m_dStartIterConsFineAllCells * m_dCellArea, 0, 14) << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion (since polygons are not constructed if no beach erosion)
      double dTmp3 = (m_bDoBeachSedimentTransport ? m_dTotalSandConsInPolygons : 0);
      double dTmp4 = (m_bDoBeachSedimentTransport ? (m_dStartIterConsSandAllCells - m_dTotalSandConsInPolygons) : 0);

      LogStream << strLeft("", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(dTmp3 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp4 * m_dCellArea, 0, 14) << "|" << strDblRight(m_dStartIterConsSandAllCells * m_dCellArea, 0, 14) << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion (since polygons are not constructed if no beach erosion)
      double dTmp5 = (m_bDoBeachSedimentTransport ? m_dTotalCoarseConsInPolygons : 0);
      double dTmp6 = (m_bDoBeachSedimentTransport ? (m_dStartIterConsCoarseAllCells - m_dTotalCoarseConsInPolygons) : 0);

      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(dTmp5 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp6 * m_dCellArea, 0, 14) << "|" << strDblRight( m_dStartIterConsCoarseAllCells * m_dCellArea, 0, 14) << "|" << endl;

      LogStream << strLeft("", 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight((dTmp1 + dTmp3 + dTmp5) * m_dCellArea, 0, 14) << "|" << strDblRight((dTmp2 + dTmp4 + dTmp6) * m_dCellArea, 0, 14) << "|" << strDblRight((m_dStartIterConsFineAllCells + m_dStartIterConsSandAllCells + m_dStartIterConsCoarseAllCells) * m_dCellArea, 0, 14) << "|" << endl;

      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;

      dFineTmp = m_dStartIterConsFineAllCells * m_dCellArea;
      dSandTmp = m_dStartIterConsSandAllCells * m_dCellArea;
      dCoarseTmp = m_dStartIterConsCoarseAllCells * m_dCellArea;

      // Shore platform erosion, consolidated sediment lost (becomes unconsolidated esediment and suspended sediment)
      LogStream << strLeft("Consolidated sediment lost via platform erosion (becomes suspended sediment)", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(-m_dThisIterActualPlatformErosionFineCons * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("Consolidated sediment lost via platform erosion (becomes unconsolidated sediment)", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(-m_dThisIterActualPlatformErosionSandCons * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(-m_dThisIterActualPlatformErosionCoarseCons * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight(-(m_dThisIterActualPlatformErosionFineCons + m_dThisIterActualPlatformErosionSandCons + m_dThisIterActualPlatformErosionCoarseCons) * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(44, '-') << "|" << endl;

      dFineTmp -= (m_dThisIterActualPlatformErosionFineCons * m_dCellArea);
      dSandTmp -= (m_dThisIterActualPlatformErosionSandCons * m_dCellArea);
      dCoarseTmp -= (m_dThisIterActualPlatformErosionCoarseCons * m_dCellArea);

      // Cliff collapse, consolidated sediment lost (becomes unconsolidated esediment and suspended sediment)
      LogStream << strLeft("Consolidated sediment lost via cliff collapse (becomes suspended sediment)", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(-m_dThisIterCliffCollapseErosionFineCons * m_dCellArea, 0, 44) << "|" << endl;
      LogStream << strLeft("Consolidated sediment lost via cliff collapse (becomes unconsolidated sediment)", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(-m_dThisIterCliffCollapseErosionSandCons * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(-m_dThisIterCliffCollapseErosionCoarseCons * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight(-(m_dThisIterCliffCollapseErosionFineCons + m_dThisIterCliffCollapseErosionSandCons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea, 0, 44)  << "|" << endl;

      dFineTmp -= (m_dThisIterCliffCollapseErosionFineCons * m_dCellArea);
      dSandTmp -= (m_dThisIterCliffCollapseErosionSandCons * m_dCellArea);
      dCoarseTmp -= (m_dThisIterCliffCollapseErosionCoarseCons * m_dCellArea);

      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("", 8) << "|" << strCentre("In Polygons", 14) << "|" << strCentre("Outside", 14) << "|" << strCentre("ALL", 14) << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion
      dTmp1 = (m_bDoBeachSedimentTransport ? dEndIterConsFineInPolygons : 0);
      dTmp2 = (m_bDoBeachSedimentTransport ? (dEndIterConsFineAllCells - dEndIterConsFineInPolygons) : 0);

      LogStream << strLeft("At end of timestep " + to_string(m_ulIter) + ", stored consolidated sediment", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(dTmp1 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp2 * m_dCellArea, 0, 14) << "|" << strDblRight(dEndIterConsFineAllCells * m_dCellArea, 0, 14)  << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion
      dTmp3 = (m_bDoBeachSedimentTransport ? dEndIterConsSandInPolygons : 0);
      dTmp4 = (m_bDoBeachSedimentTransport ? (dEndIterConsSandAllCells - dEndIterConsSandInPolygons) : 0);
      LogStream << strLeft("", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(dTmp3 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp4 * m_dCellArea, 0, 14) << "|" << strDblRight(dEndIterConsSandAllCells * m_dCellArea, 0, 14) << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion
      dTmp5 = (m_bDoBeachSedimentTransport ? dEndIterConsCoarseInPolygons : 0);
      dTmp6 = (m_bDoBeachSedimentTransport ? (dEndIterConsCoarseAllCells - dEndIterConsCoarseInPolygons) : 0);

      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(dTmp5 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp6 * m_dCellArea, 0, 14) << "|" << strDblRight(dEndIterConsCoarseAllCells * m_dCellArea, 0, 14)  << "|" << endl;

      LogStream << strLeft("", 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight((dTmp1 + dTmp3 + dTmp5) * m_dCellArea, 0, 14) << "|" << strDblRight((dTmp2 + dTmp4 + dTmp6) * m_dCellArea, 0, 14) << "|" << strDblRight((dEndIterConsFineAllCells + dEndIterConsSandAllCells + dEndIterConsCoarseAllCells) * m_dCellArea, 0, 14)  << "|" << endl;

      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;

      double dFineError = (dEndIterConsFineAllCells * m_dCellArea) - dFineTmp;
      double dSandError = (dEndIterConsSandAllCells * m_dCellArea) - dSandTmp;
      double dCoarseError = (dEndIterConsCoarseAllCells * m_dCellArea) - dCoarseTmp;

      // Mass balance check
      bool bError = false;
      string strFineErrMsg = "";
      string strSandErrMsg = "";
      string strCoarseErrMsg = "";
      string strAllErrMsg = "";

      if (! bFPIsEqual(dFineError, 0.0, MASS_BALANCE_TOLERANCE))
      {
         strFineErrMsg = MASS_BALANCE_ERROR;
         bError = true;
      }

      if (! bFPIsEqual(dSandError, 0.0, MASS_BALANCE_TOLERANCE))
      {
         strSandErrMsg = MASS_BALANCE_ERROR;
         bError = true;
      }

      if (! bFPIsEqual(dCoarseError, 0.0, MASS_BALANCE_TOLERANCE))
      {
         strCoarseErrMsg = MASS_BALANCE_ERROR;
         bError = true;
      }

      if (bError)
         strAllErrMsg = MASS_BALANCE_ERROR;

      LogStream << strLeft("Consolidated sediment mass balance check (+ve means end total > start total)", 101) << "|" << strLeft("", 8) << "|" << strLeft("", 44) << "|" << endl;
      LogStream << strLeft(strFineErrMsg, 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(dFineError, 3, 30, false) << strRightPerCent(dFineError, dEndIterConsFineAllCells, 14, 2) << "|" << endl;
      LogStream << strLeft(strSandErrMsg, 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(dSandError, 3, 30, false) << strRightPerCent(dSandError, dEndIterConsSandAllCells, 14, 2) << "|" << endl;
      LogStream << strLeft(strCoarseErrMsg, 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(dCoarseError, 3, 30, false) << strRightPerCent(dCoarseError, dEndIterConsCoarseAllCells, 14, 2)  << "|" << endl;
      LogStream << strLeft(strAllErrMsg, 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight(dFineError + dSandError + dCoarseError, 3, 30, false) << strRightPerCent(dFineError + dSandError + dCoarseError, dEndIterConsFineAllCells + dEndIterConsSandAllCells + dEndIterConsCoarseAllCells, 14, 2)  << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(44, '-') << "|" << endl;
      LogStream << endl;

      // Now look at unconsolidated sediment
      dFineTmp =
      dSandTmp =
      dCoarseTmp = 0;

      LogStream << m_ulIter << ": Unconsolidated sediment budget, all m^3." << endl;

      // Stored unconsolidated sediment, and in suspension, at start of timestep: sediment within and outside the polygons
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("", 8) << "|" << strCentre("In Polygons", 14) << "|" << strCentre("Outside", 14) << "|" << strCentre("ALL", 14) << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion (since polygons are not constructed if no beach erosion)
      dTmp1 = (m_bDoBeachSedimentTransport ? m_dStartIterSuspFineInPolygons : 0);
      dTmp2 = (m_bDoBeachSedimentTransport ? (m_dStartIterSuspFineAllCells - m_dStartIterSuspFineInPolygons) : 0);

      LogStream << strLeft("At start of timestep " + to_string(m_ulIter) + ", sediment in suspension", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(dTmp1 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp2 * m_dCellArea, 0, 14) << "|" << strDblRight(m_dStartIterSuspFineAllCells * m_dCellArea, 0, 14) << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion (since polygons are not constructed if no beach erosion)
      dTmp3 = (m_bDoBeachSedimentTransport ? m_dTotalFineUnconsInPolygons : 0);
      dTmp4 = (m_bDoBeachSedimentTransport ? (m_dStartIterUnconsFineAllCells - m_dTotalFineUnconsInPolygons) : 0);

      LogStream << strLeft("At start of timestep " + to_string(m_ulIter) + ", stored unconsolidated sediment", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(dTmp3 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp4 * m_dCellArea, 0, 14) << "|" << strDblRight(m_dStartIterUnconsFineAllCells * m_dCellArea, 0, 14) << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion (since polygons are not constructed if no beach erosion)
      dTmp5 = (m_bDoBeachSedimentTransport ? m_dTotalSandUnconsInPolygons : 0);
      dTmp6 = (m_bDoBeachSedimentTransport ? (m_dStartIterUnconsSandAllCells - m_dTotalSandUnconsInPolygons) : 0);

      LogStream << strLeft("", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(dTmp5 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp6 * m_dCellArea, 0, 14) << "|" << strDblRight(m_dStartIterUnconsSandAllCells * m_dCellArea, 0, 14) << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion (since polygons are not constructed if no beach erosion)
      double dTmp7 = (m_bDoBeachSedimentTransport ? m_dTotalCoarseUnconsInPolygons : 0);
      double dTmp8 = (m_bDoBeachSedimentTransport ? (m_dStartIterUnconsCoarseAllCells - m_dTotalCoarseUnconsInPolygons) : 0);

      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(dTmp7 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp8 * m_dCellArea, 0, 14) << "|" << strDblRight( m_dStartIterUnconsCoarseAllCells * m_dCellArea, 0, 14) << "|" << endl;

      LogStream << strLeft("", 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight((dTmp1 + dTmp3 + dTmp5 + dTmp7) * m_dCellArea, 0, 14) << "|" << strDblRight((dTmp2 + dTmp4 + dTmp6 + dTmp8) * m_dCellArea, 0, 14) << "|" << strDblRight((m_dStartIterSuspFineAllCells + m_dStartIterUnconsFineAllCells + m_dStartIterUnconsSandAllCells + m_dStartIterUnconsCoarseAllCells) * m_dCellArea, 0, 14) << "|" << endl;

      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;

      dFineTmp += ((m_dStartIterSuspFineAllCells + m_dStartIterUnconsFineAllCells) * m_dCellArea);
      dSandTmp += (m_dStartIterUnconsSandAllCells * m_dCellArea);
      dCoarseTmp += (m_dStartIterUnconsCoarseAllCells * m_dCellArea);

      // Shore platform erosion, consolidated sediment becomes unconsolidated sediment and suspended sediment
      LogStream << strLeft("Suspended sediment derived from platform erosion of consolidated sediment", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(m_dThisIterActualPlatformErosionFineCons * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("Unconsolidated sediment derived from platform erosion of consolidated sediment", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(m_dThisIterActualPlatformErosionSandCons * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(m_dThisIterActualPlatformErosionCoarseCons * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight((m_dThisIterActualPlatformErosionFineCons + m_dThisIterActualPlatformErosionSandCons + m_dThisIterActualPlatformErosionCoarseCons) * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(44, '-') << "|" << endl;

      dFineTmp += (m_dThisIterActualPlatformErosionFineCons * m_dCellArea);
      dSandTmp += (m_dThisIterActualPlatformErosionSandCons * m_dCellArea);
      dCoarseTmp += (m_dThisIterActualPlatformErosionCoarseCons * m_dCellArea);

      // Cliff collapse, consolidated sediment becomes unconsolidated sediment and suspended sediment
      LogStream << strLeft("Suspended sediment derived from cliff collapse erosion of consolidated sediment", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(m_dThisIterCliffCollapseErosionFineCons * m_dCellArea, 0, 44) << "|" << endl;
      LogStream << strLeft("Unconsolidated sediment derived from cliff collapse erosion of consolidated sediment only", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(m_dThisIterCliffCollapseErosionSandCons  * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(m_dThisIterCliffCollapseErosionCoarseCons * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight((m_dThisIterCliffCollapseErosionFineCons + m_dThisIterCliffCollapseErosionSandCons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(44, '-') << "|" << endl;

      dFineTmp += (m_dThisIterCliffCollapseErosionFineCons * m_dCellArea);
      dSandTmp += (m_dThisIterCliffCollapseErosionSandCons  * m_dCellArea);
      dCoarseTmp += (m_dThisIterCliffCollapseErosionCoarseCons * m_dCellArea);

      // Beach (unconsolidated sediment) lost from grid due to beach erosion and deposition, and to cliff collapse with talus going outside the grid
      LogStream << strLeft("Unconsolidated sediment lost from grid", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(-m_dThisIterLeftGridUnconsFine * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(-m_dThisIterLeftGridUnconsSand * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(-m_dThisIterLeftGridUnconsCoarse * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight(-(m_dThisIterLeftGridUnconsFine + m_dThisIterLeftGridUnconsSand + m_dThisIterLeftGridUnconsCoarse) * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(44, '-') << "|" << endl;

      dFineTmp -= (m_dThisIterLeftGridUnconsFine * m_dCellArea);
      dSandTmp -= (m_dThisIterLeftGridUnconsSand * m_dCellArea);
      dCoarseTmp -= (m_dThisIterLeftGridUnconsCoarse * m_dCellArea);

      // Sediment added via input events
      LogStream << strLeft("Unconsolidated sediment added via input event(s)", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(m_dThisiterUnconsFineInput * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(m_dThisiterUnconsSandInput * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(m_dThisiterUnconsCoarseInput * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(44, '-') << "|" << endl;

      dFineTmp += (m_dThisiterUnconsFineInput * m_dCellArea);
      dSandTmp += (m_dThisiterUnconsSandInput * m_dCellArea);
      dCoarseTmp += (m_dThisiterUnconsCoarseInput * m_dCellArea);

      // Any uncons sediment left over from last iter and deposited this iter?

      LogStream << strLeft("Insufficient unconsolidated sediment deposited last iteration, carried forward to this iteration", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(m_dUnconsSandNotDepositedLastIter * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(m_dUnconsCoarseNotDepositedLastIter * m_dCellArea, 0, 44)  << "|" << endl;

      dSandTmp += (m_dUnconsSandNotDepositedLastIter * m_dCellArea);
      dCoarseTmp += (m_dUnconsCoarseNotDepositedLastIter * m_dCellArea);

      // Any insufficient deposition?
      LogStream << strLeft("Insufficient unconsolidated sediment deposited this iteration", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(m_dDepositionSandDiff * m_dCellArea, 0, 44)  << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(m_dDepositionCoarseDiff * m_dCellArea, 0, 44)  << "|" << endl;

      dSandTmp -= (m_dDepositionSandDiff * m_dCellArea);
      dCoarseTmp -= (m_dDepositionCoarseDiff * m_dCellArea);

      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;
      LogStream << strLeft("", 101) << "|" << strLeft("", 8) << "|" << strCentre("In Polygons", 14) << "|" << strCentre("Outside", 14) << "|" << strCentre("ALL", 14) << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion
      dTmp1 = (m_bDoBeachSedimentTransport ? dEndIterSuspFineInPolygons : 0);
      dTmp2 = (m_bDoBeachSedimentTransport ? (dEndIterSuspFineAllCells - dEndIterSuspFineInPolygons) : 0);

      LogStream << strLeft("At end of timestep " + to_string(m_ulIter) + ", sediment in suspension", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(dTmp1 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp2 * m_dCellArea, 0, 14) << "|" << strDblRight(dEndIterSuspFineAllCells * m_dCellArea, 0, 14)  << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion
      dTmp3 = (m_bDoBeachSedimentTransport ? dEndIterUnconsFineInPolygons : 0);
      dTmp4 = (m_bDoBeachSedimentTransport ? (dEndIterUnconsFineAllCells - dEndIterUnconsFineInPolygons) : 0);

      LogStream << strLeft("At end of timestep " + to_string(m_ulIter) + ", stored unconsolidated sediment", 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(dTmp3 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp4 * m_dCellArea, 0, 14) << "|" << strDblRight(dEndIterUnconsFineAllCells * m_dCellArea, 0, 14)  << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion
      dTmp5 = (m_bDoBeachSedimentTransport ? dEndIterUnconsSandInPolygons : 0);
      dTmp6 = (m_bDoBeachSedimentTransport ? (dEndIterUnconsSandAllCells - dEndIterUnconsSandInPolygons) : 0);

      LogStream << strLeft("", 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(dTmp5 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp6 * m_dCellArea, 0, 14) << "|" << strDblRight(dEndIterUnconsSandAllCells * m_dCellArea, 0, 14) << "|" << endl;

      // Don't show inside/outside polygon values if not simulating beach erosion
      dTmp7 = (m_bDoBeachSedimentTransport ? dEndIterUnconsCoarseInPolygons : 0);
      dTmp8 = (m_bDoBeachSedimentTransport ? (dEndIterUnconsCoarseAllCells - dEndIterUnconsCoarseInPolygons) : 0);

      LogStream << strLeft("", 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(dTmp7 * m_dCellArea, 0, 14) << "|" << strDblRight(dTmp8 * m_dCellArea, 0, 14) << "|" << strDblRight(dEndIterUnconsCoarseAllCells * m_dCellArea, 0, 14)  << "|" << endl;

      LogStream << strLeft("", 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight((dTmp1 + dTmp3 + dTmp5 + dTmp7) * m_dCellArea, 0, 14) << "|" << strDblRight((dTmp2 + dTmp4 + dTmp6 + dTmp8) * m_dCellArea, 0, 14) << "|" << strDblRight((dEndIterSuspFineAllCells + dEndIterUnconsFineAllCells + dEndIterUnconsSandAllCells + dEndIterUnconsCoarseAllCells) * m_dCellArea, 0, 14)  << "|" << endl;

      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << string(14, '-') << "|" << endl;

      dFineError = ((dEndIterSuspFineAllCells + dEndIterUnconsFineAllCells) * m_dCellArea) - dFineTmp,
      dSandError = (dEndIterUnconsSandAllCells * m_dCellArea) - dSandTmp,
      dCoarseError = (dEndIterUnconsCoarseAllCells * m_dCellArea) - dCoarseTmp;

      // Mass balance check
      bError = false;
      strFineErrMsg = "";
      strSandErrMsg = "";
      strCoarseErrMsg = "";
      strAllErrMsg = "";

      if (! bFPIsEqual(dFineError, 0.0, MASS_BALANCE_TOLERANCE))
      {
         strFineErrMsg = MASS_BALANCE_ERROR;
         bError = true;
      }

      if (! bFPIsEqual(dSandError, 0.0, MASS_BALANCE_TOLERANCE))
      {
         strSandErrMsg = MASS_BALANCE_ERROR;
         bError = true;
      }

    if (! bFPIsEqual(dCoarseError, 0.0, MASS_BALANCE_TOLERANCE))
      {
         strCoarseErrMsg = MASS_BALANCE_ERROR;
         bError = true;
      }

      if (bError)
         strAllErrMsg = MASS_BALANCE_ERROR;

      LogStream << strLeft("Unconsolidated sediment mass balance check (+ve means iteration-end total > iteration-start total)", 101) << "|" << strLeft("", 8) << "|" << strLeft("", 44)  << "|" << endl;
      LogStream << strLeft(strFineErrMsg, 101) << "|" << strLeft("Fine", 8) << "|" << strDblRight(dFineError, 0, 30, false) << strRightPerCent(dFineError, (dEndIterSuspFineAllCells + dEndIterUnconsFineAllCells), 14, 2) << "|" << endl;
      LogStream << strLeft(strSandErrMsg, 101) << "|" << strLeft("Sand", 8) << "|" << strDblRight(dSandError, 0, 30, false) << strRightPerCent(dSandError, dEndIterUnconsSandAllCells, 14, 2) << "|" << endl;
      LogStream << strLeft(strCoarseErrMsg, 101) << "|" << strLeft("Coarse", 8) << "|" << strDblRight(dCoarseError, 0, 30, false) << strRightPerCent(dCoarseError, dEndIterUnconsCoarseAllCells, 14, 2) << "|" << endl;
      LogStream << strLeft(strAllErrMsg, 101) << "|" << strLeft("ALL", 8) << "|" << strDblRight(dFineError + dSandError + dCoarseError, 0, 30, false) << strRightPerCent(dFineError + dSandError + dCoarseError, dEndIterSuspFineAllCells + dEndIterUnconsFineAllCells + dEndIterUnconsSandAllCells + dEndIterUnconsCoarseAllCells, 14, 2) << "|" << endl;
      LogStream << string(101, '-') << "|" << string(8, '-') << "|" << string(44, '-') << "|" << endl;
      LogStream << endl;
   }

   // Add to grand totals: first platform erosion
   m_ldGTotPotentialPlatformErosion += m_dThisIterPotentialPlatformErosion;
   //    assert(isfinite(m_dThisIterPotentialPlatformErosion));

   m_ldGTotFineActualPlatformErosion += m_dThisIterActualPlatformErosionFineCons;
   m_ldGTotSandActualPlatformErosion += m_dThisIterActualPlatformErosionSandCons;
   m_ldGTotCoarseActualPlatformErosion += m_dThisIterActualPlatformErosionCoarseCons;

   // Erosion from cliff collapse, both consolidated and unconsolidated
   m_ldGTotCliffCollapseFine += m_dThisIterCliffCollapseErosionFineUncons;
   m_ldGTotCliffCollapseFine += m_dThisIterCliffCollapseErosionFineCons;
   m_ldGTotCliffCollapseSand += m_dThisIterCliffCollapseErosionSandUncons;
   m_ldGTotCliffCollapseSand += m_dThisIterCliffCollapseErosionSandCons;
   m_ldGTotCliffCollapseCoarse += m_dThisIterCliffCollapseErosionCoarseUncons;
   m_ldGTotCliffCollapseCoarse += m_dThisIterCliffCollapseErosionCoarseCons;

   // Deposition (with fine to suspension) of unconsolidated talus from cliff collapse
   m_ldGTotCliffTalusFineToSuspension += m_dThisIterCliffCollapseErosionFineUncons;
   m_ldGTotCliffTalusSandDeposition += m_dThisIterUnconsSandCliffDeposition;
   m_ldGTotCliffTalusCoarseDeposition += m_dThisIterUnconsCoarseCliffDeposition;

   // Erosion of unconsolidated sediment during deposition of unconsolidated cliff collapse talus
   m_ldGTotCliffCollapseFineErodedDuringDeposition += m_dThisIterCliffCollapseFineErodedDuringDeposition;
   m_ldGTotCliffCollapseSandErodedDuringDeposition += m_dThisIterCliffCollapseSandErodedDuringDeposition;
   m_ldGTotCliffCollapseCoarseErodedDuringDeposition += m_dThisIterCliffCollapseCoarseErodedDuringDeposition;

   // Beach erosion of unconsolidated sediment
   m_ldGTotPotentialBeachErosion += m_dThisIterPotentialBeachErosion;

   m_ldGTotActualFineBeachErosion += m_dThisIterBeachErosionFine;
   m_ldGTotActualSandBeachErosion += m_dThisIterBeachErosionSand;
   m_ldGTotActualCoarseBeachErosion += m_dThisIterBeachErosionCoarse;

   // Beach deposition of unconsolidated sediment
   m_ldGTotSandBeachDeposition += m_dThisIterBeachDepositionSand;
   m_ldGTotCoarseBeachDeposition += m_dThisIterBeachDepositionCoarse;

   // Unconsolidated sediment lost due to beach erosion
   m_ldGTotPotentialSedLostBeachErosion += m_dThisIterPotentialSedLostBeachErosion;

   m_ldGTotActualFineLostBeachErosion += m_dThisIterLeftGridUnconsFine;
   m_ldGTotActualSandLostBeachErosion += m_dThisIterLeftGridUnconsSand;
   m_ldGTotActualCoarseLostBeachErosion += m_dThisIterLeftGridUnconsCoarse;

   // Unconsolidated sediment input event(s)
   m_ldGTotFineSedimentInput += m_dThisiterUnconsFineInput;
   m_ldGTotSandSedimentInput += m_dThisiterUnconsSandInput;
   m_ldGTotCoarseSedimentInput += m_dThisiterUnconsCoarseInput;

   // Suspended unconsolidated sediment
   m_ldGTotSuspendedSediment += m_dThisIterFineSedimentToSuspension;

   // Shortfall in unconsolidated sediment deposition
   m_ldGTotSandDepositionDiff += m_dDepositionSandDiff;
   m_ldGTotCoarseDepositionDiff += m_dDepositionCoarseDiff;
}

