/*!
 *
 * \file simulation.cpp
 * \brief The start-of-simulation routine
 * \details TODO 001 A more detailed description of this routine.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2024
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include <assert.h>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::setprecision;

#include <cfloat>

#include <string>
using std::to_string;

#include "cme.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"

//===============================================================================================================================
//! The CSimulation constructor
//===============================================================================================================================
CSimulation::CSimulation (void)
{
   // Initialization
   m_bHaveFineSediment =
   m_bHaveSandSediment =
   m_bHaveCoarseSediment =
   m_bBasementElevSave =
   m_bSedimentTopSurfSave =
   m_bTopSurfSave =
   m_bSliceSave =
   m_bSeaDepthSave =
   m_bAvgSeaDepthSave =
   m_bWaveHeightSave =
   m_bAvgWaveHeightSave =
   m_bWaveAngleSave =
   m_bAvgWaveAngleSave =
   m_bWaveAngleAndHeightSave =
   m_bAvgWaveAngleAndHeightSave =
   m_bDeepWaterWaveAngleAndHeightSave =
   m_bBeachProtectionSave =
   m_bWaveEnergySinceCollapseSave =
   m_bMeanWaveEnergySave =
   m_bBreakingWaveHeightSave =
   m_bPotentialPlatformErosionSave =
   m_bActualPlatformErosionSave =
   m_bTotalPotentialPlatformErosionSave =
   m_bTotalActualPlatformErosionSave =
   m_bPotentialBeachErosionSave =
   m_bActualBeachErosionSave =
   m_bTotalPotentialBeachErosionSave =
   m_bTotalActualBeachErosionSave =
   m_bBeachDepositionSave =
   m_bTotalBeachDepositionSave =
   m_bLandformSave =
   m_bLocalSlopeSave =
   m_bInterventionClassSave =
   m_bInterventionHeightSave =
   m_bSuspSedSave =
   m_bAvgSuspSedSave =
   m_bFineUnconsSedSave =
   m_bSandUnconsSedSave =
   m_bCoarseUnconsSedSave =
   m_bFineConsSedSave =
   m_bSandConsSedSave =
   m_bCoarseConsSedSave =
   m_bRasterCoastlineSave =
   m_bRasterNormalSave =
   m_bActiveZoneSave =
   m_bCliffCollapseSave =
   m_bTotCliffCollapseSave =
   m_bCliffCollapseDepositionSave =
   m_bTotCliffCollapseDepositionSave =
   m_bRasterPolygonSave =
   m_bPotentialPlatformErosionMaskSave =
   m_bSeaMaskSave =
   m_bBeachMaskSave =
   m_bShadowZoneCodesSave =
   m_bSaveRegular =
   m_bCoastSave =
   m_bNormalsSave =
   m_bInvalidNormalsSave =
   m_bCoastCurvatureSave =
   m_bPolygonNodeSave =
   m_bPolygonBoundarySave =
   m_bCliffNotchSave =
   m_bShadowBoundarySave =
   m_bShadowDowndriftBoundarySave =
   m_bDeepWaterWaveAngleSave =
   m_bDeepWaterWaveHeightSave =
   m_bDeepWaterWavePeriodSave =
   m_bPolygonUnconsSedUpOrDownDriftSave =
   m_bPolygonUnconsSedGainOrLossSave =
   m_bSeaAreaTSSave =
   m_bStillWaterLevelTSSave =
   m_bActualPlatformErosionTSSave =
   m_bSuspSedTSSave =
   m_bFloodSetupSurgeTSSave =
   m_bFloodSetupSurgeRunupTSSave =
   m_bCliffCollapseDepositionTSSave =
   m_bCliffCollapseErosionTSSave =
   m_bCliffCollapseNetTSSave =
   m_bBeachErosionTSSave =
   m_bBeachDepositionTSSave =
   m_bBeachSedimentChangeNetTSSave =
   m_bSaveGISThisIter =
   m_bOutputProfileData =
   m_bOutputParallelProfileData =
   m_bOutputErosionPotentialData =
   m_bOmitSearchNorthEdge =
   m_bOmitSearchSouthEdge =
   m_bOmitSearchWestEdge =
   m_bOmitSearchEastEdge =
   m_bErodeShorePlatformAlternateDirection =
   m_bDoShorePlatformErosion =
   m_bDoCliffCollapse =
   m_bDoBeachSedimentTransport =
   m_bGDALCanWriteFloat =
   m_bGDALCanWriteInt32 =
   m_bScaleRasterOutput =
   m_bWorldFile =
   m_bSingleDeepWaterWaveValues =
   m_bHaveWaveStationData =
   m_bSedimentInput =
   m_bSedimentInputAtPoint =
   m_bSedimentInputAtCoast =
   m_bSedimentInputAlongLine =
   m_bSedimentInputThisIter =
   m_bSedimentInputEventSave =
   m_bWaveSetupSave =
   m_bStormSurgeSave =
   m_bDoRiverineFlooding =
   m_bRunUpSave =
   m_bSetupSurgeFloodMaskSave =
   m_bSetupSurgeRunupFloodMaskSave =
   m_bRasterWaveFloodLineSave =
   m_bVectorWaveFloodLineSave =
   m_bFloodLocation =
   m_bFloodSWLSetupLine =
   m_bFloodSWLSetupSurgeLine =
   m_bFloodSWLSetupSurgeRunupLine =
   m_bGISSaveDigitsSequential = false;

   m_bGDALCanCreate = true;

   m_papszGDALRasterOptions =
   m_papszGDALVectorOptions = NULL;

   m_nLayers =
   m_nCoastSmooth =
   m_nCoastSmoothWindow =
   m_nSavGolCoastPoly =
   m_nProfileSmoothWindow =
   m_nCoastNormalAvgSpacing =
   m_nCoastCurvatureInterval =
   m_nNaturalCapeNormals =
   m_nGISMaxSaveDigits =
   m_nGISSave =
   m_nUSave =
   m_nThisSave =
   m_nXGridMax =
   m_nYGridMax =
   m_nCoastMax =
   m_nCoastMin =
   // m_nNThisIterCliffCollapse =
   // m_nNTotCliffCollapse =
   m_nGlobalPolygonID =
   m_nUnconsSedimentHandlingAtGridEdges =
   m_nBeachErosionDepositionEquation =
   m_nWavePropagationModel =
   m_nSimStartSec =
   m_nSimStartMin =
   m_nSimStartHour =
   m_nSimStartDay =
   m_nSimStartMonth =
   m_nSimStartYear =
   m_nDeepWaterWaveDataNumTimeSteps =
   m_nLogFileDetail =
   m_nRunUpEquation = 
   m_nLevel =
   m_nCoastCurvatureMovingWindowSize = 0;

   // TODO 011 May wish to make this a user-supplied value
   m_nMissingValue = INT_NODATA;

   m_nXMinBoundingBox = INT_MAX;
   m_nXMaxBoundingBox = INT_MIN;
   m_nYMinBoundingBox = INT_MAX;
   m_nYMaxBoundingBox = INT_MIN;

   m_GDALWriteIntDataType =
   m_GDALWriteFloatDataType = GDT_Unknown;

   m_lGDALMaxCanWrite =
   m_lGDALMinCanWrite = 0;

   m_ulIter =
   m_ulTotTimestep =
   m_ulThisIterNumPotentialBeachErosionCells =
   m_ulThisIterNumActualBeachErosionCells =
   m_ulThisIterNumBeachDepositionCells =
   m_ulTotPotentialPlatformErosionOnProfiles =
   m_ulTotPotentialPlatformErosionBetweenProfiles =
   m_ulMissingValueBasementCells = 0;
   m_ulNumCells =
   m_ulThisIterNumSeaCells =
   m_ulThisIterNumCoastCells =
   m_ulThisIterNumPotentialPlatformErosionCells =
   m_ulThisIterNumActualPlatformErosionCells = 0;

   for (int i = 0; i < NRNG; i++)
      m_ulRandSeed[i] = 0;

   for (int i = 0; i < SAVEMAX; i++)
      m_dUSaveTime[i] = 0;

   m_dDurationUnitsMult =
   m_dNorthWestXExtCRS =
   m_dNorthWestYExtCRS =
   m_dSouthEastXExtCRS =
   m_dSouthEastYExtCRS =
   m_dExtCRSGridArea =
   m_dCellSide =
   m_dCellDiagonal =
   m_dInvCellSide =
   m_dInvCellDiagonal =
   m_dCellArea =
   m_dSimDuration =
   m_dTimeStep =
   m_dSimElapsed =
   m_dRegularSaveTime =
   m_dRegularSaveInterval =
   m_dClkLast =
   m_dCPUClock =
   m_dSeaWaterDensity =
   m_dThisIterSWL =
   m_dOrigSWL =
   m_dFinalSWL =
   m_dDeltaSWLPerTimestep =
   m_dBreakingWaveHeight =
   m_dC_0 =
   m_dL_0 =
   m_dWaveDepthRatioForWaveCalcs =
   m_dAllCellsDeepWaterWaveHeight =
   m_dAllCellsDeepWaterWaveAngle =
   m_dAllCellsDeepWaterWavePeriod =
   m_dMaxUserInputWaveHeight =
   m_dMaxUserInputWavePeriod =
   m_dR =
   m_dD50Fine =
   m_dD50Sand =
   m_dD50Coarse =
   m_dBeachSedimentDensity =
   m_dBeachSedimentPorosity =
   m_dFineErodibility =
   m_dSandErodibility =
   m_dCoarseErodibility =
   m_dFineErodibilityNormalized =
   m_dSandErodibilityNormalized =
   m_dCoarseErodibilityNormalized =
   m_dKLS =
   m_dKamphuis =
   m_dG =
   m_dInmersedToBulkVolumetric =
   m_dDepthOfClosure =
   m_dCoastNormalAvgSpacing =
   m_dCoastNormalLength =
   m_dThisIterTotSeaDepth =
   m_dThisIterPotentialSedLostBeachErosion =
   m_dThisIterLeftGridUnconsFine =                 // TODO 067
   m_dThisIterLeftGridUnconsSand =
   m_dThisIterLeftGridUnconsCoarse =
   m_dThisIterPotentialPlatformErosion =
   m_dThisIterActualPlatformErosionFineCons =
   m_dThisIterActualPlatformErosionSandCons =
   m_dThisIterActualPlatformErosionCoarseCons =
   m_dThisIterPotentialBeachErosion =
   m_dThisIterBeachErosionFine =
   m_dThisIterBeachErosionSand =
   m_dThisIterBeachErosionCoarse =
   m_dThisIterBeachDepositionSand =
   m_dThisIterBeachDepositionCoarse =
   m_dThisIterFineSedimentToSuspension =
   m_dDepositionSandDiff =
   m_dDepositionCoarseDiff =
   m_dDepthOverDBMax =
   m_dTotPotentialPlatformErosionOnProfiles =
   m_dTotPotentialPlatformErosionBetweenProfiles =
   m_dProfileMaxSlope =
   m_dMaxBeachElevAboveSWL =
   m_dCliffErosionResistance =
   m_dNotchDepthAtCollapse =
   m_dNotchBaseBelowSWL =
   m_dCliffDepositionA =
   m_dCliffDepositionPlanviewWidth =
   m_dCliffTalusMinDepositionLength =
   m_dMinCliffTalusHeightFrac =
   m_dThisIterCliffCollapseErosionFineUncons =
   m_dThisIterCliffCollapseErosionSandUncons =
   m_dThisIterCliffCollapseErosionCoarseUncons =
   m_dThisIterCliffCollapseErosionFineCons =
   m_dThisIterCliffCollapseErosionSandCons =
   m_dThisIterCliffCollapseErosionCoarseCons =
   m_dThisIterUnconsSandCliffDeposition =
   m_dThisIterUnconsCoarseCliffDeposition =
   m_dThisIterCliffCollapseFineErodedDuringDeposition =
   m_dThisIterCliffCollapseSandErodedDuringDeposition =
   m_dThisIterCliffCollapseCoarseErodedDuringDeposition =
   m_dCoastNormalRandSpacingFactor =
   m_dDeanProfileStartAboveSWL =
   m_dAccumulatedSeaLevelChange =
   m_dBreakingWaveHeightDepthRatio =
   m_dWaveDataWrapHours =
   m_dThisIterTopElevMax =
   m_dThisIterTopElevMin =
   m_dThisiterUnconsFineInput =
   m_dThisiterUnconsSandInput =
   m_dThisiterUnconsCoarseInput = 
   m_dStartIterSuspFineAllCells =
   m_dStartIterSuspFineInPolygons =
   m_dStartIterUnconsFineAllCells =
   m_dStartIterUnconsSandAllCells =
   m_dStartIterUnconsCoarseAllCells =
   m_dStartIterConsFineAllCells =
   m_dStartIterConsSandAllCells =
   m_dStartIterConsCoarseAllCells =
   m_dThisIterDiffTotWaterLevel = 
   m_dThisIterDiffWaveSetupWaterLevel = 
   m_dThisIterDiffWaveSetupSurgeWaterLevel = 
   m_dThisIterDiffWaveSetupSurgeRunupWaterLevel =
   m_dTotalFineUnconsInPolygons =
   m_dTotalSandUnconsInPolygons =
   m_dTotalCoarseUnconsInPolygons =
   m_dUnconsSandNotDepositedLastIter =
   m_dUnconsCoarseNotDepositedLastIter = 0;

   m_dMinSWL = DBL_MAX;
   m_dMaxSWL = DBL_MIN;

   for (int i = 0; i < 6; i++)
      m_dGeoTransform[i] = 0;

   // TODO 011 May wish to make this a user-supplied value
   m_dMissingValue = DBL_NODATA;

   m_ldGTotPotentialPlatformErosion =
   m_ldGTotFineActualPlatformErosion =
   m_ldGTotSandActualPlatformErosion =
   m_ldGTotCoarseActualPlatformErosion =
   m_ldGTotPotentialSedLostBeachErosion =
   m_ldGTotActualFineLostBeachErosion =
   m_ldGTotActualSandLostBeachErosion =
   m_ldGTotActualCoarseLostBeachErosion =
   m_ldGTotSandSedLostCliffCollapse =
   m_ldGTotCoarseSedLostCliffCollapse =
   m_ldGTotCliffCollapseFine =
   m_ldGTotCliffCollapseSand =
   m_ldGTotCliffCollapseCoarse =
   m_ldGTotCliffTalusFineToSuspension =
   m_ldGTotCliffTalusSandDeposition =
   m_ldGTotCliffTalusCoarseDeposition =
   m_ldGTotCliffCollapseFineErodedDuringDeposition =
   m_ldGTotCliffCollapseSandErodedDuringDeposition =
   m_ldGTotCliffCollapseCoarseErodedDuringDeposition =
   m_ldGTotPotentialBeachErosion =
   m_ldGTotActualFineBeachErosion =
   m_ldGTotActualSandBeachErosion =
   m_ldGTotActualCoarseBeachErosion =
   m_ldGTotSandBeachDeposition =
   m_ldGTotCoarseBeachDeposition =
   m_ldGTotSuspendedSediment =
   m_ldGTotSandDepositionDiff =
   m_ldGTotCoarseDepositionDiff =
   m_ldGTotFineSedimentInput =
   m_ldGTotSandSedimentInput =
   m_ldGTotCoarseSedimentInput = 0;

   for (int i = 0; i < 2; i++)
   {
      m_ulRState[i].s1 =
      m_ulRState[i].s2 =
      m_ulRState[i].s3 = 0;
   }

   m_tSysStartTime =
   m_tSysEndTime = 0;

   m_pRasterGrid = NULL;
}

//===============================================================================================================================
//! The CSimulation destructor
//===============================================================================================================================
CSimulation::~CSimulation (void)
{
   // Close output files if open
   if (LogStream && LogStream.is_open())
   {
      LogStream.flush();
      LogStream.close();
   }

   if (OutStream && OutStream.is_open())
   {
      OutStream.flush();
      OutStream.close();
   }

   if (SeaAreaTSStream && SeaAreaTSStream.is_open())
   {
      SeaAreaTSStream.flush();
      SeaAreaTSStream.close();
   }

   if (StillWaterLevelTSStream && StillWaterLevelTSStream.is_open())
   {
      StillWaterLevelTSStream.flush();
      StillWaterLevelTSStream.close();
   }

   if (PlatformErosionTSStream && PlatformErosionTSStream.is_open())
   {
      PlatformErosionTSStream.flush();
      PlatformErosionTSStream.close();
   }

   if (CliffCollapseErosionTSStream && CliffCollapseErosionTSStream.is_open())
   {
      CliffCollapseErosionTSStream.flush();
      CliffCollapseErosionTSStream.close();
   }

   if (CliffCollapseDepositionTSStream && CliffCollapseDepositionTSStream.is_open())
   {
      CliffCollapseDepositionTSStream.flush();
      CliffCollapseDepositionTSStream.close();
   }

   if (CliffCollapseNetChangeTSStream && CliffCollapseNetChangeTSStream.is_open())
   {
      CliffCollapseNetChangeTSStream.flush();
      CliffCollapseNetChangeTSStream.close();
   }

   if (FineSedSuspensionTSStream && FineSedSuspensionTSStream.is_open())
   {
      FineSedSuspensionTSStream.flush();
      FineSedSuspensionTSStream.close();
   }

   if (FloodSetupSurgeTSStream && FloodSetupSurgeTSStream.is_open())
   {
      FloodSetupSurgeTSStream.flush();
      FloodSetupSurgeTSStream.close();
   }

   if (FloodSetupSurgeRunupTSStream && FloodSetupSurgeRunupTSStream.is_open())
   {
      FloodSetupSurgeRunupTSStream.flush();
      FloodSetupSurgeRunupTSStream.close();
   }

   if (m_pRasterGrid)
      delete m_pRasterGrid;
}

//===============================================================================================================================
//! Returns the double Missing Value code
//===============================================================================================================================
double CSimulation::dGetMissingValue (void) const
{
   return m_dMissingValue;
}

//===============================================================================================================================
//! Returns the still water level (SWL)
//===============================================================================================================================
double CSimulation::dGetThisIterSWL (void) const
{
   return m_dThisIterSWL;
}

//===============================================================================================================================
//! Returns the this-iteration total water level
//===============================================================================================================================
double CSimulation::dGetThisIterTotWaterLevel (void) const
{
   return m_dThisIterDiffTotWaterLevel;
}

// //===============================================================================================================================
// //! Returns the max elevation of the beach above SWL
// //===============================================================================================================================
// double CSimulation::dGetMaxBeachElevAboveSWL (void) const
// {
//    return m_dMaxBeachElevAboveSWL;
// }

//===============================================================================================================================
// Returns the cell side length
//===============================================================================================================================
// double CSimulation::dGetCellSide(void) const
// {
//    return m_dCellSide;
// }

//===============================================================================================================================
//! Returns X grid max
//===============================================================================================================================
int CSimulation::nGetGridXMax (void) const
{
   return m_nXGridMax;
}

//===============================================================================================================================
//! Returns Y grid max
//===============================================================================================================================
int CSimulation::nGetGridYMax (void) const
{
   return m_nYGridMax;
}

//===============================================================================================================================
//! Returns D50 for fine sediment
//===============================================================================================================================
double CSimulation::dGetD50Fine (void) const
{
   return m_dD50Fine;
}

//===============================================================================================================================
//! Returns D50 for sand sediment
//===============================================================================================================================
double CSimulation::dGetD50Sand (void) const
{
   return m_dD50Sand;
}

//===============================================================================================================================
//! Returns D50 for coarse sediment
//===============================================================================================================================
double CSimulation::dGetD50Coarse (void) const
{
   return m_dD50Coarse;
}

//===============================================================================================================================
//! The nDoSimulation member function of CSimulation sets up and runs the simulation
//===============================================================================================================================
int CSimulation::nDoSimulation(int nArg, char const* pcArgv[])
{
#ifdef RANDCHECK
   CheckRand();
   return RTN_OK;
#endif

   // ================================================== initialization section ================================================
   // Hello, World!
   AnnounceStart();

   // Start the clock ticking
   StartClock();

   // Find out the folder in which the CoastalME executable sits, in order to open the .ini file (they are assumed to be in the same folder)
   if (! bFindExeDir(pcArgv[0]))
      return (RTN_ERR_CMEDIR);

   // Deal with command-line parameters
   int nRet = nHandleCommandLineParams (nArg, pcArgv);
   if (nRet != RTN_OK)
      return (nRet);

   // OK, we are off, tell the user about the licence and the start time
   AnnounceLicence();

   // Read the .ini file and get the name of the run-data file, and path for output etc.
   if (! bReadIniFile())
      return (RTN_ERR_INI);

   // We have the name of the run-data input file, so read it
   if (! bReadRunDataFile())
      return RTN_ERR_RUNDATA;

   // Check raster GIS output format
   if (! bCheckRasterGISOutputFormat())
      return (RTN_ERR_RASTER_GIS_OUT_FORMAT);

   // Check vector GIS output format
   if (! bCheckVectorGISOutputFormat())
      return (RTN_ERR_VECTOR_GIS_OUT_FORMAT);

   // Open log file
   if (! bOpenLogFile())
      return (RTN_ERR_LOGFILE);

   // Set up the time series output files
   if (! bSetUpTSFiles())
      return (RTN_ERR_TSFILE);

   // Initialize the random number generators
   InitRand0 (m_ulRandSeed[0]);
   InitRand1 (m_ulRandSeed[1]);

   // If we are doing Savitzky-Golay smoothing of the vector coastline(s), calculate the filter coefficients
   if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
      CalcSavitzkyGolayCoeffs();

   // Create the raster grid object
   m_pRasterGrid = new CGeomRasterGrid (this);

   // Read in the basement layer (must have this file), create the raster grid, then read in the basement DEM data to the array
   AnnounceReadBasementDEM();
   nRet = nReadRasterBasementDEM();
   if (nRet != RTN_OK)
      return nRet;

   // If we are simulating cliff collapse: then now that we have a value for m_dCellSide, we can check some more input parameters. Talus must be more than one cell wide, and since the number of cells must be odd, three cells is the minimum width
   if (m_bDoCliffCollapse)
   {
      int nTmp = nConvertMetresToNumCells(m_dCliffDepositionPlanviewWidth);
      if (nTmp < 3)
      {
         string strErr = ERR + "cliff deposition must have a planview width of at least three cells. The current setting of " + to_string(m_dCliffDepositionPlanviewWidth) + " m gives a planview width of " + to_string(nTmp) + " cells. Please edit " + m_strDataPathName;
         cerr << strErr << endl;
         LogStream << strErr << endl;
         OutStream << strErr << endl;
         return RTN_ERR_RUNDATA;
      }
   }

   // Do some more initialisation
   m_ulNumCells = m_nXGridMax * m_nYGridMax;

   // Mark edge cells, as defined by the basement layer
   nRet = nMarkBoundingBoxEdgeCells();
   if (nRet != RTN_OK)
      return nRet;

   //    // DEBUG CODE =================================================
   //    for (int n = 0; n < m_VEdgeCell.size(); n++)
   //    {
   //       LogStream << "[" << m_VEdgeCell[n].nGetX() << "][" << m_VEdgeCell[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VEdgeCell[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VEdgeCell[n].nGetY()) << "} " << m_VEdgeCellEdge[n] << endl;
   //    }
   //    // DEBUG CODE =================================================

   // If we are using the default cell spacing, then now that we know the size of the raster cells, we can set the size of profile spacing in m
   if (bFPIsEqual(m_dCoastNormalAvgSpacing, 0.0, TOLERANCE))
      m_dCoastNormalAvgSpacing = MIN_PROFILE_SPACING * m_dCellSide;
   else
   {
      // The user specified a profile spacing, is this too small?
      m_nCoastNormalAvgSpacing = nRound (m_dCoastNormalAvgSpacing / m_dCellSide);

      if (m_nCoastNormalAvgSpacing < MIN_PROFILE_SPACING)
      {
         cerr << ERR << "profile spacing was specified as " << m_dCoastNormalAvgSpacing << " m, which is " << m_nCoastNormalAvgSpacing << " cells. Polygon creation works poorly if profile spacing is less than " << MIN_PROFILE_SPACING << " cells, i.e. " << MIN_PROFILE_SPACING *m_dCellSide << " m" << endl;

         LogStream << ERR << "profile spacing was specified as " << m_dCoastNormalAvgSpacing << " m, which is " << m_nCoastNormalAvgSpacing << " cells. Polygon creation works poorly if profile spacing is less than " << MIN_PROFILE_SPACING << " cells, i.e. " << MIN_PROFILE_SPACING *m_dCellSide << " m" << endl;

         return RTN_ERR_PROFILESPACING;
      }
   }

   // We have at least one filename for the first layer, so add the correct number of layers. Note the the number of layers does not change during the simulation: however layers can decrease in thickness until they have zero thickness
   AnnounceAddLayers();

   for (int nX = 0; nX < m_nXGridMax; nX++)
      for (int nY = 0; nY < m_nYGridMax; nY++)
         m_pRasterGrid->m_Cell[nX][nY].AppendLayers (m_nLayers);

   // Tell the user what is happening then read in the layer files
   AnnounceReadRasterFiles();

   for (int nLayer = 0; nLayer < m_nLayers; nLayer++)
   {
      // Read in the initial fine unconsolidated sediment depth file(s)
      AnnounceReadInitialFineUnconsSedGIS (nLayer);
      nRet = nReadRasterGISFile (FINE_UNCONS_RASTER, nLayer);
      if (nRet != RTN_OK)
         return (nRet);

      // Read in the initial sand unconsolidated sediment depth file
      AnnounceReadInitialSandUnconsSedGIS (nLayer);
      nRet = nReadRasterGISFile (SAND_UNCONS_RASTER, nLayer);
      if (nRet != RTN_OK)
         return (nRet);

      // Read in the initial coarse unconsolidated sediment depth file
      AnnounceReadInitialCoarseUnconsSedGIS (nLayer);
      nRet = nReadRasterGISFile (COARSE_UNCONS_RASTER, nLayer);
      if (nRet != RTN_OK)
         return (nRet);

      // Read in the initial fine consolidated sediment depth file
      AnnounceReadInitialFineConsSedGIS (nLayer);
      nRet = nReadRasterGISFile (FINE_CONS_RASTER, nLayer);
      if (nRet != RTN_OK)
         return (nRet);

      // Read in the initial sand consolidated sediment depth file
      AnnounceReadInitialSandConsSedGIS (nLayer);
      nRet = nReadRasterGISFile (SAND_CONS_RASTER, nLayer);
      if (nRet != RTN_OK)
         return (nRet);

      // Read in the initial coarse consolidated sediment depth file
      AnnounceReadInitialCoarseConsSedGIS (nLayer);
      nRet = nReadRasterGISFile (COARSE_CONS_RASTER, nLayer);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Read in the initial suspended sediment depth file
   AnnounceReadInitialSuspSedGIS();
   nRet = nReadRasterGISFile (SUSP_SED_RASTER, 0);
   if (nRet != RTN_OK)
      return (nRet);

   // Maybe read in the landform class data, otherwise calculate this during the first timestep using identification rules
   if (! m_strInitialLandformFile.empty())
   {
      AnnounceReadLGIS();
      nRet = nReadRasterGISFile (LANDFORM_RASTER, 0);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Maybe read in intervention data
   if (! m_strInterventionClassFile.empty())
   {
      AnnounceReadICGIS();
      nRet = nReadRasterGISFile (INTERVENTION_CLASS_RASTER, 0);
      if (nRet != RTN_OK)
         return (nRet);

      AnnounceReadIHGIS();
      nRet = nReadRasterGISFile (INTERVENTION_HEIGHT_RASTER, 0);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Maybe read in the tide data
   if (! m_strTideDataFile.empty())
   {
      AnnounceReadTideData();
      nRet = nReadTideDataFile();
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Read in the erosion potential shape function data
   AnnounceReadSCAPEShapeFunctionFile();
   nRet = nReadShapeFunctionFile();
   if (nRet != RTN_OK)
      return (nRet);

   // Do we want to output the erosion potential look-up values, for checking purposes?
   if (m_bOutputErosionPotentialData)
      WriteLookUpData();

   // OK, now read in the vector files (if any)
   if (m_bHaveWaveStationData || m_bSedimentInput)
      AnnounceReadVectorFiles();

   // Maybe read in deep water wave station data
   if (m_bHaveWaveStationData)
   {
      // We are reading deep water wave height, orientation and period from a file of vector points and file time series
      AnnounceReadDeepWaterWaveValuesGIS();

      // Read in vector points
      nRet = nReadVectorGISFile (DEEP_WATER_WAVE_STATIONS_VEC);
      if (nRet != RTN_OK)
         return (nRet);

      int nWaveStations = static_cast<int> (m_VnDeepWaterWaveStationID.size());

      if (nWaveStations == 1)
         m_bSingleDeepWaterWaveValues = true;

      // Read in time series values, and initialize the vector which stores each timestep's deep water wave height, orientation and period
      nRet = nReadWaveStationTimeSeriesFile (nWaveStations);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Maybe read in sediment input event data
   if (m_bSedimentInput)
   {
      // We are reading sediment input event data
      AnnounceReadSedimentEventInputValuesGIS();

      // Read in vector points for sediment input events
      nRet = nReadVectorGISFile (SEDIMENT_INPUT_EVENT_LOCATION_VEC);
      if (nRet != RTN_OK)
         return (nRet);

      // Read in the time series values for sediment input events
      nRet = nReadSedimentInputEventTimeSeriesFile();
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Maybe read in flood input location
   if (m_bFloodLocation)
   {
      // We are reading sediment input event data
      AnnounceReadFloodLocationGIS();

      // Read in vector points for sediment input events
      nRet = nReadVectorGISFile (FLOOD_LOCATION_VEC);
      if (nRet != RTN_OK)
         return (nRet);
   }

   // Open the main output
   OutStream.open (m_strOutFile.c_str(), ios::out | ios::trunc);

   if (! OutStream)
   {
      // Error, cannot open Out file
      cerr << ERR << "cannot open " << m_strOutFile << " for output" << endl;
      return (RTN_ERR_OUTFILE);
   }

   // Write beginning-of-run information to Out and Log files
   WriteStartRunDetails();

   // Start initializing
   AnnounceInitializing();

   // Misc initialization calcs
   m_nCoastMax = COAST_LENGTH_MAX * tMax(m_nXGridMax, m_nYGridMax);           // Arbitrary but probably OK
   m_nCoastMin = tMin(m_nXGridMax, m_nYGridMax);                              // In some cases the following rule doesn't work TODO 007 Info needed
   // nRound(COAST_LENGTH_MIN_X_PROF_SPACE * m_dCoastNormalAvgSpacing / m_dCellSide);           // TODO 007 Info needed
   m_nCoastCurvatureInterval = tMax(nRound(m_dCoastNormalAvgSpacing / (m_dCellSide * 2)), 2);   // TODO 007 Info needed

   // For beach erosion/deposition, conversion from immersed weight to bulk volumetric (sand and voids) transport rate (Leo Van Rijn) TODO 007 need full reference
   m_dInmersedToBulkVolumetric = 1 / ((m_dBeachSedimentDensity - m_dSeaWaterDensity) * (1 - m_dBeachSedimentPorosity) * m_dG);

   m_bConsChangedThisIter.resize(m_nLayers, false);
   m_bUnconsChangedThisIter.resize(m_nLayers, false);

   // Normalize sediment erodibility values, so that none are > 1
   double dTmp = m_dFineErodibility + m_dSandErodibility + m_dCoarseErodibility;
   m_dFineErodibilityNormalized = m_dFineErodibility / dTmp;
   m_dSandErodibilityNormalized = m_dSandErodibility / dTmp;
   m_dCoarseErodibilityNormalized = m_dCoarseErodibility / dTmp;

   // Intialise SWL
   m_dThisIterSWL = m_dOrigSWL;

   // If SWL changes during the simulation, calculate the per-timestep increment (could be -ve)
   if (! bFPIsEqual(m_dFinalSWL, m_dOrigSWL, TOLERANCE))
   {
      m_dDeltaSWLPerTimestep = (m_dTimeStep * (m_dFinalSWL - m_dOrigSWL)) / m_dSimDuration;
      m_dAccumulatedSeaLevelChange -= m_dDeltaSWLPerTimestep;
   }

   // ===================================================== The main loop ======================================================
   // Tell the user what is happening
   AnnounceIsRunning();

   while (true)
   {
      //       // DEBUG CODE =========================================
      //       LogStream << ulGetRand0() << " " << ulGetRand1() << endl;
      //       // DEBUG CODE =========================================

      // Check that we haven't gone on too long: if not then update timestep number etc.
      if (bTimeToQuit())
         break;

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << "TIMESTEP " << m_ulIter << " " << string(154, '=') << endl;
      
      LogStream << std::fixed << setprecision(3);
      
      // Check to see if there is a new intervention in place: if so, update it on the RasterGrid array
      nRet = nUpdateIntervention();
      if (nRet != RTN_OK)
         return nRet;

      // Calculate changes due to external forcing
      nRet = nCalcExternalForcing();
      if (nRet != RTN_OK)
         return nRet;

      // Do per-timestep intialization: set up the grid cells ready for this timestep, also initialize per-timestep totals. Note that in the first timestep, all cells (including hinterland cells) are given the deep water wave values
      nRet = nInitGridAndCalcStillWaterLevel();
      if (nRet != RTN_OK)
         return nRet;

      // Next find out which cells are inundated and locate the coastline(s). This also gives to all sea cells, wave values which are the same as the deep water values. For shallow water sea cells, these wave values will be changed later, in nDoAllPropagateWaves()
      int nValidCoast = 0;
      nRet = nLocateSeaAndCoasts(nValidCoast);
      if (nRet != RTN_OK)
         return nRet;
      
      // Tell the user how the simulation is progressing
      AnnounceProgress();

      // Locate estuaries TODO 044 someday...

      // Sort out hinterland landforms
      nRet = nAssignNonCoastlineLandforms();
      if (nRet != RTN_OK)
         return nRet;

      // For each coastline, use classification rules to assign landform categories
      nRet = nAssignAllCoastalLandforms();
      if (nRet != RTN_OK)
         return nRet;

      // If we have sediment input events, then check to see whether this is time for an event to occur. If it is, then do it
      if (m_bSedimentInput)
      {
         m_bSedimentInputThisIter = false;

         nRet = nCheckForSedimentInputEvent();
         if (nRet != RTN_OK)
            return nRet;
      }

      // Create the coastline-normal profiles
      nRet = nCreateAllProfilesAndCheckForIntersection();
      if (nRet != RTN_OK)
         return nRet;

      // Tell the user how the simulation is progressing
      AnnounceProgress();
      
      // Create the coast polygons
      nRet = nCreateAllPolygons();
      if (nRet != RTN_OK)
         return nRet;

      //       // DEBUG CODE =========================================
      //       int nNODATA = 0;
      //       int nPoly0 = 0;
      //       for (int nX = 0; nX < m_nXGridMax; nX++)
      //       {
      //          for (int nY = 0; nY < m_nYGridMax; nY++)
      //          {
      //             int nTmp = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
      //             if (nTmp == INT_NODATA)
      //                nNODATA++;
      //             if (nTmp == 0)
      //                nPoly0++;
      //          }
      //       }
      //       LogStream << "Before marking polygon cells, N cells with NODATA polygon ID = " << nNODATA << endl;
      //       LogStream << "Before marking polygon cells, N cells with zero polygon ID = " << nPoly0 << endl;
      //       // DEBUG CODE =========================================

      // Mark cells of the raster grid that are within each polygon
      MarkPolygonCells();

      // Calculate the length of the shared normal between each polygon and the adjacent polygon(s)
      nRet = nDoPolygonSharedBoundaries();
      if (nRet != RTN_OK)
         return nRet;

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      //       // DEBUG CODE =========================================
      //       nNODATA = 0;
      //       nPoly0 = 0;
      //       for (int nX = 0; nX < m_nXGridMax; nX++)
      //       {
      //          for (int nY = 0; nY < m_nYGridMax; nY++)
      //          {
      //             int nTmp = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
      //             if (nTmp == INT_NODATA)
      //                nNODATA++;
      //             if (nTmp == 0)
      //                nPoly0++;
      //          }
      //       }
      //       LogStream << "After marking polygon cells, N cells with NODATA polygon ID = " << nNODATA << endl;
      //       LogStream << "After marking polygon cells, N cells with zero polygon ID = " << nPoly0 << endl;
      //       // DEBUG CODE =========================================

      // PropagateWind();

      // Give every coast point a value for deep water wave height and direction
      nRet = nSetAllCoastpointDeepWaterWaveValues();
      if (nRet != RTN_OK)
         return nRet;
      
      // Change the wave properties in all shallow water sea cells: propagate waves and define the active zone, also locate wave shadow zones
      nRet = nDoAllPropagateWaves();
      if (nRet != RTN_OK)
         return nRet;


      // Output polygon share table and pre-existing sediment table to log file
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         for (int nCoast = 0; nCoast < nValidCoast; nCoast++)
         {
            WritePolygonShareTable(nCoast);
            WritePolygonPreExistingSediment(nCoast);
         }
      }

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      //       // DEBUG CODE ===========================================
      //       string strOutFile = m_strOutPath;
      //       strOutFile += "sea_wave_height_CHECKPOINT_";
      //       strOutFile += std::to_string(m_ulIter);
      //       strOutFile += ".tif";
      //
      //       GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      //       GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
      //       pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      //       pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      //       int nn = 0;
      //       double* pdRaster = new double[m_nXGridMax * m_nYGridMax];
      //       for (int nY = 0; nY < m_nYGridMax; nY++)
      //       {
      //          for (int nX = 0; nX < m_nXGridMax; nX++)
      //          {
      //             pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
      //          }
      //       }
      //
      //       GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
      //       pBand->SetNoDataValue(m_dMissingValue);
      //       int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
      //
      //       if (nRet == CE_Failure)
      //          return RTN_ERR_GRIDCREATE;
      //
      //       GDALClose(pDataSet);
      //       delete[] pdRaster;
      //       // DEBUG CODE ===========================================
      //
      //       // DEBUG CODE ===========================================
      //       strOutFile = m_strOutPath;
      //       strOutFile += "sea_wave_angle_CHECKPOINT_";
      //       strOutFile += std::to_string(m_ulIter);
      //       strOutFile += ".tif";
      //
      //       pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      //       pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
      //       pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      //       pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      //       nn = 0;
      //       pdRaster = new double[m_nXGridMax * m_nYGridMax];
      //       for (int nY = 0; nY < m_nYGridMax; nY++)
      //       {
      //          for (int nX = 0; nX < m_nXGridMax; nX++)
      //          {
      //             pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
      //          }
      //       }
      //
      //       pBand = pDataSet->GetRasterBand(1);
      //       pBand->SetNoDataValue(m_dMissingValue);
      //       nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
      //
      //       if (nRet == CE_Failure)
      //          return RTN_ERR_GRIDCREATE;
      //
      //       GDALClose(pDataSet);
      //       delete[] pdRaster;
      //       // DEBUG CODE ===========================================

      // Save the not-deposited values, to be shown in the logfile after we've finished beach sediment movement
      m_dUnconsSandNotDepositedLastIter = m_dDepositionSandDiff;
      m_dUnconsCoarseNotDepositedLastIter = m_dDepositionCoarseDiff;

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         if (m_dDepositionSandDiff > MASS_BALANCE_TOLERANCE)
         {
            LogStream << m_ulIter << ": AT ITERATION START m_dDepositionSandDiff = " << m_dDepositionSandDiff * m_dCellArea << " m_dUnconsSandNotDepositedLastIter = " << m_dUnconsSandNotDepositedLastIter << endl;
            LogStream << m_ulIter << ": AT ITERATION START m_dDepositionCoarseDiff = " << m_dDepositionCoarseDiff * m_dCellArea << " m_dUnconsCoarseNotDepositedLastIter = " << m_dUnconsCoarseNotDepositedLastIter << endl;
         }

      if (m_bDoShorePlatformErosion)
      {
         // Calculate elevation change on the consolidated sediment which comprises the coastal platform
         nRet = nDoAllShorePlatFormErosion();
         if (nRet != RTN_OK)
            return nRet;
      }

      // Output shore platform erosion table to log file
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         for (int nCoast = 0; nCoast < nValidCoast; nCoast++)
            WritePolygonShorePlatformErosion(nCoast);
      }

      if (m_bDoCliffCollapse)
      {
         // Do all cliff collapses for this timestep (if any)
         nRet = nDoAllWaveEnergyToCoastLandforms();
         if (nRet != RTN_OK)
            return nRet;
      }

      // Output cliff collapse table to log file
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         for (int nCoast = 0; nCoast < nValidCoast; nCoast++)
            WritePolygonCliffCollapseErosion(nCoast);
      }

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      if (m_bDoBeachSedimentTransport)
      {
         // Next simulate beach erosion and deposition i.e. simulate alongshore transport of unconsolidated sediment (longshore drift) between polygons. First calculate potential sediment movement between polygons
         DoAllPotentialBeachErosion();

         // Do within-sediment redistribution of unconsolidated sediment, constraining potential sediment movement to give actual (i.e. supply-limited) sediment movement to/from each polygon in three size clases
         nRet = nDoAllActualBeachErosionAndDeposition();
         if (nRet != RTN_OK)
            return nRet;
      }

//       // Add the fine sediment that was eroded this timestep (from the shore platform, from cliff collapse, from erosion of existing fine sediment during cliff collapse talus deposition, and from beach erosion; minus the fine sediment from beach erosion that went off-grid) to the suspended sediment load
//       double dFineThisIter = m_dThisIterActualPlatformErosionFineCons + m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons + m_dThisIterCliffCollapseFineErodedDuringDeposition + m_dThisIterBeachErosionFine - m_dThisIterLeftGridUnconsFine;
//
//       m_dThisIterFineSedimentToSuspension += dFineThisIter;
      
      // Tell the user how the simulation is progressing
      AnnounceProgress();

      //       // DEBUG CODE ===========================================
      //       string strOutFile = m_strOutPath;
      //       strOutFile += "sea_wave_height_CHECKPOINT_";
      //       strOutFile += std::to_string(m_ulIter);
      //       strOutFile += ".tif";
      //
      //       GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      //       GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
      //       pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      //       pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      //       int nn = 0;
      //       double* pdRaster = new double[m_nXGridMax * m_nYGridMax];
      //       for (int nY = 0; nY < m_nYGridMax; nY++)
      //       {
      //          for (int nX = 0; nX < m_nXGridMax; nX++)
      //          {
      //             pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
      //          }
      //       }
      //
      //       GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
      //       pBand->SetNoDataValue(m_dMissingValue);
      //       int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
      //
      //       if (nRet == CE_Failure)
      //          return RTN_ERR_GRIDCREATE;
      //
      //       GDALClose(pDataSet);
      //       delete[] pdRaster;
      //       // DEBUG CODE ===========================================
      //
      //       // DEBUG CODE ===========================================
      //       strOutFile = m_strOutPath;
      //       strOutFile += "sea_wave_angle_CHECKPOINT_";
      //       strOutFile += std::to_string(m_ulIter);
      //       strOutFile += ".tif";
      //
      //       pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
      //       pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
      //       pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
      //       pDataSet->SetGeoTransform(m_dGeoTransform);
      //
      //       nn = 0;
      //       pdRaster = new double[m_nXGridMax * m_nYGridMax];
      //       for (int nY = 0; nY < m_nYGridMax; nY++)
      //       {
      //          for (int nX = 0; nX < m_nXGridMax; nX++)
      //          {
      //             pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
      //          }
      //       }
      //
      //       pBand = pDataSet->GetRasterBand(1);
      //       pBand->SetNoDataValue(m_dMissingValue);
      //       nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
      //
      //       if (nRet == CE_Failure)
      //          return RTN_ERR_GRIDCREATE;
      //
      //       GDALClose(pDataSet);
      //       delete[] pdRaster;
      //       // DEBUG CODE ===========================================

      // Do some end-of-timestep updates to the raster grid, also update per-timestep and running totals
      nRet = nUpdateGrid();
      if (nRet != RTN_OK)
         return nRet;

      // Make water level inundation on grid
      if (m_bFloodSWLSetupSurgeLine || m_bSetupSurgeFloodMaskSave)
      {
         m_nLevel = 0;
         
         nRet = nLocateFloodAndCoasts();
         if (nRet != RTN_OK)
            return nRet;
      }

      if (m_bFloodSWLSetupSurgeRunupLine || m_bSetupSurgeRunupFloodMaskSave)
      {
         // TODO 007
         m_nLevel = 1;
         
         nRet = nLocateFloodAndCoasts();
         if (nRet != RTN_OK)
            return nRet;
      }

      // Now save results, first the raster and vector GIS files if required
      m_bSaveGISThisIter = false;

      if ((m_bSaveRegular && (m_dSimElapsed >= m_dRegularSaveTime) && (m_dSimElapsed < m_dSimDuration)) || (! m_bSaveRegular && (m_dSimElapsed >= m_dUSaveTime[m_nThisSave])))
      {
         m_bSaveGISThisIter = true;

         // Save the values from the RasterGrid array into raster GIS files
         if (! bSaveAllRasterGISFiles())
            return (RTN_ERR_RASTER_FILE_WRITE);

         // Tell the user how the simulation is progressing
         AnnounceProgress();

         // Save the vector GIS files
         if (! bSaveAllVectorGISFiles())
            return (RTN_ERR_VECTOR_FILE_WRITE);

         // Tell the user how the simulation is progressing
         AnnounceProgress();
      }

      // Output per-timestep results to the .out file
      if (! bWritePerTimestepResults())
         return (RTN_ERR_TEXT_FILE_WRITE);

      // Now output time series CSV stuff
      if (! bWriteTSFiles())
         return (RTN_ERR_TIMESERIES_FILE_WRITE);

      // Tell the user how the simulation is progressing
      AnnounceProgress();

      // Update grand totals
      DoTimestepTotals();
   } // ================================================ End of main loop ======================================================

   // =================================================== post-loop tidying =====================================================
   // Tell the user what is happening
   AnnounceSimEnd();

   // Write end-of-run information to Out, Log and time-series files
   nRet = nWriteEndRunDetails();
   if (nRet != RTN_OK)
      return (nRet);

   return RTN_OK;
}
