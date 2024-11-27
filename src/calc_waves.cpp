/*!
 *
 * \file calc_waves.cpp
 * \brief Simulates wave propagation
 * \details TODO 001 A more detailed description of these routines.
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

#include <cmath>
using std::isfinite;

#include <unistd.h>

#include <string>
using std::stoi;
using std::to_string;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ios;

#include <iomanip>
using std::resetiosflags;
using std::setiosflags;
using std::setprecision;
using std::setw;

#include <algorithm>
using std::remove;
using std::reverse;
using std::sort;

#include "cme.h"
#include "coast.h"
#include "simulation.h"

#include "interpolate.h"

//===============================================================================================================================
//! Give every coast point a value for deep water wave height and direction TODO 005 This may not be realistic, maybe better to use end-of-profile value instead (how?)
//===============================================================================================================================
int CSimulation::nSetAllCoastpointDeepWaterWaveValues(void)
{
   // For each coastline, put a value for deep water wave height and direction at each coastline point
   int nNextProfile;
   int nDistFromPrevProfile = 0;
   int nDistToNextProfile = 0;

   double dPrevProfileDeepWaterWaveHeight = 0;
   double dPrevProfileDeepWaterWaveAngle = 0;
   double dPrevProfileDeepWaterWavePeriod = 0;
   double dNextProfileDeepWaterWaveHeight = 0;
   double dNextProfileDeepWaterWaveAngle = 0;
   double dNextProfileDeepWaterWavePeriod = 0;
   double dDist = 0;

   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nPoint = 0; nPoint < m_VCoast[nCoast].nGetCoastlineSize(); nPoint++)
      {
         // We are going down-coast
         if (m_VCoast[nCoast].bIsProfileStartPoint(nPoint))
         {
            // OK, a coastline-normal profile begins at this coastline point, so set the deep water wave values at this coastline point to be the values at the seaward end of the coastline normal
            int nProfile = m_VCoast[nCoast].nGetProfileNumber(nPoint);
            CGeomProfile const* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

            double dThisDeepWaterWaveHeight = pProfile->dGetProfileDeepWaterWaveHeight();
            double dThisDeepWaterWaveAngle = pProfile->dGetProfileDeepWaterWaveAngle();
            double dThisDeepWaterWavePeriod = pProfile->dGetProfileDeepWaterWavePeriod();

            m_VCoast[nCoast].SetCoastDeepWaterWaveHeight(nPoint, dThisDeepWaterWaveHeight);
            m_VCoast[nCoast].SetCoastDeepWaterWaveAngle(nPoint, dThisDeepWaterWaveAngle);
            m_VCoast[nCoast].SetCoastDeepWaterWavePeriod(nPoint, dThisDeepWaterWavePeriod);

            // Reset for next time
            nDistFromPrevProfile = 0;
            dPrevProfileDeepWaterWaveHeight = dThisDeepWaterWaveHeight;
            dPrevProfileDeepWaterWaveAngle = dThisDeepWaterWaveAngle;
            dPrevProfileDeepWaterWavePeriod = dThisDeepWaterWavePeriod;

            // Find the next profile
            nNextProfile = m_VCoast[nCoast].nGetDownCoastProfileNumber(nProfile);
            if (nNextProfile == INT_NODATA)
            {
               // We are at the end of the coast
               break;
            }

            CGeomProfile const* pNextProfile = m_VCoast[nCoast].pGetProfile(nNextProfile);

            // And the distance (in along-coast points) to the next profile
            nDistToNextProfile = pNextProfile->nGetNumCoastPoint() - nPoint;
            dDist = nDistToNextProfile;

            // And the next profile's deep water wave values
            dNextProfileDeepWaterWaveHeight = pNextProfile->dGetProfileDeepWaterWaveHeight();
            dNextProfileDeepWaterWaveAngle = pNextProfile->dGetProfileDeepWaterWaveAngle();
            dNextProfileDeepWaterWavePeriod = pNextProfile->dGetProfileDeepWaterWavePeriod();

            //             LogStream << m_ulIter << ": coast point = " << nPoint << " IS PROFILE START, dThisDeepWaterWaveHeight = " << dThisDeepWaterWaveHeight << ", dThisDeepWaterWaveAngle = " << dThisDeepWaterWaveAngle << endl;
         }

         else
         {
            // This coast point is not the start of a coastline normal, so set the deep water wave values at this coastline point to be a weighted average of those from the up-coast and down-coast profiles
            nDistFromPrevProfile++;
            nDistToNextProfile--;

            double dPrevWeight = (dDist - nDistFromPrevProfile) / dDist;
            double dNextWeight = (dDist - nDistToNextProfile) / dDist;
            double dThisDeepWaterWaveHeight = (dPrevWeight * dPrevProfileDeepWaterWaveHeight) + (dNextWeight * dNextProfileDeepWaterWaveHeight);
            double dThisDeepWaterWaveAngle = dKeepWithin360((dPrevWeight * dPrevProfileDeepWaterWaveAngle) + (dNextWeight * dNextProfileDeepWaterWaveAngle));
            double dThisDeepWaterWavePeriod = (dPrevWeight * dPrevProfileDeepWaterWavePeriod) + (dNextWeight * dNextProfileDeepWaterWavePeriod);

            m_VCoast[nCoast].SetCoastDeepWaterWaveHeight(nPoint, dThisDeepWaterWaveHeight);
            m_VCoast[nCoast].SetCoastDeepWaterWaveAngle(nPoint, dThisDeepWaterWaveAngle);
            m_VCoast[nCoast].SetCoastDeepWaterWavePeriod(nPoint, dThisDeepWaterWavePeriod);

            //             LogStream << m_ulIter << ": coast point = " << nPoint << " dThisDeepWaterWaveHeight = " << dThisDeepWaterWaveHeight << " dThisDeepWaterWaveAngle = " << dThisDeepWaterWaveAngle << endl;
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Simulates wave propagation along all coastline-normal profiles
//===============================================================================================================================
int CSimulation::nDoAllPropagateWaves(void)
{
   // Set up all-profile vectors to hold the wave attribute data at every profile point on all profiles
   vector<bool> VbBreakingAll;

   vector<double> VdXAll;
   vector<double> VdYAll;
   vector<double> VdHeightXAll;
   vector<double> VdHeightYAll;

   // Calculate wave properties for every coast
   bool bSomeNonStartOrEndOfCoastProfiles = false;
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
      int nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      // Calculate wave properties at every point along each valid profile, and for the cells under the profiles. Do this in the original (curvature-related) profile sequence
      for (int nProfile = 0; nProfile < nNumProfiles; nProfile++)
      {
         vector<bool> VbBreaking;
         vector<double> VdX;
         vector<double> VdY;
         vector<double> VdHeightX;
         vector<double> VdHeightY;

         int nRet = nCalcWavePropertiesOnProfile(nCoast, nCoastSize, nProfile, &VdX, &VdY, &VdHeightX, &VdHeightY, &VbBreaking);
         if (nRet != RTN_OK)
            return nRet;

         // Are the waves off-shore? If so, do nothing more with this profile. The wave values for cells have already been given the off-shore value
         if (VbBreaking.empty())
            continue;

         // Is this a start of coast or end of coast profile?
         if ((! m_VCoast[nCoast].pGetProfile(nProfile)->bStartOfCoast()) && (! m_VCoast[nCoast].pGetProfile(nProfile)->bEndOfCoast()))
         {
            // It is neither a start of coast or an end of coast profile, so set switch
            bSomeNonStartOrEndOfCoastProfiles = true;
         }

         // // DEBUG CODE
         // for (int nn = 0; nn < VdX.size(); nn++)
         // {
         //    LogStream << "nProfile = " << nProfile << " nn = " << nn << " VdX[nn] = " << VdX[nn] << " VdY[nn] = " << VdY[nn] << " VdHeightX[nn] = " << VdHeightX[nn] << " VdHeightY[nn] = " << VdHeightY[nn] << " VbBreaking[nn] = " << VbBreaking[nn] << endl;
         // }
         // LogStream << endl;
         // // DEBUG CODE

         // Append to the all-profile vectors
         VdXAll.insert(VdXAll.end(), VdX.begin(), VdX.end());
         VdYAll.insert(VdYAll.end(), VdY.begin(), VdY.end());
         VdHeightXAll.insert(VdHeightXAll.end(), VdHeightX.begin(), VdHeightX.end());
         VdHeightYAll.insert(VdHeightYAll.end(), VdHeightY.begin(), VdHeightY.end());
         VbBreakingAll.insert(VbBreakingAll.end(), VbBreaking.begin(), VbBreaking.end());
      }
   }

   // OK, do we have some profiles other than start of coast or end of coast profiles in the all-profile vectors? We need to check this, because GDALGridCreate() in nInterpolateWavePropertiesToWithinPolygonCells() does not work if we give it only a start-of-coast or an end-of-coast profile to work with TODO 006 Is this still true?
   if (! bSomeNonStartOrEndOfCoastProfiles)
   {
      LogStream << m_ulIter << ": waves are on-shore only for start and/or end of coast profiles" << endl;

      return RTN_OK;
   }

   // We need to also send the deepwater points from the edge of the grid to nInterpolateWavePropertiesToWithinPolygonCells(), this is necessary to prevent GDALGridCreate() leaving holes in the interpolated grid when the polygons are far from regular
   double dDeepWaterWaveX;
   double dDeepWaterWaveY;

   if (m_bSingleDeepWaterWaveValues)
   {
      // Just using the same value of deep water wave hehght and angle for all cells
      dDeepWaterWaveX = m_dAllCellsDeepWaterWaveHeight * sin(m_dAllCellsDeepWaterWaveAngle * PI / 180),
      dDeepWaterWaveY = m_dAllCellsDeepWaterWaveHeight * cos(m_dAllCellsDeepWaterWaveAngle * PI / 180);
   }

   for (int nX = 0; nX < m_nXGridMax; nX++)
   {
      if (m_pRasterGrid->m_Cell[nX][0].bIsInContiguousSea())
      {
         int nID = m_pRasterGrid->m_Cell[nX][0].nGetPolygonID();
         if (nID == INT_NODATA)
         {
            // Not in a polygon
            VdXAll.push_back(nX);
            VdYAll.push_back(0);

            if (! m_bSingleDeepWaterWaveValues)
            {
               // Not using the same value of deep water height and angle for all cells, so get this cell's deep water height and angle values
               dDeepWaterWaveX = m_pRasterGrid->m_Cell[nX][0].dGetCellDeepWaterWaveHeight() * sin(m_pRasterGrid->m_Cell[nX][0].dGetCellDeepWaterWaveAngle() * PI / 180);
               dDeepWaterWaveY = m_pRasterGrid->m_Cell[nX][0].dGetCellDeepWaterWaveHeight() * cos(m_pRasterGrid->m_Cell[nX][0].dGetCellDeepWaterWaveAngle() * PI / 180);
            }

            VdHeightXAll.push_back(dDeepWaterWaveX);
            VdHeightYAll.push_back(dDeepWaterWaveY);
         }
      }

      if (m_pRasterGrid->m_Cell[nX][m_nYGridMax - 1].bIsInContiguousSea())
      {
         int nID = m_pRasterGrid->m_Cell[nX][m_nYGridMax - 1].nGetPolygonID();
         if (nID == INT_NODATA)
         {
            // Not in a polygon
            VdXAll.push_back(nX);
            VdYAll.push_back(m_nYGridMax - 1);

            if (! m_bSingleDeepWaterWaveValues)
            {
               // Not using the same value of deep water height and angle for all cells, so get this cell's deep water height and angle values
               dDeepWaterWaveX = m_pRasterGrid->m_Cell[nX][m_nYGridMax - 1].dGetCellDeepWaterWaveHeight() * sin(m_pRasterGrid->m_Cell[nX][m_nYGridMax - 1].dGetCellDeepWaterWaveAngle() * PI / 180);
               dDeepWaterWaveY = m_pRasterGrid->m_Cell[nX][m_nYGridMax - 1].dGetCellDeepWaterWaveHeight() * cos(m_pRasterGrid->m_Cell[nX][m_nYGridMax - 1].dGetCellDeepWaterWaveAngle() * PI / 180);
            }

            VdHeightXAll.push_back(dDeepWaterWaveX);
            VdHeightYAll.push_back(dDeepWaterWaveY);
         }
      }
   }

   for (int nY = 0; nY < m_nYGridMax; nY++)
   {
      if (m_pRasterGrid->m_Cell[0][nY].bIsInContiguousSea())
      {
         int nID = m_pRasterGrid->m_Cell[0][nY].nGetPolygonID();
         if (nID == INT_NODATA)
         {
            // Not in a polygon
            VdXAll.push_back(0);
            VdYAll.push_back(nY);

            if (! m_bSingleDeepWaterWaveValues)
            {
               // Not using the same value of deep water height and angle for all cells, so get this cell's deep water height and angle values
               dDeepWaterWaveX = m_pRasterGrid->m_Cell[0][nY].dGetCellDeepWaterWaveHeight() * sin(m_pRasterGrid->m_Cell[0][nY].dGetCellDeepWaterWaveAngle() * PI / 180);
               dDeepWaterWaveY = m_pRasterGrid->m_Cell[0][nY].dGetCellDeepWaterWaveHeight() * cos(m_pRasterGrid->m_Cell[0][nY].dGetCellDeepWaterWaveAngle() * PI / 180);
            }

            VdHeightXAll.push_back(dDeepWaterWaveX);
            VdHeightYAll.push_back(dDeepWaterWaveY);
         }
      }

      if (m_pRasterGrid->m_Cell[m_nXGridMax - 1][nY].bIsInContiguousSea())
      {
         int nID = m_pRasterGrid->m_Cell[m_nXGridMax - 1][nY].nGetPolygonID();
         if (nID == INT_NODATA)
         {
            // Not in a polygon
            VdXAll.push_back(m_nXGridMax - 1);
            VdYAll.push_back(nY);

            if (! m_bSingleDeepWaterWaveValues)
            {
               // Not using the same value of deep water height and angle for all cells, so get this cell's deep water height and angle values
               dDeepWaterWaveX = m_pRasterGrid->m_Cell[m_nXGridMax - 1][nY].dGetCellDeepWaterWaveHeight() * sin(m_pRasterGrid->m_Cell[m_nXGridMax - 1][nY].dGetCellDeepWaterWaveAngle() * PI / 180);
               dDeepWaterWaveY = m_pRasterGrid->m_Cell[m_nXGridMax - 1][nY].dGetCellDeepWaterWaveHeight() * cos(m_pRasterGrid->m_Cell[m_nXGridMax - 1][nY].dGetCellDeepWaterWaveAngle() * PI / 180);
            }

            VdHeightXAll.push_back(dDeepWaterWaveX);
            VdHeightYAll.push_back(dDeepWaterWaveY);
         }
      }
   }

   //    // DEBUG CODE ============================================
   //    LogStream << "Out of loop" << endl;
   //    for (int nn = 0; nn < VdXAll.size(); nn++)
   //    {
   //       LogStream << "nn = " << nn << " VdXAll[nn] = " << VdXAll[nn] << " VdYAll[nn] = " << VdYAll[nn] << " VdHeightXAll[nn] = " << VdHeightXAll[nn] << " VdHeightYAll[nn] = " << VdHeightYAll[nn] << " VbBreakingAll[nn] = " << VbBreakingAll[nn] << endl;
   //    }
   //    LogStream << endl;
   //    // DEBUG CODE ============================================

   // Are the waves off-shore for every profile? If so, do nothing more
   if (VbBreakingAll.empty())
   {
      LogStream << m_ulIter << ": waves off-shore for all profiles" << endl;

      return RTN_OK;
   }

   // Some waves are on-shore, so interpolate the wave attributes from all profile points to all within-polygon sea cells, also update the active zone status for each cell
   int nRet = nInterpolateWavesToPolygonCells(&VdXAll, &VdYAll, &VdHeightXAll, &VdHeightYAll);
   if (nRet != RTN_OK)
      return nRet;

   //    // DEBUG CODE ===========================================
   //    string strOutFile = m_strOutPath;
   //    strOutFile += "sea_wave_height_CHECKPOINT_";
   //    strOutFile += std::to_string(m_ulIter);
   //    strOutFile += ".tif";
   //
   //    GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   //    GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
   //    pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   //    pDataSet->SetGeoTransform(m_dGeoTransform);
   //
   //    int nn = 0;
   //    double* pdRaster = new double[m_nXGridMax * m_nYGridMax];
   //    for (int nY = 0; nY < m_nYGridMax; nY++)
   //    {
   //       for (int nX = 0; nX < m_nXGridMax; nX++)
   //       {
   //          pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
   //       }
   //    }
   //
   //    GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   //    pBand->SetNoDataValue(m_dMissingValue);
   //    nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
   //
   //    if (nRet == CE_Failure)
   //       return RTN_ERR_GRIDCREATE;
   //
   //    GDALClose(pDataSet);
   //    delete[] pdRaster;
   //    // DEBUG CODE ===========================================

   //    // DEBUG CODE ===========================================
   //    strOutFile = m_strOutPath;
   //    strOutFile += "sea_wave_angle_CHECKPOINT_";
   //    strOutFile += std::to_string(m_ulIter);
   //    strOutFile += ".tif";
   //
   //    pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   //    pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
   //    pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   //    pDataSet->SetGeoTransform(m_dGeoTransform);
   //
   //    nn = 0;
   //    pdRaster = new double[m_nXGridMax * m_nYGridMax];
   //    for (int nY = 0; nY < m_nYGridMax; nY++)
   //    {
   //       for (int nX = 0; nX < m_nXGridMax; nX++)
   //       {
   //          pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
   //       }
   //    }
   //
   //    pBand = pDataSet->GetRasterBand(1);
   //    pBand->SetNoDataValue(m_dMissingValue);
   //    nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
   //
   //    if (nRet == CE_Failure)
   //       return RTN_ERR_GRIDCREATE;
   //
   //    GDALClose(pDataSet);
   //    delete[] pdRaster;
   //    // DEBUG CODE ===========================================

   // Find any shadow zones and then modify waves in and adjacent to them
   nRet = nDoAllShadowZones();
   if (nRet != RTN_OK)
      return nRet;

   // Calculate the D50 for each polygon. Also fill in any artefactual 'holes' in active zone and wave property patterns
   CalcD50AndFillWaveCalcHoles();

   //    // DEBUG CODE ===========================================
   //    string strOutFile = m_strOutPath;
   //    strOutFile += "sea_wave_height_CHECKPOINT_";
   //    strOutFile += std::to_string(m_ulIter);
   //    strOutFile += ".tif";
   //
   //    GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   //    GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
   //    pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   //    pDataSet->SetGeoTransform(m_dGeoTransform);
   //
   //    int nn = 0;
   //    double* pdRaster = new double[m_nXGridMax * m_nYGridMax];
   //    for (int nY = 0; nY < m_nYGridMax; nY++)
   //    {
   //       for (int nX = 0; nX < m_nXGridMax; nX++)
   //       {
   //          pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
   //       }
   //    }
   //
   //    GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   //    pBand->SetNoDataValue(m_dMissingValue);
   //    nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
   //
   //    if (nRet == CE_Failure)
   //       return RTN_ERR_GRIDCREATE;
   //
   //    GDALClose(pDataSet);
   //    delete[] pdRaster;
   //    // DEBUG CODE ===========================================
   //
   //    // DEBUG CODE ===========================================
   //    strOutFile = m_strOutPath;
   //    strOutFile += "sea_wave_angle_CHECKPOINT_";
   //    strOutFile += std::to_string(m_ulIter);
   //    strOutFile += ".tif";
   //
   //    pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   //    pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
   //    pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   //    pDataSet->SetGeoTransform(m_dGeoTransform);
   //
   //    nn = 0;
   //    pdRaster = new double[m_nXGridMax * m_nYGridMax];
   //    for (int nY = 0; nY < m_nYGridMax; nY++)
   //    {
   //       for (int nX = 0; nX < m_nXGridMax; nX++)
   //       {
   //          pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
   //       }
   //    }
   //
   //    pBand = pDataSet->GetRasterBand(1);
   //    pBand->SetNoDataValue(m_dMissingValue);
   //    nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
   //
   //    if (nRet == CE_Failure)
   //       return RTN_ERR_GRIDCREATE;
   //
   //    GDALClose(pDataSet);
   //    delete[] pdRaster;
   //    // DEBUG CODE ===========================================

   // Modify the wave breaking properties (wave height, wave dir, breaking depth, breaking distance) for coastline points within the shadow zone
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      for (int nProfile = 0; nProfile < nNumProfiles; nProfile++)
         ModifyBreakingWavePropertiesWithinShadowZoneToCoastline(nCoast, nProfile);
   }

   //    // DEBUG CODE ===========================================
   //    string strOutFile = m_strOutPath;
   //    strOutFile += "sea_wave_height_CHECKPOINT_";
   //    strOutFile += std::to_string(m_ulIter);
   //    strOutFile += ".tif";
   //
   //    GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   //    GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
   //    pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   //    pDataSet->SetGeoTransform(m_dGeoTransform);
   //
   //    int nn = 0;
   //    double* pdRaster = new double[m_nXGridMax * m_nYGridMax];
   //    for (int nY = 0; nY < m_nYGridMax; nY++)
   //    {
   //       for (int nX = 0; nX < m_nXGridMax; nX++)
   //       {
   //          pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();
   //       }
   //    }
   //
   //    GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   //    pBand->SetNoDataValue(m_dMissingValue);
   //    nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
   //
   //    if (nRet == CE_Failure)
   //       return RTN_ERR_GRIDCREATE;
   //
   //    GDALClose(pDataSet);
   //    delete[] pdRaster;
   //    // DEBUG CODE ===========================================
   //
   //    // DEBUG CODE ===========================================
   //    strOutFile = m_strOutPath;
   //    strOutFile += "sea_wave_angle_CHECKPOINT_";
   //    strOutFile += std::to_string(m_ulIter);
   //    strOutFile += ".tif";
   //
   //    pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   //    pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridMax, m_nYGridMax, 1, GDT_Float64, m_papszGDALRasterOptions);
   //    pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   //    pDataSet->SetGeoTransform(m_dGeoTransform);
   //
   //    nn = 0;
   //    pdRaster = new double[m_nXGridMax * m_nYGridMax];
   //    for (int nY = 0; nY < m_nYGridMax; nY++)
   //    {
   //       for (int nX = 0; nX < m_nXGridMax; nX++)
   //       {
   //          pdRaster[nn++] = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();
   //       }
   //    }
   //
   //    pBand = pDataSet->GetRasterBand(1);
   //    pBand->SetNoDataValue(m_dMissingValue);
   //    nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridMax, m_nYGridMax, pdRaster, m_nXGridMax, m_nYGridMax, GDT_Float64, 0, 0, NULL);
   //
   //    if (nRet == CE_Failure)
   //       return RTN_ERR_GRIDCREATE;
   //
   //    GDALClose(pDataSet);
   //    delete[] pdRaster;
   //    // DEBUG CODE ===========================================

   // Interpolate these wave properties for all remaining coastline points. Do this in along-coastline sequence
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
      int nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      // Interpolate these wave properties for all remaining coastline points. Do this in along-coastline sequence, but do not do this for the end-of-coastline profile (which is the final one)
      for (int n = 0; n < nNumProfiles - 1; n++)
         InterpolateWavePropertiesBetweenProfiles(nCoast, n, nNumProfiles);

      // And do the same for the coastline cells
      InterpolateWaveHeightToCoastPoints(nCoast);

      // Calculate wave energy at breaking for every point on the coastline
      for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
      {
         // Equation 4 from Walkden & Hall, 2005
         double dBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nCoastPoint);
         double dCoastPointWavePeriod = m_VCoast[nCoast].dGetCoastDeepWaterWavePeriod(nCoastPoint);
         
         // BUG 001 Why do we get -ve dBreakingWaveHeight here?
         if (bFPIsEqual(dBreakingWaveHeight, DBL_NODATA, TOLERANCE))
         {
            m_VCoast[nCoast].SetBreakingWaveHeight(nCoastPoint, 0);
         }
         else
         {
            double dErosiveWaveForce = pow(dBreakingWaveHeight, WALKDEN_HALL_PARAM_1) * pow(dCoastPointWavePeriod, WALKDEN_HALL_PARAM_2);

            // Calculate total wave energy at this coast point during this timestep
            double dWaveEnergy = dErosiveWaveForce * m_dTimeStep * 3600;
            m_VCoast[nCoast].SetWaveEnergyAtBreaking(nCoastPoint, dWaveEnergy);
         }
      }
   }
   return RTN_OK;
}

//===============================================================================================================================
//! Calculates the angle between the wave direction and a normal to the coastline tangent. If wave direction has a component which is down-coast (i.e. in the direction with increasing coast point numbers), then the angle returned is +ve. If wave direction has a component which is up-coast (i.e. in the direction with decreasing coast point numbers), then the angle returned is -ve. If waves are in an off-shore direction, DBL_NODATA is returned
//===============================================================================================================================
double CSimulation::dCalcWaveAngleToCoastNormal(double const dCoastAngle, double const dWaveAngle, int const nSeaHand)
{
   double dWaveToNormalAngle = 0;

   if (nSeaHand == LEFT_HANDED)
      // Left-handed coast
      dWaveToNormalAngle = fmod((dWaveAngle - dCoastAngle + 360), 360) - 90;
   else
      // Right-handed coast
      dWaveToNormalAngle = fmod((dWaveAngle - dCoastAngle + 360), 360) - 270;

   if ((dWaveToNormalAngle >= 90) || (dWaveToNormalAngle <= -90))
      dWaveToNormalAngle = DBL_NODATA;

   return dWaveToNormalAngle;
}

//===============================================================================================================================
//! Calculates wave properties along a coastline-normal profile using either the COVE linear wave theory approach or the external CShore model
//===============================================================================================================================
int CSimulation::nCalcWavePropertiesOnProfile(int const nCoast, int const nCoastSize, int const nProfile, vector<double>* pVdX, vector<double>* pVdY, vector<double>* pVdHeightX, vector<double>* pVdHeightY, vector<bool>* pVbBreaking)
{
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // Calculate some wave properties based on the wave period following Airy wave theory
   double dDeepWaterWavePeriod = pProfile->dGetProfileDeepWaterWavePeriod();

   m_dC_0 = (m_dG * dDeepWaterWavePeriod) / (2 * PI); // Deep water (offshore) wave celerity (m/s)
   m_dL_0 = m_dC_0 * dDeepWaterWavePeriod;            // Deep water (offshore) wave length (m)

   // Only do this for profiles without problems. Still do start- and end-of-coast profiles however
   if (! pProfile->bOKIncStartAndEndOfCoast())
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": Coast " << nCoast << ", profile " << nProfile << " has been marked invalid, will not calc wave properties on this profile" << endl;

      return RTN_OK;
   }

   int nSeaHand = m_VCoast[nCoast].nGetSeaHandedness();
   int nCoastPoint = pProfile->nGetNumCoastPoint();

   // Get the flux orientation (the orientation of a line which is tangential to the coast) at adjacent coastline points. Note special treatment for the coastline end points
   double dFluxOrientationThis = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint);
   double dFluxOrientationPrev = 0;
   double dFluxOrientationNext = 0;

   if (nCoastPoint == 0)
   {
      dFluxOrientationPrev = dFluxOrientationThis;
      dFluxOrientationNext = m_VCoast[nCoast].dGetFluxOrientation(1);
   }
   else if (nCoastPoint == nCoastSize - 1)
   {
      dFluxOrientationPrev = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint - 2);
      dFluxOrientationNext = dFluxOrientationThis;
   }
   else
   {
      dFluxOrientationPrev = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint - 1);
      dFluxOrientationNext = m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint + 1);
   }

   // Get the deep water wave orientation for this profile
   double dDeepWaterWaveAngle = pProfile->dGetProfileDeepWaterWaveAngle();

   // Calculate the angle between the deep water wave direction and a normal to the coast tangent
   double dWaveToNormalAngle = dCalcWaveAngleToCoastNormal(dFluxOrientationThis, dDeepWaterWaveAngle, nSeaHand);

   // Are the waves off-shore?
   if (bFPIsEqual(dWaveToNormalAngle, DBL_NODATA, TOLERANCE))
   {
      // They are so, do nothing (each cell under the profile has already been initialised with deep water wave height and wave direction)
      //       LogStream << m_ulIter << ": profile " << nProfile << " has sea to " << (m_VCoast[nCoast].nGetSeaHandedness() == RIGHT_HANDED ? "right" : "left") << " dWaveToNormalAngle = " << dWaveToNormalAngle << " which is off-shore" << endl;

      return RTN_OK;
   }

   //    LogStream << m_ulIter << ": profile = " << nProfile << " has sea to " << (m_VCoast[nCoast].nGetSeaHandedness() == RIGHT_HANDED ? "right" : "left") << " dWaveToNormalAngle = " << dWaveToNormalAngle << " which is " << (dWaveToNormalAngle < 0 ? "DOWN" : "UP") << "-coast" << endl;

   // Calculate the angle between the deep water wave direction and a normal to the coast tangent for the previous coast point
   double dWaveToNormalAnglePrev;
   if (nCoastPoint > 0)
   {
      // Get the deep water wave orientation for the up-coast point
      double dPrevDeepWaterWaveAngle = m_VCoast[nCoast].dGetCoastDeepWaterWaveAngle(nCoastPoint - 1);

      dWaveToNormalAnglePrev = dCalcWaveAngleToCoastNormal(dFluxOrientationPrev, dPrevDeepWaterWaveAngle, nSeaHand);
   }
   else
   {
      dWaveToNormalAnglePrev = dWaveToNormalAngle;
   }

   //    if (dWaveToNormalAnglePrev == DBL_NODATA)
   //       LogStream << "\tPrevious profile, dWaveToNormalAnglePrev = " << dWaveToNormalAnglePrev << " which is off-shore" << endl;
   //    else
   //       LogStream << "\tPrevious profile, dWaveToNormalAnglePrev = " << dWaveToNormalAnglePrev << " which is " << (dWaveToNormalAnglePrev < 0 ? "DOWN" : "UP") << "-coast" << endl;

   // Calculate the angle between the deep water wave direction and a normal to the coast tangent for the next coast point
   double dWaveToNormalAngleNext;
   if (nCoastPoint < nCoastSize - 1)
   {
      // Get the deep water wave orientation for the down-coast point
      double dNextDeepWaterWaveAngle = m_VCoast[nCoast].dGetCoastDeepWaterWaveAngle(nCoastPoint + 1);

      dWaveToNormalAngleNext = dCalcWaveAngleToCoastNormal(dFluxOrientationNext, dNextDeepWaterWaveAngle, nSeaHand);
   }
   else
   {
      dWaveToNormalAngleNext = dWaveToNormalAngle;
   }

   //    if (dWaveToNormalAngleNext == DBL_NODATA)
   //       LogStream << "\tNext profile, dWaveToNormalAngleNext = " << dWaveToNormalAngleNext << " which is off-shore" << endl;
   //    else
   //       LogStream << "\tNext profile, dWaveToNormalAngleNext = " << dWaveToNormalAngleNext << " which is " << (dWaveToNormalAngleNext < 0 ? "DOWN" : "UP") << "-coast" << endl;

   // Following Ashton and Murray (2006), if we have high-angle waves then use the flux orientation of the previous (up-coast) profile, if transitioning from diffusive to antidiffusive use flux maximizing angle (45 degrees)
   if ((dWaveToNormalAngle > 0) && (! bFPIsEqual(dWaveToNormalAnglePrev, DBL_NODATA, TOLERANCE)) && (dWaveToNormalAnglePrev > 0))
   {
      if (dWaveToNormalAngle > 45)
      {
         if (dWaveToNormalAnglePrev < 45)
         {
            dWaveToNormalAngle = 45;
            //            LogStream << "\tA1" << endl;
         }
         else
         {
            dWaveToNormalAngle = dWaveToNormalAnglePrev;
            //            LogStream << "\tA2" << endl;
         }
      }
   }
   else if ((dWaveToNormalAngle < 0) && (! bFPIsEqual(dWaveToNormalAngleNext, DBL_NODATA, TOLERANCE)) && (dWaveToNormalAngleNext < 0))
   {
      if (dWaveToNormalAngle < -45)
      {
         if (dWaveToNormalAngleNext > -45)
         {
            dWaveToNormalAngle = -45;
            //            LogStream << "\tB1" << endl;
         }
         else
         {
            dWaveToNormalAngle = dWaveToNormalAngleNext;
            //            LogStream << "\tB2" << endl;
         }
      }
   }
   else if ((dWaveToNormalAngle > 45) && (! bFPIsEqual(dWaveToNormalAnglePrev, DBL_NODATA, TOLERANCE)) && (dWaveToNormalAnglePrev > 0))
   {
      // The wave direction here has an up-coast (decreasing indices) component: so for high-angle waves use the orientation from the up-coast (previous) profile
      //      LogStream << "\tCCC" << endl;

      dWaveToNormalAngle = dFluxOrientationPrev;
   }
   else if ((dWaveToNormalAngle < -45) && (! bFPIsEqual(dWaveToNormalAngleNext, DBL_NODATA, TOLERANCE)) && (dWaveToNormalAngleNext < 0))
   {
      // The wave direction here has a down-coast (increasing indices) component: so for high-angle waves use the orientation from the down-coast (next) profile
      //      LogStream << "\tDDD" << endl;

      dWaveToNormalAngle = dFluxOrientationNext;
   }

   // Initialize the wave properties at breaking for this profile
   bool bBreaking = false;

   int nProfileSize = pProfile->nGetNumCellsInProfile();
   int nProfileBreakingDist = 0;

   double dProfileBreakingWaveHeight = DBL_NODATA;
   double dProfileBreakingWaveAngle = 0;
   double dProfileBreakingDepth = 0;
   double dProfileWaveHeight = DBL_NODATA;
   double dProfileDeepWaterWaveHeight = pProfile->dGetProfileDeepWaterWaveHeight();
   double dProfileWaveAngle = DBL_NODATA;
   double dProfileDeepWaterWaveAngle = pProfile->dGetProfileDeepWaterWaveAngle();

   vector<bool> VbWaveIsBreaking(nProfileSize, 0);

   vector<double> VdWaveHeight(nProfileSize, 0);
   vector<double> VdWaveSetupSurge(nProfileSize, 0);          // TODO 007 What is this for?
   // vector<double> VdStormSurge(nProfileSize, 0);            // TODO 007 What is this for?
   vector<double> VdWaveSetupRunUp(nProfileSize, 0);
   vector<double> VdWaveDirection(nProfileSize, 0);

   if (m_nWavePropagationModel == WAVE_MODEL_CSHORE)
   {
      // We are using CShore to propagate the waves
      double dCShoreTimeStep = 3600;   // In seconds, not important because we are not using CShore to erode the profile, just to get the hydrodynamics
      double dSurgeLevel = CSHORE_SURGE_LEVEL;

      // Set up vectors for the coastline-normal profile elevations. The length of this vector line is given by the number of cells 'under' the profile. Thus each point on the vector relates to a single cell in the grid. This assumes that all points on the profile vector are equally spaced (not quite true, depends on the orientation of the line segments which comprise the profile)
      vector<double> VdProfileZ;                   // Initial (pre-erosion) elevation of both consolidated and unconsolidated sediment for cells 'under' the profile, in CShore units
      vector<double> VdProfileDistXY;              // Along-profile distance measured from the seaward limit, in CShore units
      vector<double> VdProfileFrictionFactor;      // Along-profile friction factor from seaward limit

      // The elevation of each of these profile points is the elevation of the centroid of the cell that is 'under' the point. However we cannot always be confident that this is the 'true' elevation of the point on the vector since (unless the profile runs planview N-S or W-E) the vector does not always run exactly through the centroid of the cell
      int nRet = nGetThisProfileElevationVectorsForCShore(nCoast, nProfile, nProfileSize, &VdProfileDistXY, &VdProfileZ, &VdProfileFrictionFactor);
      if (nRet != RTN_OK)
      {
         // Could not create the profile elevation vectors
         LogStream << m_ulIter << ": could not create CShore profile elevation vectors for profile " << nProfile << endl;

         return nRet;
      }

      if (VdProfileDistXY.empty())
      {
         // The profile elevation vector was created, but was not populated
         LogStream << m_ulIter << ": could not populate CShore profile elevation vector for profile " << nProfile << endl;

         return RTN_ERR_CSHORE_EMPTY_PROFILE;
      }

      // Constrain the wave to normal angle to be between -80 and 80 degrees, this is a requirement of CShore
      dWaveToNormalAngle = tMax(dWaveToNormalAngle, -80.0);
      dWaveToNormalAngle = tMin(dWaveToNormalAngle, 80.0);

      int nProfileDistXYSize = static_cast<int>(VdProfileDistXY.size());
      vector<double>
          VdFreeSurfaceStd(nProfileDistXYSize, 0),        // This is converted to Hrms by Hrms = sqr(8)*FreeSurfaceStd
          VdSinWaveAngleRadians(nProfileDistXYSize, 0),   // This is converted to deg by asin(VdSinWaveAngleRadians)*(180/pi)
          VdFractionBreakingWaves(nProfileDistXYSize, 0); // Is 0 if no wave breaking, and 1 if all waves breaking

      // Now define the other values that CShore requires
      int
          nILine = 1,  // This is the number of cross-shore lines i.e. the number of CoastalME profiles. Only one at a time, at present
          nIProfl = 0, // 0 for fixed bottom profile, 1 for profile evolution computation
          nIPerm = 0,  // 0 for impermeable bottom, 1 for permeable bottom of stone structure
          nIOver = 0,  // 0 for no wave overtopping and overflow on crest, 1 for wave overtopping and overflow
          nIWCInt = 0, // 0 for no wave/current interaction, 1 for wave/current interaction in frequency dispersion, momentum and wave action equations
          nIRoll = 0,  // 0 for no roller effects, 1 for roller effects in governing equations
          nIWind = 0,  // 0 for no wind effects, 1 for wind shear stresses on momentum equations
          nITide = 0,  // 0 for no tidal effect on currents, 1 for longshore and cross-shore tidal currents
          nILab = 0,   // 0 for field data set, 1 for laboratory data set
          nNWave = 1,  // Number of waves at x = 0 starting from time = 0
          nNSurge = 1; // Number of water levels at x = 0 from time = 0    // TODO 007 What is this for?
      double
          dDX = VdProfileDistXY.back() / (static_cast<double>(VdProfileDistXY.size() - 1)),
          dWaveInitTime = 0,  // CShore wave start time
          dSurgeInitTime = 0; // CShore surge start time                   // TODO 007 What is this for?

#if defined CSHORE_FILE_INOUT || CSHORE_BOTH
      // Move to the CShore folder
      nRet = chdir(CSHORE_DIR.c_str());
      if (nRet != RTN_OK)
         return nRet;
#endif

#if defined CSHORE_FILE_INOUT
      // We are communicating with CShore using ASCII files, so create an input file for this profile which will be read by CShore
      nRet = nCreateCShoreInfile(nCoast, nProfile, nILine, nIProfl, nIPerm, nIOver, nIWCInt, nIRoll, nIWind, nITide, nILab, nNWave, nNSurge, dDX, dCShoreTimeStep, dWaveInitTime, dDeepWaterWavePeriod, dProfileDeepWaterWaveHeight, dWaveToNormalAngle, dSurgeInitTime, dSurgeLevel, &VdProfileDistXY, &VdProfileZ, &VdProfileFrictionFactor);
      if (nRet != RTN_OK)
         return nRet;

      // Set the error flag: this will be changed to 0 within CShore if CShore returns correctly
      nRet = -1;

      // Run CShore for this profile
      CShore(&nRet);
//      CShore(&nRet, &m_ulIter, &nProfile, &nProfile);

      // Check return code for error
      if (nRet != 0)
      {
         string strErr;

         switch (nRet)
         {
            case -1:
               strErr = to_string(m_ulIter) + ": CShore ERROR: negative depth at the first node ";
               break;

            case 2:
               strErr = to_string(m_ulIter) + ": CShore WARNING 2: negative value at end of landward marching computation ";
               break;

            case 3:
               strErr = to_string(m_ulIter) + ": CShore WARNING 3: large energy gradients at the first node: small waves with short period at sea boundary ";
               break;

            case 4:
               strErr = to_string(m_ulIter) + ": CShore WARNING 4: zero energy at the first node ";
               break;

            case 5:
               strErr = to_string(m_ulIter) + ": CShore WARNING 5: at end of landward marching computation, insufficient water depth ";
               break;

            case 7:
               strErr = to_string(m_ulIter) + ": CShore WARNING 7: did not reach convergence ";
               break;
         }

         strErr += "(coast " + to_string(nCoast) + " profile " + to_string(nProfile) + ")\n";
         
         if (nRet < -1)
         {
            // This is serious, so give up for this profile
            cerr << strErr;
            return RTN_ERR_CSHORE_ERROR;
         }

         // Not too serious, so carry on
      }

      // Fetch the CShore results by reading files written by CShore
      string strOSETUP = "OSETUP";
      string strOYVELO = "OYVELO";
      string strOPARAM = "OPARAM";

      nRet = nReadCShoreOutput(nProfile, &strOSETUP, 4, 4, &VdProfileDistXY, &VdFreeSurfaceStd);
      if (nRet != RTN_OK)
         return nRet;

      nRet = nReadCShoreOutput(nProfile, &strOYVELO, 4, 2, &VdProfileDistXY, &VdSinWaveAngleRadians);
      if (nRet != RTN_OK)
         return nRet;

      nRet = nReadCShoreOutput(nProfile, &strOPARAM, 4, 3, &VdProfileDistXY, &VdFractionBreakingWaves);
      if (nRet != RTN_OK)
         return nRet;

      // TODO 008 Read surge outputs
      // VdTSurg = {dSurgeInitTime, dCShoreTimeStep},                           // Ditto
      // VdSWLin = {dSurgeLevel, dSurgeLevel},                                  // Ditto

      // Clean up the CShore outputs
#ifdef _WIN32
      nRet = system("./clean.bat");
#else
      if (SAVE_CSHORE_OUTPUT)
      {
         string strCommand = "./save_CShore_output.sh ";
         strCommand += to_string(m_ulIter);
         strCommand += " ";
         strCommand += to_string(nCoast);
         strCommand += " ";
         strCommand += to_string(nProfile);

         nRet = system(strCommand.c_str());
         if (nRet != RTN_OK)
            return nRet;
      }

      nRet = system("./clean.sh");
#endif
      if (nRet != RTN_OK)
         return nRet;

      // And return to the CoastalME folder
      nRet = chdir(m_strCMEDir.c_str());
      if (nRet != RTN_OK)
         return nRet;
#endif

#if defined CSHORE_ARG_INOUT || CSHORE_BOTH
      // We are communicating with CShore by passing arguments
      int nOutSize = 0; // CShore will return the size of the output vectors

      // Set the error flag: this will be changed within CShore if there is a problem
      nRet = 0;

      vector<double>
          VdInitTime = {dWaveInitTime, dCShoreTimeStep},                         // Size is nNwave+1, value 1 is for the start of the CShore run, value 2 for end of CShore run
          VdTPIn = {dDeepWaterWavePeriod, dDeepWaterWavePeriod},                 // Ditto
          VdHrmsIn = {dProfileDeepWaterWaveHeight, dProfileDeepWaterWaveHeight}, // Ditto
          VdWangIn = {dWaveToNormalAngle, dWaveToNormalAngle},                   // Ditto
          VdTSurg = {dSurgeInitTime, dCShoreTimeStep},                           // Ditto
          VdSWLin = {dSurgeLevel, dSurgeLevel},                                  // Ditto
          VdFPInp = VdProfileFrictionFactor,                                     // Set the value for wave friction at every point of the normal profile
          VdXYDistFromCShoreOut(CSHOREARRAYOUTSIZE, 0),                          // Output from CShore
          VdFreeSurfaceStdOut(CSHOREARRAYOUTSIZE, 0),                            // Ditto
          VdWaveSetupSurgeOut(CSHOREARRAYOUTSIZE, 0),                            // TODO 007 Info needed
         //  VdStormSurgeOut(CSHOREARRAYOUTSIZE, 0),                             // TODO 007 Info needed
          VdWaveSetupRunUpOut(CSHOREARRAYOUTSIZE, 0),                            // TODO 007 Info needed
          VdSinWaveAngleRadiansOut(CSHOREARRAYOUTSIZE, 0),                       // Ditto
          VdFractionBreakingWavesOut(CSHOREARRAYOUTSIZE, 0);                     // Ditto

      // Call CShore using the argument-passing wrapper
      // long lIter = static_cast<long>(m_ulIter);    // Bodge to get round compiler 'invalid conversion' error
      // 
      // CShoreWrapper(&lIter, &nCoast, &nProfile, &nILine, &nIProfl, &nIPerm, &nIOver, &nIWCInt, &nIRoll, &nIWind, &nITide, &nILab, &nNWave, &nNSurge, &dDX, &m_dBreakingWaveHeightDepthRatio, &VdInitTime[0], &VdTPIn[0], &VdHrmsIn[0], &VdWangIn[0], &VdTSurg[0], &VdSWLin[0], &nProfileDistXYSize, &VdProfileDistXY[0], &VdProfileZ[0], &VdFPInp[0], &nRet, &nOutSize, &VdXYDistFromCShoreOut[0], &VdFreeSurfaceStdOut[0], &VdWaveSetupOut[0], &VdSinWaveAngleRadiansOut[0], &VdFractionBreakingWavesOut[0]);
          
      // // TEST
      // for (int nn = 0; nn < nProfileDistXYSize; nn++)
      // {
      //    LogStream << nn << "\t" << VdProfileDistXY[nn] << "\t" << VdProfileDistXY[nn] << "\t" << VdProfileZ[nn] << endl;         
      // }
          
      CShoreWrapper(&nILine,                          /* In_ILINE */
                    &nIProfl,                         /* In_IPROFL */
                    &nIPerm,                          /* In_IPERM */
                    &nIOver,                          /* In_IOVER */
                    &nIWCInt,                         /* In_IWCINT */
                    &nIRoll,                          /* In_IROLL */
                    &nIWind,                          /* In_IWIND */
                    &nITide,                          /* In_ITIDE */
                    &nILab,                           /* In_ILAB */
                    &nNWave,                          /* In_NWAVE */
                    &nNSurge,                         /* In_NSURG */                   // TODO 007 Info needed
                    &dDX,                             /* In_DX */
                    &m_dBreakingWaveHeightDepthRatio, /* In_GAMMA */
                    &VdInitTime[0],                   /* In_TWAVE */
                    &VdTPIn[0],                       /* In_TPIN */
                    &VdHrmsIn[0],                     /* In_HRMSIN */
                    &VdWangIn[0],                     /* In_WANGIN */
                    &VdTSurg[0],                      /* In_TSURG */
                    &VdSWLin[0],                      /* In_SWLIN */
                    &nProfileDistXYSize,              /* In_NBINP */
                    &VdProfileDistXY[0],              /* In_XBINP */
                    &VdProfileZ[0],                   /* In_ZBINP */
                    &VdFPInp[0],                      /* In_FBINP */
                    &nRet,                            /* Out_IError */
                    &nOutSize,                        /* Out_nOutSize */
                    &VdXYDistFromCShoreOut[0],        /* Out_XYDist */
                    &VdFreeSurfaceStdOut[0],          /* Out_FreeSurfaceStd */
                    &VdWaveSetupSurgeOut[0],          /* Out_WaveSetupSurge */         // TODO 007 Info needed
                    &VdSinWaveAngleRadiansOut[0],     /* Out_SinWaveAngleRadians */
                    &VdFractionBreakingWavesOut[0]);  /* Out_FractionBreakingWaves */
      
      // OK, now check for warnings and errors
      if (nOutSize < 2)
      {
         // CShore sometimes returns only one row of results, which contains data only for the seaward point of the profile. This happens when all other (more coastward) points give an invalid result during CShore's calculations. This is a problem. We don't want to abandon the simulation just because of this, so instead we just put some dummy data into the second row, and carry on with these two rows. The profile will get ignored later, since it is too small to be useful
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": " << WARN << "for coast " << nCoast << " profile " << nProfile << ", only " << nOutSize << " CShore output row" << endl;

         // Set dummy data in the second row
         VdXYDistFromCShoreOut[1] = 1e-5;    // Dummy data, must not be the same as VdXYDistFromCShoreOut[0] tho', or get crash in linear interpolation routine
         VdFreeSurfaceStdOut[1] = VdFreeSurfaceStdOut[0];
         VdSinWaveAngleRadiansOut[1] = VdSinWaveAngleRadiansOut[0];
         VdFractionBreakingWavesOut[1] = VdFractionBreakingWavesOut[0];
         VdWaveSetupSurgeOut[1] = VdWaveSetupSurgeOut[0];                              // TODO 007 Info needed
         // VdStormSurgeOut[1] = VdStormSurgeOut[0];                                   // TODO 007 Info needed
         // VdWaveSetupRunUpOut[1] = VdWaveSetupRunUpOut[0];

         // And increase the expected number of rows
         nOutSize = 2;
      }

      if (nRet != RTN_OK)
      {
         string strErr;

         switch (nRet)
         {
            case -1:
               strErr = to_string(m_ulIter) + ": CShore ERROR: negative depth at the first node ";
               break;

            case 2:
               strErr = to_string(m_ulIter) + ": CShore WARNING 2: negative value at end of landward marching computation ";
               break;

            case 3:
               strErr = to_string(m_ulIter) + ": CShore WARNING 3: large energy gradients at the first node: small waves with short period at sea boundary ";
               break;

            case 4:
               strErr = to_string(m_ulIter) + ": CShore WARNING 4: zero energy at the first node ";
               break;

            case 5:
               strErr = to_string(m_ulIter) + ": CShore WARNING 5: at end of landward marching computation, insufficient water depth ";
               break;

            case 7:
               strErr = to_string(m_ulIter) + ": CShore WARNING 7: did not reach convergence ";
               break;
         }

         strErr += "(coast " + to_string(nCoast) + " profile " + to_string(nProfile) + " profile length " + to_string(nOutSize) + ")\n";
         
         if (nRet < -1)
         {
            // This is serious, so give up for this profile
            cerr << strErr;
            return RTN_ERR_CSHORE_ERROR;
         }

         // Not too serious, so carry on
      }
      
      // TEST nOutsize < nOutSize occasionally
      //assert(static_cast<int>(VdFreeSurfaceStd.size()) == nOutSize);
      //LogStream << "VdFreeSurfaceStd.size() = " << nOutSize << " nOutSize = " << nOutSize << endl;
      
      // TEST
//       for (int nn = 0; nn < static_cast<int>(VdFreeSurfaceStd.size()); nn++)
//       {
//          assert(isfinite(VdProfileDistXY[nn]));
//          
//          assert(isfinite(VdXYDistFromCShoreOut[nn]));
//          assert(isfinite(VdFreeSurfaceStdOut[nn]));      
//          assert(isfinite(VdWaveSetupSurgeOut[nn]));                  // TODO 007 What is this for?
//          assert(isfinite(VdSinWaveAngleRadiansOut[nn]));
//          assert(isfinite(VdFractionBreakingWavesOut[nn]));
//          
//          assert(isfinite(VdFreeSurfaceStd[nn]));
//          assert(isfinite(VdWaveSetupSurge[nn]));                     // TODO 007 What is this for?
//          assert(isfinite(VdSinWaveAngleRadians[nn]));
//          assert(isfinite(VdFractionBreakingWaves[nn]));
//       }
      
      // Now interpolate the output
      InterpolateCShoreOutput(&VdProfileDistXY, nOutSize, &VdXYDistFromCShoreOut, &VdFreeSurfaceStdOut, &VdWaveSetupSurgeOut, &VdSinWaveAngleRadiansOut, &VdFractionBreakingWavesOut, &VdFreeSurfaceStd, &VdWaveSetupSurge, &VdSinWaveAngleRadians, &VdFractionBreakingWaves);
      
      // TEST
//       for (int nn = 0; nn < static_cast<int>(VdFreeSurfaceStd.size()); nn++)
//       {
//          assert(isfinite(VdProfileDistXY[nn]));
//          
//          assert(isfinite(VdXYDistFromCShoreOut[nn]));
//          assert(isfinite(VdFreeSurfaceStdOut[nn]));      
//          assert(isfinite(VdWaveSetupSurgeOut[nn]));                  // TODO 007 What is this for?
//          assert(isfinite(VdSinWaveAngleRadiansOut[nn]));
//          assert(isfinite(VdFractionBreakingWavesOut[nn]));
//          
//          assert(isfinite(VdFreeSurfaceStd[nn]));
//          assert(isfinite(VdWaveSetupSurge[nn]));                     // TODO 007 What is this for?
//          assert(isfinite(VdSinWaveAngleRadians[nn]));
//          assert(isfinite(VdFractionBreakingWaves[nn]));
//       }
#endif

#if defined CSHORE_BOTH
   #if ! defined _WIN32
      if (SAVE_CSHORE_OUTPUT)
      {
         string strCommand = "./save_CShore_output.sh ";
         strCommand += to_string(m_ulIter);
         strCommand += " ";
         strCommand += to_string(nCoast);
         strCommand += " ";
         strCommand += to_string(nProfile);

         nRet = system(strCommand.c_str());
         if (nRet != RTN_OK)
            return nRet;
      }
   #endif

      // Return to the CoastalME folder
      nRet = chdir(m_strCMEDir.c_str());
      if (nRet != RTN_OK)
         return nRet;
#endif

      // OK we have the CShore output, so now we must convert this to wave height and wave direction and update wave profile attributes
      for (int nProfilePoint = (nProfileSize - 1); nProfilePoint >= 0; nProfilePoint--)
      {
         int
             nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX(),
             nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();
             
         // Safety check: deal with NaN values
         if (! isfinite(VdFreeSurfaceStd[nProfilePoint]))
            VdFreeSurfaceStd[nProfilePoint] = 0;
             
         VdWaveHeight[nProfilePoint] = sqrt(8) * VdFreeSurfaceStd[nProfilePoint];
         
         // Another safety check: deal with NaN values
         if (! isfinite(VdSinWaveAngleRadians[nProfilePoint]))
         {
            VdSinWaveAngleRadians[nProfilePoint] = 0;
            VdWaveHeight[nProfilePoint] = 0;
         }
         
         // More safety checks: constrain to the interval -1 to +1 to keep asin() happy
         if (VdSinWaveAngleRadians[nProfilePoint] < -1)
            VdSinWaveAngleRadians[nProfilePoint] = -1;
         if (VdSinWaveAngleRadians[nProfilePoint] > 1)
            VdSinWaveAngleRadians[nProfilePoint] = 1;

         double dAlpha = asin(VdSinWaveAngleRadians[nProfilePoint]) * (180 / PI);
         if (nSeaHand == LEFT_HANDED)
            VdWaveDirection[nProfilePoint] = dKeepWithin360(dAlpha + 90 + dFluxOrientationThis);
         else
            VdWaveDirection[nProfilePoint] = dKeepWithin360(dAlpha + 270 + dFluxOrientationThis);
         
         // Yet another safety check: deal with NaN values
         if (! isfinite(VdFractionBreakingWaves[nProfilePoint]))
         {
            VdFractionBreakingWaves[nProfilePoint] = 0;
            VdWaveHeight[nProfilePoint] = 0;
         }

         // if ((VdFractionBreakingWaves[nProfilePoint] >= 0.10) && (! bBreaking)) // Sometimes is possible that waves break again
         if ((VdFractionBreakingWaves[nProfilePoint] >= 0.10) && (m_dDepthOfClosure >= m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth()) && (! bBreaking))
         {
            bBreaking = true;
            // assert(VdWaveHeight[nProfilePoint] >= 0);
            dProfileBreakingWaveHeight = VdWaveHeight[nProfilePoint];
            dProfileBreakingWaveAngle = VdWaveDirection[nProfilePoint];
            dProfileBreakingDepth = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(); // Water depth for the cell 'under' this point in the profile
            nProfileBreakingDist = nProfilePoint + 1;                             // At the nearest point nProfilePoint = 0, so, plus one

            //             LogStream << m_ulIter << ": CShore breaking at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} nProfile = " << nProfile << ", nProfilePoint = " << nProfilePoint << ", dBreakingWaveHeight = " << dBreakingWaveHeight << ", dBreakingWaveAngle = " << dBreakingWaveAngle << ", dProfileBreakingDepth = " << dProfileBreakingDepth << ", nProfileBreakingDist = " << nProfileBreakingDist << endl;
         }

         VbWaveIsBreaking[nProfilePoint] = bBreaking;
      }
      
      if (dProfileBreakingWaveHeight >= dProfileDeepWaterWaveHeight)
      {
         dProfileBreakingWaveHeight = DBL_NODATA; // checking poorly conditions profiles problems for cshore
      }
   }

   else if (m_nWavePropagationModel == WAVE_MODEL_COVE)
   {
      // We are using COVE's linear wave theory to propagate the waves
      double dDepthLookupMax = m_dWaveDepthRatioForWaveCalcs * dProfileDeepWaterWaveHeight;

      // Go landwards along the profile, calculating wave height and wave angle for every inundated point on the profile (don't do point zero, this is on the coastline) until the waves start to break  after breaking wave height is assumed to decrease linearly to zero at the shoreline and wave angle is equalt to wave angle at breaking
      for (int nProfilePoint = (nProfileSize - 1); nProfilePoint >= 0; nProfilePoint--)
      {
         int
             nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX(),
             nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();

         // Safety check
         if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
            continue;

         double dSeaDepth = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(); // Water depth for the cell 'under' this point in the profile
         if (dSeaDepth > dDepthLookupMax)
         {
            // Sea depth is too large relative to wave height to feel the bottom, so do nothing since each cell under the profile has already been initialised with deep water wave height and wave direction
            dProfileWaveHeight = dProfileDeepWaterWaveHeight;
            dProfileWaveAngle = dProfileDeepWaterWaveAngle;
         }
         else
         {
            if (! bBreaking)
            {
               // Start calculating wave properties using linear wave theory
               double dL = m_dL_0 * sqrt(tanh((2 * PI * dSeaDepth) / m_dL_0));                          // Wavelength (m) in intermediate-shallow waters
               double dC = m_dC_0 * tanh((2 * PI * dSeaDepth) / dL);                                    // Wave speed (m/s) set by dSeaDepth, dL and m_dC_0
               double dk = 2 * PI / dL;                                                                 // Wave number (1/m)
               double dn = ((2 * dSeaDepth * dk) / (sinh(2 * dSeaDepth * dk)) + 1) / 2;                 // Shoaling factor
               double dKs = sqrt(m_dC_0 / (dn * dC * 2));                                               // Shoaling coefficient
               double dAlpha = (180 / PI) * asin((dC / m_dC_0) * sin((PI / 180) * dWaveToNormalAngle)); // Calculate angle between wave direction and the normal to the coast tangent
               double dKr = sqrt(cos((PI / 180) * dWaveToNormalAngle) / cos((PI / 180) * dAlpha));      // Refraction coefficient
               dProfileWaveHeight = dProfileDeepWaterWaveHeight * dKs * dKr;                            // Calculate wave height, based on the previous (more seaward) wave height
               if (nSeaHand == LEFT_HANDED)
                  dProfileWaveAngle = dKeepWithin360(dAlpha + 90 + dFluxOrientationThis);
               else
                  dProfileWaveAngle = dKeepWithin360(dAlpha + 270 + dFluxOrientationThis);

               // Test to see if the wave breaks at this depth
               if (dProfileWaveHeight < (dSeaDepth * m_dBreakingWaveHeightDepthRatio))
               {
                  dProfileBreakingWaveHeight = dProfileWaveHeight;
                  dProfileBreakingWaveAngle = dProfileWaveAngle;
               }
            }
            else
            {
               // It does
               bBreaking = true;

               // Wave has already broken
               dProfileWaveAngle = dProfileBreakingWaveAngle; // Wave orientation remains equal to wave orientation at breaking

               dProfileWaveHeight = dProfileBreakingWaveHeight * (nProfilePoint / nProfileBreakingDist); // Wave height decreases linearly to zero at shoreline
               // dProfileWaveHeight = dSeaDepth * m_dBreakingWaveHeightDepthRatio; // Wave height is limited by depth
               dProfileBreakingDepth = dSeaDepth;
               nProfileBreakingDist = nProfilePoint;
            }
         }

         // Save current wave attributes
         VdWaveDirection[nProfilePoint] = dProfileWaveAngle;
         VdWaveHeight[nProfilePoint] = dProfileWaveHeight;
         VbWaveIsBreaking[nProfilePoint] = bBreaking;
      }
   }

   // Go landwards along the profile, fetching the calculated wave height and wave angle for every inundated point on this profile
   for (int nProfilePoint = (nProfileSize - 1); nProfilePoint >= 0; nProfilePoint--)
   {
      int
          nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX(),
          nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();

      // Safety check
      if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         continue;

      // Get the wave attributes calculated for this profile: wave height, wave angle, and whether is in the active zone
      double
          dWaveHeight = VdWaveHeight[nProfilePoint],
          dWaveAngle = VdWaveDirection[nProfilePoint];

      bBreaking = VbWaveIsBreaking[nProfilePoint];

      // And store the wave properties for this point in the all-profiles vectors
      pVdX->push_back(nX);
      pVdY->push_back(nY);
      pVdHeightX->push_back(dWaveHeight * sin(dWaveAngle * PI / 180));
      pVdHeightY->push_back(dWaveHeight * cos(dWaveAngle * PI / 180));
      pVbBreaking->push_back(bBreaking);
   }

   // Obtain the profile nodes near the coast
   int
       nX = pProfile->pPtiGetCellInProfile(nProfileSize - 2)->nGetX(),
       nY = pProfile->pPtiGetCellInProfile(nProfileSize - 2)->nGetY(),
       nX1 = pProfile->pPtiGetCellInProfile(nProfileSize - 1)->nGetX(),
       nY1 = pProfile->pPtiGetCellInProfile(nProfileSize - 1)->nGetY();

   double
       dXDist,
       dYDist = 0;

   // Calculate the horizontal distance between the profile points
   dXDist = tAbs(dGridCentroidXToExtCRSX(nX1) - dGridCentroidXToExtCRSX(nX));
   dYDist = tAbs(dGridCentroidYToExtCRSY(nY1) - dGridCentroidYToExtCRSY(nY));
   
   // Safety check
   if (bFPIsEqual(dXDist, 0.0, TOLERANCE))
      dXDist = 1e-3;
   
   // Safety check
   if (bFPIsEqual(dYDist, 0.0, TOLERANCE))
      dYDist = 1e-3;

   double dDiffProfileDistXY = hypot(dXDist, dYDist);

   // Compute the beach slope
   double dtanBeta = tan(tAbs(m_pRasterGrid->m_Cell[nX1][nY1].dGetSeaDepth() - m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth()) / dDiffProfileDistXY);

   // Compute the wave run-up using NIELSEN & HANSLOW (1991) & DHI (2004)
   int nValidPointsWaveHeight = 0;
   int nValidPointsWaveSetup = 0;

   for (int nPoint = 0; nPoint < static_cast<int>(VdWaveHeight.size()); nPoint++)
   {
      if (VdWaveHeight[nPoint] > 1e-4)
      {
         nValidPointsWaveHeight += 1;
      }
      else
      {
         break;
      }
   }
   nValidPointsWaveHeight -= 1;

   // TODO 007 What is this for?
   for (int nPoint = 0; nPoint < static_cast<int>(VdWaveSetupSurge.size()); nPoint++)
   {
      if (tAbs(VdWaveSetupSurge[nPoint]) < 1)    // limiting the absolute value of setup + surge if cshore run fails
      {
         nValidPointsWaveSetup += 1;
      }
      else
      {
         break;
      }
   }
   nValidPointsWaveSetup -= 1;

   // Check DBL_NODATA
   double dWaveHeight = 0;
   if (! bFPIsEqual(VdWaveHeight[nValidPointsWaveHeight], DBL_NODATA, TOLERANCE))
   {
      dWaveHeight = VdWaveHeight[nValidPointsWaveHeight];
   }

   // TODO 060 Remove these 'magic numbers'
   double dRunUp = 0;
   if (m_nRunUpEquation == 0)
   {
      // Compute the run-up using NIELSEN & HANSLOW (1991) & DHI (2004) // TODO 007
      dRunUp = 0.36 * pow(9.81, 0.5) * dtanBeta * pow(dWaveHeight, 0.5) * dDeepWaterWavePeriod;
   }
   else if (m_nRunUpEquation == 1)
   {
      // Compute the run-up using MASE 1989
      double dS0 = 2 * PI * dWaveHeight / (9.81 * dDeepWaterWavePeriod * dDeepWaterWavePeriod);
      dRunUp = 1.86 * dWaveHeight * pow(pow(dtanBeta / dS0, 0.5), 0.71);
   }
   else if (m_nRunUpEquation == 2)
   {
      // Compute the run-up using STOCKDON (2006)
      double dS0 = 2 * PI * dWaveHeight / (9.81 * dDeepWaterWavePeriod * dDeepWaterWavePeriod);
      // dRunUp = 1.1 * ((0.35 * dWaveHeight * (pow((1 / dS0) * dWaveHeight * dWaveHeight, 0.5))) + (((((1 / dS0) * dWaveHeight * dWaveHeight) * (0.563 * dWaveHeight * dWaveHeight + 0.0004)), 0.5)) / 2);

      double dH0OverL0 = (1 / dS0) * dWaveHeight;
      double dTmp1 = 0.35 * dWaveHeight * pow(dH0OverL0, 0.5);
      double dTmp2 = pow(dH0OverL0 * ((0.563 * dWaveHeight * dWaveHeight) + 0.0004), 0.5);
      dRunUp = 1.1 * (dTmp1 + (dTmp2 / 2));
   }

   if ((tAbs(dRunUp) < 1e-4) || (! isfinite(dRunUp)))
   {
      dRunUp = 0;
   }

   // Check for DBL_NODATA
   double dWaveSetupSurge = 0;               // TODO 007 Info needed
   if (! bFPIsEqual(VdWaveSetupSurge[nValidPointsWaveSetup], DBL_NODATA, TOLERANCE))
   {
      dWaveSetupSurge = VdWaveSetupSurge[nValidPointsWaveSetup];
   }
   
   if ((tAbs(dWaveSetupSurge) < 1e-4) || (! isfinite(dWaveSetupSurge)))
   {
      dWaveSetupSurge = 0;
   }

   // Update wave attributes along the coastline object. Wave height at the coast is always calculated (i.e. whether or not waves are breaking)
   // cout << "Wave Height at the coast is " << VdWaveHeight[nProfileSize - 1] << endl;
   m_VCoast[nCoast].SetCoastWaveHeight(nCoastPoint, dWaveHeight);
   m_VCoast[nCoast].SetWaveSetupSurge(nCoastPoint, dWaveSetupSurge);       // TODO 007 Info needed
   m_VCoast[nCoast].SetRunUp(nCoastPoint, dRunUp);

   if (nProfileBreakingDist > 0)
   {
      // This coast point is in the active zone, so set breaking wave height, breaking wave angle, and depth of breaking for the coast point
      m_VCoast[nCoast].SetBreakingWaveHeight(nCoastPoint, dProfileBreakingWaveHeight);
      m_VCoast[nCoast].SetBreakingWaveAngle(nCoastPoint, dProfileBreakingWaveAngle);
      m_VCoast[nCoast].SetDepthOfBreaking(nCoastPoint, dProfileBreakingDepth);
      m_VCoast[nCoast].SetBreakingDistance(nCoastPoint, nProfileBreakingDist);

      //       LogStream << m_ulIter << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " in active zone, dBreakingWaveHeight = " << dBreakingWaveHeight << endl;
   }
   else
   {
      // This coast point is not in the active zone, so breaking wave height, breaking wave angle, and depth of breaking are all meaningless
      m_VCoast[nCoast].SetBreakingWaveHeight(nCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingWaveAngle(nCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetDepthOfBreaking(nCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingDistance(nCoastPoint, INT_NODATA);

      //       LogStream << m_ulIter << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " NOT in active zone" << endl;
   }

   return RTN_OK;
}

#if defined CSHORE_FILE_INOUT
//===============================================================================================================================
//! Create and write to the CShore input file
//===============================================================================================================================
int CSimulation::nCreateCShoreInfile(int const nCoast, int const nProfile, int const nILine, int const nIProfl, int const nIPerm, int const nIOver, int const nIWcint, int const nIRoll, int const nIWind, int const nITide, int const nILab, int const nWave, int const nSurge, double const dX, double const dTimestep, double const dWaveInitTime, double const dWavePeriod, double const dHrms, double const dWaveAngle, double const dSurgeInitTime, double const dSurgeLevel, vector<double> const* pVdXdist, vector<double> const* pVdBottomElevation, vector<double> const* pVdWaveFriction)
{
   // Create the CShore input file
   ofstream CShoreOutStream;
   CShoreOutStream.open(CSHORE_INFILE.c_str(), ios::out | ios::app);
   if (CShoreOutStream.fail())
   {
      // Error, cannot open file for writing
      LogStream << m_ulIter << ": " << ERR << "cannot write to CShore input file '" << CSHORE_INFILE << "'" << endl;
      return RTN_ERR_CSHORE_FILE_INPUT;
   }

   // And write to the file
   CShoreOutStream << 3 << endl; // Number of comment lines
   CShoreOutStream << "------------------------------------------------------------" << endl;
   CShoreOutStream << "CShore input file created by CoastalME for iteration " << m_ulIter << ", coast " << nCoast << ", profile " << nProfile << endl;
   CShoreOutStream << "------------------------------------------------------------" << endl;
   CShoreOutStream << nILine << "                                         -> ILINE" << endl;
   CShoreOutStream << nIProfl << "                                         -> IPROFL" << endl;
   CShoreOutStream << nIPerm << "                                         -> IPERM" << endl;
   CShoreOutStream << nIOver << "                                         -> IOVER" << endl;
   CShoreOutStream << nIWcint << "                                         -> IWCINT" << endl;
   CShoreOutStream << nIRoll << "                                         -> IROLL" << endl;
   CShoreOutStream << nIWind << "                                         -> IWIND" << endl;
   CShoreOutStream << nITide << "                                         -> ITIDE" << endl;
   CShoreOutStream << std::fixed;
   CShoreOutStream << setw(11) << setprecision(4) << dX << "                               -> DX" << endl;
   CShoreOutStream << setw(11) << m_dBreakingWaveHeightDepthRatio << "                               -> GAMMA" << endl;
   CShoreOutStream << setw(11) << nILab << "                               -> ILAB" << endl;
   CShoreOutStream << setw(11) << nWave << "                               -> NWAVE" << endl;
   CShoreOutStream << setw(11) << nSurge << "                               -> NSURGE" << endl;

   // Line 18 of infile
   CShoreOutStream << setw(11) << setprecision(2) << dWaveInitTime;           // TWAVE(1) in CShore
   CShoreOutStream << setw(11) << setprecision(4) << dWavePeriod;             // TPIN(1) in CShore
   CShoreOutStream << setw(11) << dHrms;                                      // HRMS(1) in CShore
   CShoreOutStream << setw(11) << dWaveAngle << endl;                         // WANGIN(1) in CShore

   // Line 19 of infile
   CShoreOutStream << setw(11) << setprecision(2) << dTimestep;               // TWAVE(2) in CShore
   CShoreOutStream << setw(11) << setprecision(4) << dWavePeriod;             // TPIN(2) in CShore
   CShoreOutStream << setw(11) << dHrms;                                      // HRMS(2) in CShore
   CShoreOutStream << setw(11) << dWaveAngle << endl;                         // WANGIN(2) in CShore

   // Line 20 of infile
   CShoreOutStream << setw(11) << setprecision(2) << dSurgeInitTime;          // TSURG(1) in CShore
   CShoreOutStream << setw(11) << setprecision(4) << dSurgeLevel << endl;     // SWLIN(1) in CShore

   // Line 21 of infile
   CShoreOutStream << setw(11) << setprecision(2) << dTimestep;               // TSURG(2) in CShore
   CShoreOutStream << setw(11) << setprecision(4) << dSurgeLevel << endl;     // SWLIN(2) in CShore

   // Line 22 of infile
   CShoreOutStream << setw(8) << pVdXdist->size() << "                                  -> NBINP" << endl;

   CShoreOutStream << std::fixed << setprecision(4);
   for (unsigned int i = 0; i < pVdXdist->size(); i++)
      // These are BINP(J,1), ZBINP(J,1), FBINP(J-1,1) in CShore
      CShoreOutStream << setw(11) << pVdXdist->at(i) << setw(11) << pVdBottomElevation->at(i) << setw(11) << pVdWaveFriction->at(i) << endl;

   CShoreOutStream << endl;

   // File written, so close it
   CShoreOutStream.close();

   return RTN_OK;
}
#endif

//===============================================================================================================================
//! Get profile horizontal distance and bottom elevation vectors in CShore units
//===============================================================================================================================
int CSimulation::nGetThisProfileElevationVectorsForCShore(int const nCoast, int const nProfile, int const nProfSize, vector<double> *VdDistXY, vector<double> *VdVZ, vector<double> *VdFricF)
{
   bool bIsBehindIntervention = false;

   int nX1 = 0;
   int nY1 = 0;

   double dXDist;
   double dYDist;
   double dProfileDistXY = 0;
   double dProfileFricFact = 0;
   double dPrevDist = -1;

   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   for (int i = nProfSize - 1; i >= 0; i--)
   {
      int nX = pProfile->pPtiVGetCellsInProfile()->at(i).nGetX();
      int nY = pProfile->pPtiVGetCellsInProfile()->at(i).nGetY();

      // Calculate the horizontal distance relative to the most seaward point
      if (i == nProfSize - 1)
         dProfileDistXY = 0;
      else
      {
         dXDist = dGridCentroidXToExtCRSX(nX1) - dGridCentroidXToExtCRSX(nX),
         dYDist = dGridCentroidYToExtCRSY(nY1) - dGridCentroidYToExtCRSY(nY),
         dProfileDistXY = dProfileDistXY + hypot(dXDist, dYDist);
      }

      // Before we store the X-Y distance, must check that it is not the same as the previously-stored distance (if it is, we get zero-divide errors in CShore). If they are the same, skip this one
      if (bFPIsEqual(dProfileDistXY, dPrevDist, TOLERANCE))
         continue;

      // Update the cell indexes, the initial cell is now the previous one
      nX1 = nX;
      nY1 = nY;

      // Get the number of the highest layer with non-zero thickness
      int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

      // Safety check
      if (nTopLayer == INT_NODATA)
         return RTN_ERR_NO_TOP_LAYER;

      if (nTopLayer == NO_NONZERO_THICKNESS_LAYERS)
         // TODO 009 We are down to basement, decide what to do
         return RTN_OK;

      // Get the elevation for both consolidated and unconsolidated sediment on this cell
      double VdProfileZ = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() + m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight() - m_dThisIterSWL;

      // // TODO 072 There is a problem here, is it because of -ve VdProfileZ values?
      // if (VdProfileZ < 0)
      // {
      //    VdProfileZ = 0;
      //    LogStream << m_ulIter << ": " << WARN << " elevation is -ve, changing for CShore profile elevation vector " << nProfile << endl;
      // }

      // Check that landward location is +ve
      if (i == 0)
      {
         if (VdProfileZ < 0)
         {
            VdVZ->push_back(0.5);      // TODO 053 Set it to a small +ve elevation. However there must be a better way of doing this
            
            // Could not create the profile elevation vectors
            LogStream << m_ulIter << ": " << WARN << " landward location is -ve, changing for CShore profile elevation vector " << nProfile << endl;
         }
         else
         {
            VdVZ->push_back(VdProfileZ);
         }
      }
      else
      {
         VdVZ->push_back(VdProfileZ);
      }

      // Now store the X-Y plane distance from the start of the profile
      VdDistXY->push_back(dProfileDistXY);

      // Get the landform type at each point along the profile
      double dInterventionHeight = m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight();

      // Modify default friction factor if a structural intervention is found, otherwise use the default
      if (dInterventionHeight > 0 || bIsBehindIntervention)
      {
         // Use an arbitrarily high value if structure is present
         dProfileFricFact = 100 * CSHORE_FRICTION_FACTOR;
         bIsBehindIntervention = true;
      }
      else
         dProfileFricFact = CSHORE_FRICTION_FACTOR;

      // Store the friction factor
      VdFricF->push_back(dProfileFricFact);

      // For next time round
      dPrevDist = dProfileDistXY;
   }

   return RTN_OK;
}

#if defined CSHORE_FILE_INOUT
//===============================================================================================================================
//! Reads a CShore output file and creates a vector holding interpolated values
//===============================================================================================================================
int CSimulation::nReadCShoreOutput(int const nProfile, string const *strCShoreFilename, int const nExpectedColumns, int const nCShorecolumn, vector<double> const* pVdProfileDistXYCME, vector<double>* pVdInterpolatedValues)
{
   // Read in the first column (contains XY distance relative to seaward limit) and CShore column from the CShore output file
   ifstream InStream;
   InStream.open(strCShoreFilename->c_str(), ios::in);

   // Did it open OK?
   if (! InStream.is_open())
   {
      // Error: cannot open CShore file for input
      LogStream << m_ulIter << ": " << ERR << "for profile " << nProfile << ", cannot open " << *strCShoreFilename << " for input" << endl;

      return RTN_ERR_READING_CSHORE_FILE_OUTPUT;
   }

   // Opened OK, so set up the vectors to hold the CShore output data
   vector<double> VdXYDistCShore;
   vector<double> VdValuesCShore;

   // And read in the data
   int n = -1;
   int nExpectedRows = 0;
   string strLineIn;
   while (getline(InStream, strLineIn))
   {
      n++;

      if (n == 0)
      {
         // Read in the header line
         vector<string> VstrItems = VstrSplit(&strLineIn, SPACE);

         if (! bIsStringValidInt(VstrItems[1]))
         {
            string strErr = ERR + "invalid integer for number of expected rows '" + VstrItems[1] + "' in " + *strCShoreFilename + "\n";
            cerr << strErr;
            LogStream << strErr;

            return RTN_ERR_READING_CSHORE_FILE_OUTPUT;
         }

         // Get the number of expected rows
         nExpectedRows = stoi(VstrItems[1]);
      }
      else
      {
         // Read in a data line
         vector<string> VstrItems = VstrSplit(&strLineIn, SPACE);

         int nCols = static_cast<int>(VstrItems.size());
         if (nCols != nExpectedColumns)
         {
            // Error: did not read the expected number of CShore output columns
            LogStream << m_ulIter << ": " << ERR << "for profile " << nProfile << ", expected " << nExpectedColumns << " CShore output columns but read " << nCols << " columns from header section of file " << *strCShoreFilename << endl;

            return RTN_ERR_READING_CSHORE_FILE_OUTPUT;
         }

         // Number of columns is OK
         VdXYDistCShore.push_back(strtod(VstrItems[0].c_str(), NULL));
         VdValuesCShore.push_back(strtod(VstrItems[nCShorecolumn - 1].c_str(), NULL));
      }
   }

   // Check that we have read nExpectedRows from the file
   int nReadRows = static_cast<int>(VdXYDistCShore.size());
   if (nReadRows != nExpectedRows)
   {
      // Error: did not get nExpectedRows CShore output rows
      LogStream << m_ulIter << ": " << ERR << "for profile " << nProfile << ", expected " << nExpectedRows << " CShore output rows, but read " << nReadRows << " rows from file " << *strCShoreFilename << endl;

      return RTN_ERR_READING_CSHORE_FILE_OUTPUT;
   }

   if (nReadRows < 2)
   {
      // CShore sometimes returns only one row, which contains data for the seaward point of the profile. This happens when all other (more coastward) points give an invalid result during CShore's calculations. This is a problem. We don't want to abandon the simulation just because of this, so instead we just duplicate the row, so that the profile will later get marked as invalid
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << m_ulIter << ": " << WARN << "for profile " << nProfile << ", only " << nReadRows << " CShore output rows in file " << *strCShoreFilename << endl;

      // Duplicate the data
      VdXYDistCShore.push_back(VdXYDistCShore[0]);
      VdValuesCShore.push_back(VdValuesCShore[0]);

      // And increase the expected number of rows
      nReadRows++;
   }

   // The output is OK, so change the origin of the across-shore distance from the CShore convention to the one used here (i.e. with the origin at the shoreline)
   vector<double> VdXYDistCShoreTmp(nReadRows, 0);
   for (int i = 0; i < nReadRows; i++)
      VdXYDistCShoreTmp[i] = VdXYDistCShore[nReadRows - 1] - VdXYDistCShore[i];

   // Reverse the XY-distance and value vectors (i.e. first point is at the shoreline and must be in strictly ascending order)
   reverse(VdXYDistCShoreTmp.begin(), VdXYDistCShoreTmp.end());

   // Similarly, reverse the CShore output
   reverse(VdValuesCShore.begin(), VdValuesCShore.end());

   // Using a simple linear interpolation approach
   vector<double> VdDistXYCopy(pVdProfileDistXYCME->begin(), pVdProfileDistXYCME->end());

   // assertVdXYDistCShoreTmp.size() == VdValuesCShore.size());
   *pVdInterpolatedValues = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, &VdValuesCShore, &VdDistXYCopy);

   return RTN_OK;
}
#endif

#if defined CSHORE_ARG_INOUT || CSHORE_BOTH
//===============================================================================================================================
//! Interpolates CShore output
//===============================================================================================================================
void CSimulation::InterpolateCShoreOutput(vector<double> const* pVdProfileDistXYCME, int const nOutSize, vector<double> const* pVdXYDistFromCShoreOut, vector<double> const* pVdFreeSurfaceStdCShore, vector<double> const* pVdWaveSetupSurgeCShore, vector<double> const* pVdSinWaveAngleRadiansCShore, vector<double> const* pVdFractionBreakingWavesCShore, vector<double>* pVdFreeSurfaceStdCME, vector<double>* pVdWaveSetupSurgeCME, vector<double>* pVdSinWaveAngleRadiansCME, vector<double>* pVdFractionBreakingWavesCME)
{
   // The CShore cross-shore distance has its origin at the seaward end, so create a copy of the valid part of this which is in the CME convention (i.e. with the origin at the shoreline)
   vector<double> VdXYDistCShoreTmp(nOutSize, 0);
   for (int i = 0; i < nOutSize; i++)
      VdXYDistCShoreTmp[i] = pVdXYDistFromCShoreOut->at(nOutSize - 1) - pVdXYDistFromCShoreOut->at(i);

   vector<double> VdFreeSurfaceStdCShoreTmp(pVdFreeSurfaceStdCShore->begin(), pVdFreeSurfaceStdCShore->begin() + nOutSize);
   vector<double> VdWaveSetupSurgeCShoreTmp(pVdWaveSetupSurgeCShore->begin(), pVdWaveSetupSurgeCShore->begin() + nOutSize);   // TODO 057 Sometimes get -ve values here, is this OK?
   vector<double> VdSinWaveAngleRadiansCShoreTmp(pVdSinWaveAngleRadiansCShore->begin(), pVdSinWaveAngleRadiansCShore->begin() + nOutSize);
   vector<double> VdFractionBreakingWavesCShoreTmp(pVdFractionBreakingWavesCShore->begin(), pVdFractionBreakingWavesCShore->begin() + nOutSize);
   
   // And reverse these to fit the CME convention, with the origin at the shoreline
   reverse(VdFreeSurfaceStdCShoreTmp.begin(), VdFreeSurfaceStdCShoreTmp.end());
   reverse(VdWaveSetupSurgeCShoreTmp.begin(), VdWaveSetupSurgeCShoreTmp.end());
   // reverse(VdStormSurgeCShoreTmp.begin(), VdStormSurgeCShoreTmp.end());
   // reverse(VdWaveSetupRunUpCShoreTmp.begin(), VdWaveSetupRunUpCShoreTmp.end());
   reverse(VdSinWaveAngleRadiansCShoreTmp.begin(), VdSinWaveAngleRadiansCShoreTmp.end());
   reverse(VdFractionBreakingWavesCShoreTmp.begin(), VdFractionBreakingWavesCShoreTmp.end());

   // Using the simple linear approach
   vector<double> VdDistXYCopy(pVdProfileDistXYCME->begin(), pVdProfileDistXYCME->end());

   *pVdFreeSurfaceStdCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdFreeSurfaceStdCShore, &VdDistXYCopy);
   *pVdWaveSetupSurgeCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdWaveSetupSurgeCShore, &VdDistXYCopy);
   // *pVdStormSurgeCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdStormSurgeCShore, &VdDistXYCopy);
   // *pVdWaveSetupRunUpCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdWaveSetupRunUpCShore, &VdDistXYCopy);
   *pVdSinWaveAngleRadiansCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdSinWaveAngleRadiansCShore, &VdDistXYCopy);
   *pVdFractionBreakingWavesCME = VdInterpolateCShoreProfileOutput(&VdXYDistCShoreTmp, pVdFractionBreakingWavesCShore, &VdDistXYCopy);
}
#endif

//===============================================================================================================================
//! Modifies the wave breaking properties at coastline points of profiles within the shadow zone
//===============================================================================================================================
void CSimulation::ModifyBreakingWavePropertiesWithinShadowZoneToCoastline(int const nCoast, int const nProfIndex)
{
   int nProfile = m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(nProfIndex);
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // Only do this for profiles without problems, including the start and end-of-coast profile
   if (!pProfile->bOKIncStartAndEndOfCoast())
      return;

   int
       nThisCoastPoint = pProfile->nGetNumCoastPoint(),
       nProfileSize = pProfile->nGetNumCellsInProfile();

   int nThisBreakingDist = m_VCoast[nCoast].nGetBreakingDistance(nThisCoastPoint);
   double
       dThisBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nThisCoastPoint), // This could be DBL_NODATA
       dThisBreakingWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nThisCoastPoint),
       dThisBreakingDepth = m_VCoast[nCoast].dGetDepthOfBreaking(nThisCoastPoint);
   bool
       bModfiedWaveHeightisBreaking = false,
       bProfileIsinShadowZone = false;

   // Traverse the profile landwards, checking if any profile cell is within the shadow zone
   for (int nProfilePoint = (nProfileSize - 1); nProfilePoint >= 0; nProfilePoint--)
   {
      int
          nX = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetX(),
          nY = pProfile->pPtiGetCellInProfile(nProfilePoint)->nGetY();

      // If there is any cell profile  within the shadow zone and waves are breaking then modify wave breaking properties otherwise continue
      if (m_pRasterGrid->m_Cell[nX][nY].bIsinAnyShadowZone())
      {
         bProfileIsinShadowZone = true;

         // Check if the new wave height is breaking
         double
             // dSeaDepth = m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(),
             dWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight(),
             dWaveAngle = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle();

         // Check that wave height at the given point is lower than maximum real wave height. If breaking wave height is expected that no good wave height are obtained, so, do not take it
         if (dWaveHeight > (m_dDepthOfClosure * m_dBreakingWaveHeightDepthRatio) && (! bModfiedWaveHeightisBreaking) && (! bFPIsEqual(dThisBreakingWaveHeight, DBL_NODATA, TOLERANCE)))
         {
            // It is breaking
            bModfiedWaveHeightisBreaking = true;

            dThisBreakingWaveHeight = m_dDepthOfClosure * m_dBreakingWaveHeightDepthRatio;
            dThisBreakingWaveAngle = dWaveAngle;
            dThisBreakingDepth = m_dDepthOfClosure;
            nThisBreakingDist = nProfilePoint;
         }
      }
   }

   // Update breaking wave properties along coastal line object (Wave height, dir, distance). TODO 010 Update the active zone cells
   if (bProfileIsinShadowZone && bModfiedWaveHeightisBreaking) // Modified wave height is still breaking
   {
      // This coast point is in the active zone, so set breaking wave height, breaking wave angle, and depth of breaking for the coast point TODO 007 Where does the 0.78 come from? TODO 011 Should it be an input variable or a named constant?
      if (dThisBreakingWaveHeight > dThisBreakingDepth * 0.78)
      {
         dThisBreakingWaveHeight = dThisBreakingDepth * 0.78; // Likely CShore output wave height is not adequately reproduced due to input profile and wave properties. TODO 007 Info needed. Does something need to be changed then?
      }

      // assert(dThisBreakingWaveHeight >= 0);
      m_VCoast[nCoast].SetBreakingWaveHeight(nThisCoastPoint, dThisBreakingWaveHeight);
      m_VCoast[nCoast].SetBreakingWaveAngle(nThisCoastPoint, dThisBreakingWaveAngle);
      m_VCoast[nCoast].SetDepthOfBreaking(nThisCoastPoint, dThisBreakingDepth);
      m_VCoast[nCoast].SetBreakingDistance(nThisCoastPoint, nThisBreakingDist);

      //       LogStream << m_ulIter << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " in active zone, dBreakingWaveHeight = " << dBreakingWaveHeight << endl;
   }
   else if (bProfileIsinShadowZone && (! bModfiedWaveHeightisBreaking))
   {
      // This coast point is no longer in the active zone
      m_VCoast[nCoast].SetBreakingWaveHeight(nThisCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingWaveAngle(nThisCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetDepthOfBreaking(nThisCoastPoint, DBL_NODATA);
      m_VCoast[nCoast].SetBreakingDistance(nThisCoastPoint, INT_NODATA);

      //       LogStream << m_ulIter << ": nProfile = " << nProfile << ", nCoastPoint = " << nCoastPoint << " NOT in active zone" << endl;
   }

   return;
}

//===============================================================================================================================
//! Interpolates wave properties from profiles to the coastline points between two profiles. Do this by weighting the wave properties so that they change smoothly between the two profiles
//===============================================================================================================================
void CSimulation::InterpolateWavePropertiesBetweenProfiles(int const nCoast, int const nProfIndex, int const nNumProfiles)
{
   int nProfile = m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(nProfIndex);
   CGeomProfile const* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // Only do this for profiles without problems, including the start-of-coast profile (but not the end-of-coast profile)
   // if (!pProfile->bOKIncStartOfCoast())
   //    return;

   int nThisCoastPoint = pProfile->nGetNumCoastPoint();

   // For the breaking wave stuff, to go into the in-between coastline points
   int const nThisBreakingDist = m_VCoast[nCoast].nGetBreakingDistance(nThisCoastPoint);
   double
       dThisBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nThisCoastPoint), // This could be DBL_NODATA
       dThisBreakingWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nThisCoastPoint),
       dThisBreakingDepth = m_VCoast[nCoast].dGetDepthOfBreaking(nThisCoastPoint),
       dThisWaveSetupSurge = m_VCoast[nCoast].dGetWaveSetupSurge(nThisCoastPoint),
       //  dThisStormSurge = m_VCoast[nCoast].dGetStormSurge(nThisCoastPoint),
       dThisRunUp = m_VCoast[nCoast].dGetRunUp(nThisCoastPoint);

   // Get the number of the next profile along the coast. If this next profile has a problem, go to the one after that, etc
   int nNextProfile = -1;
   for (int nNextProfIndex = nProfIndex + 1; nNextProfIndex < nNumProfiles; nNextProfIndex++)
   {
      nNextProfile = m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(nNextProfIndex);

      if (m_VCoast[nCoast].pGetProfile(nNextProfile)->bOKIncStartAndEndOfCoast())
         break;
   }

   // If last profile has a problem, take the value of this profile
   int const
       nNextCoastPoint = m_VCoast[nCoast].pGetProfile(nNextProfile)->nGetNumCoastPoint(),
       nDistBetween = nNextCoastPoint - nThisCoastPoint;

   if (nDistBetween <= 0)
      // Nothing to do
      return;

   int nNextBreakingDist = m_VCoast[nCoast].nGetBreakingDistance(nNextCoastPoint);
   double
       dNextBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nNextCoastPoint), // This could be DBL_NODATA
       dNextBreakingWaveAngle = m_VCoast[nCoast].dGetBreakingWaveAngle(nNextCoastPoint),
       dNextBreakingDepth = m_VCoast[nCoast].dGetDepthOfBreaking(nNextCoastPoint),
       dNextWaveSetupSurge = m_VCoast[nCoast].dGetWaveSetupSurge(nNextCoastPoint),
       dNextRunUp = m_VCoast[nCoast].dGetRunUp(nNextCoastPoint);

   // OK, fill coast point between profiles for setupsurge and runup
   for (int n = nThisCoastPoint; n <= nNextCoastPoint; n++)
   {
      // Fill first wave setup and surge
      int nDist = n - nThisCoastPoint;
      double
          dThisWeight = (nDistBetween - nDist) / static_cast<double>(nDistBetween),
          dNextWeight = 1 - dThisWeight,
          dWaveSetupSurge = 0,
          dRunUp = 0;

      dWaveSetupSurge = (dThisWeight * dThisWaveSetupSurge) + (dNextWeight * dNextWaveSetupSurge);
      m_VCoast[nCoast].SetWaveSetupSurge(n, dWaveSetupSurge);

      dRunUp = (dThisWeight * dThisRunUp) + (dNextWeight * dNextRunUp);
      m_VCoast[nCoast].SetRunUp(n, dRunUp);
   }

   // If both this profile and the next profile are not in the active zone, then do no more
   if ((bFPIsEqual(dThisBreakingWaveHeight, DBL_NODATA, TOLERANCE)) && (bFPIsEqual(dNextBreakingWaveHeight, DBL_NODATA, TOLERANCE)))
   {
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ": both profile " << nProfile << " at coast point " << nThisCoastPoint << ", and profile " << nNextProfile << " at coast point " << nNextCoastPoint << ", are not in the active zone" << endl;
      
      // Set the breaking wave height, breaking wave angle, and depth of breaking to DBL_NODATA
      for (int n = nThisCoastPoint; n < nNextCoastPoint; n++)
      {
         m_VCoast[nCoast].SetBreakingWaveHeight(nThisCoastPoint, DBL_NODATA);
         m_VCoast[nCoast].SetBreakingWaveAngle(nThisCoastPoint, DBL_NODATA);
         m_VCoast[nCoast].SetDepthOfBreaking(nThisCoastPoint, DBL_NODATA);
         m_VCoast[nCoast].SetBreakingDistance(nThisCoastPoint, DBL_NODATA);
      }
      return;
   }

   // OK, at least one of the two profiles is in the active zone
   if (bFPIsEqual(dThisBreakingWaveHeight, DBL_NODATA, TOLERANCE)) 
   {
      // The next profile must be in the active zone, so use values from the next profile
      for (int n = nThisCoastPoint; n < nNextCoastPoint; n++)
      {
         // Set the breaking wave height, breaking wave angle, and depth of breaking for this coast point TODO 056 This assert sometimes fails: why?
         // assert(dNextBreakingWaveHeight >= 0);
         m_VCoast[nCoast].SetBreakingWaveHeight(n, dNextBreakingWaveHeight);
         m_VCoast[nCoast].SetBreakingWaveAngle(n, dNextBreakingWaveAngle);
         m_VCoast[nCoast].SetDepthOfBreaking(n, dNextBreakingDepth);
         m_VCoast[nCoast].SetBreakingDistance(n, nNextBreakingDist);
      }
      return;
   }

   if (bFPIsEqual(dNextBreakingWaveHeight, DBL_NODATA, TOLERANCE))
   {
      // This profile must be in the active zone, so use values from this profile
      for (int n = nThisCoastPoint + 1; n <= nNextCoastPoint; n++)
      {
         // Set the breaking wave height, breaking wave angle, and depth of breaking for this coast point TODO 056 This assert sometimes fails: why?
         // assert(dThisBreakingWaveHeight >= 0);
         m_VCoast[nCoast].SetBreakingWaveHeight(n, dThisBreakingWaveHeight);
         m_VCoast[nCoast].SetBreakingWaveAngle(n, dThisBreakingWaveAngle);
         m_VCoast[nCoast].SetDepthOfBreaking(n, dThisBreakingDepth);
         m_VCoast[nCoast].SetBreakingDistance(n, nThisBreakingDist);
      }

      return;
   }

   // The values for both this profile point and the next profile point are fine, so do a weighted interpolation between this profile and the next profile
   for (int n = nThisCoastPoint + 1; n < nNextCoastPoint; n++)
   {
      int nDist = n - nThisCoastPoint;

      double
          dBreakingWaveHeight = DBL_NODATA,
          dBreakingWaveAngle = DBL_NODATA,
          dBreakingDepth = DBL_NODATA,
          dBreakingDist = DBL_NODATA;

      if ((dNextBreakingDepth > 0) && (dThisBreakingDepth > 0))
      {
         double
             dThisWeight = (nDistBetween - nDist) / static_cast<double>(nDistBetween),
             dNextWeight = 1 - dThisWeight;

         dBreakingWaveHeight = (dThisWeight * dThisBreakingWaveHeight) + (dNextWeight * dNextBreakingWaveHeight),
         dBreakingWaveAngle = (dThisWeight * dThisBreakingWaveAngle) + (dNextWeight * dNextBreakingWaveAngle),
         dBreakingDepth = (dThisWeight * dThisBreakingDepth) + (dNextWeight * dNextBreakingDepth),
         dBreakingDist = (dThisWeight * nThisBreakingDist) + (dNextWeight * nNextBreakingDist);
      }
      else if (dNextBreakingDepth > 0)
      {
         dBreakingWaveHeight = dNextBreakingWaveHeight,
         dBreakingWaveAngle = dNextBreakingWaveAngle,
         dBreakingDepth = dNextBreakingDepth,
         dBreakingDist = nNextBreakingDist;
      }
      else if (dThisBreakingDepth > 0)
      {
         dBreakingWaveHeight = dThisBreakingWaveHeight,
         dBreakingWaveAngle = dThisBreakingWaveAngle,
         dBreakingDepth = dThisBreakingDepth,
         dBreakingDist = nThisBreakingDist;
      }

      // Set the breaking wave height, breaking wave angle, and depth of breaking for this coast point TODO 056 This assert sometimes fails: why?
      // assert(dBreakingWaveHeight >= 0);
      m_VCoast[nCoast].SetBreakingWaveHeight(n, dBreakingWaveHeight);
      m_VCoast[nCoast].SetBreakingWaveAngle(n, dBreakingWaveAngle);
      m_VCoast[nCoast].SetDepthOfBreaking(n, dBreakingDepth);
      m_VCoast[nCoast].SetBreakingDistance(n, nRound(dBreakingDist));
   }
}

//===============================================================================================================================
//! Linearly interpolates wave properties from profiles to the coastline cells between two profiles
//===============================================================================================================================
void CSimulation::InterpolateWaveHeightToCoastPoints(int const nCoast)
{
   int nCoastPoints = m_VCoast[nCoast].nGetCoastlineSize();

   // Initialize all vectors pairs (x,y) for each variable
   vector<int> nVCoastWaveHeightX;
   vector<double> dVCoastWaveHeightY;

   // Search all coast points for non NaN values and store them into temporary variables for later interporlation
   for (int n = 0; n < nCoastPoints; n++)
   {
      double dCoastWaveHeight = m_VCoast[nCoast].dGetCoastWaveHeight(n);

      if (! bFPIsEqual(dCoastWaveHeight, DBL_NODATA, TOLERANCE))
      {
         nVCoastWaveHeightX.push_back(n);
         dVCoastWaveHeightY.push_back(dCoastWaveHeight);
      }
   }

   // Interpolate all coast points. Check first that x,y have more than 3 points and are both of equal size
   if ((nVCoastWaveHeightX.size() >= 3) && (nVCoastWaveHeightX.size() == dVCoastWaveHeightY.size()))
   {
      for (int n = 0; n < nCoastPoints; n++)
      {
         double dInterpCoastWaveHeight = dGetInterpolatedValue(&nVCoastWaveHeightX, &dVCoastWaveHeightY, n, false);
         m_VCoast[nCoast].SetCoastWaveHeight(n, dInterpCoastWaveHeight);
      }
   }
   else
   {
      for (int n = 0; n < nCoastPoints; n++)
      {
         m_VCoast[nCoast].SetCoastWaveHeight(n, 0);
      }
   }
}

//===============================================================================================================================
//! Calculates tangents to a coastline: the tangent is assumed to be the orientation of energy/sediment flux along a coast. The tangent is specified as an angle (in degrees) measured clockwise from north. Based on a routine by Martin Hurst
//===============================================================================================================================
void CSimulation::CalcCoastTangents(int const nCoast)
{
   int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
   {
      double
          dXDiff,
          dYDiff;

      if (nCoastPoint == 0)
      {
         // For the point at the start of the coastline: use the straight line from 'this' point to the next point
         CGeom2DPoint PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);      // In external CRS
         CGeom2DPoint PtAfter = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint + 1); // In external CRS

         dXDiff = PtAfter.dGetX() - PtThis.dGetX();
         dYDiff = PtAfter.dGetY() - PtThis.dGetY();
      }
      else if (nCoastPoint == nCoastSize - 1)
      {
         // For the point at the end of the coastline: use the straight line from the point before to 'this' point
         CGeom2DPoint PtBefore = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint - 1); // In external CRS
         CGeom2DPoint PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint);       // In external CRS

         dXDiff = PtThis.dGetX() - PtBefore.dGetX();
         dYDiff = PtThis.dGetY() - PtBefore.dGetY();
      }
      else
      {
         // For coastline points not at the start or end of the coast: start with a straight line which links the coastline points before and after 'this' coastline point
         CGeom2DPoint PtBefore = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint - 1); // In external CRS
         CGeom2DPoint PtAfter = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint + 1);  // In external CRS

         dXDiff = PtAfter.dGetX() - PtBefore.dGetX();
         dYDiff = PtAfter.dGetY() - PtBefore.dGetY();
      }

      // Calculate angle between line and north point, measured clockwise (the azimuth)
      if (bFPIsEqual(dYDiff, 0.0, TOLERANCE))
      {
         // The linking line runs either W-E or E-W
         if (dXDiff > 0)
            // It runs W to E
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 90);
         else
            // It runs E to W
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 270);
      }
      else if (bFPIsEqual(dXDiff, 0.0, TOLERANCE))
      {
         // The linking line runs N-S or S-N
         if (dYDiff > 0)
            // It runs S to N
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 0);
         else
            // It runs N to S
            m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, 180);
      }
      else
      {
         // The linking line runs neither W-E nor N-S so we have to work a bit harder to find the angle between it and the azimuth
         double dAzimuthAngle;
         if (dXDiff > 0)
            dAzimuthAngle = (180 / PI) * (PI * 0.5 - atan(dYDiff / dXDiff));
         else
            dAzimuthAngle = (180 / PI) * (PI * 1.5 - atan(dYDiff / dXDiff));

         m_VCoast[nCoast].SetFluxOrientation(nCoastPoint, dAzimuthAngle);
      }

      //      LogStream << m_ulIter << ": nCoastPoint = " << nCoastPoint << " FluxOrientation = " << m_VCoast[nCoast].dGetFluxOrientation(nCoastPoint) << endl;
   }
}

//===============================================================================================================================
//! Calculates an average d50 for each polygon. Also fills in 'holes' in active zone and wave calcs i.e. orphan cells which should have been included in the active zone but which have been omitted because of rounding problems
//===============================================================================================================================
void CSimulation::CalcD50AndFillWaveCalcHoles(void)
{
   vector<int> VnPolygonD50Count(m_nGlobalPolygonID + 1, 0);
   vector<double> VdPolygonD50(m_nGlobalPolygonID + 1, 0);

   for (int nX = 0; nX < m_nXGridMax; nX++)
   {
      for (int nY = 0; nY < m_nYGridMax; nY++)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
         {
            // This is a sea cell, first get polygon ID
            int nID = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();

            // Is it in the active zone?
            bool bActive = m_pRasterGrid->m_Cell[nX][nY].bIsInActiveZone();

            if (bActive)
            {
               // It is in the active zone. Does it have unconsolidated sediment on it? Test this using the UnconD50 value: if dGetUnconsD50() returns DBL_NODATA, there is no unconsolidated sediment
               double dTmpd50 = m_pRasterGrid->m_Cell[nX][nY].dGetUnconsD50();
               if (! bFPIsEqual(dTmpd50, DBL_NODATA, TOLERANCE))
               {
                  // It does have unconsolidated sediment, so which polygon is this cell in?
                  if (nID != INT_NODATA)
                  {
                     VnPolygonD50Count[nID]++;
                     VdPolygonD50[nID] += dTmpd50;
                  }
               }
            }

            //             // Now fill in wave calc holes
            //             if (m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight() == 0)
            //             {
            //                if (nID == INT_NODATA)
            //                   m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(m_dAllCellsDeepWaterWaveHeight);
            //             }
            //
            //             if (m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle() == 0)
            //             {
            //                if (nID == INT_NODATA)
            //                   m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(m_dAllCellsDeepWaterWaveAngle);
            //             }

            // Next look at the cell's N-S and W-E neighbours
            int
                nXTmp,
                nYTmp,
                nActive = 0,
                nShadow = 0,
                nShadowNum = 0,
                nDownDrift = 0,
                nDownDriftNum = 0,
                nCoast = 0,
                nRead = 0;
            double
                dWaveHeight = 0,
                dWaveAngle = 0;

            // North
            nXTmp = nX;
            nYTmp = nY - 1;
            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveAngle());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();
                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }

                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();
                  if (nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }

            // East
            nXTmp = nX + 1;
            nYTmp = nY;
            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveAngle());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();
                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }

                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();
                  if (nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }

            // South
            nXTmp = nX;
            nYTmp = nY + 1;
            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveAngle());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();
                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }

                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();
                  if (nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }

            // West
            nXTmp = nX - 1;
            nYTmp = nY;
            if (bIsWithinValidGrid(nXTmp, nYTmp))
            {
               if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInContiguousSea())
               {
                  nRead++;
                  dWaveHeight += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveHeight();
                  dWaveAngle += (m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetWaveAngle());

                  if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsInActiveZone())
                     nActive++;

                  int nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetShadowZoneNumber();
                  if (nTmp != 0)
                  {
                     nShadow++;
                     nShadowNum = nTmp;
                  }

                  nTmp = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetDownDriftZoneNumber();
                  if (nTmp > 0)
                  {
                     nDownDrift++;
                     nDownDriftNum = nTmp;
                  }
               }
               else
                  nCoast++;
            }

            if (nRead > 0)
            {
               // Calculate the average of neighbours
               dWaveHeight /= nRead;
               dWaveAngle /= nRead;
               dWaveAngle = dKeepWithin360(dWaveAngle);

               // If this sea cell has four active-zone neighbours, then it must also be in the active zone: give it wave height and orientation which is the average of its neighbours
               if (nActive == 4)
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetInActiveZone(true);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(dWaveAngle);
               }

               // If this sea cell has a wave height which is the same as its deep-water wave height, but its neighbours have a different average wave height, then give it the average of its neighbours
               double dDeepWaterWaveHeight = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight();
               if ((bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight(), dDeepWaterWaveHeight, TOLERANCE)) && (! bFPIsEqual(dDeepWaterWaveHeight, dWaveHeight, TOLERANCE)))
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
               }

               // If this sea cell has a wave orientation which is the same as its deep-water wave orientation, but its neighbours have a different average wave orientation, then give it the average of its neighbours
               double dDeepWaterWaveAngle = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveAngle();
               if ((bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle(), dDeepWaterWaveAngle, TOLERANCE)) && (! bFPIsEqual(dDeepWaterWaveAngle, dWaveAngle, TOLERANCE)))
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(dWaveAngle);
               }

               // Is this sea cell is not already marked as in a shadow zone (note could be marked as in a shadow zone but not yet processed: a -ve number)?
               int nShadowZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetShadowZoneNumber();
               if (nShadowZoneCode <= 0)
               {
                  // If the cell has four neighbours which are all in a shadow zone, or four neighbours some of which are shadow zone and the remainder downdrift zone, or four neighbours some of which are shadow zone and the remainder coast; then it should also be in the shadow zone: give it the average of its neighbours
                  if ((nShadow == 4) || (nShadow + nDownDrift == 4) || (nShadow + nCoast == 4))
                  {
                     m_pRasterGrid->m_Cell[nX][nY].SetShadowZoneNumber(nShadowNum);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                     m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(dWaveAngle);
                  }
               }

               // If this sea cell is not marked as in a downdrift zone but has four neighbours which are in a downdrift zone, then it should also be in the downdrift zone: give it the average of its neighbours
               int nDownDriftZoneCode = m_pRasterGrid->m_Cell[nX][nY].nGetDownDriftZoneNumber();
               if ((nDownDriftZoneCode == 0) && (nDownDrift == 4))
               {
                  m_pRasterGrid->m_Cell[nX][nY].SetDownDriftZoneNumber(nDownDriftNum);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveHeight(dWaveHeight);
                  m_pRasterGrid->m_Cell[nX][nY].SetWaveAngle(dWaveAngle);
               }
            }
         }
      }
   }

   // Calculate the average d50 for every polygon
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
      {
         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
         int nID = pPolygon->nGetGlobalID();

         if (VnPolygonD50Count[nID] > 0)
            VdPolygonD50[nID] /= VnPolygonD50Count[nID];

         pPolygon->SetAvgUnconsD50(VdPolygonD50[nID]);
      }
   }
}
