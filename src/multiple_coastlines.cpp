/*!

   \file multiple_coastlines.cpp
   \brief Routines relating to ukltiple coastlines
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License

*/

/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include <assert.h>

#include <iostream>
using std::endl;

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "2di_point.h"
#include "multi_line.h"

class CRWCoast;      // Forward declaration

//===============================================================================================================================
//! Checks all profiles on all coasts for intersections between profiles belonging to different coasts
//===============================================================================================================================
int CSimulation::nDoMultipleCoastlines(void)
{
   // Check all coastlines
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // Check all profiles
      for (int n = 0; n < m_VCoast[nCoast].nGetNumProfiles(); n++)
      {
         CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfileWithDownCoastSeq(n);
         int const nProfile = pProfile->nGetProfileID();
         bool const bCoastStart = pProfile->bStartOfCoast();
         bool const bCoastEnd = pProfile->bEndOfCoast();

         // Check every cell that is 'under' this profile
         for (int nCell = 0; nCell < pProfile->nGetNumCellsInProfile(); nCell++)
         {
            CGeom2DIPoint const* pCell = pProfile->pPtiGetCellInProfile(nCell);
            int const nX = pCell->nGetX();
            int const nY = pCell->nGetY();

            // Have we hit a cell which is 'under' a coast-normal profile belonging to another coast? NOTE Is a problem if get more than two coast normals passing through this cell
            int nHitProfileCoast = m_pRasterGrid->m_Cell[nX][nY].nGetProfileCoastID();
            if ((nHitProfileCoast != INT_NODATA) && (nHitProfileCoast != nCoast))
            {
               // Yes, we have hit a profile which belongs to a different coast
               int const nHitProfile = m_pRasterGrid->m_Cell[nX][nY].nGetProfileID();

               // Truncate both profiles
               int const nRtn = nTruncateProfilesDifferentCoasts(nCoast, nProfile, nCell, nHitProfileCoast, nHitProfile, nX, nY, bCoastStart, bCoastEnd);
               if (nRtn != RTN_OK)
                  return nRtn;
            }
            else
            {
               // Now try an adjacent point, does not matter wehich one NOTE Is a problem if get more than two coast normals passing through this cell
               nHitProfileCoast = m_pRasterGrid->m_Cell[nX][nY+1].nGetProfileCoastID();
               if ((nHitProfileCoast != INT_NODATA) && (nHitProfileCoast != nCoast))
               {
                  // Yes, we have hit a profile which belongs to a different coast
                  int const nHitProfile = m_pRasterGrid->m_Cell[nX][nY+1].nGetProfileID();

                  // Truncate both profiles
                  int const nRtn = nTruncateProfilesDifferentCoasts(nCoast, nProfile, nCell, nHitProfileCoast, nHitProfile, nX, nY+1, bCoastStart, bCoastEnd);
                  if (nRtn != RTN_OK)
                     return nRtn;
               }
            }

            // Have we hit a cell which is 'under' another coast?
            if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               // Yes this is a coastline cell, as well as a coast-normal profile cell
               int const nHitCoast = m_pRasterGrid->m_Cell[nX][nY].nGetCoastline();
               if (nHitCoast != nCoast)
               {
                  // We have hit a different coastline, so truncate this profile
                  int const nRtn = nTruncateProfileHitDifferentCoast(nCoast, nProfile, nCell, nX, nY, bCoastStart, bCoastEnd);
                  if (nRtn != RTN_OK)
                     return nRtn;
               }
            }
            else
            {
               // We also need to check an adjacent cell, doesn't matter which one
               if (m_pRasterGrid->m_Cell[nX][nY+1].bIsCoastline())
               {
                  // Yes this is a coastline cell, as well as a coast-normal profile cell
                  int const nHitCoast = m_pRasterGrid->m_Cell[nX][nY+1].nGetCoastline();
                  if (nHitCoast != nCoast)
                  {
                  // We have hit a different coastline, so truncate this profile
                  int const nRtn = nTruncateProfileHitDifferentCoast(nCoast, nProfile, nCell, nX, nY, bCoastStart, bCoastEnd);
                  if (nRtn != RTN_OK)
                     return nRtn;
                  }
               }
            }
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Truncates two intersecting coast-normal profile belonging to different coasts
//===============================================================================================================================
int CSimulation::nTruncateProfilesDifferentCoasts(int const nCoast, int const nProfile, int nCell, int const nHitProfileCoast, int const nHitProfile, int nX, int nY, bool const bStartCoastEdgeProfile, bool const bEndCoastEdgeProfile)
{
   if ((nProfile == INT_NODATA) || (nHitProfile == INT_NODATA))
   {
      // Should never happen
      return RTN_ERR_CELL_MARKED_PROFILE_COAST_BUT_NOT_PROFILE;
   }

   // OK, get pointers to 'this' profile and to the hit profile
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
   CGeomProfile* pHitProfile = m_VCoast[nHitProfileCoast].pGetProfile(nHitProfile);

   // Get details of the cells in this profile
   vector<CGeom2DIPoint>* pVProfileCells = pProfile->pPtiVGetCellsInProfile();
   int const nProfileLen = static_cast<int>(pVProfileCells->size());

   // Are either of the profiles grid-edge profiles?
   bool bProfileGridEdge = false;
   bool bHitProfileGridEdge = false;
   if (pProfile->bStartOfCoast() || pProfile->bEndOfCoast())
      bProfileGridEdge = true;
   if (pHitProfile->bStartOfCoast() || pHitProfile->bEndOfCoast())
      bHitProfileGridEdge = true;

   // Are both profiles grid-edge profiles?
   if (bProfileGridEdge && bHitProfileGridEdge)
   {
      // Yes, so we need to treat these profiles differently. Get the start point of each profile
      CGeom2DIPoint const* pPtiProfileStart = pProfile->pPtiGetStartPoint();
      CGeom2DIPoint const* pPtiHitProfileStart = pHitProfile->pPtiGetStartPoint();

      // And get the distance between these
      double const dDist = dGetDistanceBetween(pPtiProfileStart, pPtiHitProfileStart);
      nCell = static_cast<int>(dDist - GAP_BETWEEN_DIFFERENT_COAST_PROFILES);

      // Safety check
      nCell = tMax(nCell, 0);

      // Finally get the grid CRS location of the new profile endpoint
      nX = pVProfileCells->at(nCell).nGetX();
      nY = pVProfileCells->at(nCell).nGetY();
   }

   // Calculate the truncated length of the list of cells in 'this' profile
   int const nProfileNewLen = tMax(nCell - GAP_BETWEEN_DIFFERENT_COAST_PROFILES, MIN_PROFILE_SIZE);

   // Next unmark the cells that will no longer be 'under' the profile
   for (int nn = nProfileNewLen-1; nn < nProfileLen; nn++)
   {
      int const nXThis = pVProfileCells->at(nn).nGetX();
      int const nYThis = pVProfileCells->at(nn).nGetY();

      if ((m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileID() == nProfile) && (m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileCoastID() == nCoast))
         m_pRasterGrid->m_Cell[nXThis][nYThis].SetCoastAndProfileID(INT_NODATA, INT_NODATA);
   }

   // And truncate the list of cells in the profile
   pVProfileCells->resize(nProfileNewLen);
   pProfile->SetCellsInProfile(pVProfileCells);

   // Set the profile's end point (grid CRS)
   CGeom2DIPoint PtiLast = pVProfileCells->back();
   pProfile->SetEndPoint(&PtiLast);

   // Now truncate the profile's CGeomMultiLine (external CRS)
   int nRtn = nTruncateProfileMultiLineDifferentCoasts(pProfile, nX, nY);
   if (nRtn != RTN_OK)
      return nRtn;

   // And flag as truncated
   pProfile->SetTruncatedDifferentCoast(true);

   if (m_nLogFileDetail >= LOG_FILE_ALL)
   {
      string strTmp;
      if (bStartCoastEdgeProfile || bEndCoastEdgeProfile)
      {
         strTmp += " grid-edge ";

         if (bStartCoastEdgeProfile)
            strTmp += "start";

         if (bEndCoastEdgeProfile)
            strTmp += "end";
      }

      LogStream << m_ulIter << ": coast " << nCoast << strTmp << " profile " << nProfile << " hit profile belonging to another coast at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}. Profile truncated, length of profile " << nProfile << " was " << nProfileLen << " cells, is now " << nProfileNewLen << " cells." << endl;
   }

   // Get details of the cells in the hit profile
   vector<CGeom2DIPoint>* pVHitProfileCells = pHitProfile->pPtiVGetCellsInProfile();
   int const nHitProfileLen = static_cast<int>(pVHitProfileCells->size());

   int const nHitCell = pHitProfile->nGetIndexOfCellInProfile(nX, nY);
   if (nHitCell == INT_NODATA)
      return RTN_ERR_CELL_NOT_FOUND_IN_HIT_PROFILE_DIFFERENT_COASTS;

   // Calculate the truncated length of the list of cells in the hit profile
   int const nHitProfileNewLen = tMax(nHitCell - GAP_BETWEEN_DIFFERENT_COAST_PROFILES, MIN_PROFILE_SIZE);

   // Next unmark the cells that will no longer be 'under' the hit profile
   for (int nn = nHitProfileNewLen-1; nn < nHitProfileLen; nn++)
   {
      int const nXThis = pVHitProfileCells->at(nn).nGetX();
      int const nYThis = pVHitProfileCells->at(nn).nGetY();

      if ((m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileID() == nHitProfile) && (m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileCoastID() == nHitProfileCoast))
         m_pRasterGrid->m_Cell[nXThis][nYThis].SetCoastAndProfileID(INT_NODATA, INT_NODATA);
   }

   // And truncate the hit profile
   pVHitProfileCells->resize(nHitProfileNewLen);
   pHitProfile->SetCellsInProfile(pVHitProfileCells);

   // Set the hit profile's end point (grid CRS)
   PtiLast = pVHitProfileCells->back();
   pHitProfile->SetEndPoint(&PtiLast);

   // Now truncate the profile's CGeomMultiLine (external CRS)
   nRtn = nTruncateProfileMultiLineDifferentCoasts(pProfile, nX, nY);
   if (nRtn != RTN_OK)
      return nRtn;

   // And flag as truncated
   pHitProfile->SetTruncatedDifferentCoast(true);

   if (m_nLogFileDetail >= LOG_FILE_ALL)
   {
      string strTmp;
      if (bStartCoastEdgeProfile || bEndCoastEdgeProfile)
      {
         strTmp += " grid-edge ";

         if (bStartCoastEdgeProfile)
            strTmp += "start";

         if (bEndCoastEdgeProfile)
            strTmp += "end";
      }

      LogStream << m_ulIter << ": coast " << nHitProfileCoast << strTmp << " profile " << nHitProfile << " also truncated, length of profile " << nHitProfile << " was " << nHitProfileLen << " cells, is now " << nHitProfileNewLen << " cells." << endl;
   }

   // // DEBUG CODE ================
   // m_nGISSave++;
   // if (! bWriteVectorGISFile(VECTOR_PLOT_COAST, &VECTOR_PLOT_COAST_TITLE))
   //    return false;
   // if (! bWriteVectorGISFile(VECTOR_PLOT_NORMALS, &VECTOR_PLOT_NORMALS_TITLE))
   //    return false;
   // if (! bWriteVectorGISFile(VECTOR_PLOT_INVALID_NORMALS, &VECTOR_PLOT_INVALID_NORMALS_TITLE))
   //    return false;
   // if (! bWriteRasterGISFile(RASTER_PLOT_NORMAL_PROFILE, &RASTER_PLOT_NORMAL_PROFILE_TITLE))
   //    return false;
   // if (! bWriteRasterGISFile(RASTER_PLOT_POLYGON, &RASTER_PLOT_POLYGON_TITLE))
   //    return false;
   // // DEBUG CODE ================

   return RTN_OK;
}

//===============================================================================================================================
//! Truncates a profile which has hit a different coast
//===============================================================================================================================
int CSimulation::nTruncateProfileHitDifferentCoast(int const nCoast, int const nProfile, int const nCell, int const nX, int const nY, bool const bStartCoastEdgeProfile, bool const bEndCoastEdgeProfile)
{
   // OK, get a pointer to 'this' profile
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // Now get details of the cells in this profile
   vector<CGeom2DIPoint>* pVProfileCells = pProfile->pPtiVGetCellsInProfile();
   int nProfileLen = static_cast<int>(pVProfileCells->size());

   // Calculate the truncated length of the list of cells in 'this' profile
   int nProfileNewLen = tMax(nCell - GAP_BETWEEN_DIFFERENT_COAST_PROFILES, MIN_PROFILE_SIZE);

   // Get the lengths of the adjacent, this-coast, profiles
   CGeomProfile const* pUpCoastProfile = pProfile->pGetUpCoastAdjacentProfile();
   CGeomProfile const* pDownCoastProfile = pProfile->pGetDownCoastAdjacentProfile();
   int nUpCoastProfileLen = INT_NODATA;
   int nDownCoastProfileLen = INT_NODATA;
   if (pUpCoastProfile != NULL)
      nUpCoastProfileLen = pUpCoastProfile->nGetNumCellsInProfile();
   if (pDownCoastProfile != NULL)
      nDownCoastProfileLen = pDownCoastProfile->nGetNumCellsInProfile();

   // And calculate the average length of adjacent profiles
   int nAvgAdjacentProfileLen;
   if (pUpCoastProfile == NULL)
      nAvgAdjacentProfileLen = nDownCoastProfileLen;
   else if (pDownCoastProfile == NULL)
      nAvgAdjacentProfileLen = nUpCoastProfileLen;
   else
      nAvgAdjacentProfileLen = (nUpCoastProfileLen + nDownCoastProfileLen) / 2;

   // Use the average length of adjacent profiles to further truncated this profile, if necessary
   nProfileNewLen = tMin(nProfileNewLen, nAvgAdjacentProfileLen);

   // We have the new profile length so next unmark the cells that will no longer be 'under' the profile
   for (int nn = nProfileNewLen-1; nn < nProfileLen; nn++)
   {
      int const nXThis = pVProfileCells->at(nn).nGetX();
      int const nYThis = pVProfileCells->at(nn).nGetY();

      if ((m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileID() == nProfile) && (m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileCoastID() == nCoast))
         m_pRasterGrid->m_Cell[nXThis][nYThis].SetCoastAndProfileID(INT_NODATA, INT_NODATA);
   }

   // And truncate the list of cells in the profile
   pVProfileCells->resize(nProfileNewLen);
   pProfile->SetCellsInProfile(pVProfileCells);

   // Set the profile's end point (grid CRS)
   CGeom2DIPoint const PtiLast = pVProfileCells->back();
   pProfile->SetEndPoint(&PtiLast);

   // Now truncate the profile's CGeomMultiLine (external CRS)
   int const nRtn = nTruncateProfileMultiLineDifferentCoasts(pProfile, nX, nY);
   if (nRtn != RTN_OK)
      return nRtn;

   // And flag as truncated
   pProfile->SetTruncatedDifferentCoast(true);

   if (m_nLogFileDetail >= LOG_FILE_ALL)
   {
      string strTmp;
      if (bStartCoastEdgeProfile || bEndCoastEdgeProfile)
      {
         strTmp += " grid-edge ";

         if (bStartCoastEdgeProfile)
            strTmp += "start";

         if (bEndCoastEdgeProfile)
            strTmp += "end";
      }

      LogStream << m_ulIter << ": coast " << nCoast << strTmp << " profile " << nProfile << " hit another coast at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}. Profile truncated, length of profile " << nProfile << " was " << nProfileLen << " cells, is now " << nProfileNewLen << " cells." << endl;
   }


   // DEBUG CODE ================
   m_nGISSave++;
   if (! bWriteVectorGISFile(VECTOR_PLOT_COAST, &VECTOR_PLOT_COAST_TITLE))
      return false;
   if (! bWriteVectorGISFile(VECTOR_PLOT_NORMALS, &VECTOR_PLOT_NORMALS_TITLE))
      return false;
   if (! bWriteVectorGISFile(VECTOR_PLOT_INVALID_NORMALS, &VECTOR_PLOT_INVALID_NORMALS_TITLE))
      return false;
   if (! bWriteRasterGISFile(RASTER_PLOT_NORMAL_PROFILE, &RASTER_PLOT_NORMAL_PROFILE_TITLE))
      return false;
   if (! bWriteRasterGISFile(RASTER_PLOT_POLYGON, &RASTER_PLOT_POLYGON_TITLE))
      return false;
   // DEBUG CODE ================

   return RTN_OK;
}

//===============================================================================================================================
//! Truncates the CGeomMultiLine (external CRS) of a profile which has hit a different coast, or a cost-normal profile belonging to a different coast
//===============================================================================================================================
int CSimulation::nTruncateProfileMultiLineDifferentCoasts(CGeomProfile* pProfile, int const nX, int const nY)
{
   // Find which multiline line segment 'contains' this extCRS end point, then truncate at this point
   bool bFound = false;

   int const nProfileLineSegments = pProfile->CGeomMultiLine::nGetNumLineSegments();
   for (int nSeg = 0; nSeg < nProfileLineSegments; nSeg++)
   {
      // Search each segment
      int const nNumCoinc = pProfile->CGeomMultiLine::nGetNumCoincidentProfilesInLineSegment(nSeg);
      for (int nCoinc = 0; nCoinc < nNumCoinc; nCoinc++)
      {
         vector<CGeom2DPoint>& pVPt = pProfile->CGeomMultiLine::pGetPoints();

         for (int nLin = 0; nLin < static_cast<int>(pVPt.size())-1; nLin++)
         {
            double const dX1 = pVPt[nLin].dGetX();
            double const dY1 = pVPt[nLin].dGetY();
            double const dX2 = pVPt[nLin+1].dGetX();
            double const dY2 = pVPt[nLin+1].dGetY();

            double const dXMin = tMin(dX1, dX2);
            double const dXMax = tMax(dX1, dX2);
            double const dYMin = tMin(dY1, dY2);
            double const dYMax = tMax(dY1, dY2);

            double const dX = dGridXToExtCRSX(nX);
            double const dY = dGridYToExtCRSY(nY);

            if ((dX >= dXMin) && (dX <= dXMax) && (dY >= dYMin) && (dY <= dYMax))
            {
               bFound = true;

               pProfile->TruncateProfile(nSeg+1);
               pVPt.push_back(CGeom2DPoint(dX, dY));
               pProfile->CGeomMultiLine::SetPoints(pVPt);

               break;
            }
         }
         if (bFound)
             break;
      }
      if (bFound)
         break;
   }

   if (! bFound)
      return RTN_ERR_POINT_NOT_FOUND_IN_MULTILINE_DIFFERENT_COASTS;

   return RTN_OK;
}
