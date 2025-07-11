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

#include "simulation.h"
#include "coast.h"

class CRWCoast;      // Forward declaration

//===============================================================================================================================
//! Checks all profiles on all coasts for collisions
//===============================================================================================================================
int CSimulation::nDoMultipleCoastlines(void)
{
   // Check all coastlines
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // Check all profiles
      for (int nProfile = 0; nProfile < m_VCoast[nCoast].nGetNumProfiles(); nProfile++)
      {
         CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
         bool bCoastStart = pProfile->bStartOfCoast();
         bool bCoastEnd = pProfile->bEndOfCoast();

         // Check every cell that is 'under' the profile
         for (int nCell = 0; nCell < pProfile->nGetNumCellsInProfile(); nCell++)
         {
            CGeom2DIPoint* pCell = pProfile->pPtiGetCellInProfile(nCell);
            int nX = pCell->nGetX();
            int nY = pCell->nGetY();

            // Have we hit a cell which is 'under' a coast-normal profile belonging to another coast? NOTE Is a problem if get more than two coast normals passing through this cell
            int nHitProfileCoast = m_pRasterGrid->m_Cell[nX][nY].nGetProfileCoastID();
            if ((nHitProfileCoast != INT_NODATA) && (nHitProfileCoast != nCoast))
            {
               // Yes, we have hit a profile which belongs to a different coast
               int nHitProfile = m_pRasterGrid->m_Cell[nX][nY].nGetProfileID();

               // TEST
               if ((nCoast == 0) && (nProfile == 16))
                  LogStream << endl;

               // Truncate this profile (provided it has not already been truncated)
               int nRtn = nHitProfileTruncate(nCoast, nProfile, nHitProfileCoast, nHitProfile, nX, nY, bCoastStart, bCoastEnd);
               if (nRtn != RTN_OK)
                  return nRtn;

               // And also truncate this profile (provided it has not already been truncated)
               nRtn = nThisProfileTruncate(nCoast, nProfile, nX, nY, bCoastStart, bCoastEnd);
               if (nRtn != RTN_OK)
                  return nRtn;
            }

            // Have we hit a cell which is 'under' another coast?
            if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               // Yes this is a coastline cell, as well as a coast-normal profile cell
               int nHitCoast = m_pRasterGrid->m_Cell[nX][nY].nGetCoastline();
               if (nHitCoast != nCoast)
               {
                  // We have hit a different coastline
                  // TODO
                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                  {
                     LogStream << m_ulIter << ": coast " << nCoast << " profile " << nProfile << " HIT ANOTHER COAST " << nHitCoast << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}." << endl;
                  }
               }
            }
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Truncates a coast-normal profile belonging to another coast, due to collision between normals belonging to different coasts
//===============================================================================================================================
int CSimulation::nHitProfileTruncate(int const nCoast, int const nProfile, int const nHitProfileCoast, int const nHitProfile, int const nX, int const nY, bool const bStartCoastEdgeProfile, bool const bEndCoastEdgeProfile)
{
   if (nHitProfile == INT_NODATA)
   {
      // Should never happen
      return RTN_ERR_CELL_MARKED_PROFILE_COAST_BUT_NOT_PROFILE;
   }

   // OK, get a pointer to the hit profile
   CGeomProfile* pHitProfile = m_VCoast[nHitProfileCoast].pGetProfile(nHitProfile);

   // Has the hit profile already been truncated because of hitting a different coast?
   if (pHitProfile->bTruncatedDifferentCoast())
      // Yes it has, so do not truncate it again
      return RTN_OK;

   // Now get details of the cells in the hit profile
   vector<CGeom2DIPoint>* pVHitProfileCells = pHitProfile->pPtiVGetCellsInProfile();
   int nHitProfileLen = static_cast<int>(pVHitProfileCells->size());

   // Calculate the truncated length of the hit profile
   int nHitProfileNewLen = tMax((nHitProfileLen - GAP_BETWEEN_DIFFERENT_COAST_PROFILES) / 2, 2);

   // Next unmark the cells that will no longer be 'under' the hit profile
   for (int nn = nHitProfileNewLen-1; nn < nHitProfileLen; nn++)
   {
      int nXThis = pVHitProfileCells->at(nn).nGetX();
      int nYThis = pVHitProfileCells->at(nn).nGetY();

      if ((m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileID() == nHitProfile) && (m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileCoastID() == nHitProfileCoast))
         m_pRasterGrid->m_Cell[nXThis][nYThis].SetCoastAndProfileID(INT_NODATA, INT_NODATA);
   }

   // And truncate the hit profile
   pVHitProfileCells->resize(nHitProfileNewLen);
   pHitProfile->SetCellsInProfile(pVHitProfileCells);

   // Set the cells end point (grid CRS)
   CGeom2DIPoint PtiLast = pVHitProfileCells->back();
   pHitProfile->SetEndPoint(&PtiLast);

   // Next get the vector of points in the hit profile
   vector<CGeom2DPoint>* pVHitProfilePoints = pHitProfile->pPtVGetCellsInProfileExtCRS();

   // Truncate the vector of points
   pVHitProfilePoints->resize(nHitProfileNewLen);
   pHitProfile->SetPointsInProfile(pVHitProfilePoints);

   // And flag as truncated
   // pHitProfile->SetTruncatedDifferentCoast(true);

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

      LogStream << m_ulIter << ": coast " << nCoast << strTmp << " profile " << nProfile << " hit profile " << nHitProfile << " belonging to another coast " << nHitProfileCoast << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}. Hit profile truncated, length of profile " << nHitProfile << " was " << nHitProfileLen << " cells, is now " << nHitProfileNewLen << " cells." << endl;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Truncates a coast-normal profile, due to collision between normals belonging to different coasts
//===============================================================================================================================
int CSimulation::nThisProfileTruncate(int const nCoast, int const nProfile, int const nX, int const nY, bool const bStartCoastEdgeProfile, bool const bEndCoastEdgeProfile)
{
   if (nCoast == INT_NODATA)
   {
      // Should never happen
      return RTN_ERR_CELL_MARKED_PROFILE_COAST_BUT_NOT_PROFILE;
   }

   // OK, get a pointer to the profile
   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // Has the hit profile already been truncated because of hitting a different coast?
   if (pProfile->bTruncatedDifferentCoast())
      // Yes it has so do not truncate it again
      return RTN_OK;

   // Now get details of the cells in this profile
   vector<CGeom2DIPoint>* pVProfileCells = pProfile->pPtiVGetCellsInProfile();
   int nProfileLen = static_cast<int>(pVProfileCells->size());

   // Calculate the truncated length of the profile
   int nProfileNewLen = tMax((nProfileLen - GAP_BETWEEN_DIFFERENT_COAST_PROFILES) / 2, 2);

   // Next unmark the cells that will no longer be 'under' the hit profile
   for (int nn = nProfileNewLen-1; nn < nProfileLen; nn++)
   {
      int nXThis = pVProfileCells->at(nn).nGetX();
      int nYThis = pVProfileCells->at(nn).nGetY();

      if ((m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileID() == nProfile) && (m_pRasterGrid->m_Cell[nXThis][nYThis].nGetProfileCoastID() == nCoast))
         m_pRasterGrid->m_Cell[nXThis][nYThis].SetCoastAndProfileID(INT_NODATA, INT_NODATA);
   }

   // And truncate the profile
   pVProfileCells->resize(nProfileNewLen);
   pProfile->SetCellsInProfile(pVProfileCells);

   // Set the cells end point (grid CRS)
   CGeom2DIPoint PtiLast = pVProfileCells->back();
   pProfile->SetEndPoint(&PtiLast);

   // Next get the vector of points in the profile
   vector<CGeom2DPoint>* pVProfilePoints = pProfile->pPtVGetCellsInProfileExtCRS();

   // Truncate the vector of points
   pVProfilePoints->resize(nProfileNewLen);
   pProfile->SetPointsInProfile(pVProfilePoints);

   // And flag as truncated
   // pProfile->SetTruncatedDifferentCoast(true);

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

   return RTN_OK;
}
