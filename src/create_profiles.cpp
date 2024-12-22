/*!
 *
 * \file create_profiles.cpp
 * \brief Creates profiles which are approximately normal to the coastline, these will become inter-polygon boundaries
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
#include <cfloat>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::setiosflags;
using std::setprecision;

#include <algorithm>
using std::find;
using std::sort;

#include <utility>
using std::make_pair;
using std::pair;

#include "cme.h"
#include "simulation.h"
#include "coast.h"

//===============================================================================================================================
//! Function used to sort coastline curvature values when locating start points of normal profiles
//===============================================================================================================================
bool bCurvaturePairCompareDescending(const pair<int, double> &prLeft, const pair<int, double> &prRight)
{
   // Sort in descending order (i.e. most concave first)
   return prLeft.second > prRight.second;
}

//===============================================================================================================================
//! Create coastline-normal profiles for all coastlines at locations of greatest concave curvature of the vector coastline. Finally, grid-edge normal profiles are created
//===============================================================================================================================
int CSimulation::nCreateAllProfiles(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Creating profiles" << endl;
   
   for (unsigned int nCoast = 0; nCoast < m_VCoast.size(); nCoast++)
   {
      int nProfile = 0;
      int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

      // Create a bool vector to mark coast points which have been searched
      vector<bool> bVCoastPointDone(nCoastSize, false);

      // Now create a vector of pairs: the first value of the pair is the coastline point, the second is the coastline's curvature at that point
      vector<pair<int, double>> prVCurvature;
      for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
      {
         double dCurvature;
         if (m_VCoast[nCoast].pGetCoastLandform(nCoastPoint)->nGetLandFormCategory() != LF_CAT_INTERVENTION)
         {
            // Not an intervention coast point, so just store the smoothed curvature
            dCurvature = m_VCoast[nCoast].dGetSmoothCurvature(nCoastPoint);
         }
         else
         {
            // This is an intervention coast point, which is likely to have some sharp angles. So store the average of the smooth and detailed curvature
            dCurvature = (m_VCoast[nCoast].dGetSmoothCurvature(nCoastPoint) + m_VCoast[nCoast].dGetDetailedCurvature(nCoastPoint)) / 2;
         }
         prVCurvature.push_back(make_pair(nCoastPoint, dCurvature));
      }

      // Sort this pair vector in descending order, so that the most concave smoothed-curvature points are first
      sort(prVCurvature.begin(), prVCurvature.end(), bCurvaturePairCompareDescending);

      // // DEBUG CODE =======================================================
      // for (int n = 0; n < prVCurvature.size(); n++)
      // {
      //    LogStream << prVCurvature[n].first << "\t" << prVCurvature[n].second << endl;
      // }
      // LogStream << endl << endl;
      // // DEBUG CODE =======================================================

      // And mark points at and near the start and end of the coastline so that they don't get searched (will be creating 'special' start- and end-of-coast profiles at these end points later)
      for (int n = 0; n < m_nCoastNormalAvgSpacing; n++)
      {
         if (n < nCoastSize)
            bVCoastPointDone[n] = true;

         int m = nCoastSize - n - 1;
         if (m >= 0)
            bVCoastPointDone[m] = true;
      }

      // Now locate all normal profiles at points of maximum smoothed convexity
      LocateAllProfiles(nCoast, nProfile, m_nCoastNormalAvgSpacing, &bVCoastPointDone, &prVCurvature);

      // Did we fail to create any normal profiles? If so, quit
      if (nProfile < 0)
      {
         string strErr = ERR + "timestep " + strDblToStr(m_ulIter) + ": could not create profiles for coastline " + strDblToStr(nCoast);
         if (m_ulIter == 1)
            strErr += ". Check the SWL";
         strErr += "\n";

         cerr << strErr;
         LogStream << strErr;

         return RTN_ERR_NO_PROFILES_1;
      }

      // Locate a 'special' profile at the beginning of the coastline, and put this onto the raster grid
      int nRet = nLocateGridEdgeProfile(true, nCoast, nProfile);
      if (nRet != RTN_OK)
         return nRet;

      // Locate a 'special' profile at the end of the coastline, and put this onto the raster grid
      nRet = nLocateGridEdgeProfile(false, nCoast, ++nProfile);
      if (nRet != RTN_OK)
         return nRet;

      // Create an index to the profiles in along-coast sequence
      m_VCoast[nCoast].CreateAlongCoastProfileIndex();

      // // DEBUG CODE =======================================================
      // int nProf = m_VCoast[nCoast].nGetNumProfiles();
      // for (int n = 0; n < nProf; n++)
      // {
      //    LogStream << n << "\t" << m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(n);
      //
      //    if (m_VCoast[nCoast].pGetCoastLandform(n)->nGetLandFormCategory() == LF_CAT_INTERVENTION)
      //       LogStream << "\tI";
      //
      //    LogStream << endl;
      // }
      // // DEBUG CODE =======================================================

      // // DEBUG CODE =======================================================
      // for (int n = 0; n < nCoastSize; n++)
      // {
      //    LogStream << n << "\t";
      //
      //    LogStream << m_VCoast[nCoast].dGetDetailedCurvature(n) << "\t";
      //
      //    LogStream << m_VCoast[nCoast].dGetSmoothCurvature(n) << "\t";
      //
      //    if (m_VCoast[nCoast].pGetCoastLandform(n)->nGetLandFormCategory() == LF_CAT_INTERVENTION)
      //       LogStream << "I\t";
      //
      //    if (m_VCoast[nCoast].bIsProfileStartPoint(n))
      //       LogStream << "profile " << m_VCoast[nCoast].nGetProfileNumber(n);
      //
      //    LogStream << endl;
      // }
      // LogStream << endl;
      // // DEBUG CODE =======================================================
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Create profiles normal to the coastline, modify these if they intersect, then put valid profiles onto the raster grid
//===============================================================================================================================
int CSimulation::nCreateAllProfilesAndCheckValidity(void)
{
   // Create all coastline-normal profiles, in coastline-concave-curvature sequence i.e. the first profiles are created 'around' the most concave bits of coast. An index is also created which allows profiles to be accessed in along-coast sequence. Also create 'special' profiles at the start and end of the coast, and put these onto the raster griid now
   int nRet = nCreateAllProfiles();
   if (nRet != RTN_OK)
      return nRet;

   // Check to see which coastline-normal profiles intersect. Then modify intersecting profiles so that the sections of each profile seaward of the point of intersection are 'shared' i.e. are multi-lines. This creates the boundaries of the triangular polygons
   nRet = nModifyAllIntersectingProfiles();
   if (nRet != RTN_OK)
      return nRet;

   // Again check the normal profiles for insufficient length: is the water depth at the end point less than the depth of closure? We do this again because some profiles may have been shortened as a result of intersection. Do once for every coastline object
   for (unsigned int nCoast = 0; nCoast < m_VCoast.size(); nCoast++)
   {
      for (int nProf = 0; nProf < m_VCoast[nCoast].nGetNumProfiles(); nProf++)
      {
         CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProf);

         if (pProfile->bProfileOK())
         {
            int nSize = pProfile->nGetProfileSize();

            // Safety check
            if (nSize == 0)
            {
               pProfile->SetTooShort(true);
               continue;
            }

            CGeom2DPoint const* pPtEnd = pProfile->pPtGetPointInProfile(nSize - 1);
            CGeom2DIPoint PtiEnd = PtiExtCRSToGrid(pPtEnd);
            int nXEnd = PtiEnd.nGetX();
            int nYEnd = PtiEnd.nGetY();

            // Safety check: the point may be outside the grid because of rounding, so keep it within the grid
            nXEnd = tMax(nXEnd, 0);
            nXEnd = tMin(nXEnd, m_nXGridMax - 1);
            nYEnd = tMax(nYEnd, 0);
            nYEnd = tMin(nYEnd, m_nYGridMax - 1);

            if (m_pRasterGrid->m_Cell[nXEnd][nYEnd].dGetSeaDepth() < m_dDepthOfClosure)
            {
               // if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               //    LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProf << " is invalid, is too short for depth of closure " << m_dDepthOfClosure << " at end point [" << nXEnd << "][" << nYEnd << "] = {" << pPtEnd->dGetX() << ", " << pPtEnd->dGetY() << "}, flagging as too short" << endl;

               pProfile->SetTooShort(true);
            }
         }
      }
   }

   // Put all valid coastline-normal profiles (apart from the profiles at the start and end of the coast, since they have already been done) onto the raster grid. But if the profile is not long enough, or crosses a coastline or hits dry land, then mark the profile as invalid
   nRet = nPutAllProfilesOntoGrid();
   if (nRet != RTN_OK)
      return nRet;

   return RTN_OK;
}

//===============================================================================================================================
//! Locate the normal profiles on the coastline
//===============================================================================================================================
void CSimulation::LocateAllProfiles(int const nCoast, int &nProfile, int const nProfileHalfAvgSpacing, vector<bool>* pbVCoastPointDone, vector<pair<int, double>> const* prVCurvature)
{
   int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   // Work along the vector of curvature pairs starting at the convex end
   for (int n = nCoastSize - 1; n >= 0; n--)
   {
      // Have we searched all the coastline points?
      int nStillToSearch = 0;
      for (int m = 0; m < nCoastSize; m++)
         if (! pbVCoastPointDone->at(m))
            nStillToSearch++;

      if (nStillToSearch == 0)
         // Nope, we are done here
         return;

      // This convex point on the coastline is a potential location for a normal
      int nNormalPoint = prVCurvature->at(n).first;

      // Ignore each end of the coastline
      if ((nNormalPoint == 0) || (nNormalPoint == nCoastSize - 1))
         continue;

      if (! pbVCoastPointDone->at(nNormalPoint))
      {
         // We have not already searched this coast point. Is it an intervention coast point?
         bool bIntervention = false;
         if (m_VCoast[nCoast].pGetCoastLandform(nNormalPoint)->nGetLandFormCategory() == LF_CAT_INTERVENTION)
            bIntervention = true;

         // Try locating a profile here
         int nRet = nLocateProfile(nCoast, nNormalPoint, nProfile, bIntervention);

         // Mark this coast point as searched
         pbVCoastPointDone->at(nNormalPoint) = true;

         if (nRet != RTN_OK)
            // This potential profile is no good (has hit coast, or hit dry land, etc.) so forget about it
            continue;

         // This profile is fine
         // CGeom2DIPoint PtiThis = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nNormalPoint);
         // CGeom2DPoint PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nNormalPoint);

         // if (m_nLogFileDetail >= LOG_FILE_ALL)
         //    LogStream << m_ulIter << ": Coast " << nCoast << " profile " << nProfile << " profile created at coast point " << nNormalPoint << " [" << PtiThis.nGetX() << "][" << PtiThis.nGetY() << "] = {" << PtThis.dGetX() << ", " << PtThis.dGetY() << "} (smoothed curvature = " << m_VCoast[nCoast].dGetSmoothCurvature(nNormalPoint) << ", detailed curvature = " << m_VCoast[nCoast].dGetDetailedCurvature(nNormalPoint) << ")" << endl;

         nProfile++;

         // We need to mark points on either side of this profile so that we don't get profiles which are too close together. This is straightforward for non-intervention profiles, but best-placed profiles on narrow intervention structures may need to be quite close. So guess in a value for number of points to mark on interventions TODO 011
         double dNumToMark = nProfileHalfAvgSpacing;
         if (bIntervention)
            dNumToMark = 2;

         // Mark points on either side of the profile
         for (int m = 1; m < dNumToMark; m++)
         {
            int nTmpPoint = nNormalPoint + m;
            if (nTmpPoint < nCoastSize)
               pbVCoastPointDone->at(nTmpPoint)= true;

            nTmpPoint = nNormalPoint - m;
            if (nTmpPoint >= 0)
               pbVCoastPointDone->at(nTmpPoint) = true;
         }
      }
   }
}

//===============================================================================================================================
//! Locates a single coastline-normal profile (which may be an intervention profile)
//===============================================================================================================================
int CSimulation::nLocateProfile(int const nCoast, int const nProfileStartPoint, int &nProfile, bool const bIntervention)
{
   // OK, we have flagged the start point of this new coastline-normal profile, so create it. Make the start of the profile the centroid of the actual cell that is marked as coast (not the cell under the smoothed vector coast, they may well be different)
   int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();

   CGeom2DPoint PtStart; // In external CRS
   PtStart.SetX(dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nProfileStartPoint)->nGetX()));
   PtStart.SetY(dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nProfileStartPoint)->nGetY()));

   CGeom2DPoint PtEnd;   // In external CRS
   CGeom2DIPoint PtiEnd; // In grid CRS
   int nRet = nGetCoastNormalEndPoint(nCoast, nProfileStartPoint, nCoastSize, &PtStart, m_dCoastNormalLength, &PtEnd, &PtiEnd, bIntervention);
   if (nRet == RTN_ERR_NO_SOLUTION_FOR_ENDPOINT)
   {
      // Could not solve end-point equation, so forget about this profile
      return nRet;      
   }

   int nXEnd = PtiEnd.nGetX();
   int nYEnd = PtiEnd.nGetY();

   // Safety check: is the end point in the contiguous sea?
   if (! m_pRasterGrid->m_Cell[nXEnd][nYEnd].bIsInContiguousSea())
   {
      // if (m_nLogFileDetail >= LOG_FILE_ALL)
      //    LogStream << m_ulIter << ": coast " << nCoast << ", possible profile with start point " << nProfileStartPoint << " has inland end point at [" << nXEnd << "][" << nYEnd << "] = {" << dGridCentroidXToExtCRSX(nXEnd) << ", " << dGridCentroidYToExtCRSY(nYEnd) << "}, ignoring" << endl;

      return RTN_ERR_PROFILE_ENDPOINT_IS_INLAND;
   }

   // Safety check: is the water depth at the end point less than the depth of closure?
   if (m_pRasterGrid->m_Cell[nXEnd][nYEnd].dGetSeaDepth() < m_dDepthOfClosure)
   {
      // if (m_nLogFileDetail >= LOG_FILE_ALL)
      //    LogStream << m_ulIter << ": coast " << nCoast << ", possible profile with start point " << nProfileStartPoint << " is too short for depth of closure " << m_dDepthOfClosure << " at end point [" << nXEnd << "][" << nYEnd << "] = {" << dGridCentroidXToExtCRSX(nXEnd) << ", " << dGridCentroidYToExtCRSY(nYEnd) << "}, ignoring" << endl;

      return RTN_ERR_PROFILE_END_INSUFFICIENT_DEPTH;
   }

   // No problems, so create the new profile
   m_VCoast[nCoast].AppendProfile(nProfileStartPoint, nProfile);

   // And create the profile's coastline-normal vector (start and end points are in external CRS)
   vector<CGeom2DPoint> VNormal;
   VNormal.push_back(PtStart);
   VNormal.push_back(PtEnd);

   CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);
   pProfile->SetAllPointsInProfile(&VNormal); // Set only two points: the start and end points

   // Create the profile's CGeomMultiLine then set nProfile as the only co-incident profile of the only line segment
   pProfile->AppendLineSegment();
   pProfile->AppendCoincidentProfileToLineSegments(make_pair(nProfile, 0));

   // LogStream << m_ulIter << ": in nLocateProfile() coast " << nCoast << " profile " << nProfile << " created at coast point (" << nProfileStartPoint << ") from [" << PtStart.dGetX() << "][" << PtStart.dGetY() << "] to [" << PtEnd.dGetX() << "][" << PtEnd.dGetY() << "]" << endl;

   return RTN_OK;
}

//===============================================================================================================================
//! Creates a 'special' profile at each end of a coastline, at the edge of the raster grid. This profile is not necessarily normal to the coastline since it goes along the grid's edge
//===============================================================================================================================
int CSimulation::nLocateGridEdgeProfile(bool const bCoastStart, int const nCoast, int &nProfile)
{
   int nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
   int nHandedness = m_VCoast[nCoast].nGetSeaHandedness();
   int nProfileLen = nRound(m_dCoastNormalLength / m_dCellSide);  // Profile length in grid CRS
   int nProfileStartEdge;

   CGeom2DIPoint PtiProfileStart; // In grid CRS
   CGeomProfile* pProfile = NULL;
   vector<CGeom2DIPoint> VPtiNormalPoints; // In grid CRS

   if (bCoastStart)
   {
      // At start of coast
      PtiProfileStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(0); // Grid CRS
      nProfileStartEdge = m_VCoast[nCoast].nGetStartEdge();
   }
   else
   {
      // At end of coast
      PtiProfileStart = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastSize - 1); // Grid CRS
      nProfileStartEdge = m_VCoast[nCoast].nGetEndEdge();
   }

   VPtiNormalPoints.push_back(PtiProfileStart);

   // Find the start cell in the list of edge cells
   auto it = find(m_VEdgeCell.begin(), m_VEdgeCell.end(), PtiProfileStart);
   if (it == m_VEdgeCell.end())
   {
      // Not found. This can happen because of rounding problems, i.e. the cell which was stored as the first cell of the raster coastline
      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
         LogStream << m_ulIter << ": " << ERR << " when constructing start-of-coast profile, [" << PtiProfileStart.nGetX() << "][" << PtiProfileStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiProfileStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiProfileStart.nGetY()) << "} not found in list of edge cells" << endl;

      return RTN_ERR_COAST_CANT_FIND_EDGE_CELL;
   }

   // Found
   int nPos = static_cast<int>(it - m_VEdgeCell.begin());

   // Now construct the edge profile, searching for edge cells
   for (int n = 0; n < nProfileLen; n++)
   {
      if (bCoastStart)
      {
         // At start of coast
         if (nHandedness == LEFT_HANDED)
         {
            // The list of edge cells is in clockwise sequence, go in this direction
            nPos++;

            if (nPos >= static_cast<int>(m_VEdgeCell.size()))
            {
               // We've reached the end of the list of edge cells before the profile is long enough. OK, we can live with this
               break;
            }
         }
         else // Right-handed
         {
            // The list of edge cells is in clockwise sequence, go in the opposite direction
            nPos--;

            if (nPos < 0)
            {
               // We've reached the beginning of the list of edge cells before the profile is long enough. OK, we can live with this
               break;
            }
         }
      }
      else
      {
         // At end of coast
         if (nHandedness == LEFT_HANDED)
         {
            // The list of edge cells is in clockwise sequence, go in the opposite direction
            nPos--;

            if (nPos < 0)
            {
               // We've reached the beginning of the list of edge cells before the profile is long enough. OK, we can live with this
               break;
            }
         }
         else // Right-handed
         {
            // The list of edge cells is in clockwise sequence, go in this direction
            nPos++;

            if (nPos >= static_cast<int>(m_VEdgeCell.size()))
            {
               // We've reached the end of the list of edge cells before the profile is long enough. OK, we can live with this
               break;
            }
         }
      }

      if (m_VEdgeCellEdge[nPos] != nProfileStartEdge)
      {
         // We've reached the end of a grid side before the profile is long enough. OK, we can live with this
         break;
      }

      // Append this grid-edge cell, making sure that there is no gap between this and the previously-appended cell (if there is, will get problems with flood fill)
      AppendEnsureNoGap(&VPtiNormalPoints, &m_VEdgeCell[nPos]);
   }

   // Create the new start-of-coast profile
   if (bCoastStart)
   {
      // Start of coast
      m_VCoast[nCoast].AppendProfile(0, nProfile);
      pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

      // Mark this as a start-of-coast profile
      pProfile->SetStartOfCoast(true);
   }
   else
   {
      // End of coast
      m_VCoast[nCoast].AppendProfile(nCoastSize - 1, nProfile);
      pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

      // Mark this as an end-of-coast profile
      pProfile->SetEndOfCoast(true);
   }

   for (unsigned int n = 0; n < VPtiNormalPoints.size(); n++)
   {
      int nX = VPtiNormalPoints[n].nGetX();
      int nY = VPtiNormalPoints[n].nGetY();

      // Mark each cell in the raster grid
      m_pRasterGrid->m_Cell[nX][nY].SetProfile(nProfile);

      // Store the raster-grid coordinates in the profile object
      pProfile->AppendCellInProfile(nX, nY);

      CGeom2DPoint Pt(dGridCentroidXToExtCRSX(nX), dGridCentroidYToExtCRSY(nY)); // In external CRS

      // Store the external coordinates in the profile object. Note that for this 'special' profile, the coordinates of the cells and the coordinates of points on the profile itself are identical, this is not the case for ordinary profiles
      pProfile->AppendCellInProfileExtCRS(&Pt);

      // Also store the same external ccordinatespoint in this 'special' profile's coastline-normal vector (note that this will contain many points, whereas the coastline-normal vector of ordinary profiles contains only the two start and end points
      pProfile->AppendPointInProfile(&Pt);
   }

   // Get the deep water wave height and orientation values at the end of the profile
   int nEndX = VPtiNormalPoints.back().nGetX();
   int nEndY = VPtiNormalPoints.back().nGetY();

   double dDeepWaterWaveHeight = m_pRasterGrid->m_Cell[nEndX][nEndY].dGetCellDeepWaterWaveHeight();
   double dDeepWaterWaveAngle = m_pRasterGrid->m_Cell[nEndX][nEndY].dGetCellDeepWaterWaveAngle();
   double dDeepWaterWavePeriod = m_pRasterGrid->m_Cell[nEndX][nEndY].dGetCellDeepWaterWavePeriod();

   // And store them for this profile
   pProfile->SetProfileDeepWaterWaveHeight(dDeepWaterWaveHeight);
   pProfile->SetProfileDeepWaterWaveAngle(dDeepWaterWaveAngle);
   pProfile->SetProfileDeepWaterWavePeriod(dDeepWaterWavePeriod);

   // Create the profile's CGeomMultiLine then set nProfile as the only co-incident profile of the only line segment
   pProfile->AppendLineSegment();
   pProfile->AppendCoincidentProfileToLineSegments(make_pair(nProfile, 0));

   // if (m_nLogFileDetail >= LOG_FILE_ALL)
   //    LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " created at coast " << (bCoastStart ? "start" : "end") << " point " << (bCoastStart ? 0 : nCoastSize - 1) << ", from [" << PtiProfileStart.nGetX() << "][" << PtiProfileStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiProfileStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiProfileStart.nGetY()) << "} to [" << VPtiNormalPoints.back().nGetX() << "][" << VPtiNormalPoints.back().nGetY() << "] = {" << dGridCentroidXToExtCRSX(VPtiNormalPoints.back().nGetX()) << ", " << dGridCentroidYToExtCRSY(VPtiNormalPoints.back().nGetY()) << "}" << endl;

   return RTN_OK;
}

//===============================================================================================================================
//! Finds the end point of a coastline-normal line, given the start point on the vector coastline. All co-ordinates are in the external CRS
//===============================================================================================================================
int CSimulation::nGetCoastNormalEndPoint(int const nCoast, int const nStartCoastPoint, int const nCoastSize, CGeom2DPoint const* pPtStart, double const dLineLength, CGeom2DPoint* pPtEnd, CGeom2DIPoint* pPtiEnd, bool const bIntervention)
{
   int const AVGSIZE = 21;    // TODO 011 This should be a user input

   CGeom2DPoint PtBefore;
   CGeom2DPoint PtAfter;

   if (bIntervention)
   {
      // This is an intervention profile, so just use one point on either side (coordinates in external CRS). Note this this assumes that this intervention profile is not at the start or end of the coastline
      PtBefore = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nStartCoastPoint-1);
      PtAfter = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nStartCoastPoint+1);
   }
   else
   {
      // This is not an intervention, so put a maximum of AVGSIZE points before the start point into a vector
      vector<CGeom2DPoint> PtBeforeToAverage;
      for (int n = 1; n <= AVGSIZE; n++)
      {
         int nPoint = nStartCoastPoint - n;
         if (nPoint < 0)
            break;

         PtBeforeToAverage.push_back(*m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nPoint));
      }

      // Put a maximum of AVGSIZE points after the start point into a vector
      vector<CGeom2DPoint> PtAfterToAverage;
      for (int n = 1; n <= AVGSIZE; n++)
      {
         int nPoint = nStartCoastPoint + n;
         if (nPoint > nCoastSize-1)
            break;

         PtAfterToAverage.push_back(*m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nPoint));
      }

      // Now average each of these vectors of points: results are in PtBefore and PtAfter (coordinates in external CRS)
      PtBefore = PtAverage(&PtBeforeToAverage);
      PtAfter = PtAverage(&PtAfterToAverage);
   }

   // Get the y = a * x + b equation of the straight line linking the coastline points before and after 'this' coastline point. For this linking line, slope a = (y2 - y1) / (x2 - x1)
   double dYDiff = PtAfter.dGetY() - PtBefore.dGetY();
   double dXDiff = PtAfter.dGetX() - PtBefore.dGetX();

   double dXEnd1 = 0, dXEnd2 = 0, dYEnd1 = 0, dYEnd2 = 0;

   if (bFPIsEqual(dYDiff, 0.0, TOLERANCE))
   {
      // The linking line runs W-E or E-W, so a straight line at right angles to this runs N-S or S-N. Calculate the two possible end points for this coastline-normal profile
      dXEnd1 = dXEnd2 = pPtStart->dGetX();
      dYEnd1 = pPtStart->dGetY() + dLineLength;
      dYEnd2 = pPtStart->dGetY() - dLineLength;
   }
   else if (bFPIsEqual(dXDiff, 0.0, TOLERANCE))
   {
      // The linking line runs N-S or S-N, so a straight line at right angles to this runs W-E or E-W. Calculate the two possible end points for this coastline-normal profile
      dYEnd1 = dYEnd2 = pPtStart->dGetY();
      dXEnd1 = pPtStart->dGetX() + dLineLength;
      dXEnd2 = pPtStart->dGetX() - dLineLength;
   }
   else
   {
      // The linking line runs neither W-E nor N-S so we have to work a bit harder to find the end-point of the coastline-normal profile
      double dA = dYDiff / dXDiff;

      // Now calculate the equation of the straight line which is perpendicular to this linking line
      double dAPerp = -1 / dA;
      double dBPerp = pPtStart->dGetY() - (dAPerp * pPtStart->dGetX());

      // Calculate the end point of the profile: first do some substitution then rearrange as a quadratic equation i.e. in the form Ax^2 + Bx + C = 0 (see http://math.stackexchange.com/questions/228841/how-do-i-calculate-the-intersections-of-a-straight-line-and-a-circle)
      double dQuadA = 1 + (dAPerp * dAPerp);
      double dQuadB = 2 * ((dBPerp * dAPerp) - (dAPerp * pPtStart->dGetY()) - pPtStart->dGetX());
      double dQuadC = ((pPtStart->dGetX() * pPtStart->dGetX()) + (pPtStart->dGetY() * pPtStart->dGetY()) + (dBPerp * dBPerp) - (2 * pPtStart->dGetY() * dBPerp) - (dLineLength * dLineLength));

      // Solve for x and y using the quadratic formula x = (−B ± sqrt(B^2 − 4AC)) / 2A
      double dDiscriminant = (dQuadB * dQuadB) - (4 * dQuadA * dQuadC);
      if (dDiscriminant < 0)
      {
         LogStream << ERR << "timestep " << m_ulIter << ": discriminant < 0 when finding profile end point on coastline " << nCoast << ", from coastline point " << nStartCoastPoint << "), ignored" << endl;
         return RTN_ERR_NO_SOLUTION_FOR_ENDPOINT;
      }

      dXEnd1 = (-dQuadB + sqrt(dDiscriminant)) / (2 * dQuadA);
      dYEnd1 = (dAPerp * dXEnd1) + dBPerp;
      dXEnd2 = (-dQuadB - sqrt(dDiscriminant)) / (2 * dQuadA);
      dYEnd2 = (dAPerp * dXEnd2) + dBPerp;
   }

   // We have two possible solutions, so decide which of the two endpoints to use then create the profile end-point (coordinates in external CRS)
   int nSeaHand = m_VCoast[nCoast].nGetSeaHandedness(); // Assumes handedness is either 0 or 1 (i.e. not -1)
   *pPtEnd = PtChooseEndPoint(nSeaHand, &PtBefore, &PtAfter, dXEnd1, dYEnd1, dXEnd2, dYEnd2);

   // Check that pPtiEnd is not off the grid. Note that pPtiEnd is NOT (necessarily) a cell centroid
   pPtiEnd->SetXY(nRound(dExtCRSXToGridX(pPtEnd->dGetX())), nRound(dExtCRSYToGridY(pPtEnd->dGetY())));
   if (! bIsWithinValidGrid(pPtiEnd))
   {
      // The end point is off the grid, so constrain it to be within the valid grid
      CGeom2DIPoint PtiStart(nRound(dExtCRSXToGridX(pPtStart->dGetX())), nRound(dExtCRSYToGridY(pPtStart->dGetY())));
      KeepWithinValidGrid(&PtiStart, pPtiEnd);

      pPtEnd->SetX(dGridCentroidXToExtCRSX(pPtiEnd->nGetX()));
      pPtEnd->SetY(dGridCentroidYToExtCRSY(pPtiEnd->nGetY()));

      // LogStream << m_ulIter << ": changed endpoint for profile, is now [" << pPtiEnd->nGetX() << "][" << pPtiEnd->nGetY() << "] = {" << pPtEnd->dGetX() << ", " << pPtEnd->dGetY() << "}. The profile starts at coastline point " << nStartCoastPoint << " = {" << pPtStart->dGetX() << ", " << pPtStart->dGetY() << "}" << endl;
      
      return RTN_ERR_PROFILE_ENDPOINT_AT_GRID_EDGE;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Choose which end point to use for the coastline-normal profile
//===============================================================================================================================
CGeom2DPoint CSimulation::PtChooseEndPoint(int const nHand, CGeom2DPoint const *PtBefore, CGeom2DPoint const *PtAfter, double const dXEnd1, double const dYEnd1, double const dXEnd2, double const dYEnd2)
{
   CGeom2DPoint PtChosen;

   // All co-ordinates here are in the external CRS, so the origin of the grid is the bottom left
   if (nHand == RIGHT_HANDED)
   {
      // The sea is to the right of the linking line. So which way is the linking line oriented? First check the N-S component
      if (PtAfter->dGetY() > PtBefore->dGetY())
      {
         // We are going S to N and the sea is to the right: the normal endpoint is to the E. We want the larger of the two x values
         if (dXEnd1 > dXEnd2)
         {
            PtChosen.SetX(dXEnd1);
            PtChosen.SetY(dYEnd1);
         }
         else
         {
            PtChosen.SetX(dXEnd2);
            PtChosen.SetY(dYEnd2);
         }
      }
      else if (PtAfter->dGetY() < PtBefore->dGetY())
      {
         // We are going N to S and the sea is to the right: the normal endpoint is to the W. We want the smaller of the two x values
         if (dXEnd1 < dXEnd2)
         {
            PtChosen.SetX(dXEnd1);
            PtChosen.SetY(dYEnd1);
         }
         else
         {
            PtChosen.SetX(dXEnd2);
            PtChosen.SetY(dYEnd2);
         }
      }
      else
      {
         // No N-S component i.e. the linking line is exactly W-E. So check the W-E component
         if (PtAfter->dGetX() > PtBefore->dGetX())
         {
            // We are going W to E and the sea is to the right: the normal endpoint is to the s. We want the smaller of the two y values
            if (dYEnd1 < dYEnd2)
            {
               PtChosen.SetX(dXEnd1);
               PtChosen.SetY(dYEnd1);
            }
            else
            {
               PtChosen.SetX(dXEnd2);
               PtChosen.SetY(dYEnd2);
            }
         }
         else // Do not check for (PtAfter->dGetX() == PtBefore->dGetX()), since this would mean the two points are co-incident
         {
            // We are going E to W and the sea is to the right: the normal endpoint is to the N. We want the larger of the two y values
            if (dYEnd1 > dYEnd2)
            {
               PtChosen.SetX(dXEnd1);
               PtChosen.SetY(dYEnd1);
            }
            else
            {
               PtChosen.SetX(dXEnd2);
               PtChosen.SetY(dYEnd2);
            }
         }
      }
   }
   else // nHand == LEFT_HANDED
   {
      // The sea is to the left of the linking line. So which way is the linking line oriented? First check the N-S component
      if (PtAfter->dGetY() > PtBefore->dGetY())
      {
         // We are going S to N and the sea is to the left: the normal endpoint is to the W. We want the smaller of the two x values
         if (dXEnd1 < dXEnd2)
         {
            PtChosen.SetX(dXEnd1);
            PtChosen.SetY(dYEnd1);
         }
         else
         {
            PtChosen.SetX(dXEnd2);
            PtChosen.SetY(dYEnd2);
         }
      }
      else if (PtAfter->dGetY() < PtBefore->dGetY())
      {
         // We are going N to S and the sea is to the left: the normal endpoint is to the E. We want the larger of the two x values
         if (dXEnd1 > dXEnd2)
         {
            PtChosen.SetX(dXEnd1);
            PtChosen.SetY(dYEnd1);
         }
         else
         {
            PtChosen.SetX(dXEnd2);
            PtChosen.SetY(dYEnd2);
         }
      }
      else
      {
         // No N-S component i.e. the linking line is exactly W-E. So check the W-E component
         if (PtAfter->dGetX() > PtBefore->dGetX())
         {
            // We are going W to E and the sea is to the left: the normal endpoint is to the N. We want the larger of the two y values
            if (dYEnd1 > dYEnd2)
            {
               PtChosen.SetX(dXEnd1);
               PtChosen.SetY(dYEnd1);
            }
            else
            {
               PtChosen.SetX(dXEnd2);
               PtChosen.SetY(dYEnd2);
            }
         }
         else // Do not check for (PtAfter->dGetX() == PtBefore->dGetX()), since this would mean the two points are co-incident
         {
            // We are going E to W and the sea is to the left: the normal endpoint is to the S. We want the smaller of the two y values
            if (dYEnd1 < dYEnd2)
            {
               PtChosen.SetX(dXEnd1);
               PtChosen.SetY(dYEnd1);
            }
            else
            {
               PtChosen.SetX(dXEnd2);
               PtChosen.SetY(dYEnd2);
            }
         }
      }
   }

   return PtChosen;
}

//===============================================================================================================================
//! Checks all coastline-normal profiles for intersection, and modifies those that intersect
//===============================================================================================================================
int CSimulation::nModifyAllIntersectingProfiles(void)
{
   // Do once for every coastline object
   int nCoastLines = static_cast<int>(m_VCoast.size());
   for (int nCoast = 0; nCoast < nCoastLines; nCoast++)
   {
      int nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      // Go along the coast, looking at profiles which are increasingly distant from the first profile
      int nMaxDist = nNumProfiles / 2; // Arbitrary
      for (int nDist = 1; nDist < nMaxDist; nDist++)
      {
         // Do once for every profile
         for (int nFirst = 0; nFirst < nNumProfiles; nFirst++)
         {
            // In coastline curvature sequence
            int nFirstProfile = m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(nFirst);
            if (nFirstProfile < 0)
               // Not found
               return RTN_ERR_BAD_INDEX;

            //            LogStream << m_ulIter << ": nFirst = " << nFirst << " nFirstProfile = " << nFirstProfile << endl;

            // Don't modify the start- or end-of coastline normals
            CGeomProfile* pFirstProfile = m_VCoast[nCoast].pGetProfile(nFirstProfile);
            if ((pFirstProfile->bStartOfCoast()) || (pFirstProfile->bEndOfCoast()))
               continue;

            // Do this in alternate directions: first down-coast (in the direction of increasing coast point numbers) then up-coast
            for (int nDirection = DIRECTION_DOWNCOAST; nDirection <= DIRECTION_UPCOAST; nDirection++)
            {
               int nSecond;
               if (nDirection == DIRECTION_DOWNCOAST)
                  nSecond = nFirst + nDist;
               else
                  nSecond = nFirst - nDist;

               if ((nSecond < 0) || (nSecond > nNumProfiles - 1))
                  // Make sure that we don't go beyond the ends of the along-coast index
                  continue;

               int nSecondProfile = m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(nSecond);

               //               LogStream << m_ulIter << ": " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast search, nFirst = " << nFirst << " nSecond = " << nSecond << " (profiles " << nFirstProfile << " and " << nSecondProfile << ")" << endl;

               // Only check these profiles for intersection if both are problem-free
               CGeomProfile* pSecondProfile = m_VCoast[nCoast].pGetProfile(nSecondProfile);
               if (!(pFirstProfile->bProfileOK()) || (!pSecondProfile->bProfileOK()))
                  continue;

               // Only check these two profiles for intersection if they are are not co-incident in the final line segment of both profiles (i.e. the profiles have not already intersected)
               if ((pFirstProfile->bFindProfileInCoincidentProfilesOfLastLineSegment(nSecondProfile)) || (pSecondProfile->bFindProfileInCoincidentProfilesOfLastLineSegment(nFirstProfile)))
                  continue;

               // OK go for it
               int nProf1LineSeg = 0;
               int nProf2LineSeg = 0;
               double dIntersectX = 0;
               double dIntersectY = 0;
               double dAvgEndX = 0;
               double dAvgEndY = 0;

               if (bCheckForIntersection(pFirstProfile, pSecondProfile, nProf1LineSeg, nProf2LineSeg, dIntersectX, dIntersectY, dAvgEndX, dAvgEndY))
               {
                  // The profiles intersect. Is the point of intersection already present in the first profile (i.e. because there has already been an intersection at this point between the first profile and some other profile)?
                  int nPoint = -1;
                  if (pFirstProfile->bIsPointInProfile(dIntersectX, dIntersectY, nPoint))
                  {
                     //                      LogStream << m_ulIter << ": ^^^^ profiles " << nFirstProfile << " and " << nSecondProfile << " intersect, but point {" << dIntersectX << ", " << dIntersectY << "} is already present in profile " << nFirstProfile << " as point " << nPoint << endl;

                     // Truncate the second profile and merge it with the first profile
                     TruncateOneProfileRetainOtherProfile(nCoast, nSecondProfile, nFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg, true);
                  }

                  // Is the point of intersection already present in the second profile?
                  else if (pSecondProfile->bIsPointInProfile(dIntersectX, dIntersectY, nPoint))
                  {
                     //                      LogStream << m_ulIter << ": ^^^^ profiles " << nFirstProfile << " and " << nSecondProfile << " intersect, but point {" << dIntersectX << ", " << dIntersectY << "} is already present in profile " << nSecondProfile << " as point " << nPoint << endl;

                     // Truncate the first profile and merge it with the second profile
                     TruncateOneProfileRetainOtherProfile(nCoast, nFirstProfile, nSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg, true);
                  }

                  else
                  {
                     // The point of intersection is not already present in either profile, so get the number of line segments of each profile
                     int nFirstProfileLineSegments = pFirstProfile->nGetNumLineSegments();
                     int nSecondProfileLineSegments = pSecondProfile->nGetNumLineSegments();

                     //                   assert(nProf1LineSeg < nFirstProfileLineSegments);
                     //                   assert(nProf2LineSeg < nSecondProfileLineSegments);

                     //  Next check whether the point of intersection is on the final line segment of both profiles
                     if ((nProf1LineSeg == (nFirstProfileLineSegments - 1)) && (nProf2LineSeg == (nSecondProfileLineSegments - 1)))
                     {
                        // Yes, the point of intersection is on the final line segment of both profiles, so merge the profiles seaward of the point of intersection
                        MergeProfilesAtFinalLineSegments(nCoast, nFirstProfile, nSecondProfile, nFirstProfileLineSegments, nSecondProfileLineSegments, dIntersectX, dIntersectY, dAvgEndX, dAvgEndY);

                        //                         LogStream << m_ulIter << ": " << ((nDirection == DIRECTION_DOWNCOAST) ? "down" : "up") << "-coast search, end-segment intersection between profiles {" << nFirstProfile << "} and {" << nSecondProfile << "} at [" << dIntersectX << ", " << dIntersectY << "] in line segment [" << nProf1LineSeg << "] of " << nFirstProfileLineSegments << ", and line segment [" << nProf2LineSeg << "] of " << nSecondProfileLineSegments << ", respectively" << endl;

                        //                         // DEBUG CODE =============================
                        //                         int nSizeTmp = pFirstProfile->nGetProfileSize();
                        //                         CGeom2DPoint PtEndTmp = *pFirstProfile->pPtGetPointInProfile(nSizeTmp-1);
                        //
                        //                         LogStream << m_ulIter << ": end of first profile (" << nFirstProfile << ") is point " << nSizeTmp-1 << " at [" << dExtCRSXToGridX(PtEndTmp.dGetX()) << "][" << dExtCRSYToGridY(PtEndTmp.dGetY()) << "} = {" << PtEndTmp.dGetX() << ", " << PtEndTmp.dGetY() << "}" << endl;
                        //
                        //                         nSizeTmp = pSecondProfile->nGetProfileSize();
                        //                         PtEndTmp = *pSecondProfile->pPtGetPointInProfile(nSizeTmp-1);
                        //
                        //                         LogStream << m_ulIter << ": end of second profile (" << nSecondProfile << ") is point " << nSizeTmp-1 << " at [" << dExtCRSXToGridX(PtEndTmp.dGetX()) << "][" << dExtCRSYToGridY(PtEndTmp.dGetY()) << "} = {" << PtEndTmp.dGetX() << ", " << PtEndTmp.dGetY() << "}" << endl;
                        //                         // DEBUG CODE =============================
                     }
                     else
                     {
                        // The profiles intersect, but the point of intersection is not on the final line segment of both profiles. One of the profiles will be truncated, the other profile will be retained
                        //                         LogStream << m_ulIter << ": " << ((nDirection == DIRECTION_DOWNCOAST) ? "down" : "up") << "-coast search, intersection (NOT both end segments) between profiles {" << nFirstProfile << "} and {" << nSecondProfile << "} at [" << dIntersectX << ", " << dIntersectY << "] in line segment [" << nProf1LineSeg << "] of " << nFirstProfileLineSegments << ", and line segment [" << nProf2LineSeg << "] of " << nSecondProfileLineSegments << ", respectively" << endl;

                        // Decide which profile to truncate, and which to retain
                        if (nFirstProfileLineSegments > nSecondProfileLineSegments)
                           // Truncate the second profile, since it has a smaller number of line segments
                           TruncateOneProfileRetainOtherProfile(nCoast, nSecondProfile, nFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg, false);

                        else if (nFirstProfileLineSegments < nSecondProfileLineSegments)
                           // Truncate the first profile, since it has a smaller number of line segments
                           TruncateOneProfileRetainOtherProfile(nCoast, nFirstProfile, nSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg, false);

                        else
                        {
                           // Both profiles have the same number of line segments, so choose randomnly
                           if (dGetRand0d1() >= 0.5)
                              TruncateOneProfileRetainOtherProfile(nCoast, nFirstProfile, nSecondProfile, dIntersectX, dIntersectY, nProf1LineSeg, nProf2LineSeg, false);
                           else
                              TruncateOneProfileRetainOtherProfile(nCoast, nSecondProfile, nFirstProfile, dIntersectX, dIntersectY, nProf2LineSeg, nProf1LineSeg, false);
                        }
                     }
                  }

                  //                   int
                  //                      nProfile1NumSegments = pFirstProfile->nGetNumLineSegments(),
                  //                      nProfile2NumSegments = pSecondProfile->nGetNumLineSegments(),
                  //                      nProfile1Size = pFirstProfile->nGetProfileSize(),
                  //                      nProfile2Size = pSecondProfile->nGetProfileSize();
                  //
                  //                   assert(pFirstProfile->nGetNumLineSegments() > 0);
                  //                   assert(pSecondProfile->nGetNumLineSegments() > 0);
                  //                   assert(nProfile1Size == nProfile1NumSegments+1);
                  //                   assert(nProfile2Size == nProfile2NumSegments+1);
               }
            }
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Checks all line segments of a pair of coastline-normal profiles for intersection. If the lines intersect, returns true with numbers of the line segments at which intersection occurs in nProfile1LineSegment and nProfile1LineSegment, the intersection point in dXIntersect and dYIntersect, and the 'average' seaward endpoint of the two intersecting profiles at dXAvgEnd and dYAvgEnd
//===============================================================================================================================
bool CSimulation::bCheckForIntersection(CGeomProfile *const pVProfile1, CGeomProfile *const pVProfile2, int &nProfile1LineSegment, int& nProfile2LineSegment, double& dXIntersect, double& dYIntersect, double& dXAvgEnd, double& dYAvgEnd)
{
   // For both profiles, look at all line segments
   int nProfile1NumSegments = pVProfile1->nGetNumLineSegments();
   int nProfile2NumSegments = pVProfile2->nGetNumLineSegments();
   //       nProfile1Size = pVProfile1->nGetProfileSize(),
   //       nProfile2Size = pVProfile2->nGetProfileSize();

   //    assert(nProfile1Size == nProfile1NumSegments+1);
   //    assert(nProfile2Size == nProfile2NumSegments+1);

   for (int i = 0; i < nProfile1NumSegments; i++)
   {
      for (int j = 0; j < nProfile2NumSegments; j++)
      {
         // In external coordinates
         double dX1 = pVProfile1->pPtVGetPoints()->at(i).dGetX();
         double dY1 = pVProfile1->pPtVGetPoints()->at(i).dGetY();
         double dX2 = pVProfile1->pPtVGetPoints()->at(i + 1).dGetX();
         double dY2 = pVProfile1->pPtVGetPoints()->at(i + 1).dGetY();

         double dX3 = pVProfile2->pPtVGetPoints()->at(j).dGetX();
         double dY3 = pVProfile2->pPtVGetPoints()->at(j).dGetY();
         double dX4 = pVProfile2->pPtVGetPoints()->at(j + 1).dGetX();
         double dY4 = pVProfile2->pPtVGetPoints()->at(j + 1).dGetY();

         // Uses Cramer's Rule to solve the equations. Modified from code at http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect (in turn based on Andre LeMothe's "Tricks of the Windows Game Programming Gurus")
         double dDiffX1 = dX2 - dX1;
         double dDiffY1 = dY2 - dY1;
         double dDiffX2 = dX4 - dX3;
         double dDiffY2 = dY4 - dY3;

         double dS = -999;
         double dT = -999;
         double dTmp = 0;

         dTmp = -dDiffX2 * dDiffY1 + dDiffX1 * dDiffY2;
         if (! bFPIsEqual(dTmp, 0.0, TOLERANCE))
            dS = (-dDiffY1 * (dX1 - dX3) + dDiffX1 * (dY1 - dY3)) / dTmp;

         dTmp = -dDiffX2 * dDiffY1 + dDiffX1 * dDiffY2;
         if (! bFPIsEqual(dTmp, 0.0, TOLERANCE))
            dT = (dDiffX2 * (dY1 - dY3) - dDiffY2 * (dX1 - dX3)) / dTmp;

         if (dS >= 0 && dS <= 1 && dT >= 0 && dT <= 1)
         {
            // Collision detected, calculate intersection co-ordinates
            dXIntersect = dX1 + (dT * dDiffX1);
            dYIntersect = dY1 + (dT * dDiffY1);

            // And calc the average end-point co-ordinates
            dXAvgEnd = (dX2 + dX4) / 2;
            dYAvgEnd = (dY2 + dY4) / 2;

            // Get the line segments at which intersection occurred
            nProfile1LineSegment = i;
            nProfile2LineSegment = j;

            //             LogStream << "\t" << "INTERSECTION dX2 = " << dX2 << " dX4 = " << dX4 << " dY2 = " << dY2 << " dY4 = " << dY4 << endl;
            return true;
         }
      }
   }

   // No collision
   return false;
}

//===============================================================================================================================
//! Puts the coastline-normal profiles (apart from the 'special' ones at the start and end of the coast) onto the raster grid, i.e. rasterizes multi-line vector objects onto the raster grid. Note that this doesn't work if the vector has already been interpolated to fit on the grid i.e. if distances between vector points are just one cell apart
//===============================================================================================================================
int CSimulation::nPutAllProfilesOntoGrid(void)
{
   int nValidProfiles = 0;

   // Do once for every coastline
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // How many profiles on this coast?
      int nProfiles = m_VCoast[nCoast].nGetNumProfiles();
      if (nProfiles == 0)
      {
         // This can happen if the coastline is very short, so just give a warning and carry on with the next coastline
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << WARN << m_ulIter << ": coast " << nCoast << " has no profiles" << endl;
         continue;
      }

      // Now do this loop for every profile in the list, apart from the last two (i.e. the profiles at the start and end of the coast)
      for (int nProfile = 0; nProfile < nProfiles - 2; nProfile++)
      {
         CGeomProfile *const pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

         // If this profile has a problem, then forget about it. Still do start- and end-of-coast profiles however
         if (!pProfile->bOKIncStartAndEndOfCoast())
            continue;

         int nPoints = pProfile->nGetProfileSize();
         if (nPoints < 2)
         {
            // Need at least two points in the profile, so this profile is invalid: mark it
            pProfile->SetTooShort(true);

            if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
               LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " is invalid, has only " << nPoints << " points" << endl;

            continue;
         }

         // OK, go for it: set up temporary vectors to hold the x-y coords (in grid CRS) of the cells which we will mark
         vector<CGeom2DIPoint> VCellsToMark;
         vector<bool> bVShared;
         bool bTooShort = false;
         bool bTruncated = false;
         bool bHitCoast = false;
         bool bHitLand = false;
         bool bHitAnotherProfile = false;

         RasterizeProfile(nCoast, nProfile, &VCellsToMark, &bVShared, bTooShort, bTruncated, bHitCoast, bHitLand, bHitAnotherProfile);

         if ((bTruncated && (! ACCEPT_SHORT_PROFILES)) || bTooShort || bHitCoast || bHitLand || bHitAnotherProfile)
            continue;

         // This profile is fine
         nValidProfiles++;

         for (unsigned int k = 0; k < VCellsToMark.size(); k++)
         {
            // So mark each cell in the raster grid
            m_pRasterGrid->m_Cell[VCellsToMark[k].nGetX()][VCellsToMark[k].nGetY()].SetProfile(nProfile);

            // Store the raster-grid coordinates in the profile object
            pProfile->AppendCellInProfile(VCellsToMark[k].nGetX(), VCellsToMark[k].nGetY());

            // And also store the external coordinates in the profile object
            pProfile->AppendCellInProfileExtCRS(dGridCentroidXToExtCRSX(VCellsToMark[k].nGetX()), dGridCentroidYToExtCRSY(VCellsToMark[k].nGetY()));

            // Mark the shared (i.e. multi-line) parts of the profile (if any)
            //             if (bVShared[k])
            //             {
            //                pProfile->AppendPointShared(true);
            // //                  LogStream << m_ulIter << ": profile " << j << " point " << k << " marked as shared" << endl;
            //             }
            //             else
            //             {
            //                pProfile->AppendPointShared(false);
            // //                  LogStream << m_ulIter << ": profile " << nProfile << " point " << k << " marked as NOT shared" << endl;
            //             }
         }

         // Get the deep water wave height and orientation values at the end of the profile
         double dDeepWaterWaveHeight = m_pRasterGrid->m_Cell[VCellsToMark.back().nGetX()][VCellsToMark.back().nGetY()].dGetCellDeepWaterWaveHeight();
         double dDeepWaterWaveAngle = m_pRasterGrid->m_Cell[VCellsToMark.back().nGetX()][VCellsToMark.back().nGetY()].dGetCellDeepWaterWaveAngle();
         double dDeepWaterWavePeriod = m_pRasterGrid->m_Cell[VCellsToMark.back().nGetX()][VCellsToMark.back().nGetY()].dGetCellDeepWaterWavePeriod();

         // And store them for this profile
         pProfile->SetProfileDeepWaterWaveHeight(dDeepWaterWaveHeight);
         pProfile->SetProfileDeepWaterWaveAngle(dDeepWaterWaveAngle);
         pProfile->SetProfileDeepWaterWavePeriod(dDeepWaterWavePeriod);
      }
   }

   if (nValidProfiles == 0)
   {
      // Problem! No valid profiles, so quit
      cerr << m_ulIter << ": " << ERR << "no coastline-normal profiles created" << endl;
      return RTN_ERR_NO_PROFILES_2;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Given a pointer to a coastline-normal profile, returns an output vector of cells which are 'under' every line segment of the profile. If there is a problem with the profile (e.g. a rasterized cell is dry land or coast, or the profile has to be truncated) then we pass this back as an error code
//===============================================================================================================================
void CSimulation::RasterizeProfile(int const nCoast, int const nProfile, vector<CGeom2DIPoint>* pVIPointsOut, vector<bool>* pbVShared, bool& bTooShort, bool& bTruncated, bool& bHitCoast, bool& bHitLand, bool& bHitAnotherProfile)
{
   CGeomProfile* const pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   pVIPointsOut->clear();
   int nSeg = 0;
   int nSegments = pProfile->nGetNumLineSegments();
   int nProfiles = m_VCoast[nCoast].nGetNumProfiles(); // TODO 015 This is a bodge, needed if we hit a profile which belongs to a different coast object

   // This contains a list of 'other' profiles which intersect this profile, but which are OK (i.e. they are coincident profiles of this profile). Added for efficiency
   vector<int> VnOtherProfileOK;

   // Do for every segment of this profile
   for (nSeg = 0; nSeg < nSegments; nSeg++)
   {
      // Do once for every line segment
      vector<CGeom2DPoint> PtVSegment;
      PtVSegment.push_back(*pProfile->pPtGetPointInProfile(nSeg));
      PtVSegment.push_back(*pProfile->pPtGetPointInProfile(nSeg + 1)); // This is OK

      bool bShared = false;
      if (pProfile->nGetNumCoincidentProfilesInLineSegment(nSeg) > 1)
         bShared = true;

      // The start point of the normal, must convert from the external CRS to grid CRS. If this is the first line segment of the profile, then the start point is the centroid of a coastline cell
      double dXStart = dExtCRSXToGridX(PtVSegment[0].dGetX());
      double dYStart = dExtCRSYToGridY(PtVSegment[0].dGetY());

      // The end point of the normal, again convert from the external CRS to grid CRS. Note too that it could be off the grid
      double dXEnd = dExtCRSXToGridX(PtVSegment[1].dGetX());
      double dYEnd = dExtCRSYToGridY(PtVSegment[1].dGetY());
          
      // Safety check
      if (bFPIsEqual(dXStart, dXEnd, TOLERANCE) && bFPIsEqual(dYStart, dYEnd, TOLERANCE))
         continue;

      // Interpolate between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
      double dXInc = dXEnd - dXStart;
      double dYInc = dYEnd - dYStart;
      double dLength = tMax(tAbs(dXInc), tAbs(dYInc));

      dXInc /= dLength;
      dYInc /= dLength;

      double dX = dXStart;
      double dY = dYStart;

      // Process each interpolated point
      for (int m = 0; m <= nRound(dLength); m++)
      {
         // int nX = static_cast<int>(dX);
         // int nY = static_cast<int>(dY);
         int nX = nRound(dX);
         int nY = nRound(dY);

         // Do some checking of this interpolated point, but only if this is not a grid-edge profile (these profiles are always valid)
         if ((!pProfile->bStartOfCoast()) && (!pProfile->bEndOfCoast()))
         {
            // Is the interpolated point within the valid raster grid?
            if (! bIsWithinValidGrid(nX, nY))
            {
               // It is outside the valid grid, so mark this profile and quit the loop
               bTruncated = true;

               if (! ACCEPT_SHORT_PROFILES)
                  pProfile->SetTruncated(true);

               if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                  LogStream << m_ulIter << ": profile " << nProfile << " is invalid, truncated at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

               break;
            }

            // If this is the first line segment of the profile, then once we are clear of the coastline (say, when m > 1), check if this profile hits land at this interpolated point. NOTE Get problems here since if the coastline vector has been heavily smoothed, this can result is 'false positives' profiles marked as invalid which are not actually invalid, because the profile hits land when m = 0 or m = 1. This results in some cells being flagged as profile cells which are actually inland
            if (((nSeg == 0) && (m > 1)) || (nSeg > 0))
            {
               if (m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
               {
                  // We've hit a coastline so set a switch and mark the profile, then quit
                  bHitCoast = true;
                  pProfile->SetHitCoast(true);

                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ": Profile " << nProfile << " is invalid, hit coast at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

                  return;
               }

               // Two diagonal(ish) raster lines can cross each other without any intersection, so must also test an adjacent cell for intersection (does not matter which adjacent cell)
               if (bIsWithinValidGrid(nX, nY + 1) && m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline())
               {
                  // We've hit a coastline so set a switch and mark the profile, then quit
                  bHitCoast = true;
                  pProfile->SetHitCoast(true);

                  if (m_nLogFileDetail >= LOG_FILE_ALL)
                     LogStream << m_ulIter << ": Profile " << nProfile << " is invalid, hit coast at [" << nX << "][" << nY + 1 << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY + 1) << "}" << endl;

                  return;
               }

               //                if (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea())
               //                {
               //                   // We've hit dry land so set a switch and mark the profile
               //                   bHitLand = true;
               //                   pProfile->SetHitLand(true);
               //
               //                   LogStream << m_ulIter << ": profile " << nProfile << " HIT LAND at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, elevation = " << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() << ", SWL = " << m_dThisIterSWL << endl;
               //
               //                   LogStream << "On [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, sea depth = " << m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth() << ", bIsInContiguousSea = " << (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() ? "true" : "false") << ", landform = " << (m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea() ? "sea" : "not sea") << endl;
               //
               //                   LogStream << "On [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, elevation of consolidated sediment = " << m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopForLayerAboveBasement(m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement()) << ", total cliff collapse = " << m_pRasterGrid->m_Cell[nX][nY].dGetTotCliffCollapse() << ", total beach deposition = " << m_pRasterGrid->m_Cell[nX][nY].dGetTotBeachDeposition() << endl;
               //                   return;

               //                }
            }

            // Now check to see if we hit another profile which is not a coincident normal to this normal
            if (m_pRasterGrid->m_Cell[nX][nY].bIsProfile())
            {
               // We've hit a raster cell which is already marked as 'under' a normal profile. Get the number of the profile which marked this cell
               int nHitProfile = m_pRasterGrid->m_Cell[nX][nY].nGetProfile();

               // Search for the 'other' profile in the list of OK profiles
               if (std::find(VnOtherProfileOK.begin(), VnOtherProfileOK.end(), nHitProfile) == VnOtherProfileOK.end())
               {
                  // Not found, so it must be the first time we've encountered this 'other' profile
                  //                   LogStream << "VVVV Profile " << nProfile << " touches profile " << nHitProfile << endl;
                  VnOtherProfileOK.push_back(nHitProfile);

                  // TODO 015 Bodge in case we hit a profile which belongs to a different coast
                  if (nHitProfile > nProfiles - 1)
                  {
                     bHitAnotherProfile = true;
                     pProfile->SetHitAnotherProfile(true);

                     if (m_nLogFileDetail >= LOG_FILE_ALL)
                        LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " is invalid, crosses another profile A1 (" << nHitProfile << ") at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

                     return;
                  }

                  // Only set the flag if this isn't a coincident normal to one or other of the profiles
                  // // DEBUG CODE ==============================
                  //                   bool
                  //                      bOtherConcidentToThis = pProfile->bFindProfileInCoincidentProfiles(nHitProfile),
                  //                      bThisCoincidentToOther = m_VCoast[nCoast].pGetProfile(nHitProfile)->bFindProfileInCoincidentProfiles(nProfile);
                  //
                  //                   if (bOtherConcidentToThis)
                  //                      LogStream << "Profile " << nHitProfile << " is coincident with " << nProfile << endl;
                  //
                  //                   if (bThisCoincidentToOther)
                  //                      LogStream << "Profile " << nProfile << " is coincident with " << nHitProfile << endl;
                  //
                  //                   if ((! bOtherConcidentToThis) && (! bThisCoincidentToOther))
                  //                   // DEBUG CODE ==============================

                  // Note only need to test whether the 'other' profile is concident with his since the concidence relationship is mutual
                  // if ((! pProfile->bFindProfileInCoincidentProfiles(nHitProfile)) && (! m_VCoast[nCoast].pGetProfile(nHitProfile)->bFindProfileInCoincidentProfiles(nProfile)))
                  if (!pProfile->bFindProfileInCoincidentProfiles(nHitProfile))
                  {
                     bHitAnotherProfile = true;
                     pProfile->SetHitAnotherProfile(true);

                     if (m_nLogFileDetail >= LOG_FILE_ALL)
                        LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " is invalid, crosses another profile B1 (" << nHitProfile << ") at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

                     return;
                  }
               }
            }

            // Try diagonal too
            int
                nXTmp = nX + 1,
                nYTmp = nY;

            if (bIsWithinValidGrid(nXTmp, nYTmp) && m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsProfile())
            {
               // We've hit a raster cell which is already marked as 'under' a normal profile. Get the number of the profile which marked this cell
               int nHitProfile = m_pRasterGrid->m_Cell[nXTmp][nYTmp].nGetProfile();

               // Search for the 'other' profile in the list of OK profiles
               if (std::find(VnOtherProfileOK.begin(), VnOtherProfileOK.end(), nHitProfile) == VnOtherProfileOK.end())
               {
                  // Not found, so it must be the first time we've encountered this 'other' profile
                  //                   LogStream << "VVVV Profile " << nProfile << " touches profile " << nHitProfile << " diagonally" << endl;
                  VnOtherProfileOK.push_back(nHitProfile);

                  // TODO 015 Bodge in case we hit a profile which belongs to a different coast
                  if (nHitProfile > nProfiles - 1)
                  {
                     bHitAnotherProfile = true;
                     pProfile->SetHitAnotherProfile(true);

                     // if (m_nLogFileDetail >= LOG_FILE_ALL)
                     //    LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " is invalid, crosses another profile diagonally A2 (" << nHitProfile << ") at [" << nXTmp << "][" << nYTmp << "] = {" << dGridCentroidXToExtCRSX(nXTmp) << ", " << dGridCentroidYToExtCRSY(nYTmp) << "}" << endl;

                     return;
                  }

                  // Only set the flag if this isn't a coincident normal to one or other of the profiles
                  //                   // DEBUG CODE ==============================
                  //                   bool
                  //                      bOtherConcidentToThis = pProfile->bFindProfileInCoincidentProfiles(nHitProfile),
                  //                      bThisCoincidentToOther = m_VCoast[nCoast].pGetProfile(nHitProfile)->bFindProfileInCoincidentProfiles(nProfile);
                  //
                  //                   if (bOtherConcidentToThis)
                  //                      LogStream << "Profile " << nHitProfile << " is coincident with " << nProfile << endl;
                  //
                  //                   if (bThisCoincidentToOther)
                  //                      LogStream << "Profile " << nProfile << " is coincident with " << nHitProfile << endl;
                  //
                  //                   if ((! bOtherConcidentToThis) && (! bThisCoincidentToOther))
                  //                   // DEBUG CODE ==============================

                  // Note only need to test whether the 'other' profile is concident with his since the concidence relationship is mutual
                  // if ((! pProfile->bFindProfileInCoincidentProfiles(nHitProfile)) && (! m_VCoast[nCoast].pGetProfile(nHitProfile)->bFindProfileInCoincidentProfiles(nProfile)))
                  if (! pProfile->bFindProfileInCoincidentProfiles(nHitProfile))
                  {
                     bHitAnotherProfile = true;
                     pProfile->SetHitAnotherProfile(true);

                     // if (m_nLogFileDetail >= LOG_FILE_ALL)
                     //    LogStream << m_ulIter << ": coast " << nCoast << ", profile " << nProfile << " is invalid, crosses another profile diagonally B2 (" << nHitProfile << ") at [" << nXTmp << "][" << nYTmp << "] = {" << dGridCentroidXToExtCRSX(nXTmp) << ", " << dGridCentroidYToExtCRSY(nYTmp) << "}" << endl;

                     return;
                  }
               }
            }
         }

         // Append this point to the output vector
         pVIPointsOut->push_back(CGeom2DIPoint(nX, nY)); // Is in raster-grid co-ordinates
         pbVShared->push_back(bShared);

         // And increment for next time
         dX += dXInc;
         dY += dYInc;
      }

      if (bTruncated)
         break;
   }

   if (bTruncated)
   {
      if (nSeg < (nSegments - 1))
         // We are truncating the profile, so remove any line segments after this one
         pProfile->TruncateLineSegments(nSeg);

      // Shorten the vector input. Ignore CPPCheck errors here, since we know that pVIPointsOut is not empty
      int
          nLastX = pVIPointsOut->at(pVIPointsOut->size() - 1).nGetX(),
          nLastY = pVIPointsOut->at(pVIPointsOut->size() - 1).nGetY();

      pProfile->pPtGetPointInProfile(nSeg + 1)->SetX(dGridCentroidXToExtCRSX(nLastX));
      pProfile->pPtGetPointInProfile(nSeg + 1)->SetY(dGridCentroidYToExtCRSY(nLastY));
   }

   if (pVIPointsOut->size() < 3)
   {
      // For shore platform coastline-normal profiles, we cannot have very short normal profiles with less than 3 cells, since we cannot calculate along-profile slope properly for such short profiles
      bTooShort = true;
      pProfile->SetTooShort(true);

      if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      {
         // Ignore CPPCheck errors here, since we know that pVIPointsOut is not empty
         LogStream << m_ulIter << ": profile " << nProfile << " is invalid, is too short, only " << pVIPointsOut->size() << " points, HitLand?" << bHitLand << ". From [" << pVIPointsOut->at(0).nGetX() << "][" << pVIPointsOut->at(0).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pVIPointsOut->at(0).nGetX()) << ", " << dGridCentroidYToExtCRSY(pVIPointsOut->at(0).nGetY()) << "} to [" << pVIPointsOut->at(pVIPointsOut->size() - 1).nGetX() << "][" << pVIPointsOut->at(pVIPointsOut->size() - 1).nGetY() << "] = {" << dGridCentroidXToExtCRSX(pVIPointsOut->at(pVIPointsOut->size() - 1).nGetX()) << ", " << dGridCentroidYToExtCRSY(pVIPointsOut->at(pVIPointsOut->size() - 1).nGetY()) << "}" << endl;
      }
   }
}

//===============================================================================================================================
//! Merges two profiles which intersect at their final (most seaward) line segments, seaward of their point of intersection
//===============================================================================================================================
void CSimulation::MergeProfilesAtFinalLineSegments(int const nCoast, int const nFirstProfile, int const nSecondProfile, int const nFirstProfileLineSegments, int const nSecondProfileLineSegments, double const dIntersectX, double const dIntersectY, double const dAvgEndX, double const dAvgEndY)
{
   // The point of intersection is on the final (most seaward) line segment of both profiles. Put together a vector of coincident profile numbers (with no duplicates) for both profiles
   int nCombinedLastSeg = 0;
   vector<pair<int, int>> prVCombinedProfilesCoincidentProfilesLastSeg;
   CGeomProfile* pFirstProfile = m_VCoast[nCoast].pGetProfile(nFirstProfile);
   CGeomProfile* pSecondProfile = m_VCoast[nCoast].pGetProfile(nSecondProfile);

   for (unsigned int n = 0; n < pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(nFirstProfileLineSegments - 1)->size(); n++)
   {
      pair<int, int> prTmp;
      prTmp.first = pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(nFirstProfileLineSegments - 1)->at(n).first;
      prTmp.second = pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(nFirstProfileLineSegments - 1)->at(n).second;

      bool bFound = false;
      for (unsigned int m = 0; m < prVCombinedProfilesCoincidentProfilesLastSeg.size(); m++)
      {
         if (prVCombinedProfilesCoincidentProfilesLastSeg[m].first == prTmp.first)
         {
            bFound = true;
            break;
         }
      }

      if (! bFound)
      {
         prVCombinedProfilesCoincidentProfilesLastSeg.push_back(prTmp);
         nCombinedLastSeg++;
      }
   }

   for (unsigned int n = 0; n < pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(nSecondProfileLineSegments - 1)->size(); n++)
   {
      pair<int, int> prTmp;
      prTmp.first = pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(nSecondProfileLineSegments - 1)->at(n).first;
      prTmp.second = pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(nSecondProfileLineSegments - 1)->at(n).second;

      bool bFound = false;
      for (unsigned int m = 0; m < prVCombinedProfilesCoincidentProfilesLastSeg.size(); m++)
      {
         if (prVCombinedProfilesCoincidentProfilesLastSeg[m].first == prTmp.first)
         {
            bFound = true;
            break;
         }
      }

      if (! bFound)
      {
         prVCombinedProfilesCoincidentProfilesLastSeg.push_back(prTmp);
         nCombinedLastSeg++;
      }
   }

   // Increment the number of each line segment
   for (int m = 0; m < nCombinedLastSeg; m++)
      prVCombinedProfilesCoincidentProfilesLastSeg[m].second++;

   vector<pair<int, int>> prVFirstProfileCoincidentProfilesLastSeg = *pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(nFirstProfileLineSegments - 1);
   vector<pair<int, int>> prVSecondProfileCoincidentProfilesLastSeg = *pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(nSecondProfileLineSegments - 1);
   int nNumFirstProfileCoincidentProfilesLastSeg = static_cast<int>(prVFirstProfileCoincidentProfilesLastSeg.size());
   int nNumSecondProfileCoincidentProfilesLastSeg = static_cast<int>(prVSecondProfileCoincidentProfilesLastSeg.size());

   //   LogStream << m_ulIter << ": END-SEGMENT INTERSECTION between profiles {" << nFirstProfile << "} and {" << nSecondProfile << "} at line segment " << nFirstProfileLineSegments-1 << "/" << nFirstProfileLineSegments-1 << ", and line segment " << nSecondProfileLineSegments-1 << "/" << nSecondProfileLineSegments-1 << ", respectively. Both truncated at [" << dIntersectX << ", " << dIntersectY << "] then profiles {" << nFirstProfile << "} and {" << nSecondProfile << "} extended to [" << dAvgEndX << ", " << dAvgEndY << "]" << endl;

   // Truncate the first profile, and all co-incident profiles, at the point of intersection
   for (int n = 0; n < nNumFirstProfileCoincidentProfilesLastSeg; n++)
   {
      int nThisProfile = prVFirstProfileCoincidentProfilesLastSeg[n].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);
      int nProfileLength = pThisProfile->nGetProfileSize();

      // This is the final line segment of the first 'main' profile. We are assuming that it is also the final line segment of all co-incident profiles. This is fine, altho' each profile may well have a different number of line segments landwards i.e. the number of the line segment may be different for each co-incident profile
      pThisProfile->SetPointInProfile(nProfileLength - 1, dIntersectX, dIntersectY);
   }

   // Truncate the second profile, and all co-incident profiles, at the point of intersection
   for (int n = 0; n < nNumSecondProfileCoincidentProfilesLastSeg; n++)
   {
      int nThisProfile = prVSecondProfileCoincidentProfilesLastSeg[n].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);
      int nProfileLength = pThisProfile->nGetProfileSize();

      // This is the final line segment of the second 'main' profile. We are assuming that it is also the final line segment of all co-incident profiles. This is fine, altho' each profile may well have a different number of line segments landwards i.e. the number of the line segment may be different for each co-incident profile
      pThisProfile->SetPointInProfile(nProfileLength - 1, dIntersectX, dIntersectY);
   }

   // Append a new straight line segment to the existing line segment(s) of the first profile, and to all co-incident profiles
   for (int nThisLineSeg = 0; nThisLineSeg < nNumFirstProfileCoincidentProfilesLastSeg; nThisLineSeg++)
   {
      int nThisProfile = prVFirstProfileCoincidentProfilesLastSeg[nThisLineSeg].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Update this profile
      pThisProfile->AppendPointInProfile(dAvgEndX, dAvgEndY);

      // Append details of the combined profiles
      pThisProfile->AppendLineSegment();
      for (int m = 0; m < nCombinedLastSeg; m++)
         pThisProfile->AppendCoincidentProfileToLineSegments(prVCombinedProfilesCoincidentProfilesLastSeg[m]);
   }

   // Append a new straight line segment to the existing line segment(s) of the second profile, and to all co-incident profiles
   for (int nThisLineSeg = 0; nThisLineSeg < nNumSecondProfileCoincidentProfilesLastSeg; nThisLineSeg++)
   {
      int nThisProfile = prVSecondProfileCoincidentProfilesLastSeg[nThisLineSeg].first;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Update this profile
      pThisProfile->AppendPointInProfile(dAvgEndX, dAvgEndY);

      // Append details of the combined profiles
      pThisProfile->AppendLineSegment();
      for (int m = 0; m < nCombinedLastSeg; m++)
         pThisProfile->AppendCoincidentProfileToLineSegments(prVCombinedProfilesCoincidentProfilesLastSeg[m]);
   }

   // START: FOR CHECKING PURPOSES ****************************************************************
   //    int nFirstProfileLineSeg= pFirstProfile->nGetNumLineSegments();
   //    int nSecondProfileLineSeg = pSecondProfile->nGetNumLineSegments();
   //
   //    LogStream << "\tProfile {" << nFirstProfile << "} now has " << nFirstProfileLineSeg << " line segments" << endl;
   //    for (int m = 0; m < nFirstProfileLineSeg; m++)
   //    {
   //       vector<pair<int, int> > prVCoincidentProfiles = *pFirstProfile->pprVGetPairedCoincidentProfilesForLineSegment(m);
   //       LogStream << "\tCo-incident profiles and line segments for line segment " << m << " of profile {" << nFirstProfile << "} are {";
   //       for (int nn = 0; nn < prVCoincidentProfiles.size(); nn++)
   //          LogStream << " " << prVCoincidentProfiles[nn].first << "[" << prVCoincidentProfiles[nn].second << "] ";
   //       LogStream << " }" << endl;
   //    }
   //    LogStream << "\tProfile {" << nSecondProfile << "} now has " << nSecondProfileLineSeg << " line segments" << endl;
   //    for (int m = 0; m < nSecondProfileLineSeg; m++)
   //    {
   //       vector<pair<int, int> > prVCoincidentProfiles = *pSecondProfile->pprVGetPairedCoincidentProfilesForLineSegment(m);
   //       LogStream << "\tCo-incident profiles and line segments for line segment " << m << " of profile {" << nSecondProfile << "} are {";
   //       for (int nn = 0; nn < prVCoincidentProfiles.size(); nn++)
   //          LogStream << " " << prVCoincidentProfiles[nn].first << "[" << prVCoincidentProfiles[nn].second << "] ";
   //       LogStream << " }" << endl;
   //    }
   // END: FOR CHECKING PURPOSES ******************************************************************
}

//===============================================================================================================================
//! Truncates one intersecting profile at the point of intersection, and retains the other profile
//===============================================================================================================================
void CSimulation::TruncateOneProfileRetainOtherProfile(int const nCoast, int const nProfileToTruncate, int const nProfileToRetain, double const dIntersectX, double const dIntersectY, int const nProfileToTruncateIntersectLineSeg, int const nProfileToRetainIntersectLineSeg, bool const bAlreadyPresent)
{
   // Insert the intersection point into the main retain-profile if it is not already in the profile, and do the same for all co-incident profiles of the main retain-profile. Also add details of the to-truncate profile (and all its coincident profiles) to every line segment of the main to-retain profile which is seaward of the point of intersection
   int nRet = nInsertPointIntoProfilesIfNeededThenUpdate(nCoast, nProfileToRetain, dIntersectX, dIntersectY, nProfileToRetainIntersectLineSeg, nProfileToTruncate, nProfileToTruncateIntersectLineSeg, bAlreadyPresent);
   if (nRet != RTN_OK)
   {
      LogStream << m_ulIter << ": error in nInsertPointIntoProfilesIfNeededThenUpdate()" << endl;
      return;
   }

   // Get all profile points of the main retain-profile seawards from the intersection point, and do the same for the corresponding line segments (including coincident profiles). This also includes details of the main to-truncate profile (and all its coincident profiles)
   CGeomProfile* pProfileToRetain = m_VCoast[nCoast].pGetProfile(nProfileToRetain);
   vector<CGeom2DPoint> PtVProfileLastPart;
   vector<vector<pair<int, int>>> prVLineSegLastPart;
   if (bAlreadyPresent)
   {
      PtVProfileLastPart = pProfileToRetain->PtVGetThisPointAndAllAfter(nProfileToRetainIntersectLineSeg);
      prVLineSegLastPart = pProfileToRetain->prVVGetAllLineSegAfter(nProfileToRetainIntersectLineSeg);
   }
   else
   {
      PtVProfileLastPart = pProfileToRetain->PtVGetThisPointAndAllAfter(nProfileToRetainIntersectLineSeg + 1);
      prVLineSegLastPart = pProfileToRetain->prVVGetAllLineSegAfter(nProfileToRetainIntersectLineSeg + 1);
   }

   //    assert(PtVProfileLastPart.size() > 1);
   //    assert(prVLineSegLastPart.size() > 0);

   // Truncate the truncate-profile at the point of intersection, and do the same for all its co-incident profiles. Then append the profile points of the main to-retain profile seaward from the intersection point, and do the same for the corresponding line segments (including coincident profiles)
   TruncateProfileAndAppendNew(nCoast, nProfileToTruncate, nProfileToTruncateIntersectLineSeg, &PtVProfileLastPart, &prVLineSegLastPart);

   //    assert(m_VCoast[nCoast].pGetProfile(nProfileToTruncate)->nGetProfileSize() > 1);
   //    assert(pProfileToRetain->nGetNumLineSegments() > 0);
   //    assert(m_VCoast[nCoast].pGetProfile(nProfileToTruncate)->nGetNumLineSegments() > 0);
}

//===============================================================================================================================
//! Inserts an intersection point into the profile that is to be retained, if that point is not already present in the profile, then does the same for all co-incident profiles. Finally adds the numbers of the to-truncate profile (and all its coincident profiles) to the seaward line segments of the to-retain profile and all its coincident profiles
//===============================================================================================================================
int CSimulation::nInsertPointIntoProfilesIfNeededThenUpdate(int const nCoast, int const nMainProfile, double const dIntersectX, double const dIntersectY, int const nMainProfileIntersectLineSeg, int const nProfileToTruncate, int const nProfileToTruncateIntersectLineSeg, bool const bAlreadyPresent)
{
   // START: FOR CHECKING PURPOSES ****************************************************************
   // Get the index numbers of all coincident profiles for the 'main' to-retain profile for the line segment in which intersection occurred
   //    vector<pair<int, int> > prVRetainCoincidentProfilesCHECK1 = *m_VCoast[nCoast].pGetProfile(nMainProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nMainProfileIntersectLineSeg);
   //    int nNumRetainCoincidentCHECK1 = prVRetainCoincidentProfilesCHECK1.size();
   //    for (int nn = 0; nn < nNumRetainCoincidentCHECK1; nn++)
   //    {
   //       int nThisProfile = prVRetainCoincidentProfilesCHECK1[nn].first;
   //       LogStream << "\tBEFORE nInsertPointIntoProfilesIfNeededThenUpdate(): " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-retain profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points ";
   //       for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   //          LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   //       LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles and their line segments are ";
   //       for (int mm = 0; mm < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); mm++)
   //       {
   //          vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(mm);
   //          LogStream << "{ ";
   //          for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(mm); nn++)
   //             LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   //          LogStream << "} ";
   //       }
   //       LogStream << endl;
   //
   //       for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   //       {
   //          CGeom2DPoint
   //             Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   //             Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   //          if (Pt1 == Pt2)
   //             LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   //       }
   //    }
   //    // Get the index numbers of all coincident profiles for the 'main' to-truncate profile for the line segment in which intersection occurred
   //    vector<pair<int, int> > prVTruncateCoincidentProfilesCHECK1 = *m_VCoast[nCoast].pGetProfile(nProfileToTruncate)->pprVGetPairedCoincidentProfilesForLineSegment(nProfileToTruncateIntersectLineSeg);
   //    int nNumTruncateCoincidentCHECK1 = prVTruncateCoincidentProfilesCHECK1.size();
   //    for (int nn = 0; nn < nNumTruncateCoincidentCHECK1; nn++)
   //    {
   //       int nThisProfile = prVTruncateCoincidentProfilesCHECK1[nn].first;
   //       LogStream << "\tBEFORE nInsertPointIntoProfilesIfNeededThenUpdate(): " << (nThisProfile == nProfileToTruncate ? "MAIN" : "COINCIDENT") << " to-truncate profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points ";
   //       for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   //          LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   //       LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles and their line segments are ";
   //       for (int mm = 0; mm < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); mm++)
   //       {
   //          vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(mm);
   //          LogStream << "{ ";
   //          for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(mm); nn++)
   //             LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   //          LogStream << "} ";
   //       }
   //       LogStream << endl;
   //
   //       for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   //       {
   //          CGeom2DPoint
   //             Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   //             Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   //          if (Pt1 == Pt2)
   //             LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   //       }
   //    }
   // END: FOR CHECKING PURPOSES ******************************************************************

   // Get the index numbers of all coincident profiles for the 'main' to-retain profile for the line segment in which intersection occurs
   vector<pair<int, int>> prVCoincidentProfiles = *m_VCoast[nCoast].pGetProfile(nMainProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nMainProfileIntersectLineSeg);
   int nNumCoincident = static_cast<int>(prVCoincidentProfiles.size());
   vector<int> nLineSegAfterIntersect(nNumCoincident, -1); // The line segment after the point of intersection, for each co-incident profile

   // Do this for the main profile and all profiles which are co-incident for this line segment
   for (int nn = 0; nn < nNumCoincident; nn++)
   {
      int nThisProfile = prVCoincidentProfiles[nn].first;  // The number of this profile
      int nThisLineSeg = prVCoincidentProfiles[nn].second; // The line segment of this profile
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Is the intersection point already present in the to-retain profile?
      if (! bAlreadyPresent)
      {
         // It is not already present, so insert it and also update the associated multi-line
         if (!pThisProfile->bInsertIntersection(dIntersectX, dIntersectY, nThisLineSeg))
         {
            // Error
            LogStream << WARN << m_ulIter << ": cannot insert a line segment after the final line segment (" << nThisLineSeg << ") for " << (nThisProfile == nMainProfile ? "main" : "co-incident") << " profile (" << nThisProfile << "), abandoning" << endl;

            return RTN_ERR_CANNOT_INSERT_POINT;
         }

         //          LogStream << "\tIntersection point NOT already in " << (nThisProfile == nMainProfile ? "main" : "co-incident") << " profile {" << nThisProfile << "}, inserted it as point " << nThisLineSeg+1 << endl;
      }

      // Get the line segment after intersection
      nLineSegAfterIntersect[nn] = nThisLineSeg + 1;
   }

   //    for (int nn = 0; nn < nNumCoincident; nn++)
   //       LogStream << "\tFor profile {" << prVCoincidentProfiles[nn].first << "} line segment [" << nLineSegAfterIntersect[nn] << "] is immediately after the intersection point" << endl;

   // Get the coincident profiles for the to-truncate profile, at the line segment where intersection occurs
   vector<pair<int, int>> prVToTruncateCoincidentProfiles = *m_VCoast[nCoast].pGetProfile(nProfileToTruncate)->pprVGetPairedCoincidentProfilesForLineSegment(nProfileToTruncateIntersectLineSeg);
   int nNumToTruncateCoincident = static_cast<int>(prVToTruncateCoincidentProfiles.size());

   // Now add the number of the to-truncate profile, and all its coincident profiles, to all line segments which are seaward of the point of intersection. Do this for the main profile and all profiles which are co-incident for this line segment
   for (int nn = 0; nn < nNumCoincident; nn++)
   {
      int nThisProfile = prVCoincidentProfiles[nn].first; // The number of this profile
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      // Get the number of line segments for this to-retain profile (will have just increased, if we just inserted a point)
      int nNumLineSegs = pThisProfile->nGetNumLineSegments();

      // Do for all line segments seaward of the point of intersection
      for (int nLineSeg = nLineSegAfterIntersect[nn], nIncr = 0; nLineSeg < nNumLineSegs; nLineSeg++, nIncr++)
      {
         //         LogStream << "\tnThisProfile = " << nThisProfile << " nThisLineSeg = " << nThisLineSeg << " nLineSeg = " << nLineSeg << " nNumLineSegs = " << nNumLineSegs << endl;

         // Have no idea how this can happen, but check for the moment
         //          if (nThisProfile == nProfileToTruncate)
         //          {
         //             LogStream << "\t*** ERROR nThisProfile = " << nThisProfile << " nProfileToTruncate = " << nProfileToTruncate << ", ignoring" << endl;
         //             continue;
         //          }

         // Add the number of the to-truncate profile, and all its coincident profiles, to this line segment
         for (int m = 0; m < nNumToTruncateCoincident; m++)
         {
            int nProfileToAdd = prVToTruncateCoincidentProfiles[m].first;
            int nProfileToAddLineSeg = prVToTruncateCoincidentProfiles[m].second;

            //            LogStream << "\tAdding " << (nProfileToAdd == nProfileToTruncate ? "main" : "co-incident") << " truncate-profile number {" << nProfileToAdd << "}, line segment [" << nProfileToAddLineSeg + nIncr << "] to line segment " << nLineSeg << " of " << (nThisProfile == nMainProfile ? "main" : "co-incident") << " to-retain profile {" << nThisProfile << "}" << endl;

            pThisProfile->AddCoincidentProfileToExistingLineSegment(nLineSeg, nProfileToAdd, nProfileToAddLineSeg + nIncr);
         }
      }
   }

   // START: FOR CHECKING PURPOSES ****************************************************************
   // Get the index numbers of all coincident profiles for the 'main' profile for the line segment in which intersection occurred
   //    vector<pair<int, int> > prVCoincidentProfilesCHECK2 = *m_VCoast[nCoast].pGetProfile(nMainProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nMainProfileIntersectLineSeg);
   //    int nNumCoincidentCHECK2 = prVCoincidentProfilesCHECK2.size();
   //    for (int nn = 0; nn < nNumCoincidentCHECK2; nn++)
   //    {
   //       int nThisProfile = prVCoincidentProfilesCHECK2[nn].first;
   //       LogStream << "\tAFTER nInsertPointIntoProfilesIfNeededThenUpdate(): " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-retain profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points ";
   //       for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   //          LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   //       LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles and their line segments are ";
   //       for (int nLineSeg = 0; nLineSeg < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); nLineSeg++)
   //       {
   //          vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nLineSeg);
   //          LogStream << "{ ";
   //          for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(nLineSeg); nn++)
   //             LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   //          LogStream << "} ";
   //       }
   //       LogStream << endl;
   //
   //       for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   //       {
   //          CGeom2DPoint
   //             Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   //             Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   //          if (Pt1 == Pt2)
   //             LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   //       }
   //    }
   // END: FOR CHECKING PURPOSES ******************************************************************

   return RTN_OK;
}

//===============================================================================================================================
//! Truncate a profile at the point of intersection, and do the same for all its co-incident profiles
//===============================================================================================================================
void CSimulation::TruncateProfileAndAppendNew(int const nCoast, int const nMainProfile, int const nMainProfileIntersectLineSeg, vector<CGeom2DPoint> const* pPtVProfileLastPart, vector<vector<pair<int, int>>> const* pprVLineSegLastPart)
{
   // START: FOR CHECKING PURPOSES ****************************************************************
   // Get the index numbers of all coincident profiles for the 'main' profile for the line segment in which intersection occurred
   //    vector<pair<int, int> > prVCoincidentProfilesCHECK1 = *m_VCoast[nCoast].pGetProfile(nMainProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nMainProfileIntersectLineSeg);
   //    int nNumCoincidentCHECK1 = prVCoincidentProfilesCHECK1.size();
   //
   //    LogStream << "\tTruncating profile {" << nMainProfile << "}, intersection is at [" << dIntersectX << ", " << dIntersectY << "] in line segment " << nMainProfileIntersectLineSeg << endl;
   //    for (int nn = 0; nn < nNumCoincidentCHECK1; nn++)
   //    {
   //       int nThisProfile = prVCoincidentProfilesCHECK1[nn].first;
   //       LogStream << "\tBEFORE TruncateProfileAndAppendNew(): " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-truncate profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points (";
   //       for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   //          LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   //       LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles are ";
   //       for (int nLineSeg = 0; nLineSeg < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); nLineSeg++)
   //       {
   //          vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nLineSeg);
   //          LogStream << "{ ";
   //          for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(nLineSeg); nn++)
   //             LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   //          LogStream << "} ";
   //       }
   //       LogStream << endl;
   //
   //       for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   //       {
   //          CGeom2DPoint
   //             Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   //             Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   //          if (Pt1 == Pt2)
   //             LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   //       }
   //    }
   //    LogStream << "\tPart-profile to append is ";
   //    for (int mm = 0; mm < pPtVProfileLastPart->size(); mm++)
   //       LogStream << "[" << pPtVProfileLastPart->at(mm).dGetX() << ", " << pPtVProfileLastPart->at(mm).dGetY() << "] ";
   //    LogStream << endl;
   //    LogStream << "\tPart line-segment to append is ";
   //    for (int mm = 0; mm < pprVLineSegLastPart->size(); mm++)
   //    {
   //       vector<pair<int, int> > prVTmp = pprVLineSegLastPart->at(mm);
   //       LogStream << "{ ";
   //       for (int nn = 0; nn < prVTmp.size(); nn++)
   //          LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   //       LogStream << "} ";
   //    }
   //    LogStream << endl;
   // END: FOR CHECKING PURPOSES ******************************************************************

   // Get the index numbers of all coincident profiles for the 'main' profile for the line segment in which intersection occurs
   vector<pair<int, int>> prVCoincidentProfiles = *m_VCoast[nCoast].pGetProfile(nMainProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nMainProfileIntersectLineSeg);
   int nNumCoincident = static_cast<int>(prVCoincidentProfiles.size());

   for (int nn = 0; nn < nNumCoincident; nn++)
   {
      // Do this for the main to-truncate profile, and do the same for all its co-incident profiles
      int nThisProfile = prVCoincidentProfiles[nn].first;
      int nThisProfileLineSeg = prVCoincidentProfiles[nn].second;
      CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

      //       if (nThisProfile == nMainProfile)
      //          assert(nThisProfileLineSeg == nMainProfileIntersectLineSeg);

      // Truncate the profile
      //      LogStream << "\tTruncating " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-truncate profile {" << nThisProfile << "} at line segment " << nThisProfileLineSeg+1 << endl;
      pThisProfile->TruncateProfile(nThisProfileLineSeg + 1);

      // Reduce the number of line segments for this profile
      pThisProfile->TruncateLineSegments(nThisProfileLineSeg + 1);

      // Append the profile points from the last part of the retain-profile
      for (unsigned int mm = 0; mm < pPtVProfileLastPart->size(); mm++)
      {
         CGeom2DPoint Pt = pPtVProfileLastPart->at(mm);
         pThisProfile->AppendPointInProfile(&Pt);
      }

      // Append the line segments, and their co-incident profile numbers, from the last part of the retain-profile
      for (unsigned int mm = 0; mm < pprVLineSegLastPart->size(); mm++)
      {
         vector<pair<int, int>> prVTmp = pprVLineSegLastPart->at(mm);

         pThisProfile->AppendLineSegment(&prVTmp);
      }

      // Fix the line seg numbers for this profile
      vector<int> nVProf;
      vector<int> nVProfsLineSeg;
      for (int nSeg = 0; nSeg < pThisProfile->nGetNumLineSegments(); nSeg++)
      {
         for (int nCoinc = 0; nCoinc < pThisProfile->nGetNumCoincidentProfilesInLineSegment(nSeg); nCoinc++)
         {
            int nProf = pThisProfile->nGetProf(nSeg, nCoinc);
            int nProfsLineSeg = pThisProfile->nGetProfsLineSeg(nSeg, nCoinc);

            auto it = find(nVProf.begin(), nVProf.end(), nProf);
            if (it == nVProf.end())
            {
               // Not found
               nVProf.push_back(nProf);
               nVProfsLineSeg.push_back(nProfsLineSeg);
            }
            else
            {
               // Found
               int nPos = static_cast<int>(it - nVProf.begin());
               int nNewProfsLineSeg = nVProfsLineSeg[nPos];
               nNewProfsLineSeg++;

               nVProfsLineSeg[nPos] = nNewProfsLineSeg;
               pThisProfile->SetProfsLineSeg(nSeg, nCoinc, nNewProfsLineSeg);
            }
         }
      }

      //       assert(pThisProfile->nGetProfileSize() > 1);
   }

   // START: FOR CHECKING PURPOSES ****************************************************************
   // Get the index numbers of all coincident profiles for the 'main' to-truncate profile for the line segment in which intersection occurred
   //    vector<pair<int, int> > prVToTruncateCoincidentProfilesCHECK2 = *m_VCoast[nCoast].pGetProfile(nMainProfile)->pprVGetPairedCoincidentProfilesForLineSegment(nMainProfileIntersectLineSeg);
   //    int nNumToTruncateCoincidentCHECK2 = prVToTruncateCoincidentProfilesCHECK2.size();
   //    for (int nn = 0; nn < nNumToTruncateCoincidentCHECK2; nn++)
   //    {
   //       int nThisProfile = prVToTruncateCoincidentProfilesCHECK2[nn].first;
   //       LogStream << "\tAFTER TruncateProfileAndAppendNew(): " << (nThisProfile == nMainProfile ? "MAIN" : "COINCIDENT") << " to-truncate profile {" << nThisProfile << "} has " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize() << " points (";
   //       for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize(); nn++)
   //          LogStream << "[" << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetX() << ", " << m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nn)->dGetY() << "] ";
   //       LogStream << "), and " << m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments() << " line segments, co-incident profiles are ";
   //       for (int mm = 0; mm < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumLineSegments(); mm++)
   //       {
   //          vector<pair<int, int> > prVTmp = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pprVGetPairedCoincidentProfilesForLineSegment(mm);
   //          LogStream << "{ ";
   //          for (int nn = 0; nn < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetNumCoincidentProfilesInLineSegment(mm); nn++)
   //             LogStream << prVTmp[nn].first << "[" << prVTmp[nn].second << "] ";
   //          LogStream << "} ";
   //       }
   //       LogStream << endl;
   //
   //       for (int nPoint = 0; nPoint < m_VCoast[nCoast].pGetProfile(nThisProfile)->nGetProfileSize()-1; nPoint++)
   //       {
   //          CGeom2DPoint
   //             Pt1 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint),
   //             Pt2 = *m_VCoast[nCoast].pGetProfile(nThisProfile)->pPtGetPointInProfile(nPoint+1);
   //
   //          if (Pt1 == Pt2)
   //             LogStream << m_ulIter << ": IDENTICAL POINTS before changes, in profile {" << nThisProfile << "} points = " << nPoint << " and " << nPoint+1 << endl;
   //       }
   //    }
   // END: FOR CHECKING PURPOSES ******************************************************************
}
