/*!
 *
 * \file do_shore_platform_erosion.cpp
 * \brief Erodes the consolidated sediment of the shore platform, extrapolating from erosion calculated on coastline-normal profiles. Eroded sediment from the shore platform becomes unconsolidated sediment stored in coastal polygons
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

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::setiosflags;

#include "cme.h"
#include "hermite_cubic.h"
#include "interpolate.h"
#include "simulation.h"
#include "coast.h"

//===============================================================================================================================
//! Does platform erosion on all coastlines by first calculating platform erosion on coastline-normal profiles, then extrapolating this to cells between the profiles
//===============================================================================================================================
int CSimulation::nDoAllShorePlatFormErosion(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Calculating shore platform erosion" << endl;
   
   // TODO 023 Only do potential erosion if cell is in a polygon

   // Set direction
   static bool bForward = true;

   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int const nNumProfiles = m_VCoast[nCoast].nGetNumProfiles();

      // Calculate potential erosion on every coastline-normal profile. Can do this in the original, curvature-related, sequence of profiles
      int nProfile = 0;
      for ((bForward ? nProfile = 0 : nProfile = (nNumProfiles - 1)); (bForward ? nProfile < nNumProfiles : nProfile >= 0); (bForward ? nProfile++ : nProfile--))
      {
         // Calculate potential platform erosion along the length of this profile
         int const nRet = nCalcPotentialPlatformErosionOnProfile(nCoast, nProfile);
         if (nRet != RTN_OK)
            return nRet;
      }

      // Calculate potential platform erosion between the coastline-normal profiles. Do this in along-coastline sequence
      int nProfIndex = 0;
      for ((bForward ? nProfIndex = 0 : nProfIndex = (nNumProfiles - 1)); (bForward ? nProfIndex < nNumProfiles : nProfIndex >= 0); (bForward ? nProfIndex++ : nProfIndex--))
      {
         // Calculate potential erosion for sea cells between this profile and the next profile (or up to the edge of the grid) on these cells
         int nRet = nCalcPotentialPlatformErosionBetweenProfiles(nCoast, nProfIndex, DIRECTION_DOWNCOAST);
         if (nRet != RTN_OK)
            return nRet;

         nRet = nCalcPotentialPlatformErosionBetweenProfiles(nCoast, nProfIndex, DIRECTION_UPCOAST);
         if (nRet != RTN_OK)
            return nRet;
      }
   }

   // If desired, swap direction for next timestep
   if (m_bErodeShorePlatformAlternateDirection)
      bForward = ! bForward;

   // Fills in 'holes' in the potential platform erosion i.e. orphan cells which get omitted because of rounding problems
   FillPotentialPlatformErosionHoles();

   // Do the same for beach protection
   FillInBeachProtectionHoles();

   // Finally calculate actual platform erosion on all sea cells (both on profiles, and between profiles)
   for (int nX = 0; nX < m_nXGridMax; nX++)
   {
      for (int nY = 0; nY < m_nYGridMax; nY++)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bPotentialPlatformErosion())
            // Calculate actual (supply-limited) shore platform erosion on each cell that has potential platform erosion, also add the eroded sand/coarse sediment to that cells's polygon, ready to be redistributed within the polygon during beach erosion/deposition
            DoActualPlatformErosionOnCell(nX, nY);
      }
   }

   if (m_nLogFileDetail >= LOG_FILE_ALL)
   {
      LogStream << m_ulIter << ": total potential shore platform erosion (m^3) = " << m_dThisIterPotentialPlatformErosion * m_dCellArea << " (on profiles = " << m_dTotPotentialPlatformErosionOnProfiles * m_dCellArea << ", between profiles = " << m_dTotPotentialPlatformErosionBetweenProfiles * m_dCellArea << ")" << endl;

      LogStream << m_ulIter << ": total actual shore platform erosion (m^3) = " << (m_dThisIterActualPlatformErosionFineCons + m_dThisIterActualPlatformErosionSandCons + m_dThisIterActualPlatformErosionCoarseCons) * m_dCellArea << " (fine = " << m_dThisIterActualPlatformErosionFineCons * m_dCellArea << ", sand = " << m_dThisIterActualPlatformErosionSandCons * m_dCellArea << ", coarse = " << m_dThisIterActualPlatformErosionCoarseCons * m_dCellArea << ")" << endl << endl;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Calculates potential (i.e. unconstrained by available sediment) erosional lowering of the shore platform for a single coastline-normal profile, due to wave action. This routine uses a behavioural rule to modify the original surface elevation profile geometry, in which erosion rate/slope = f(d/Db) based on Walkden & Hall (2005). Originally coded in Matlab by Andres Payo
//===============================================================================================================================
int CSimulation::nCalcPotentialPlatformErosionOnProfile(int const nCoast, int const nProfile)
{
   CGeomProfile *const pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // Only work on this profile if it is problem-free TODO 024 Or if it has just hit dry land?
   if (!pProfile->bOKIncStartAndEndOfCoast()) //  || (pProfile->nGetProblemCode() == PROFILE_DRYLAND))
      return RTN_OK;

   // Get the length of the profile (in cells) and the index of the coast point at which this profile starts
   int const
       nProfSize = pProfile->nGetNumCellsInProfile(),
       nCoastPoint = pProfile->nGetNumCoastPoint();

   // Get the breaking depth for this profile from the coastline point
   double dDepthOfBreaking = m_VCoast[nCoast].dGetDepthOfBreaking(nCoastPoint);
   if (bFPIsEqual(dDepthOfBreaking, DBL_NODATA, TOLERANCE))
      // This profile is not in the active zone, so no platform erosion here
      return RTN_OK;
   
   if (bFPIsEqual(dDepthOfBreaking, 0.0, TOLERANCE))
   {
      // Safety check, altho' this shouldn't happen
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ": depth of breaking is zero on coast " << nCoast << " at profile " << nProfile << ", continuing however." << endl;
      
      return RTN_OK;   
   }

   //   LogStream << m_ulIter << ": ON PROFILE nProfile = " << nProfile << " dDepthOfBreaking = " << dDepthOfBreaking << endl;

   // Get the height of the associated breaking wave from the coast point: this height is used in beach protection calcs
   double const dBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nCoastPoint);

   // Calculate the length of the profile in external CRS units
   int const nSegments = pProfile->nGetProfileSize() - 1;
   double dProfileLenXY = 0;
   for (int nSeg = 0; nSeg < nSegments; nSeg++)
   {
      // Do once for every line segment
      double const
          dSegStartX = pProfile->pPtGetPointInProfile(nSeg)->dGetX(),
          dSegStartY = pProfile->pPtGetPointInProfile(nSeg)->dGetY(),
          dSegEndX = pProfile->pPtGetPointInProfile(nSeg + 1)->dGetX(), // Is OK
          dSegEndY = pProfile->pPtGetPointInProfile(nSeg + 1)->dGetY();

      double const dSegLen = hypot(dSegStartX - dSegEndX, dSegStartY - dSegEndY);
      dProfileLenXY += dSegLen;
   }

   // Next calculate the average distance between profile points, again in external CRS units. Assume that the sample points are equally spaced along the profile (not quite true)
   double const dSpacingXY = dProfileLenXY / (nProfSize - 1);

   // Set up vectors for the coastline-normal profile elevations. The length of this vector line is given by the number of cells 'under' the profile. Thus each point on the vector relates to a single cell in the grid. This assumes that all points on the profile vector are equally spaced (not quite true, depends on the orientation of the line segments which comprise the profile)
   // The elevation of each of these profile points is the elevation of the centroid of the cell that is 'under' the point. However we cannot always be confident that this is the 'true' elevation of the point on the vector since (unless the profile runs planview N-S or W-E) the vector does not always run exactly through the centroid of the cell
   vector<double>
       VdProfileZ(nProfSize, 0),      // Initial (pre-erosion) elevation of both consolidated and unconsolidated sediment for cells 'under' the profile
       VdProfileDistXY(nProfSize, 0), // Along-profile distance measured from the coast, in external CRS units
       dVConsProfileZ(nProfSize, 0),  // Initial (pre-erosion) elevation of consolidated sediment only for cells 'under' the profile
       dVConsZDiff(nProfSize, 0),
       dVConsSlope(nProfSize, 0);

   for (int i = 0; i < nProfSize; i++)
   {
      int const
          nX = pProfile->pPtiVGetCellsInProfile()->at(i).nGetX(),
          nY = pProfile->pPtiVGetCellsInProfile()->at(i).nGetY();

      // Get the number of the highest layer with non-zero thickness
      int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

      // Safety check
      if (nTopLayer == INT_NODATA)
         return RTN_ERR_NO_TOP_LAYER;

      if (nTopLayer == NO_NONZERO_THICKNESS_LAYERS)
         // TODO 025 We are down to basement
         return RTN_OK;

      // Get the elevation for consolidated sediment only on this cell
      dVConsProfileZ[i] = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopForLayerAboveBasement(nTopLayer);

      // Get the elevation for both consolidated and unconsolidated sediment on this cell
      VdProfileZ[i] = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();

      // And store the X-Y plane distance from the start of the profile
      VdProfileDistXY[i] = i * dSpacingXY;
   }

   for (int i = 0; i < nProfSize - 1; i++)
   {
      // For the consolidated-only profile, get the Z differences (already in external CRS units)
      dVConsZDiff[i] = dVConsProfileZ[i] - dVConsProfileZ[i + 1];

      // Calculate dZ/dXY, the Z slope (i.e. rise over run) in the XY direction. Note that we use the elevation difference on the seaward side of 'this' point
      dVConsSlope[i] = dVConsZDiff[i] / dSpacingXY;
   }

   // Sort out the final value
   dVConsSlope[nProfSize - 1] = dVConsSlope[nProfSize - 2];

   if (m_nProfileSmoothWindow > 0)
   {
      // Smooth the vector of slopes for the consolidated-only profile
      dVConsSlope = dVSmoothProfileSlope(&dVConsSlope);
   }

   vector<double>
       dVProfileDepthOverDB(nProfSize, 0),      // Depth over wave breaking depth at the coastline-normal sample points
       dVProfileErosionPotential(nProfSize, 0); // Erosion potential at the coastline-normal sample points

   // Calculate the erosion potential along this profile using the shape function
   double dTotalErosionPotential = 0;
   for (int i = 0; i < nProfSize; i++)
   {
      // Use the actual depth of water here (i.e. the depth to the top of the unconsolidated sediment, including the thickness of consolidated sediment beneath it)
      dVProfileDepthOverDB[i] = m_dThisIterSWL - VdProfileZ[i];
      dVProfileDepthOverDB[i] /= dDepthOfBreaking;

      // Constrain dDepthOverDB[i] to be between 0 (can get small -ve values due to rounding errors) and m_dDepthOverDBMax
      dVProfileDepthOverDB[i] = tMax(dVProfileDepthOverDB[i], 0.0);
      dVProfileDepthOverDB[i] = tMin(dVProfileDepthOverDB[i], m_dDepthOverDBMax);

      // And then use the look-up table to find the value of erosion potential at this point on the profile
      dVProfileErosionPotential[i] = dLookUpErosionPotential(dVProfileDepthOverDB[i]);

      // If erosion potential (a -ve value) is tiny, set it to zero
      if (dVProfileErosionPotential[i] > -SEDIMENT_ELEV_TOLERANCE)
         dVProfileErosionPotential[i] = 0;

      // Keep track of the total erosion potential for this profile
      dTotalErosionPotential += dVProfileErosionPotential[i];
   }

   // Constrain erosion potential at every point on the profile, so that the integral of erosion potential on the whole profile is unity (Walkden and Hall 2005). Note that here, erosion potential is -ve so we must constrain to -1
   for (int i = 0; i < nProfSize; i++)
   {
      if (dTotalErosionPotential < 0)
         dVProfileErosionPotential[i] /= (-dTotalErosionPotential);
   }

   vector<double>
       dVRecessionXY(nProfSize, 0),
       dVSCAPEXY(nProfSize, 0);

   // Calculate recession at every point on the coastline-normal profile
   for (int i = 0; i < nProfSize; i++)
   {
      // dRecession = dForce * (dBeachProtection / dR) * dErosionPotential * dSlope * dTime
      // where:
      //     dVRecession [m] is the landward migration distance defined in the profile relative (XY) CRS
      //     dForce is given by Equation 4 in Walkden & Hall, 2005
      //     dVBeachProtection [1] is beach protection factor [1, 0] = [no protection, fully protected]. (This is calculated later, see dCalcBeachProtectionFactor())
      //     dVR  [m^(9/4)s^(2/3)] is the material strength and some hydrodynamic constant
      //     dVProfileErosionPotential [?] is the erosion potential at each point along the profile
      //     dVSlope [1] is the along-profile slope
      //     m_dTimeStep * 3600 [s] is the time interval in seconds
      //
      // dRecession is horizontal recession (along the XY direction):
      //
      // dVRecessionXY[i] = (dForce * dVBeachProtection[i] * dVErosionPotentialFunc[i] * dVSlope[i] * m_dTimeStep * 3600) / dVR[i]
      //
      // XY recession must be -ve or zero. If it is +ve then it represents accretion not erosion, which must be described by a different set of equations. So we also need to constrain XY recession to be <= 0
      dVRecessionXY[i] = tMin(m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nCoastPoint) * dVProfileErosionPotential[i] * dVConsSlope[i] / m_dR, 0.0);
      dVSCAPEXY[i] = VdProfileDistXY[i] - dVRecessionXY[i];
   }

   vector<double> dVChangeElevZ(nProfSize, 0);

   // We have calculated the XY-plane recession at every point on the profile, so now convert this to a change in Z-plane elevation at every inundated point on the profile (not the coast point). Again we use the elevation difference on the seaward side of 'this' point
   for (int i = 1; i < nProfSize - 1; i++)
   {
      // Vertical lowering dZ = dXY * tan(a), where tan(a) is the slope of the SCAPE profile in the XY direction
      double
          dSCAPEHorizDist = dVSCAPEXY[i + 1] - dVSCAPEXY[i],
          dSCAPEVertDist = dVConsProfileZ[i] - dVConsProfileZ[i + 1],
          dSCAPESlope = dSCAPEVertDist / dSCAPEHorizDist,
          dDeltaZ = dVRecessionXY[i] * dSCAPESlope;
      // Safety check: if thickness model has some jumps, dVConsProfileZ might be very high, limiting dSCAPESlope to 0 because all time erode a high fix quantity
      if (dSCAPESlope > 1)
      {
         dDeltaZ = 0;
      }

      int const
          nX = pProfile->pPtiVGetCellsInProfile()->at(i).nGetX(),
          nY = pProfile->pPtiVGetCellsInProfile()->at(i).nGetY();

      // Store the local slope of the consolidated sediment, this is just for output display purposes
      m_pRasterGrid->m_Cell[nX][nY].SetLocalConsSlope(dVConsSlope[i]);

      // dDeltaZ is zero or -ve: if dDeltaZ is zero then do nothing, if -ve then remove some sediment from this cell
      if (dDeltaZ < 0)
      {
         // If there has already been potential erosion on this cell, then it must be a shared line segment (i.e. has co-incident profiles)
         double dPrevPotentialErosion = -m_pRasterGrid->m_Cell[nX][nY].dGetPotentialPlatformErosion();
         if (dPrevPotentialErosion < 0)
         {
            // Average the two values
            //               LogStream << m_ulIter << ": [" << nX << "][" << nY << "] under profile " << nProfile << " has previous potential platform erosion = " << dPrevPotentialErosion << endl;
            dDeltaZ = ((dDeltaZ + dPrevPotentialErosion) / 2);
         }

         // Constrain the lowering so we don't get negative slopes or +ve erosion amounts (dDeltaZ must be -ve), this is implicit in SCAPE
         dDeltaZ = tMax(dDeltaZ, -dVConsZDiff[i]);
         dDeltaZ = tMin(dDeltaZ, 0.0);
         dVChangeElevZ[i] = dDeltaZ;

         // Set the potential (unconstrained) erosion for this cell, is a +ve value
         m_pRasterGrid->m_Cell[nX][nY].SetPotentialPlatformErosion(-dDeltaZ);

         // Update this-timestep totals
         m_ulThisIterNumPotentialPlatformErosionCells++;
         m_dThisIterPotentialPlatformErosion -= dDeltaZ; // Since dDeltaZ is a -ve value
//          assert(isfinite(m_dThisIterPotentialPlatformErosion));
//          assert(m_dThisIterPotentialPlatformErosion >= 0);

         // Increment the check values
         m_ulTotPotentialPlatformErosionOnProfiles++;
         m_dTotPotentialPlatformErosionOnProfiles -= dDeltaZ;
      }

      // Finally, calculate the beach protection factor, this will be used in estimating actual (supply-limited) erosion
      double dBeachProtectionFactor = dCalcBeachProtectionFactor(nX, nY, dBreakingWaveHeight);
      m_pRasterGrid->m_Cell[nX][nY].SetBeachProtectionFactor(dBeachProtectionFactor);
   }

   // If desired, save this coastline-normal profile data for checking purposes
   if (m_bOutputProfileData)
   {
      int nRet = nSaveProfile(nProfile, nCoast, nProfSize, &VdProfileDistXY, &dVConsProfileZ, &dVProfileDepthOverDB, &dVProfileErosionPotential, &dVConsSlope, &dVRecessionXY, &dVChangeElevZ, pProfile->pPtiVGetCellsInProfile(), &dVSCAPEXY);
      if (nRet != RTN_OK)
         return nRet;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Calculates potential platform erosion on cells to one side of a given coastline-normal profile, up to the next profile
//===============================================================================================================================
int CSimulation::nCalcPotentialPlatformErosionBetweenProfiles(int const nCoast, int const nProfIndex, int const nDirection)
{
   // Get the number of the coastline-normal profile
   int const nProfile = m_VCoast[nCoast].nGetProfileFromAlongCoastProfileIndex(nProfIndex);
   CGeomProfile *const pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

   // Only work on this profile if it is problem-free
   if (! pProfile->bOKIncStartAndEndOfCoast())
      return RTN_OK;

   int const
       nProfSize = pProfile->nGetNumCellsInProfile(),
       nCoastProfileStart = pProfile->nGetNumCoastPoint(),
       nProfileStartX = pProfile->pPtiVGetCellsInProfile()->at(0).nGetX(),
       nProfileStartY = pProfile->pPtiVGetCellsInProfile()->at(0).nGetY(),
       nCoastMax = m_VCoast[nCoast].nGetCoastlineSize();
   int
       nDistFromProfile = 0,
       nParCoastXLast = nProfileStartX,
       nParCoastYLast = nProfileStartY;

   // Start at the coast end of this coastline-normal profile, then move one cell forward along the coast, then construct a parallel profile from this new coastline start cell. Calculate erosion along this parallel profile in the same way as above. Move another cell forward along the coastline, do the same. Keep going until hit another profile
   while (true)
   {
      // Increment the distance from the start coast-normal profile
      nDistFromProfile++;

      // Start the parallel profile nDistFromProfile cells along the coastline from the coastline-normal profile, direction depending on nDirection
      int nThisPointOnCoast = nCoastProfileStart;
      if (nDirection == DIRECTION_DOWNCOAST)
         nThisPointOnCoast += nDistFromProfile;
      else
         nThisPointOnCoast -= nDistFromProfile;

      // Have we reached the beginning of the coast?
      if ((nDirection == DIRECTION_UPCOAST) && (nThisPointOnCoast < 0))
      {
//         LogStream << m_ulIter << ": LEAVING LOOP since hit nThisPointOnCoast = " << nThisPointOnCoast << " while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << nProfile << ", dist from profile = " <<  nDistFromProfile << endl;

         break;
      }

      // Have we reached the end of the coast?
      if ((nDirection == DIRECTION_DOWNCOAST) && (nThisPointOnCoast >= nCoastMax))
      {
         //         LogStream << m_ulIter << ": LEAVING LOOP since hit nThisPointOnCoast = " << nThisPointOnCoast << " while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << nProfile << ", dist from profile = " <<  nDistFromProfile << endl;

         break;
      }

      //       LogStream << m_ulIter << ": from profile " << nProfile << ", doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast, dist from profile = " <<  nDistFromProfile << endl;

      double dDepthOfBreaking = m_VCoast[nCoast].dGetDepthOfBreaking(nThisPointOnCoast);
      if (bFPIsEqual(dDepthOfBreaking, DBL_NODATA, TOLERANCE))
      {
         // This parallel profile is not in the active zone, so no platform erosion here. Move on to the next point along the coastline in this direction
         if (m_nLogFileDetail == LOG_FILE_ALL)
            LogStream << m_ulIter << ": not in active zone at coastline " << nCoast << " coast point " << nThisPointOnCoast << " when constructing parallel profile for potential platform erosion. Working from profile " << nProfile << ", " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast, dist from profile = " << nDistFromProfile << endl;

         continue;
      }

//      LogStream << m_ulIter << ": BETWEEN PROFILES " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile " << nProfile << " nThisPointOnCoast = " << nThisPointOnCoast << " dDepthOfBreaking = " << dDepthOfBreaking << " nParProfSize = " << nParProfSize << endl;

      // All is OK, so get the grid co-ordinates of this point, which is the coastline start point for the parallel profile
      int const nParCoastX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPointOnCoast)->nGetX();
      int const nParCoastY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPointOnCoast)->nGetY();

      if ((nParCoastX == nParCoastXLast) && (nParCoastY == nParCoastYLast))
      {
         // Should not happen, but could do due to rounding errors
         if (m_nLogFileDetail >= LOG_FILE_ALL)
            LogStream << WARN << m_ulIter << ": rounding problem on coast " << nCoast << " profile " << nProfile << " at [" << nParCoastX << "][" << nParCoastY << "]" << endl;

         // So move on to the next point along the coastline in this direction
         continue;
      }

      // Is this coastline start point the start point of an adjacent coastline-normal vector?
      if (m_pRasterGrid->m_Cell[nParCoastX][nParCoastY].bIsProfile())
      {
//         LogStream << m_ulIter << ": coast " << nCoast << ", LEAVING LOOP since hit another profile at nThisPointOnCoast = " << nThisPointOnCoast << " while doing potential platform erosion " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile = " << nProfile << ", dist from profile = " <<  nDistFromProfile << endl;
         break;
      }

      // Get the height of the associated breaking wave from the coast point: this height is used in beach protection calcs. Note that it will be DBL_NODATA if not in active zone
      double const dBreakingWaveHeight = m_VCoast[nCoast].dGetBreakingWaveHeight(nThisPointOnCoast);

      // OK, now construct a parallel profile
      vector<CGeom2DIPoint> PtiVGridParProfile; // Integer coords (grid CRS) of cells under the parallel profile
      vector<CGeom2DPoint> PtVExtCRSParProfile; // Co-ordinates (external CRS) of cells under the parallel profile

      ConstructParallelProfile(nProfileStartX, nProfileStartY, nParCoastX, nParCoastY, nProfSize, pProfile->pPtiVGetCellsInProfile(), &PtiVGridParProfile, &PtVExtCRSParProfile);

      int const nParProfSize = static_cast<int>(PtiVGridParProfile.size());
      // We have a parallel profile which starts at the coast, but is it long enough to be useful? May have been cut short because it extended outside the grid, or we hit an adjacent profile
      if (nParProfSize < 3)
      {
         // We cannot use this parallel profile, it is too short to calculate along-profile slope, so abandon it and move on to the next parallel profile in this direction
         nParCoastXLast = nParCoastX;
         nParCoastYLast = nParCoastY;
         //         LogStream << m_ulIter << ": parallel profile abandoned since too short, starts at [" << nParCoastX << "][" << nParCoastY << "] coastline point " << nThisPointOnCoast << ", length = " << nParProfSize << endl;
         continue;
      }

      // This parallel profile is OK, so calculate potential erosion along it. First calculate the length of the parallel profile in external CRS units
      double const dParProfileLenXY = dGetDistanceBetween(&PtVExtCRSParProfile[0], &PtVExtCRSParProfile[nParProfSize - 1]);

      // Next calculate the distance between profile points, again in external CRS units. Assume that the sample points are equally spaced along the parallel profile (not quite true)
      double const dParSpacingXY = dParProfileLenXY / (nParProfSize - 1);
      //      LogStream << "dParSpacingXY = " << dParSpacingXY << endl;

      vector<double> dVParProfileZ(nParProfSize, 0);      // Initial (pre-erosion) elevation of both consolidated and unconsolidated sediment for cells 'under' the parallel profile
      vector<double> dVParProfileDistXY(nParProfSize, 0); // Along-profile distance measured from the coast, in external CRS units
      vector<double> dVParConsProfileZ(nParProfSize, 0);  // Initial (pre-erosion) elevation of consolidated sediment only for cells 'under' the parallel profile
      vector<double> dVParConsZDiff(nParProfSize, 0);
      vector<double> dVParConsSlope(nParProfSize, 0);

      for (int i = 0; i < nParProfSize; i++)
      {
         int const nXPar = PtiVGridParProfile[i].nGetX();
         int const nYPar = PtiVGridParProfile[i].nGetY();

         // Is this a sea cell?
         if (! m_pRasterGrid->m_Cell[nXPar][nYPar].bIsInundated())
         {
            // It isn't so move along, nothing to do here
            //             LogStream << m_ulIter << " : [" << nXPar << "][" << nYPar << "] is not inundated" << endl;
            continue;
         }

         // Is this cell in a polygon?
         int nPolyID = m_pRasterGrid->m_Cell[nXPar][nYPar].nGetPolygonID();
         if (nPolyID == INT_NODATA)
         {
            // It isn't. This can happen at the seaward end of polygons TODO 026 Is it a problem?
            //             LogStream << m_ulIter << " : [" << nXPar << "][" << nYPar << "] = {" << dGridCentroidXToExtCRSX(nXPar) << ", " << dGridCentroidYToExtCRSY(nYPar) << "} is not in a polygon" << endl;
            continue;
         }

         // Get the number of the highest layer with non-zero thickness
         int const nTopLayer = m_pRasterGrid->m_Cell[nXPar][nYPar].nGetTopNonZeroLayerAboveBasement();

         // Safety check
         if (nTopLayer == INT_NODATA)
            return RTN_ERR_NO_TOP_LAYER;

         if (nTopLayer == NO_NONZERO_THICKNESS_LAYERS)
            // TODO 025 We are down to basement
            return RTN_OK;

         // Get the elevation for consolidated sediment only on this cell
         dVParConsProfileZ[i] = m_pRasterGrid->m_Cell[nXPar][nYPar].dGetConsSedTopForLayerAboveBasement(nTopLayer);

         // Get the elevation for both consolidated and unconsolidated sediment on this cell
         dVParProfileZ[i] = m_pRasterGrid->m_Cell[nXPar][nYPar].dGetSedimentTopElev();

         // And store the X-Y plane distance from the start of the profile
         dVParProfileDistXY[i] = i * dParSpacingXY;
      }

      for (int i = 0; i < nParProfSize - 1; i++)
      {
         // For the consolidated-only profile, get the Z differences (already in external CRS units)
         dVParConsZDiff[i] = dVParConsProfileZ[i] - dVParConsProfileZ[i + 1];

         // Calculate dZ/dXY, the Z slope (i.e. rise over run) in the XY direction. Note that we use the elevation difference on the seaward side of 'this' point
         dVParConsSlope[i] = dVParConsZDiff[i] / dParSpacingXY;
      }

      // Sort out the final slope value
      dVParConsSlope[nParProfSize - 1] = dVParConsSlope[nParProfSize - 2];

      if (m_nProfileSmoothWindow > 0)
      {
         // Smooth the vector of slopes for the consolidated-only profile
         dVParConsSlope = dVSmoothProfileSlope(&dVParConsSlope);
      }

      // Initialize the parallel profile vector with depth / m_dWaveBreakingDepth
      vector<double> dVParProfileDepthOverDB(nParProfSize, 0);      // Depth / wave breaking depth at the parallel profile sample points
      vector<double> dVParProfileErosionPotential(nParProfSize, 0); // Erosion potential at the parallel profile sample points

      // Calculate the erosion potential along this profile using the shape function which we read in earlier
      double dTotalErosionPotential = 0;

      // Safety check TODO 061 Is this safety check to depth of breaking a reasonable thing to do?
      if (dDepthOfBreaking <= 0.0)
         dDepthOfBreaking = 1e-10;

      for (int i = 0; i < nParProfSize; i++)
      {
         // Use the actual depth of water here (i.e. the depth to the top of the unconsolidated sediment, including the thickness of consolidated sediment beneath it)
         dVParProfileDepthOverDB[i] = m_dThisIterSWL - dVParProfileZ[i];
         dVParProfileDepthOverDB[i] /= dDepthOfBreaking;

         // Constrain dDepthOverDB[i] to be between 0 (can get small -ve due to rounding errors) and m_dDepthOverDBMax
         dVParProfileDepthOverDB[i] = tMax(dVParProfileDepthOverDB[i], 0.0);
         dVParProfileDepthOverDB[i] = tMin(dVParProfileDepthOverDB[i], m_dDepthOverDBMax);

         // And then use the look-up table to find the value of erosion potential at this point on the profile
         dVParProfileErosionPotential[i] = dLookUpErosionPotential(dVParProfileDepthOverDB[i]);

         // If erosion potential (a -ve value) is tiny, set it to zero
         if (dVParProfileErosionPotential[i] > -SEDIMENT_ELEV_TOLERANCE)
            dVParProfileErosionPotential[i] = 0;

         // Keep track of the total erosion potential for this profile
         dTotalErosionPotential += dVParProfileErosionPotential[i];
      }

      // Constrain erosion potential at every point on the profile, so that the integral of erosion potential on the whole profile is unity (Walkden and Hall 2005). Note that here, erosion potential is -ve so we must constrain to -1
      for (int i = 0; i < nParProfSize; i++)
      {
         if (dTotalErosionPotential < 0)
            dVParProfileErosionPotential[i] /= (-dTotalErosionPotential);
      }

      vector<double>
          dVParRecessionXY(nParProfSize, 0),
          dVParSCAPEXY(nParProfSize, 0);

      // Calculate recession at every point on the parallel profile
      for (int i = 0; i < nParProfSize; i++)
      {
         // dRecession = dForce * (dBeachProtection / dR) * dErosionPotential * dSlope * dTime
         // where:
         //     dVRecession [m] is the landward migration distance defined in the profile relative (XY) CRS
         //     dForce is given by Equation 4 in Walkden & Hall, 2005
         //     dVBeachProtection [1] is beach protection factor [1, 0] = [no protection, fully protected] (This is calculated later, see dCalcBeachProtectionFactor())
         //     dVR  [m^(9/4)s^(2/3)] is the material strength and some hydrodynamic constant
         //     dVProfileErosionPotential [?] is the erosion potential at each point along the profile
         //     dVSlope [1] is the along-profile slope
         //     m_dTimeStep * 3600 [s] is the time interval in seconds
         //
         // dRecession is horizontal recession (along the XY direction):
         //
         // dVRecessionXY[i] = (dForce * dVBeachProtection[i] * dVErosionPotentialFunc[i] * dVSlope[i] * m_dTimeStep * 3600) / dVR[i]
         //
         //
         // XY recession must be -ve or zero. If it is +ve then it represents accretion not erosion, which must be described by a different set of equations. So we also need to constrain XY recession to be <= 0
         dVParRecessionXY[i] = tMin(m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nThisPointOnCoast) * dVParProfileErosionPotential[i] * dVParConsSlope[i] / m_dR, 0.0);
         dVParSCAPEXY[i] = dVParProfileDistXY[i] - dVParRecessionXY[i];

         //         LogStream << m_ulIter << ": [" << nXPar << "][" << nYPar << "] = {" << dGridCentroidXToExtCRSX(nXPar) << ", " <<  dGridCentroidYToExtCRSY(nYPar) << "} wave energy = " << m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nThisPointOnCoast) << " erosion potential = " << dVParProfileErosionPotential[i] << " slope = " << dVParProfileSlope[i] << " dVParZDiff[i] = " << dVParZDiff[i] << " nParProfSize = " << nParProfSize << endl;
      }

      vector<double> dVParDeltaZ(nParProfSize, 0);

      // We have calculated the XY-plane recession at every point on the profile, so now convert this to a change in Z-plane elevation at every inundated point on the profile (not the coast point). Again we use the elevation difference on the seaward side of 'this' point
      for (int i = 1; i < nParProfSize - 1; i++)
      {
         // Vertical lowering dZ = dXY * tan(a), where tan(a) is the slope of the SCAPE profile in the XY direction
         double dSCAPEHorizDist = dVParSCAPEXY[i + 1] - dVParSCAPEXY[i];
         
         // Safety check
         if (bFPIsEqual(dSCAPEHorizDist, 0.0, TOLERANCE))
            continue;
             
         double
             dSCAPEVertDist = dVParConsProfileZ[i] - dVParConsProfileZ[i + 1],
             dSCAPESlope = dSCAPEVertDist / dSCAPEHorizDist,
             dDeltaZ = dVParRecessionXY[i] * dSCAPESlope;
             // Safety check: if thickness model has some jumps, dVConsProfileZ might be very high, limiting dSCAPESlope to 0 because all time erode a high fix quantity
         if (dSCAPESlope > 1)
         {
            dDeltaZ = 0;
         }

         int const
             nXPar = PtiVGridParProfile[i].nGetX(),
             nYPar = PtiVGridParProfile[i].nGetY();

         // Store the local slope of the consolidated sediment, this is just for output display purposes
         m_pRasterGrid->m_Cell[nXPar][nYPar].SetLocalConsSlope(dVParConsSlope[i]);

         // dDeltaZ is zero or -ve: if dDeltaZ is zero then do nothing, if -ve then remove some sediment from this cell
         if (dDeltaZ < 0)
         {
            // Has this cell already been eroded during this timestep?
            if (m_pRasterGrid->m_Cell[nXPar][nYPar].bPotentialPlatformErosion())
            {
               // It has
               double const dPrevPotentialErosion = -m_pRasterGrid->m_Cell[nXPar][nYPar].dGetPotentialPlatformErosion();

               //                  LogStream << m_ulIter << ": [" << nXPar << "][" << nYPar << "] parallel profile " << nDistFromProfile << " coast points " << (nDirection == DIRECTION_DOWNCOAST ? "down" : "up") << "-coast from profile " << nProfile << " has previous potential platform erosion = " << dPrevPotentialErosion << ", current potential platform erosion = " << dDeltaZ << ", max value = " << tMin(dPrevPotentialErosion, dDeltaZ) << endl;

               // Use the larger of the two -ve values
               dDeltaZ = tMin(dPrevPotentialErosion, dDeltaZ);

               // Adjust this-timestep totals, since this cell has already been eroded
               m_ulThisIterNumPotentialPlatformErosionCells--;
               m_dThisIterPotentialPlatformErosion += dPrevPotentialErosion; // Since dPrevPotentialErosion is +ve
                                                                             //                assert(isfinite(m_dThisIterPotentialPlatformErosion));
                                                                             //                assert(m_dThisIterPotentialPlatformErosion >= 0);

               // And also adjust the check values
               m_ulTotPotentialPlatformErosionBetweenProfiles--;
               m_dTotPotentialPlatformErosionBetweenProfiles += dPrevPotentialErosion; // Since -ve
            }

            // Constrain the lowering so we don't get negative slopes or +ve erosion amounts (dDeltaZ must be -ve), this is implicit in SCAPE
            dDeltaZ = tMax(dDeltaZ, -dVParConsZDiff[i]);
            dDeltaZ = tMin(dDeltaZ, 0.0);
            dVParDeltaZ[i] = dDeltaZ;

            // Set the potential (unconstrained) erosion for this cell, it is a +ve value
            m_pRasterGrid->m_Cell[nXPar][nYPar].SetPotentialPlatformErosion(-dDeltaZ);
            //               LogStream << "[" << nXPar << "][" << nYPar << "] = {" << dGridCentroidXToExtCRSX(nXPar) << ", " <<  dGridCentroidYToExtCRSY(nYPar) << "} has potential platform erosion = " << -dDeltaZ << endl;

            // Update this-timestep totals
            m_ulThisIterNumPotentialPlatformErosionCells++;
            m_dThisIterPotentialPlatformErosion -= dDeltaZ; // Since dDeltaZ is a -ve value
                                                            //             assert(isfinite(m_dThisIterPotentialPlatformErosion));
                                                            //             assert(m_dThisIterPotentialPlatformErosion >= 0);

            // Increment the check values
            m_ulTotPotentialPlatformErosionBetweenProfiles++;
            m_dTotPotentialPlatformErosionBetweenProfiles -= dDeltaZ; // Since -ve
         }

         // Finally, calculate the beach protection factor, this will be used in estimating actual (supply-limited) erosion
         double const dBeachProtectionFactor = dCalcBeachProtectionFactor(nXPar, nYPar, dBreakingWaveHeight);
         m_pRasterGrid->m_Cell[nXPar][nYPar].SetBeachProtectionFactor(dBeachProtectionFactor);
      }

      // If desired, save this parallel coastline-normal profile for checking purposes
      if (m_bOutputParallelProfileData)
      {
         int const nRet = nSaveParProfile(nCoast, nProfile, nParProfSize, nDirection, nDistFromProfile, &dVParProfileDistXY, &dVParConsProfileZ, &dVParProfileDepthOverDB, &dVParProfileErosionPotential, &dVParConsSlope, &dVParRecessionXY, &dVParDeltaZ, pProfile->pPtiVGetCellsInProfile(), &dVParSCAPEXY);
         if (nRet != RTN_OK)
            return nRet;
      }

      // Update for next time round the loop
      nParCoastXLast = nParCoastX;
      nParCoastYLast = nParCoastY;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Calculates actual (constrained by available sediment) erosion of the consolidated shore platform on a single sea cell
//===============================================================================================================================
void CSimulation::DoActualPlatformErosionOnCell(int const nX, int const nY)
{
//    LogStream << m_ulIter << ": doing platform erosion on cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

   // Get the beach protection factor, which quantifies the extent to which unconsolidated sediment on the shore platform (beach) protects the shore platform
   double const dBeachProtectionFactor = m_pRasterGrid->m_Cell[nX][nY].dGetBeachProtectionFactor();
   if (bFPIsEqual(dBeachProtectionFactor, 0.0, TOLERANCE))
      // The beach is sufficiently thick to prevent any platform erosion (or we are down to basement)
      return;

   // Get the potential depth of potential erosion, considering beach protection
   double dThisPotentialErosion = m_pRasterGrid->m_Cell[nX][nY].dGetPotentialPlatformErosion() * dBeachProtectionFactor;

   // We will be eroding the topmost layer that has non-zero thickness
   int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

   // Safety check
   if (nThisLayer == INT_NODATA)
   {
      cerr << ERR << "no sediment layer in DoActualPlatformErosionOnCell()" << endl;
      return;
   }

   if (nThisLayer == NO_NONZERO_THICKNESS_LAYERS)
      // No layer with non-zero thickness left, we are down to basement
      return;

   // OK, we have a layer that can be eroded so find out how much consolidated sediment we have available on this cell
   double
       dExistingAvailableFine = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetFineDepth(),
       dExistingAvailableSand = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetSandDepth(),
       dExistingAvailableCoarse = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetCoarseDepth();

   // Now partition the total lowering for this cell between the three size fractions: do this by relative erodibility
   int
       nFineWeight = (dExistingAvailableFine > 0 ? 1 : 0),
       nSandWeight = (dExistingAvailableSand > 0 ? 1 : 0),
       nCoarseWeight = (dExistingAvailableCoarse > 0 ? 1 : 0);

   double
       dTotErodibility = (nFineWeight * m_dFineErodibilityNormalized) + (nSandWeight * m_dSandErodibilityNormalized) + (nCoarseWeight * m_dCoarseErodibilityNormalized),
       dTotActualErosion = 0,
       dSandEroded = 0,
       dCoarseEroded = 0;

   if (nFineWeight)
   {
      // Erode some fine-sized consolidated sediment
      double dFineLowering = (m_dFineErodibilityNormalized * dThisPotentialErosion) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the cell
      double dFineEroded = tMin(dExistingAvailableFine, dFineLowering);
      double dRemaining = dExistingAvailableFine - dFineEroded;

      dTotActualErosion += dFineEroded;

      // Set the value for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->SetFineDepth(dRemaining);

      // And set the changed-this-timestep switch
      m_bConsChangedThisIter[nThisLayer] = true;

      // And increment the per-timestep total, also add to the suspended sediment load
      m_dThisIterActualPlatformErosionFineCons += dFineEroded;
      m_dThisIterFineSedimentToSuspension += dFineEroded;
   }

   if (nSandWeight)
   {
      // Erode some sand-sized consolidated sediment
      double dSandLowering = (m_dSandErodibilityNormalized * dThisPotentialErosion) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dSandEroded = tMin(dExistingAvailableSand, dSandLowering);
      double dRemaining = dExistingAvailableSand - dSandEroded;

      dTotActualErosion += dSandEroded;

      // Set the new value of sand consolidated sediment depth for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->SetSandDepth(dRemaining);
      
      // And add this to the depth of sand unconsolidated sediment for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->AddSandDepth(dSandEroded);      

      // Set the changed-this-timestep switch
      m_bConsChangedThisIter[nThisLayer] = true;

      // And increment the per-timestep total
      m_dThisIterActualPlatformErosionSandCons += dSandEroded;
   }

   if (nCoarseWeight)
   {
      // Erode some coarse-sized consolidated sediment
      double dCoarseLowering = (m_dCoarseErodibilityNormalized * dThisPotentialErosion) / dTotErodibility;

      // Make sure we don't get -ve amounts left on the source cell
      dCoarseEroded = tMin(dExistingAvailableCoarse, dCoarseLowering);
      double dRemaining = dExistingAvailableCoarse - dCoarseEroded;

      dTotActualErosion += dCoarseEroded;

      // Set the new value of coarse consolidated sediment depth for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->SetCoarseDepth(dRemaining);
      
      // And add this to the depth of coarse unconsolidated sediment for this layer
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->AddCoarseDepth(dCoarseEroded);

      // Set the changed-this-timestep switch
      m_bConsChangedThisIter[nThisLayer] = true;

      // And increment the per-timestep total
      m_dThisIterActualPlatformErosionCoarseCons += dCoarseEroded;
   }

   // Did we erode anything?
   if (dTotActualErosion > 0)
   {
      // We did, so set the actual erosion value for this cell
      m_pRasterGrid->m_Cell[nX][nY].SetActualPlatformErosion(dTotActualErosion);

      // Recalculate the elevation of every layer
      m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

      // And update the cell's sea depth
      m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

      // Update per-timestep totals
      m_ulThisIterNumActualPlatformErosionCells++;

      // Add eroded sand/coarse sediment for this cell to the polygon that contains the cell, ready for redistribution during beach erosion/deposition (fine sediment has already been dealt with)
      int nPolyID = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
      if (nPolyID == INT_NODATA)
      {
         // Can get occasional problems with polygon rasterization near the coastline, so also search the eight adjacent cells
         int nDirection[] = {NORTH, NORTH_EAST, EAST, SOUTH_EAST, SOUTH, SOUTH_WEST, WEST, NORTH_WEST};
         Rand1Shuffle(nDirection, 8);

         for (int n = 0; n < 8; n++)
         {
            if (nDirection[n] == NORTH)
            {
               int
                   nXAdj = nX,
                   nYAdj = nY - 1;
               if (nYAdj >= 0)
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();
                  if (nPolyID != INT_NODATA)
                     break;
               }
            }
            else if (nDirection[n] == NORTH_EAST)
            {
               int
                   nXAdj = nX + 1,
                   nYAdj = nY - 1;
               if ((nXAdj < m_nXGridMax) && (nYAdj >= 0))
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();
                  if (nPolyID != INT_NODATA)
                     break;
               }
            }
            else if (nDirection[n] == EAST)
            {
               int
                   nXAdj = nX + 1,
                   nYAdj = nY;
               if (nXAdj < m_nXGridMax)
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();
                  if (nPolyID != INT_NODATA)
                     break;
               }
            }
            else if (nDirection[n] == SOUTH_EAST)
            {
               int
                   nXAdj = nX + 1,
                   nYAdj = nY + 1;
               if ((nXAdj < m_nXGridMax) && (nYAdj < m_nYGridMax))
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();
                  if (nPolyID != INT_NODATA)
                     break;
               }
            }
            else if (nDirection[n] == SOUTH)
            {
               int
                   nXAdj = nX,
                   nYAdj = nY + 1;
               if (nYAdj < m_nYGridMax)
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();
                  if (nPolyID != INT_NODATA)
                     break;
               }
            }
            else if (nDirection[n] == SOUTH_WEST)
            {
               int
                   nXAdj = nX - 1,
                   nYAdj = nY + 1;
               if ((nXAdj >= 0) && (nXAdj < m_nXGridMax))
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();
                  if (nPolyID != INT_NODATA)
                     break;
               }
            }
            else if (nDirection[n] == WEST)
            {
               int
                   nXAdj = nX - 1,
                   nYAdj = nY;
               if (nXAdj >= 0)
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();
                  if (nPolyID != INT_NODATA)
                     break;
               }
            }
            else if (nDirection[n] == NORTH_WEST)
            {
               int
                   nXAdj = nX - 1,
                   nYAdj = nY - 1;
               if ((nXAdj >= 0) && (nYAdj >= 0))
               {
                  nPolyID = m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetPolygonID();
                  if (nPolyID != INT_NODATA)
                     break;
               }
            }
         }
      }

      if (nPolyID != INT_NODATA)
      {
         // Add this to the polygon's total of unconsolidated sand/coarse sediment, to be deposited or moved later. These values are +ve (deposition)
         m_pVCoastPolygon[nPolyID]->AddPlatformErosionUnconsSand(dSandEroded);
         m_pVCoastPolygon[nPolyID]->AddPlatformErosionUnconsCoarse(dCoarseEroded);
      }
      else
      {
         // Uh-oh, we have a problem
         if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
            LogStream << m_ulIter << ": " << WARN << "platform erosion on cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} but this is not in a polygon" << endl;

         // m_dDepositionSandDiff and m_dDepositionCoarseDiff are both +ve
         m_dDepositionSandDiff += dSandEroded;
         m_dDepositionCoarseDiff += dCoarseEroded;
      }
   }
}

//===============================================================================================================================
//! Creates a look-up table for erosion potential, given depth over DB
//===============================================================================================================================
bool CSimulation::bCreateErosionPotentialLookUp(vector<double> *VdDepthOverDBIn, vector<double> *VdErosionPotentialIn, vector<double> *VdErosionPotentialFirstDerivIn)
{
   // Set up a temporary vector to hold the incremental DepthOverDB values
   vector<double> VdDepthOverDB;
   double dTempDOverDB = 0;

   while (dTempDOverDB <= 1.1) // Arbitrary max value, we will adjust this later
   {
      VdDepthOverDB.push_back(dTempDOverDB); // These are the incremental sample values of DepthOverDB
      dTempDOverDB += DEPTH_OVER_DB_INCREMENT;

      m_VdErosionPotential.push_back(0); // This will hold the corresponding value of erosion potential for each sample point
   }

   int nSize = static_cast<int>(VdDepthOverDB.size());
   vector<double>
       VdDeriv(nSize, 0),   // First derivative at the sample points: calculated by the spline function but not subsequently used
       VdDeriv2(nSize, 0.), // Second derivative at the sample points, ditto
       VdDeriv3(nSize, 0.); // Third derivative at the sample points, ditto

   // Calculate the value of erosion potential (is a -ve value) for each of the sample values of DepthOverDB, and store it for use in the look-up function
   hermite_cubic_spline_value(static_cast<int>(VdDepthOverDBIn->size()), &(VdDepthOverDBIn->at(0)), &(VdErosionPotentialIn->at(0)), &(VdErosionPotentialFirstDerivIn->at(0)), nSize, &(VdDepthOverDB[0]), &(m_VdErosionPotential[0]), &(VdDeriv[0]), &(VdDeriv2[0]), &(VdDeriv3[0]));

   // Tidy the erosion potential look-up data: cut off values (after the first) for which erosion potential is no longer -ve
   int nLastVal = -1;

   for (int n = 1; n < nSize - 1; n++)
      if (m_VdErosionPotential[n] > 0)
      {
         nLastVal = n;
         break;
      }

   if (nLastVal > 0)
   {
      // Erosion potential is no longer -ve at this value of DepthOverDB, so set the maximum value of DepthOverDB that will be used in the simulation (any DepthOverDB value greater than this produces zero erosion potential)
      m_dDepthOverDBMax = VdDepthOverDB[nLastVal];
      m_VdErosionPotential.erase(m_VdErosionPotential.begin() + nLastVal + 1, m_VdErosionPotential.end());
      m_VdErosionPotential.back() = 0;
   }

   else
      // Erosion potential is unbounded, i.e. it is still -ve when we have reached the end of the look-up vector
      return false;

   // All OK
   return true;
}

//===============================================================================================================================
//! The erosion potential lookup: it returns a value for erosion potential given a value of Depth Over DB
//===============================================================================================================================
double CSimulation::dLookUpErosionPotential(double const dDepthOverDB) const
{
   // If dDepthOverDB exceeds the maximum value which we calculated earlier, erosion potential is zero
   if (dDepthOverDB > m_dDepthOverDBMax)
      return 0;

   // OK, dDepthOverDB is less than the maximum so look up a corresponding value for erosion potential. The look-up index is dDepthOverDB divided by (the Depth Over DB increment used when creating the look-up vector). But since this look-up index may not be an integer, split the look-up index into integer and fractional parts and deal with each separately
   double dErosionPotential = dGetInterpolatedValue(&m_VdDepthOverDB, &m_VdErosionPotential, dDepthOverDB, false);

   return dErosionPotential;
}

//===============================================================================================================================
//! Calculates the (inverse) beach protection factor as in SCAPE: 0 is fully protected, 1 = no protection
//===============================================================================================================================
double CSimulation::dCalcBeachProtectionFactor(int const nX, int const nY, double const dBreakingWaveHeight)
{
   // Safety check
   if (bFPIsEqual(dBreakingWaveHeight, DBL_NODATA, TOLERANCE))
      return 0;

   // We are considering the unconsolidated sediment (beach) of the topmost layer that has non-zero thickness
   int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

   // Safety check
   if (nThisLayer == INT_NODATA)
   {
      cerr << ERR << "no sediment layer in dCalcBeachProtectionFactor()" << endl;
      return 0;
   }

   if (nThisLayer == NO_NONZERO_THICKNESS_LAYERS)
      // There are no layers with non-zero thickness left (i.e. we are down to basement) so no beach protection
      return 0;

   // In SCAPE, 0.23 * the significant breaking wave height is assumed to be the maximum depth of beach that waves can penetrate to erode a platform. For depths less than this, the beach protective ability is assumed to vary linearly
   double
       dBeachDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->dGetUnconsolidatedThickness(),
       dMaxPenetrationDepth = BEACH_PROTECTION_HB_RATIO * dBreakingWaveHeight,
       dFactor = 0;

   if (dMaxPenetrationDepth > 0)
      dFactor = tMax(1 - (dBeachDepth / dMaxPenetrationDepth), 0.0);

   //    LogStream << m_ulIter << ": cell[" << nX << "][" << nY << "] has beach protection factor = " << dFactor << endl;

   return dFactor;
}

//===============================================================================================================================
//! Fills in 'holes' in the beach protection i.e. orphan cells which get omitted because of rounding problems
//===============================================================================================================================
void CSimulation::FillInBeachProtectionHoles(void)
{
   for (int nX = 0; nX < m_nXGridMax; nX++)
   {
      for (int nY = 0; nY < m_nYGridMax; nY++)
      {
         if ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
         {
            // This is a sea cell, and it has an initialized beach protection value. So look at its N-S and W-E neighbours
            int
                nXTmp,
                nYTmp,
                nAdjacent = 0;
            double
                dBeachProtection = 0;

            // North
            nXTmp = nX;
            nYTmp = nY - 1;
            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // East
            nXTmp = nX + 1;
            nYTmp = nY;
            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // South
            nXTmp = nX;
            nYTmp = nY + 1;
            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // West
            nXTmp = nX - 1;
            nYTmp = nY;
            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor(), DBL_NODATA, TOLERANCE)))
            {
               nAdjacent++;
               dBeachProtection += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetBeachProtectionFactor();
            }

            // If this sea cell has four neighbours with initialized beach protection values, then assume that it should not have an uninitialized beach protection value. Set it to the average of its neighbours
            if (nAdjacent == 4)
            {
               m_pRasterGrid->m_Cell[nX][nY].SetBeachProtectionFactor(dBeachProtection / 4);
            }
         }
      }
   }
}

//===============================================================================================================================
//! Fills in 'holes' in the potential platform erosion i.e. orphan cells which get omitted because of rounding problems
//===============================================================================================================================
void CSimulation::FillPotentialPlatformErosionHoles(void)
{
   for (int nX = 0; nX < m_nXGridMax; nX++)
   {
      for (int nY = 0; nY < m_nYGridMax; nY++)
      {
         if ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
         {
            // This is a sea cell, it has a zero potential platform erosion value. So look at its N-S and W-E neighbours
            int
                nXTmp,
                nYTmp,
                nAdjacent = 0;
            double
                dPotentialPlatformErosion = 0;

            // North
            nXTmp = nX;
            nYTmp = nY - 1;
            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
            {
               nAdjacent++;
               dPotentialPlatformErosion += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion();
            }

            // East
            nXTmp = nX + 1;
            nYTmp = nY;
            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
            {
               nAdjacent++;
               dPotentialPlatformErosion += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion();
            }

            // South
            nXTmp = nX;
            nYTmp = nY + 1;
            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
            {
               nAdjacent++;
               dPotentialPlatformErosion += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion();
            }

            // West
            nXTmp = nX - 1;
            nYTmp = nY;
            if ((bIsWithinValidGrid(nXTmp, nYTmp)) && (! bFPIsEqual(m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion(), 0.0, TOLERANCE)))
            {
               nAdjacent++;
               dPotentialPlatformErosion += m_pRasterGrid->m_Cell[nXTmp][nYTmp].dGetPotentialPlatformErosion();
            }

            // If this sea cell has four neighbours with non-zero potential platform erosion values, then assume that it should not have a zero potential platform erosion value. Set it to the average of its neighbours
            if (nAdjacent == 4)
            {
               double dThisPotentialPlatformErosion = dPotentialPlatformErosion / 4;

               m_pRasterGrid->m_Cell[nX][nY].SetPotentialPlatformErosion(dThisPotentialPlatformErosion);

               // Update this-timestep totals
               m_ulThisIterNumPotentialPlatformErosionCells++;
               m_dThisIterPotentialPlatformErosion += dThisPotentialPlatformErosion;
               //                assert(isfinite(m_dThisIterPotentialPlatformErosion));

               // Increment the check values
               m_ulTotPotentialPlatformErosionBetweenProfiles++;
               m_dTotPotentialPlatformErosionBetweenProfiles += dThisPotentialPlatformErosion;
            }
         }
      }
   }
}

//===============================================================================================================================
//! Constructs a parallel coastline-normal profile
//===============================================================================================================================
void CSimulation::ConstructParallelProfile(int const nProfileStartX, int const nProfileStartY, int const nParCoastX, int const nParCoastY, int const nProfSize, vector<CGeom2DIPoint>* const pPtViGridProfile, vector<CGeom2DIPoint>* pPtiVGridParProfile, vector<CGeom2DPoint>* pPtVExtCRSParProfile)
{
   // OK, we have the coastline start point for the parallel profile. Now construct a temporary profile, parallel to the coastline-normal profile, starting from this point
   int
       nXOffset = nParCoastX - nProfileStartX,
       nYOffset = nParCoastY - nProfileStartY;

   // Append co-ord values for the temporary profile
   for (int nProfileStartPoint = 0; nProfileStartPoint < nProfSize; nProfileStartPoint++)
   {
      // Get the grid co-ordinates of the cell which is this distance seaward, from the coastline-normal profile
      int
          nXProf = pPtViGridProfile->at(nProfileStartPoint).nGetX(),
          nYProf = pPtViGridProfile->at(nProfileStartPoint).nGetY();

      // Now calculate the grid co-ordinates of this cell, which is potentially in the parallel profile
      int
          nXPar = nXProf + nXOffset,
          nYPar = nYProf + nYOffset;

      // Is this cell within the grid? If not, cut short the profile
      if (! bIsWithinValidGrid(nXPar, nYPar))
      {
         //         LogStream << "NOT WITHIN GRID [" << nXPar << "][" << nYPar << "]" << endl;
         return;
      }

      // Have we hit an adjacent coastline-normal profile? If so, cut short
      if (m_pRasterGrid->m_Cell[nXPar][nYPar].bIsProfile())
      {
         //         LogStream << "HIT PROFILE " << m_pRasterGrid->m_Cell[nXPar][nYPar].nGetProfile() << " at [" << nXPar << "][" << nYPar << "] = {" << dGridCentroidXToExtCRSX(nXPar) << ", " <<  dGridCentroidYToExtCRSY(nYPar) << "}" << endl;
         return;
      }

      // OK, append the cell details
      pPtiVGridParProfile->push_back(CGeom2DIPoint(nXPar, nYPar));
      pPtVExtCRSParProfile->push_back(CGeom2DPoint(dGridCentroidXToExtCRSX(nXPar), dGridCentroidYToExtCRSY(nYPar)));
   }
}
