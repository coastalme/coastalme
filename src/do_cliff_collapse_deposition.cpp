/*!
   \file do_cliff_collapse_deposition.cpp
   \brief Following cliff collapse, distributes both consolidated and unconsolidated sediment from the collapse onto the shore polygons as unconsolidated talus
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License
*/

/* ==============================================================================================================================
   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
==============================================================================================================================*/
#include <assert.h>

#include <cmath>

#include <iostream>
using std::cerr;
using std::endl;
using std::ios;

#include "cme.h"
#include "simulation.h"
#include "cliff.h"
#include "coast_landform.h"
#include "2d_point.h"

//===============================================================================================================================
//! Deposit the unconsolidated sediment from cliff collapse as talus on the cell on which collapse occurred.
//===============================================================================================================================
int CSimulation::nDoCliffCollapseTalusDeposition(int const nCoast, CRWCliff const* pCliff, double const dSandFromCollapse, double const dCoarseFromCollapse, double const dPreCollapseCellElev, double const dPostCollapseCellElevNoTalus)
{
   // Check: is there some sand- or coarse-sized sediment to deposit?
   if ((dSandFromCollapse + dCoarseFromCollapse) < SEDIMENT_ELEV_TOLERANCE)
      return RTN_OK;

   LogStream << m_ulIter << "; in nDoCliffCollapseTalusDeposition() dSandFromCollapse = " << dSandFromCollapse << " dCoarseFromCollapse = " << dCoarseFromCollapse << endl;

   // Get the cliff cell's grid coords
   int const nX = pCliff->pPtiGetCellMarkedAsCliff()->nGetX();
   int const nY = pCliff->pPtiGetCellMarkedAsCliff()->nGetY();

   // Get the number of the highest layer with non-zero thickness
   int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();

   // Safety checks
   if (nTopLayer == INT_NODATA)
      return RTN_ERR_NO_TOP_LAYER;

   // Get a pointer to this layer
   CRWCellLayer* pLayer = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer);

   // And get a pointer to the cell layer's talus object
   CRWCellTalus* pTalus = pLayer->pGetTalus();















   // // OK, we have some sand- and/or coarse-sized sediment to deposit
   // int const nStartPoint = pCliff->nGetPointOnCoast();
   // int const nCoastSize = m_VCoast[nCoast].nGetCoastlineSize();
   //
   //
   // // // And get the cliff cell's ext crs coords
   // // double const dXCliff = dGridXToExtCRSX(nXCliff);
   // // double const dYCliff = dGridYToExtCRSY(nYCliff);
   //
   // // Get this cell's polygon
   // int const nPoly = m_pRasterGrid->m_Cell[nXCliff][nYCliff].nGetPolygonID();
   // if (nPoly == INT_NODATA)
   // {
   //    // This cell isn't in a polygon
   //    LogStream << m_ulIter << " : in nDoCliffCollapse(), [" << nXCliff << "][" << nYCliff << "] = {" << dGridCentroidXToExtCRSX(nXCliff) << ", " << dGridCentroidYToExtCRSY(nYCliff) << "} is not in a polygon" << endl;
   //
   //    return RTN_ERR_CLIFF_NOT_IN_POLYGON;
   // }
   //
   // CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
   //
   // // OK, now set up the planview sequence talus deposition. First we deposit a volume of talus which fits "under" a Dean profile starting from the coast (i.e. the cliff collapse cell). Then we deposit along the two profiles which start from the coast cells on either side of the cliff collapse cell, and then on two more profiles starting on either side of the last two coast cells, and so on. This holds the valid coast start points for each cliff collapse Dean profile
   // vector<int> VnTalusProfileCoastStartPoint(nCoastSize);
   //
   // // This is the coast start point for the first Dean profile
   // VnTalusProfileCoastStartPoint[0] = nStartPoint;
   // int nn = 1;
   // int nSigned = 1;
   // int nCount = 1;
   //
   // do
   // {
   //    int nTmpPoint;
   //
   //    if ((nCount % 2) != 0)
   //       nTmpPoint = nStartPoint + nSigned;
   //    else
   //    {
   //       nTmpPoint = nStartPoint - nSigned;
   //       nSigned++;
   //    }
   //
   //    // Is this coast start point valid?
   //    if ((nTmpPoint < 0) || (nTmpPoint > (nCoastSize - 1)))
   //    {
   //       // No, it is outside the grid, so find another
   //       nCount++;
   //       continue;
   //    }
   //    else
   //    {
   //       // It is valid
   //       VnTalusProfileCoastStartPoint[nn] = nTmpPoint;
   //       nn++;
   //       nCount++;
   //    }
   // } while (nn < nCoastSize);
   //
   // bool bHitFirstCoastPoint = false;
   // bool bHitLastCoastPoint = false;
   //
   // double dTotSandToDepositAllProfiles = dSandFromCollapse;     // Note that this can increase, if we get erosion because Dean profile is lower than profile at a point
   // double dTotCoarseToDepositAllProfiles = dCoarseFromCollapse; // Note that this can increase, if we get erosion because Dean profile is lower than profile at a point
   // double dTotSandDepositedAllProfiles = 0;
   // double dTotCoarseDepositedAllProfiles = 0;
   //
   // // Process each cliff collapse deposition profile
   // for (int nDepProfile = 0; nDepProfile < nCoastSize; nDepProfile++)
   // {
   //    bool bDoSandDepositionOnThisProfile = false;
   //    bool bSandDepositionCompletedOnThisProfile = false;
   //    bool bDoCoarseDepositionOnThisProfile = false;
   //    bool bCoarseDepositionCompletedOnThisProfile = false;
   //
   //    int const nRemainingProfiles = tMax(1, m_nDefaultTalusWidthInCells - nDepProfile);
   //
   //    // This is the minimum planview length (in cells) of the Dean profile. The initial length will be increased if we can't deposit sufficient talus
   //    int nTalusProfileLenInCells = m_nTalusProfileMinLenInCells;
   //
   //    // Calculate the target amount to be deposited on each talus profile, assuming the "preferred" talus width
   //    double dTargetSandToDepositOnThisProfile = (dTotSandToDepositAllProfiles - dTotSandDepositedAllProfiles) / nRemainingProfiles;
   //    double dTargetCoarseToDepositOnThisProfile = (dTotCoarseToDepositAllProfiles - dTotCoarseDepositedAllProfiles) / nRemainingProfiles;
   //
   //    // Use this to make sure that, if we have both sand and coarse to deposit, we don't drop all the sand on a cell and then be unable to deposit any coarse
   //    double dSandProp = 0.5;
   //    double dCoarseProp = 1 - dSandProp;
   //
   //    if (dTargetSandToDepositOnThisProfile + dTargetCoarseToDepositOnThisProfile > 0)
   //    {
   //       dSandProp = dTargetSandToDepositOnThisProfile / (dTargetSandToDepositOnThisProfile + dTargetCoarseToDepositOnThisProfile);
   //       dCoarseProp = 1 - dSandProp;
   //    }
   //
   //    double dSandDepositedOnThisProfile = 0;
   //    double dCoarseDepositedOnThisProfile = 0;
   //
   //    if (dTotSandToDepositAllProfiles > 0)
   //    {
   //       bDoSandDepositionOnThisProfile = true;
   //
   //       if (bFPIsEqual(dTotSandToDepositAllProfiles - dTotSandDepositedAllProfiles, 0.0, MASS_BALANCE_TOLERANCE))
   //          bSandDepositionCompletedOnThisProfile = true;
   //    }
   //    else
   //       bSandDepositionCompletedOnThisProfile = true;
   //
   //    if (dTotCoarseToDepositAllProfiles > 0)
   //    {
   //       bDoCoarseDepositionOnThisProfile = true;
   //
   //       if (bFPIsEqual(dTotCoarseToDepositAllProfiles - dTotCoarseDepositedAllProfiles, 0.0, MASS_BALANCE_TOLERANCE))
   //          bCoarseDepositionCompletedOnThisProfile = true;
   //    }
   //    else
   //       bCoarseDepositionCompletedOnThisProfile = true;
   //
   //    if (bSandDepositionCompletedOnThisProfile && bCoarseDepositionCompletedOnThisProfile)
   //    {
   //       LogStream << m_ulIter << ": break 2 for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nDepProfile = " << nDepProfile << endl << "\tdTargetSandToDepositOnThisProfile = " << dTargetSandToDepositOnThisProfile << " dSandDepositedOnThisProfile = " << dSandDepositedOnThisProfile << " dTotSandToDepositAllProfiles = " << dTotSandToDepositAllProfiles << " dTotSandDepositedAllProfiles = " << dTotSandDepositedAllProfiles << endl << "\tdTargetCoarseToDepositOnThisProfile = " << dTargetCoarseToDepositOnThisProfile << " dCoarseDepositedOnThisProfile = " << dCoarseDepositedOnThisProfile << " dTotCoarseToDepositAllProfiles = " <<  dTotCoarseToDepositAllProfiles << " dTotCoarseDepositedAllProfiles = " << dTotCoarseDepositedAllProfiles << endl;
   //
   //       break;
   //    }
   //
   //    // Get the start point of this cliff collapse deposition profile
   //    int const nThisPoint = VnTalusProfileCoastStartPoint[nDepProfile];
   //
   //    // Are we at the start or end of the coast?
   //    if (nThisPoint == 0)
   //       bHitFirstCoastPoint = true;
   //
   //    if (nThisPoint == nCoastSize - 1)
   //       bHitLastCoastPoint = true;
   //
   //    CGeom2DPoint PtStart;
   //    CGeom2DPoint PtEnd;
   //
   //    // Make the start of the deposition profile the cliff cell that is marked as coast (not the cell under the smoothed vector coast, they may well be different)
   //    PtStart.SetX(dGridCentroidXToExtCRSX(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPoint)->nGetX()));
   //    PtStart.SetY(dGridCentroidYToExtCRSY(m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nThisPoint)->nGetY()));
   //
   //    // Set the initial fraction of cliff height, this will be increased (max is one) if we can't deposit sufficient talus
   //    double dCliffHeightFraction = m_dMinCliffTalusHeightFrac;
   //
   //    bool bJustDepositWhatWeCan = false;
   //
   //    // The initial seaward offset, in cells. This will be increased if we can't deposit sufficient talus
   //    int nSeawardOffset = 1;
   //
   //    // Process this profile
   //    do
   //    {
   //       if (bJustDepositWhatWeCan || (bSandDepositionCompletedOnThisProfile && bCoarseDepositionCompletedOnThisProfile))
   //       {
   //          LogStream << m_ulIter << ": break 3 for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tnSeawardOffset = " << nSeawardOffset << " dCliffHeightFraction = " << dCliffHeightFraction << " nDepProfile = " << nDepProfile << endl << "\tbJustDepositWhatWeCan = " << bJustDepositWhatWeCan << "\tdTargetSandToDepositOnThisProfile = " << dTargetSandToDepositOnThisProfile << " dSandDepositedOnThisProfile = " << dSandDepositedOnThisProfile << " dTotSandToDepositAllProfiles = " << dTotSandToDepositAllProfiles << " dTotSandDepositedAllProfiles = " << dTotSandDepositedAllProfiles << endl << "\tdTargetCoarseToDepositOnThisProfile = " << dTargetCoarseToDepositOnThisProfile << " dCoarseDepositedOnThisProfile = " << dCoarseDepositedOnThisProfile << " dTotCoarseToDepositAllProfiles = " <<  dTotCoarseToDepositAllProfiles << " dTotCoarseDepositedAllProfiles = " << dTotCoarseDepositedAllProfiles << endl;
   //
   //          break;
   //       }
   //
   //       // Has the seaward offset reached the arbitrary limit?
   //       if (nSeawardOffset >= MAX_SEAWARD_OFFSET_FOR_CLIFF_TALUS)
   //       {
   //          // It has, so constrain it as a safety check
   //          nSeawardOffset = MAX_SEAWARD_OFFSET_FOR_CLIFF_TALUS;
   //
   //          // Has the cliff height fraction reached the max value?
   //          if (dCliffHeightFraction >= 0.99)
   //          {
   //             // It has, so constrain it as a safety check
   //             dCliffHeightFraction = 1;
   //
   //             // Is the talus length also at its arbitrary limit?
   //             if (nTalusProfileLenInCells >= MAX_CLIFF_TALUS_LENGTH)
   //             {
   //                // The talus length is also at this limit. So there is nothing more we can increase on this profile. Just deposit what we can and move on to the next profile
   //                bJustDepositWhatWeCan = true;
   //
   //                // LogStream << m_ulIter << ": \tjust deposit what we can: talus length = " << nTalusProfileLenInCells << " cells, seaward offset = " << nSeawardOffset << ", cliff height fraction = " <<  dCliffHeightFraction << endl;
   //             }
   //             else
   //             {
   //                // The talus length is not at its arbitrary limit, so extend it
   //                nTalusProfileLenInCells += CLIFF_COLLAPSE_LENGTH_INCREMENT;
   //
   //                // LogStream << m_ulIter << ": \ttalus length increased: talus length = " << nTalusProfileLenInCells << " cells, seaward offset = " << nSeawardOffset << ", cliff height fraction = " <<  dCliffHeightFraction << endl;
   //
   //                continue;
   //             }
   //          }
   //          else
   //          {
   //             // The cliff height fraction is not at its limit, so increment it
   //             dCliffHeightFraction += CLIFF_COLLAPSE_HEIGHT_INCREMENT;
   //
   //             // LogStream << "Cliff height fraction increased: talus length = " << nTalusProfileLenInCells << " cells, seaward offset = " << nSeawardOffset << ", cliff height fraction = " <<  dCliffHeightFraction << endl;
   //
   //             continue;
   //          }
   //       }
   //       else
   //       {
   //          // The seaward offset is not at its maximum, so extend it
   //          nSeawardOffset++;
   //
   //          // LogStream << m_ulIter << ": \tseaward offset increased: talus length = " << nTalusProfileLenInCells << " cells, seaward offset = " << nSeawardOffset << ", cliff height fraction = " <<  dCliffHeightFraction << endl;
   //
   //          continue;
   //       }
   //
   //       if (bHitFirstCoastPoint && bHitLastCoastPoint)
   //       {
   //          // Uh-oh, we've reached both ends of the coast (!) and we can't increase anything any more
   //          LogStream << m_ulIter << ": unable to deposit enough unconsolidated sediment (talus) from cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tnSeawardOffset = " << nSeawardOffset << " dCliffHeightFraction = " << dCliffHeightFraction << " nDepProfile = " << nDepProfile << endl << "\tdTargetSandToDepositOnThisProfile = " << dTargetSandToDepositOnThisProfile << " dSandDepositedOnThisProfile = " << dSandDepositedOnThisProfile << " dTotSandToDepositAllProfiles = " << dTotSandToDepositAllProfiles << " dTotSandDepositedAllProfiles = " << dTotSandDepositedAllProfiles << endl << "\tdTargetCoarseToDepositOnThisProfile = " << dTargetCoarseToDepositOnThisProfile << " dCoarseDepositedOnThisProfile = " << dCoarseDepositedOnThisProfile << " dTotCoarseToDepositAllProfiles = " << dTotCoarseToDepositAllProfiles << " dTotCoarseDepositedAllProfiles = " << dTotCoarseDepositedAllProfiles << endl;
   //
   //          return RTN_ERR_CLIFF_CANNOT_DEPOSIT_ALL;
   //       }
   //
   //       // Now construct a deposition collapse profile from the start point, it is one cell longer than the specified length because it includes the cliff point in the profile. Calculate its length in external CRS units, the way it is done here is approximate but probably OK
   //       double const dThisProfileLength = (nTalusProfileLenInCells + nSeawardOffset + 1) * m_dCellSide;
   //
   //       // Get the end point of this coastline-normal line
   //       CGeom2DIPoint PtiEnd; // In grid CRS
   //       int const nRtn = nGetCoastNormalEndPoint(nCoast, nThisPoint, nCoastSize, &PtStart, dThisProfileLength, &PtEnd, &PtiEnd, false);
   //
   //       // Safety check
   //       if (nRtn == RTN_ERR_NO_SOLUTION_FOR_ENDPOINT)
   //       {
   //          LogStream << m_ulIter << ": could not find a solution for the end point of the Dean profile at {" << PtStart.dGetX() << ", " << PtStart.dGetY() << "}" << endl;
   //
   //          return nRtn;
   //       }
   //
   //       // Safety check
   //       if (PtStart == PtEnd)
   //       {
   //          // This would give a zero-length profile, and a zero-divide error during rasterization. So just move on to the next profile
   //          LogStream << "Zero length" << endl;
   //          break;
   //       }
   //
   //       // OK, both the start and end points of this deposition profile are within the grid
   //       // LogStream << m_ulIter << ": nWidthDistSigned = " << nWidthDistSigned << " cliff collapse profile from " << PtStart.dGetX() << ", " << PtStart.dGetY() << " to " << PtEnd.dGetX() << ", " << PtEnd.dGetY() << " with length (inc. cliff point) = " << dThisProfileLength << endl;
   //
   //       vector<CGeom2DPoint> VTmpProfile;
   //       VTmpProfile.push_back(PtStart);
   //       VTmpProfile.push_back(PtEnd);
   //       vector<CGeom2DIPoint> VCellsUnderProfile;
   //
   //       // Now get the raster cells under this profile
   //       RasterizeCliffCollapseProfile(&VTmpProfile, &VCellsUnderProfile);
   //
   //       int const nRasterProfileLength = static_cast<int>(VCellsUnderProfile.size());
   //
   //       // Check now, for the case where the profile is very short
   //       if (nRasterProfileLength - nSeawardOffset < 3)
   //       {
   //          // Can't do anything with this very short profile, since (nRasterProfileLength - nSeawardOffset - 2) later will give zero or -ve dInc. So just move on to the next profile
   //          break;
   //       }
   //       else if (nRasterProfileLength - nSeawardOffset == 3)
   //       {
   //          // Can't increase offset any more, or get zero divide with (nRasterProfileLength - nSeawardOffset - 2) later. So just deposit what we can and then move on to the next profile
   //          bJustDepositWhatWeCan = true;
   //       }
   //
   //       vector<double> dVProfileNow(nRasterProfileLength, 0);
   //       vector<bool> bVProfileValid(nRasterProfileLength, true);
   //
   //       LogStream << m_ulIter << ": for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nDepProfile = " << nDepProfile << endl << "\tnRasterProfileLength = " << nRasterProfileLength << " nSeawardOffset = " << nSeawardOffset << " nRasterProfileLength - nSeawardOffset - 2 = " << nRasterProfileLength - nSeawardOffset - 2 << endl;
   //
   //       // Calculate the existing elevation for all points along the deposition profile
   //       for (int n = 0; n < nRasterProfileLength; n++)
   //       {
   //          int const nX = VCellsUnderProfile[n].nGetX();
   //          int const nY = VCellsUnderProfile[n].nGetY();
   //
   //          dVProfileNow[n] = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();
   //
   //          // Don't allow cliff collapse talus onto intervention cells TODO 078 Is this realistic? Should it change with different types on intervention?
   //          if (bIsInterventionCell(nX, nY))
   //             bVProfileValid[n] = false;
   //       }
   //
   //       // Now calculate the elevation of the talus top at the shoreline
   //       double const dCliffHeight = dPreCollapseCellElev - dPostCollapseCellElevNoTalus;
   //       double const dTalusTopElev = dPostCollapseCellElevNoTalus + (dCliffHeight * dCliffHeightFraction);
   //
   //       // LogStream << "Elevations: cliff top = " << dPreCollapseCellElev << " cliff base = " << dPostCollapseCellElevNoTalus << " talus top = " << dTalusTopElev << endl;
   //
   //       // if (dPreCollapseCellElev < dPostCollapseCellElevNoTalus)
   //       // LogStream << "*** ERROR, cliff top is lower than cliff base" << endl;
   //
   //       // Next calculate the talus slope length in external CRS units, this is approximate but probably OK
   //       double const dTalusSlopeLength = dThisProfileLength - ((nSeawardOffset - 1) * m_dCellSide);
   //
   //       // If user has not supplied a value for m_dCliffDepositionA, then solve for dA so that the elevations at end of the existing profile, and at the end of the Dean equilibrium profile, are the same
   //       double dA = 0;
   //
   //       if (bFPIsEqual(m_dCliffDepositionA, 0.0, TOLERANCE))
   //          dA = m_dCliffDepositionA;
   //       else
   //          dA = (dTalusTopElev - dVProfileNow[nRasterProfileLength - 1]) / pow(dTalusSlopeLength, DEAN_POWER);
   //
   //       // assert((nRasterProfileLength - nSeawardOffset - 2) > 0);
   //       double const dInc = dTalusSlopeLength / (nRasterProfileLength - nSeawardOffset - 2);
   //       vector<double> dVDeanProfile(nRasterProfileLength);
   //
   //       // Calculate the Dean equilibrium profile of the talus h(y) = A * y^(2/3) where h(y) is the distance below the talus-top elevation (the highest point in the Dean profile) at a distance y from the cliff (the landward start of the profile)
   //       CalcDeanProfile(&dVDeanProfile, dInc, dTalusTopElev, dA, true, nSeawardOffset, dTalusTopElev);
   //
   //       // Get the total difference in elevation between the two profiles (Dean profile - present profile). Since we want the Dean profile to be higher than the present profile, a good result is a +ve number
   //       double const dTotElevDiff = dSubtractProfiles(&dVDeanProfile, &dVProfileNow, &bVProfileValid);
   //
   //       //          // DEBUG CODE -----------------------------------------------------
   //       // LogStream << endl;
   //       // LogStream << "dTalusSlopeLength = " << dTalusSlopeLength << " dA = " << dA << endl;
   //       // LogStream << "dDistFromTalusStart - dInc = " << dDistFromTalusStart - dInc << " dThisProfileLength - nSeawardOffset - 2 = " << dThisProfileLength - nSeawardOffset - 2 << endl;
   //       // LogStream << "Profile now (inc. cliff cell) = ";
   //       // for (int n = 0; n < nRasterProfileLength; n++)
   //       // {
   //       // int
   //       // nX = VCellsUnderProfile[n].nGetX(),
   //       // nY = VCellsUnderProfile[n].nGetY();
   //       // dVProfileNow[n] = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev();
   //       // LogStream << dVProfileNow[n] << " ";
   //       // }
   //       // LogStream << endl;
   //       // LogStream << "Dean equilibrium profile (inc. cliff cell) = ";
   //       // for (int n = 0; n < nRasterProfileLength; n++)
   //       // {
   //       // LogStream << dVDeanProfile[n] << " ";
   //       // }
   //       // LogStream << endl;
   //       // LogStream << "Difference (inc. cliff cell) = ";
   //       // for (int n = 0; n < nRasterProfileLength; n++)
   //       // {
   //       // LogStream << dVDeanProfile[n] - dVProfileNow[n] << " ";
   //       // }
   //       // LogStream << endl;
   //       //          // DEBUG CODE -----------------------------------------------------
   //
   //       // If we are not in a "just deposit what we can" situation, then for this planview profile, does the Dean equilibrium profile allow us to deposit all the talus sediment which we need to get rid of?
   //       if (! bJustDepositWhatWeCan && (dTotElevDiff < (dTargetSandToDepositOnThisProfile + dTargetCoarseToDepositOnThisProfile)))
   //       {
   //          // No it doesn't, so try again with a larger seaward offset and/or a longer Dean profile length
   //          LogStream << m_ulIter << ": bJustDepositWhatWeCan = " << bJustDepositWhatWeCan << " nSeawardOffset = " << nSeawardOffset << " dTotElevDiff = " << dTotElevDiff << endl;
   //
   //          continue;
   //       }
   //
   //       // OK, now process all cells in this profile, including the first one (which is where the cliff collapse occurred)
   //       for (int n = 0; n < nRasterProfileLength; n++)
   //       {
   //          // Are we depositing sand talus sediment on this profile?
   //          if (bDoSandDepositionOnThisProfile)
   //          {
   //             if (bFPIsEqual(dTargetSandToDepositOnThisProfile - dSandDepositedOnThisProfile, 0.0, MASS_BALANCE_TOLERANCE))
   //                bSandDepositionCompletedOnThisProfile = true;
   //             else
   //                bSandDepositionCompletedOnThisProfile = false;
   //          }
   //
   //          // Are we depositing coarse talus sediment on this profile?
   //          if (bDoCoarseDepositionOnThisProfile)
   //          {
   //             if (bFPIsEqual(dTargetCoarseToDepositOnThisProfile - dCoarseDepositedOnThisProfile, 0.0, MASS_BALANCE_TOLERANCE))
   //                bCoarseDepositionCompletedOnThisProfile = true;
   //             else
   //                bCoarseDepositionCompletedOnThisProfile = false;
   //          }
   //
   //          // If we have deposited enough, then break out of the loop
   //          if (bSandDepositionCompletedOnThisProfile && bCoarseDepositionCompletedOnThisProfile)
   //          {
   //             LogStream << m_ulIter << ": break 1 for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tbJustDepositWhatWeCan = " << bJustDepositWhatWeCan << " nSeawardOffset = " << nSeawardOffset << " dCliffHeightFraction = " << dCliffHeightFraction << " nDepProfile = " << nDepProfile << endl << "\tdTargetSandToDepositOnThisProfile = " << dTargetSandToDepositOnThisProfile << " dSandDepositedOnThisProfile = " << dSandDepositedOnThisProfile << " dTotSandToDepositAllProfiles = " << dTotSandToDepositAllProfiles << " dTotSandDepositedAllProfiles = " << dTotSandDepositedAllProfiles << endl << "\tdTargetCoarseToDepositOnThisProfile = " << dTargetCoarseToDepositOnThisProfile << " dCoarseDepositedOnThisProfile = " << dCoarseDepositedOnThisProfile << " dTotCoarseToDepositAllProfiles = " <<  dTotCoarseToDepositAllProfiles << " dTotCoarseDepositedAllProfiles = " << dTotCoarseDepositedAllProfiles << endl;
   //
   //             break;
   //          }
   //
   //          // Nope, we still have some talus left to deposit
   //          int const nX = VCellsUnderProfile[n].nGetX();
   //          int const nY = VCellsUnderProfile[n].nGetY();
   //
   //          // Don't allow cliff collapse talus onto intervention cells TODO 078 Is this realistic? Should it change with different types on intervention?
   //          if (bIsInterventionCell(nX, nY))
   //             continue;
   //
   //          int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();
   //
   //          // Safety check
   //          if (nTopLayer == INT_NODATA)
   //             return RTN_ERR_NO_TOP_LAYER;
   //
   //          if (nTopLayer == NO_NONZERO_THICKNESS_LAYERS)
   //          {
   //             // TODO 021 Improve this
   //             cerr << m_ulIter << ": all layers have zero thickness" << endl;
   //             return RTN_ERR_CLIFF_CANNOT_DEPOSIT_ALL;
   //          }
   //
   //          // Only do deposition on this cell if its elevation is below the Dean elevation
   //          if (dVDeanProfile[n] > dVProfileNow[n])
   //          {
   //             // At this point along the profile, the Dean profile is higher than the present profile. So we can deposit some sediment on this cell
   //             double dSandToDeposit = 0;
   //
   //             if (bDoSandDepositionOnThisProfile)
   //             {
   //                dSandToDeposit = (dVDeanProfile[n] - dVProfileNow[n]) * dSandProp;
   //                dSandToDeposit = tMin(dSandToDeposit, (dTargetSandToDepositOnThisProfile - dSandDepositedOnThisProfile), (dTotSandToDepositAllProfiles - dTotSandDepositedAllProfiles));
   //
   //                m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSandDepth(dSandToDeposit);
   //
   //                // Set the changed-this-timestep switch
   //                m_bUnconsChangedThisIter[nTopLayer] = true;
   //
   //                dSandDepositedOnThisProfile += dSandToDeposit;
   //                dTotSandDepositedAllProfiles += dSandToDeposit;
   //             }
   //
   //             double dCoarseToDeposit = 0;
   //
   //             if (bDoCoarseDepositionOnThisProfile)
   //             {
   //                dCoarseToDeposit = (dVDeanProfile[n] - dVProfileNow[n]) * dCoarseProp;
   //                dCoarseToDeposit = tMin(dCoarseToDeposit, (dTargetCoarseToDepositOnThisProfile - dCoarseDepositedOnThisProfile), (dTotCoarseToDepositAllProfiles - dTotCoarseDepositedAllProfiles));
   //
   //                m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarseDepth(dCoarseToDeposit);
   //
   //                // Set the changed-this-timestep switch
   //                m_bUnconsChangedThisIter[nTopLayer] = true;
   //
   //                dCoarseDepositedOnThisProfile += dCoarseToDeposit;
   //                dTotCoarseDepositedAllProfiles += dCoarseToDeposit;
   //             }
   //
   //             // Now update the cell's layer elevations
   //             m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();
   //
   //             // Update the cell's sea depth
   //             m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();
   //
   //             // Update the cell's talus deposition, and total talus deposition, values
   //             m_pRasterGrid->m_Cell[nX][nY].AddSandTalusDeposition(dSandToDeposit);
   //             m_pRasterGrid->m_Cell[nX][nY].AddCoarseTalusDeposition(dCoarseToDeposit);
   //
   //             // And set the landform category
   //             CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetLandform();
   //             int const nCat = pLandform->nGetLFCategory();
   //
   //             if ((nCat != LF_CAT_SEDIMENT_INPUT) && (nCat != LF_CAT_SEDIMENT_INPUT_SUBMERGED) && (nCat != LF_CAT_SEDIMENT_INPUT_NOT_SUBMERGED))
   //                pLandform->SetLFSubCategory(LF_SUBCAT_DRIFT_TALUS);
   //          }
   //          else if (dVDeanProfile[n] < dVProfileNow[n])
   //          {
   //             // Here, the Dean profile is lower than the existing profile, so we must remove some sediment from this cell  TODO 075 What if bedrock sticks above Dean profile?
   //             double const dThisLowering = dVProfileNow[n] - dVDeanProfile[n];
   //
   //             // Find out how much sediment we have available on this cell
   //             double const dExistingAvailableFine = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetFineDepth();
   //             double const dExistingAvailableSand = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
   //             double const dExistingAvailableCoarse = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();
   //
   //             // Now partition the total lowering for this cell between the three size fractions: do this by relative erodibility
   //             int const nFineWeight = (dExistingAvailableFine > 0 ? 1 : 0);
   //             int const nSandWeight = (dExistingAvailableSand > 0 ? 1 : 0);
   //             int const nCoarseWeight = (dExistingAvailableCoarse > 0 ? 1 : 0);
   //
   //             double const dTotErodibility = (nFineWeight * m_dFineErodibilityNormalized) + (nSandWeight * m_dSandErodibilityNormalized) + (nCoarseWeight * m_dCoarseErodibilityNormalized);
   //
   //             if (nFineWeight)
   //             {
   //                // Erode some fine-sized sediment
   //                double const dFineLowering = (m_dFineErodibilityNormalized * dThisLowering) / dTotErodibility;
   //
   //                // Make sure we don't get -ve amounts left on the cell
   //                double const dFine = tMin(dExistingAvailableFine, dFineLowering);
   //                double const dRemaining = dExistingAvailableFine - dFine;
   //
   //                // Set the value for this layer
   //                m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetFineDepth(dRemaining);
   //
   //                // And set the changed-this-timestep switch
   //                m_bUnconsChangedThisIter[nTopLayer] = true;
   //
   //                // And increment the per-timestep total for fine sediment eroded during cliff collapse deposition (note that this gets added in to the suspended load elsewhere, so no need to do it here)
   //                m_dThisIterCliffCollapseFineErodedDuringDeposition += dFine;
   //
   //                // LogStream << m_ulIter << ": FINE erosion during cliff collapse talus deposition = " << dFine * m_dCellArea << endl;
   //
   //                // Also add to the suspended load
   //                m_dThisIterFineSedimentToSuspension += dFine;
   //             }
   //
   //             if (nSandWeight)
   //             {
   //                // Erode some sand-sized sediment
   //                double const dSandLowering = (m_dSandErodibilityNormalized * dThisLowering) / dTotErodibility;
   //
   //                // Make sure we don't get -ve amounts left on the source cell
   //                double const dSandToErode = tMin(dExistingAvailableSand, dSandLowering);
   //                double const dRemaining = dExistingAvailableSand - dSandToErode;
   //
   //                // Set the value for this layer
   //                m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dRemaining);
   //
   //                // Set the changed-this-timestep switch
   //                m_bUnconsChangedThisIter[nTopLayer] = true;
   //
   //                // And increment the per-timestep total for sand sediment eroded during cliff collapse deposition
   //                m_dThisIterCliffCollapseSandErodedDuringDeposition += dSandToErode;
   //
   //                // Increase the all-profiles and this-profile sand deposition targets
   //                dTargetSandToDepositOnThisProfile += dSandToErode;
   //                dTotSandToDepositAllProfiles += dSandToErode;
   //
   //                // LogStream << m_ulIter << ": SAND erosion during cliff collapse talus deposition = " << dSandToErode * m_dCellArea << endl;
   //
   //                // Store the depth of sand sediment eroded during Dean profile deposition of sand cliff collapse talus
   //                pPolygon->AddCliffCollapseSandErodedDeanProfile(dSandToErode);
   //             }
   //
   //             if (nCoarseWeight)
   //             {
   //                // Erode some coarse-sized sediment
   //                double const dCoarseLowering = (m_dCoarseErodibilityNormalized * dThisLowering) / dTotErodibility;
   //
   //                // Make sure we don't get -ve amounts left on the source cell
   //                double const dCoarseToErode = tMin(dExistingAvailableCoarse, dCoarseLowering);
   //                double const dRemaining = dExistingAvailableCoarse - dCoarseToErode;
   //
   //                // Set the value for this layer
   //                m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dRemaining);
   //
   //                // Set the changed-this-timestep switch
   //                m_bUnconsChangedThisIter[nTopLayer] = true;
   //
   //                // And increment the per-timestep total for coarse sediment eroded during cliff collapse deposition
   //                m_dThisIterCliffCollapseCoarseErodedDuringDeposition += dCoarseToErode;
   //
   //                // Increase the all-profiles and this-profile coarse deposition targets
   //                dTargetCoarseToDepositOnThisProfile += dCoarseToErode;
   //                dTotCoarseToDepositAllProfiles += dCoarseToErode;
   //
   //                // LogStream << m_ulIter << ": COARSE erosion during cliff collapse talus deposition = " << dCoarseToErode * m_dCellArea << endl;
   //
   //                // Store the depth of coarse sediment eroded during Dean profile deposition of coarse cliff collapse talus
   //                pPolygon->AddCliffCollapseCoarseErodedDeanProfile(dCoarseToErode);
   //             }
   //
   //             // for this cell, recalculate the elevation of every layer
   //             m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();
   //
   //             // And update the cell's sea depth
   //             m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();
   //          }
   //       } // All cells in this profile
   //
   //       // OK we have either processed all cells in this profile, or we have deposited enough talus sediment on this profile
   //       break;
   //    } while (true); // The seaward offset etc. loop
   //
   //    // LogStream << m_ulIter << ": \tleft seaward offset loop for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tnSeawardOffset = " << nSeawardOffset << " dCliffHeightFraction = " << dCliffHeightFraction << " nDepProfile = " << nDepProfile << " bJustDepositWhatWeCan = " << bJustDepositWhatWeCan << endl << "\tdTargetSandToDepositOnThisProfile = " << dTargetSandToDepositOnThisProfile << " dSandDepositedOnThisProfile = " << dSandDepositedOnThisProfile << " dTotSandToDepositAllProfiles = " << dTotSandToDepositAllProfiles << " dTotSandDepositedAllProfiles = " << dTotSandDepositedAllProfiles << endl << "\tdTargetCoarseToDepositOnThisProfile = " << dTargetCoarseToDepositOnThisProfile << " dCoarseDepositedOnThisProfile = " << dCoarseDepositedOnThisProfile << " dTotCoarseToDepositAllProfiles = " <<  dTotCoarseToDepositAllProfiles << " dTotCoarseDepositedAllProfiles = " << dTotCoarseDepositedAllProfiles << endl;
   //
   //    // Have we deposited enough for this cliff collapse?
   //    if (bDoSandDepositionOnThisProfile)
   //    {
   //       if (bFPIsEqual(dTotSandToDepositAllProfiles - dTotSandDepositedAllProfiles, 0.0, MASS_BALANCE_TOLERANCE))
   //          bSandDepositionCompletedOnThisProfile = true;
   //    }
   //
   //    if (bDoCoarseDepositionOnThisProfile)
   //    {
   //       if (bFPIsEqual(dTotCoarseToDepositAllProfiles - dTotCoarseDepositedAllProfiles, 0.0, MASS_BALANCE_TOLERANCE))
   //          bCoarseDepositionCompletedOnThisProfile = true;
   //    }
   //
   //    if (bSandDepositionCompletedOnThisProfile && bCoarseDepositionCompletedOnThisProfile)
   //    {
   //       LogStream << m_ulIter << ": bSandDepositionCompletedOnThisProfile && bCoarseDepositionCompletedOnThisProfile for cliff collapse at [" << nXCliff << "][" << nYCliff << "] nStartPoint = " << nStartPoint << " nThisPoint = " << nThisPoint << endl << "\tnbJustDepositWhatWeCan = " << bJustDepositWhatWeCan << " SeawardOffset = " << nSeawardOffset << " dCliffHeightFraction = " << dCliffHeightFraction << " nDepProfile = " << nDepProfile << endl << "\tdTargetSandToDepositOnThisProfile = " << dTargetSandToDepositOnThisProfile << " dSandDepositedOnThisProfile = " << dSandDepositedOnThisProfile << " dTotSandToDepositAllProfiles = " << dTotSandToDepositAllProfiles << " dTotSandDepositedAllProfiles = " << dTotSandDepositedAllProfiles << endl << "\tdTargetCoarseToDepositOnThisProfile = " << dTargetCoarseToDepositOnThisProfile << " dCoarseDepositedOnThisProfile = " << dCoarseDepositedOnThisProfile << " dTotCoarseToDepositAllProfiles = " <<  dTotCoarseToDepositAllProfiles << " dTotCoarseDepositedAllProfiles = " << dTotCoarseDepositedAllProfiles << endl;
   //    }
   // } // Process each deposition profile
   //
   // // Safety check for sand sediment
   // if (! bFPIsEqual((dTotSandToDepositAllProfiles - dTotSandDepositedAllProfiles), 0.0, MASS_BALANCE_TOLERANCE))
   //    LogStream << ERR << m_ulIter << ": non-zero dTotSandToDepositAllProfiles = " << dTotSandToDepositAllProfiles << endl;
   //
   // // Ditto for coarse sediment
   // if (! bFPIsEqual((dTotCoarseToDepositAllProfiles - dTotCoarseDepositedAllProfiles), 0.0, MASS_BALANCE_TOLERANCE))
   //    LogStream << ERR << m_ulIter << ": non-zero dTotCoarseToDepositAllProfiles = " << dTotCoarseToDepositAllProfiles << endl;
   //
   // // Store the total depths of cliff collapse deposition for this polygon
   // pPolygon->AddCliffCollapseUnconsSandDeposition(dTotSandDepositedAllProfiles);
   // pPolygon->AddCliffCollapseUnconsCoarseDeposition(dTotCoarseDepositedAllProfiles);
   //
   // // Increment this-timestep totals for cliff collapse deposition
   // m_dThisIterUnconsSandCliffDeposition += dTotSandDepositedAllProfiles;
   // m_dThisIterUnconsCoarseCliffDeposition += dTotCoarseDepositedAllProfiles;
   //
   // // LogStream << endl;
   // // LogStream << "\tdTotSandToDepositAllProfiles = " << dTotSandToDepositAllProfiles << " dTotSandDepositedAllProfiles = " << dTotSandDepositedAllProfiles << endl;
   // // LogStream << "\tdTotCoarseToDepositAllProfiles = " << dTotCoarseToDepositAllProfiles << " dTotCoarseDepositedAllProfiles = " << dTotCoarseDepositedAllProfiles << endl;
   // // LogStream << endl << "****************************************" << endl << endl;

   return RTN_OK;
}

//===============================================================================================================================
//! Given the start and end points of a cliff-collapse normal profile, returns an output vector of cells which are 'under' the vector line
//===============================================================================================================================
void CSimulation::RasterizeCliffCollapseProfile(vector<CGeom2DPoint> const* pVPointsIn, vector<CGeom2DIPoint>* pVIPointsOut) const
{
   pVIPointsOut->clear();

   // The start point of the normal is the centroid of a coastline cell. Convert from the external CRS to grid CRS
   double const dXStart = dExtCRSXToGridX(pVPointsIn->at(0).dGetX());
   double const dYStart = dExtCRSYToGridY(pVPointsIn->at(0).dGetY());

   // The end point of the normal, again convert from the external CRS to grid CRS. Note too that it could be off the grid
   double const dXEnd = dExtCRSXToGridX(pVPointsIn->at(1).dGetX());
   double const dYEnd = dExtCRSYToGridY(pVPointsIn->at(1).dGetY());

   // Interpolate between cells by a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
   double dXInc = dXEnd - dXStart;
   double dYInc = dYEnd - dYStart;
   double const dLength = tMax(tAbs(dXInc), tAbs(dYInc));

   dXInc /= dLength;
   dYInc /= dLength;

   double dX = dXStart;
   double dY = dYStart;

   // Process each interpolated point
   int const nLength = nRound(dLength);

   for (int m = 0; m <= nLength; m++)
   {
      int nX = nRound(dX);
      int nY = nRound(dY);

      // Make sure the interpolated point is within the raster grid (can get this kind of problem due to rounding)
      if (! bIsWithinValidGrid(nX, nY))
         KeepWithinValidGrid(nRound(dXStart), nRound(dYStart), nX, nY);

      // This point is fine, so append it to the output vector
      pVIPointsOut->push_back(CGeom2DIPoint(nX, nY)); // In raster grid coordinates

      // And increment for next time
      dX += dXInc;
      dY += dYInc;
   }
}

