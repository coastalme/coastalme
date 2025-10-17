/*!
   \file assign_landforms.cpp
   \brief Assigns landform categories to coastlines and coastal cells, and to all other dryland cells
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "cliff.h"
#include "drift.h"
#include "intervention.h"

//===============================================================================================================================
//! Each timestep, classify coastal landforms and assign a coastal landform object to every point on every coastline. If, for a given cell, the coastal landform class has not changed then it inherits values from the previous timestep
//===============================================================================================================================
int CSimulation::nAssignLandformsForAllCoasts(void)
{
   // For each coastline, put a coastal landform at every point along the coastline
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nCoastPoint = 0; nCoastPoint < m_VCoast[nCoast].nGetCoastlineSize(); nCoastPoint++)
      {
         // Get the coords of the grid cell marked as coastline for the coastal landform object
         int const nX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetX();
         int const nY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetY();

         // Store the coastline number and the number of the coastline point in the cell so we can get these quickly later
         m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCoast(nCoast);
         m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetPointOnCoast(nCoastPoint);

         // OK, start assigning coastal landforms. First, is there an intervention here?
         if (bIsInterventionCell(nX, nY) || m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight() > 0)
         {
            // There is, so create an intervention object on the vector coastline with these attributes
            CACoastLandform* pIntervention = new CRWIntervention(&m_VCoast[nCoast], nCoast, nCoastPoint);
            m_VCoast[nCoast].AppendCoastLandform(pIntervention);

#if _DEBUG
            LogStream << nCoastPoint << " [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} " << m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() << " " << m_pRasterGrid->m_Cell[nX][nY].dGetInterventionHeight() << endl;
#endif

            continue;
         }

         // OK this landform is something other than an intervention. Is there talus on the cell?
         if (m_pRasterGrid->m_Cell[nX][nY].dGetTalusDepth() > 0)
         {
            // This is a talus cell
            m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFCategory(LF_CAT_DRIFT);
            m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_TALUS);

            CACoastLandform* pDrift = new CRWDrift(&m_VCoast[nCoast], nCoast, nCoastPoint);
            m_VCoast[nCoast].AppendCoastLandform(pDrift);

            continue;
         }

         // Not an intervention, not talus: do some safety checks. Note that layer 0 is the first layer above basement
         int const nLayer = m_pRasterGrid->m_Cell[nX][nY].nGetLayerAtElev(m_dThisIterSWL);
         if (nLayer == ELEV_IN_BASEMENT)
         {
            // Should never happen
            LogStream << m_ulIter << ": SWL (" << m_dThisIterSWL << ") is in basement on cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, cannot assign coastal landform for coastline " << nCoast << endl;

            return RTN_ERR_CANNOT_ASSIGN_COASTAL_LANDFORM;
         }

         if (nLayer == ELEV_ABOVE_SEDIMENT_TOP)
         {
            // Again, should never happen
            LogStream << m_ulIter << ": SWL (" << m_dThisIterSWL << ") is above sediment-top elevation (" << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevOmitTalus() << ") on cell [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, cannot assign coastal landform for coastline " << nCoast << endl;

            return RTN_ERR_CANNOT_ASSIGN_COASTAL_LANDFORM;
         }

         // Next, check what we have at SWL on this cell: is it unconsolidated or consolidated sediment?
         double const dConsSedTop = m_pRasterGrid->m_Cell[nX][nY].dGetConsSedTopElevForLayerAboveBasement(nLayer);
         if (dConsSedTop >= m_dThisIterSWL)
         {
            // We have consolidated sediment at or above SWL on this cell. Are we considering cliff collapse?
            if (m_bDoCliffCollapse)
            {
               // OK we are considering cliff collapse, and we have consolidated sediment at SWL, so this is a cliff cell. Get the existing landform category for this cell
               if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() == LF_CAT_CLIFF)
               {
                  // This cell was a cliff in the previous timestep, so get the pre-existing data stored in the cell
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_CLIFF_ON_COASTLINE);
                  double dAccumWaveEnergy = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetAccumWaveEnergy();
                  double dNotchIncision = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetCliffNotchIncisionDepth();
                  double dNotchApexElev = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetCliffNotchApexElev();

                  // Create a cliff object on the vector coastline with these attributes
                  CACoastLandform* pCliff = new CRWCliff(&m_VCoast[nCoast], nCoast, nCoastPoint, m_dCellSide, dNotchIncision, dNotchApexElev, dAccumWaveEnergy);
                  m_VCoast[nCoast].AppendCoastLandform(pCliff);

#if _DEBUG
                  double const dSedTopElevIncTalus = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevIncTalus();
                  double const dSedTopElevNoTalus = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevOmitTalus();

                  assert(dNotchApexElev < dSedTopElevNoTalus);

                  LogStream << m_ulIter << ": \tcontinues to be a cliff at [" << nX << "][" << nY << "] dAccumWaveEnergy = " << dAccumWaveEnergy << " dNotchApexElev = " << dNotchApexElev << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << " dSedTopElevIncTalus = " << dSedTopElevIncTalus << " dNotchIncision = " << dNotchIncision << endl;
#endif
               }
               else
               {
                  // This was not a cliff in the previous timestep, but it is now
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFCategory(LF_CAT_CLIFF);
                  m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_CLIFF_ON_COASTLINE);

                  // Get the pre-existing wave energy stored in the cell
                  double dAccumWaveEnergy = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetAccumWaveEnergy();

                  // The DBL_NODATA Values indicate that the cliff object does not have an incised notch
                  double dNotchIncision = DBL_NODATA;
                  double dNotchApexElev = DBL_NODATA;

                  // Would the new notch apex elevation be below the top of the cell's consolidated sediment?
                  double const dSedTopElevNoTalus = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevOmitTalus();
                  if (m_dThisIterNewNotchApexElev < dSedTopElevNoTalus)
                  {
                     // Yes it would, so this cliff object has a notch
                     dNotchIncision = 0;
                     dNotchApexElev = m_dThisIterNewNotchApexElev;
                  }

                  // Create a cliff object on the vector coastline with these attributes
                  CACoastLandform* pCliff = new CRWCliff(&m_VCoast[nCoast], nCoast, nCoastPoint, m_dCellSide, dNotchIncision, dNotchApexElev, dAccumWaveEnergy);
                  m_VCoast[nCoast].AppendCoastLandform(pCliff);

#if _DEBUG
                  double const dSedTopElevIncTalus = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevIncTalus();

                  if (bFPIsEqual(dNotchIncision, 0.0, TOLERANCE))
                     LogStream << m_ulIter << ": \tcliff created at [" << nX << "][" << nY << "] dAccumWaveEnergy = " << dAccumWaveEnergy << " dNotchApexElev = " << dNotchApexElev << " dNotchTopElev = " << DBL_NODATA << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << " dSedTopElevIncTalus = " << dSedTopElevIncTalus << " dNotchIncision = " << dNotchIncision << endl;
                  else
                     LogStream << m_ulIter << ": \tNO NOTCH cliff created at [" << nX << "][" << nY << "] dAccumWaveEnergy = " << dAccumWaveEnergy << " dSedTopElevNoTalus = " << dSedTopElevNoTalus << " dSedTopElevIncTalus = " << dSedTopElevIncTalus << endl;
#endif
               }
            }
            else
            {
               // We have consolidated sediment at SWL but we are not considering cliff collapse. Get the pre-existing wave energy stored in the cell
               double dAccumWaveEnergy = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->dGetAccumWaveEnergy();

               // Create a cliff object on the vector coastline
               CACoastLandform* pCliff = new CRWCliff(&m_VCoast[nCoast], nCoast, nCoastPoint, m_dCellSide, DBL_NODATA, DBL_NODATA, dAccumWaveEnergy);
               m_VCoast[nCoast].AppendCoastLandform(pCliff);

#if _DEBUG
               LogStream << m_ulIter << ": \tcliff created (cliff collapse not considered) at " << nX << "][" << nY << "] dAccumWaveEnergy = " << dAccumWaveEnergy << " dSedTopElevNoTalus = " << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevOmitTalus() << " dSedTopElevIncTalus = " << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevIncTalus() << endl;
#endif
            }
         }
         else
         {
            // We have unconsolidated sediment at SWL, so this is a drift cell: create a drift object on the vector coastline with these attributes
            CACoastLandform* pDrift = new CRWDrift(&m_VCoast[nCoast], nCoast, nCoastPoint);
            m_VCoast[nCoast].AppendCoastLandform(pDrift);

            // Safety check
            if (m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory() != LF_CAT_DRIFT)
               m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_DRIFT_MIXED);

#if _DEBUG
            LogStream << m_ulIter << ": drift created at [" << nX << "][" << nY << "]" << endl;
#endif
         }
      }
   }

   // // DEBUG CODE ============================================================================================================================================
   // for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
   // {
   //    for (int j = 0; j < m_VCoast[i].nGetCoastlineSize(); j++)
   //    {
   //       int nX = m_VCoast[i].pPtiGetCellMarkedAsCoastline(j)->nGetX();
   //       int nY = m_VCoast[i].pPtiGetCellMarkedAsCoastline(j)->nGetY();
   //
   //       LogStream << m_ulIter << ": coast cell " << j << " at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} has landform category = ";
   //
   //       int nCat = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFCategory();
   //       int nSubCat = m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->nGetLFSubCategory();
   //
   //       switch (nCat)
   //       {
   //          case LF_CAT_HINTERLAND:
   //             LogStream << "hinterland";
   //             break;
   //
   //          case LF_CAT_SEA:
   //             LogStream << "sea";
   //             break;
   //
   //          case LF_CAT_CLIFF:
   //             LogStream << "cliff";
   //             break;
   //
   //          case LF_CAT_DRIFT:
   //             LogStream << "drift";
   //             break;
   //
   //          case LF_CAT_INTERVENTION:
   //             LogStream << "intervention";
   //             break;
   //
   //          case LF_NONE:
   //             LogStream << "none";
   //             break;
   //
   //          default:
   //             LogStream << "NONE";
   //             break;
   //       }
   //
   //       LogStream << " and landform subcategory = ";
   //
   //       switch (nSubCat)
   //       {
   //          case LF_SUBCAT_CLIFF_ON_COASTLINE:
   //             LogStream << "cliff on coastline";
   //             break;
   //
   //          case LF_SUBCAT_CLIFF_INLAND:
   //             LogStream << "cliff inland";
   //             break;
   //
   //          case LF_SUBCAT_DRIFT_MIXED:
   //             LogStream << "mixed drift";
   //             break;
   //
   //          case LF_SUBCAT_DRIFT_TALUS:
   //             LogStream << "talus";
   //             break;
   //
   //          case LF_SUBCAT_DRIFT_BEACH:
   //             LogStream << "beach";
   //             break;
   //
   //          case LF_SUBCAT_DRIFT_DUNES:
   //             LogStream << "dunes";
   //             break;
   //
   //          case LF_SUBCAT_INTERVENTION_STRUCT:
   //             LogStream << "structural intervention";
   //             break;
   //
   //          case LF_SUBCAT_INTERVENTION_NON_STRUCT:
   //             LogStream << "non-structural intervention";
   //             break;
   //
   //          case LF_NONE:
   //             LogStream << "none";
   //             break;
   //
   //          default:
   //             LogStream << "NONE";
   //             break;
   //       }
   //       LogStream << endl;
   //    }
   // }
   // // DEBUG CODE ============================================================================================================================================

   return RTN_OK;
}

//===============================================================================================================================
//! At the end of each timestep, this routine stores the attributes from a single coastal landform object in the grid cell 'under' the object, ready for the next timestep
//===============================================================================================================================
int CSimulation::nLandformToGrid(int const nCoast, int const nPoint)
{
   // What is the coastal landform here?
   CACoastLandform* pCoastLandform = m_VCoast[nCoast].pGetCoastLandform(nPoint);
   int const nCategory = pCoastLandform->nGetLandFormCategory();
   if (nCategory == LF_CAT_CLIFF)
   {
      // It's a cliff
      CRWCliff const* pCliff = reinterpret_cast<CRWCliff*>(pCoastLandform);

      int const nX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nPoint)->nGetX();
      int const nY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nPoint)->nGetY();

      // if (! pCliff->bHasCollapsed())
      // {
      //    // The cliff has not collapsed. Get attribute values from the cliff object
      //    double const dNotchBaseElev = pCliff->dGetNotchApexElev();
      //    double const dNotchIncision = pCliff->dGetNotchIncision();
      //
      //    // And store some attribute values in the cliff cell
      //    m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetLFSubCategory(LF_SUBCAT_CLIFF_ON_COASTLINE);
      //    m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchApexElev(dNotchBaseElev);
      //    m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchIncisionDepth(dNotchIncision);
      // }
      // else
      // {
      //    // // The cliff has collapsed: all sediment above the base of the erosional notch is gone from this cliff object via cliff collapse, so this cell is no longer a cliff
      //    // m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();
      //    //
      //    // // Check the x-y extremities of the contiguous sea for the bounding box (used later in wave propagation)
      //    // if (nX < m_nXMinBoundingBox)
      //    //    m_nXMinBoundingBox = nX;
      //    //
      //    // if (nX > m_nXMaxBoundingBox)
      //    //    m_nXMaxBoundingBox = nX;
      //    //
      //    // if (nY < m_nYMinBoundingBox)
      //    //    m_nYMinBoundingBox = nY;
      //    //
      //    // if (nY > m_nYMaxBoundingBox)
      //    //    m_nYMaxBoundingBox = nY;
      //
      //    int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetNumOfTopLayerAboveBasement();
      //
      //    // Safety check
      //    if (nTopLayer == INT_NODATA)
      //       return RTN_ERR_NO_TOP_LAYER;
      //
      //    // Update the cell's layer elevations
      //    m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();
      //
      //    // And update the cell's sea depth
      //    m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();
      // }

      // Always accumulate wave energy
      m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetAccumWaveEnergy(pCliff->dGetTotAccumWaveEnergy());
   }
   else if (nCategory == LF_CAT_DRIFT)
   {
      // It's drift, so calculate D50 TODO 002 Why might we need this?
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Each timestep, classify landforms for cells that are not on the coastline
//===============================================================================================================================
int CSimulation::nAssignLandformsForAllCells(void)
{
   // First pass: collect information about cells that need to be changed. This avoids race conditions from reading neighbour cells while writing to current cells
   vector<vector<int>> vCellUpdates(m_nXGridSize, vector<int>(m_nYGridSize, -1));

   // Read-only phase: determine what changes need to be made
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif

   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         // Get this cell's landform category
         CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetLandform();
         int const nCat = pLandform->nGetLFCategory();

         // Store what action to take (to avoid writing during read phase)
         int nAction = -1; // -1 = no change, others defined below

         if (m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue())
         {
            nAction = 0; // Set to unknown landform
         }
         else if ((nCat == LF_CAT_SEDIMENT_INPUT) || (nCat == LF_CAT_SEDIMENT_INPUT_SUBMERGED) || (nCat == LF_CAT_SEDIMENT_INPUT_NOT_SUBMERGED))
         {
            if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
               nAction = 1; // Set to submerged sediment input
            else
               nAction = 2; // Set to not submerged sediment input
         }
         else if (nCat == LF_CAT_SEA)
         {
            // This is a sea cell. Is it surrounded by drift cells, or drift and cliff?
            if (bSurroundedByDriftCells(nX, nY))
               nAction = 3; // Set to beach
            else if (m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevOmitTalus() > m_dThisIterSWL)
               nAction = 4; // Set to island
            // else keep as sea (no action needed)
         }
         else if (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
         {
            // Is not coastline
            if (nCat == LF_CAT_CLIFF)
               nAction = 5; // Set to former cliff
            else if ((nCat != LF_CAT_DRIFT) && (nCat != LF_CAT_INTERVENTION))
               nAction = 6; // Set to hinterland
         }

         vCellUpdates[nX][nY] = nAction;
      }
   }

   // Write phase: apply the changes
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         int const nAction = vCellUpdates[nX][nY];

         if (nAction == -1)
            continue; // No change

         CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetLandform();

         switch (nAction)
         {
         case 0: // Set to unknown landform
            pLandform->SetLFCategory(LF_NONE);
            break;

         case 1: // Set to submerged sediment input
            pLandform->SetLFCategory(LF_CAT_SEDIMENT_INPUT_SUBMERGED);
            m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();
            break;

         case 2: // Set to not submerged sediment input
            pLandform->SetLFCategory(LF_CAT_SEDIMENT_INPUT_NOT_SUBMERGED);
            break;

         case 3: // Set to beach
            pLandform->SetLFSubCategory(LF_SUBCAT_DRIFT_BEACH);
            break;

         case 4: // Set to island
            pLandform->SetLFCategory(LF_CAT_ISLAND);
            break;

         case 5: // Set to former cliff
            pLandform->SetLFSubCategory(LF_SUBCAT_CLIFF_INLAND);
            break;

         case 6: // Set to hinterland
            pLandform->SetLFCategory(LF_CAT_HINTERLAND);
            break;
         }
      }
   }

   return RTN_OK;
}
//===============================================================================================================================
//! Returns true if this cell has four drift cells surrounding it
//===============================================================================================================================
bool CSimulation::bSurroundedByDriftCells(int const nX, int const nY)
{
   int nXTmp;
   int nYTmp;
   int nAdjacent = 0;

   // North
   nXTmp = nX;
   nYTmp = nY - 1;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLFCategory();

      if ((nCat == LF_CAT_DRIFT) || (nCat == LF_CAT_CLIFF))
         nAdjacent++;
   }

   // East
   nXTmp = nX + 1;
   nYTmp = nY;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLFCategory();

      if ((nCat == LF_CAT_DRIFT) || (nCat == LF_CAT_CLIFF))
         nAdjacent++;
   }

   // South
   nXTmp = nX;
   nYTmp = nY + 1;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLFCategory();

      if ((nCat == LF_CAT_DRIFT) || (nCat == LF_CAT_CLIFF))
         nAdjacent++;
   }

   // West
   nXTmp = nX - 1;
   nYTmp = nY;

   if (bIsWithinValidGrid(nXTmp, nYTmp))
   {
      CRWCellLandform const* pLandform = m_pRasterGrid->m_Cell[nXTmp][nYTmp].pGetLandform();
      int const nCat = pLandform->nGetLFCategory();

      if ((nCat == LF_CAT_DRIFT) || (nCat == LF_CAT_CLIFF))
         nAdjacent++;
   }

   if (nAdjacent == 4)
   {
      // This cell has four LF_CAT_DRIFT neighbours
      return true;
   }

   return false;
}
