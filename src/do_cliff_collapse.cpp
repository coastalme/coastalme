/*!
   \file do_cliff_collapse.cpp
   \brief Collapses cliffs if a critical notch depth is exceeded
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

//===============================================================================================================================
//! Update accumulated wave energy in coastal landform objects. If the object is a cliff, then deepen the incised notch. If the notch is sufficiently deep, cliff collapse occurs.
//! CoastalME's representation of notch incision is based on Trenhaile, A.S. (2015). Coastal notches: Their morphology, formation, and function. Earth-Science Reviews 150, 285-304
//===============================================================================================================================
int CSimulation::nDoAllWaveEnergyToCoastLandforms(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Calculating cliff collapse" << endl;

   int nRet;

   // First go along each coastline and at each point on the coastline, update the total wave energy which it has experienced TODO Note that currently, only cliff objects respond to accumulated wave energy
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      for (int nCoastPoint = 0; nCoastPoint < m_VCoast[nCoast].nGetCoastlineSize(); nCoastPoint++)
      {
         CACoastLandform* pCoastLandform = m_VCoast[nCoast].pGetCoastLandform(nCoastPoint);

         // First get wave energy for the coastal landform object
         double const dWaveHeightAtCoast = m_VCoast[nCoast].dGetCoastWaveHeight(nCoastPoint);

         // If the waves at this point are off-shore, then do nothing, just move to next coast point
         if (bFPIsEqual(dWaveHeightAtCoast, DBL_NODATA, TOLERANCE))
            continue;

         // OK we have on-shore waves so get the previously-calculated wave energy
         double const dWaveEnergy = m_VCoast[nCoast].dGetWaveEnergyAtBreaking(nCoastPoint);
         // assert(isfinite(dWaveEnergy));

         // And save the accumulated value
         pCoastLandform->IncTotAccumWaveEnergy(dWaveEnergy);

         // int const nX = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetX();
         // int const nY = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint)->nGetY();
         int const nCat = pCoastLandform->nGetLandFormCategory();
         // double const dTopElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevOmitTalus();

         // Is this a cliff?
         if ((nCat == LF_CAT_CLIFF) || (nCat == LF_SUBCAT_CLIFF_ON_COASTLINE) || (nCat == LF_SUBCAT_CLIFF_INLAND))
         {
            // It is, so get the cliff object
            CRWCliff* pCliff = reinterpret_cast<CRWCliff*>(pCoastLandform);

            // And do the notch incision, if any. Note that we consider sediment eroded due to notch incision to be still in place until cliff collapse, i.e. the sediment which filled the notch, pre-incision, is assumed to remain there. If the notch is eventually incised sufficiently to cause cliff collapse, then the sediment from the notch volume is included with the above-notch talus
            IncreaseCliffNotchIncision(pCliff, dWaveEnergy);

            // OK, is the notch now extended enough to cause collapse (either because the overhang is greater than the threshold overhang, or because there is no sediment remaining)?
            if (pCliff->bReadyToCollapse(m_dNotchIncisionAtCollapse))
            {
               // // DEBUG CODE ============================================================================================================================================
               // // Get total depths of sand consolidated and unconsolidated for every cell
               // if (m_ulIter == 5)
               // {
               // double dTmpSandCons = 0;
               // double dTmpSandUncons = 0;
               // for (int nX1 = 0; nX1 < m_nXGridSize; nX1++)
               // {
               // for (int nY1 = 0; nY1 < m_nYGridSize; nY1++)
               // {
               // dTmpSandCons += m_pRasterGrid->m_Cell[nX1][nY1].dGetConsSandDepthAllLayers();
               //
               // dTmpSandUncons += m_pRasterGrid->m_Cell[nX1][nY1].dGetUnconsSandDepthAllLayers();
               // }
               // }
               //
               //    // Get the cliff cell's grid coords
               // int nXCliff = pCliff->pPtiGetCellMarkedAsCliff()->nGetX();
               // int nYCliff = pCliff->pPtiGetCellMarkedAsCliff()->nGetY();
               //
               //    // Get this cell's polygon
               // int nPoly = m_pRasterGrid->m_Cell[nXCliff][nYCliff].nGetPolygonID();
               //
               // LogStream << endl;
               // LogStream << "*****************************" << endl;
               // LogStream << m_ulIter << ": before cliff collapse on nPoly = " << nPoly << " total consolidated sand = " << dTmpSandCons * m_dCellArea << " total unconsolidated sand = " << dTmpSandUncons * m_dCellArea << endl;
               // }
               // // DEBUG CODE ============================================================================================================================================

               // It is ready to collapse
               int nNotchLayer;
               double dCliffElevPreCollapse = 0;
               double dCliffElevPostCollapse = 0;
               double dFineCollapse = 0;
               double dSandCollapse = 0;
               double dCoarseCollapse = 0;

               // So do the cliff collapse
               nRet = nDoCliffCollapse(nCoast, pCliff, dFineCollapse, dSandCollapse, dCoarseCollapse, nNotchLayer, dCliffElevPreCollapse, dCliffElevPostCollapse);
               if (nRet != RTN_OK)
               {
                  if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
                     LogStream << m_ulIter << ": " << WARN << "problem with cliff collapse, continuing however" << endl;
               }

               // Deposit all sand and/or coarse sediment derived from this cliff collapse as talus, on the cell on which collapse occurred
               nRet = nDoCliffCollapseTalusDeposition(nCoast, pCliff, dSandCollapse, dCoarseCollapse, nNotchLayer);
               if (nRet != RTN_OK)
                  return nRet;

               // // DEBUG CODE ============================================================================================================================================
               // // Get total depths of sand consolidated and unconsolidated for every cell
               // if (m_ulIter == 5)
               // {
               // double dTmpSandCons = 0;
               // double dTmpSandUncons = 0;
               // for (int nX1 = 0; nX1 < m_nXGridSize; nX1++)
               // {
               // for (int nY1 = 0; nY1 < m_nYGridSize; nY1++)
               // {
               // dTmpSandCons += m_pRasterGrid->m_Cell[nX1][nY1].dGetConsSandDepthAllLayers();
               //
               // dTmpSandUncons += m_pRasterGrid->m_Cell[nX1][nY1].dGetUnconsSandDepthAllLayers();
               // }
               // }
               //
               //    // Get the cliff cell's grid coords
               // int nXCliff = pCliff->pPtiGetCellMarkedAsCliff()->nGetX();
               // int nYCliff = pCliff->pPtiGetCellMarkedAsCliff()->nGetY();
               //
               //    // Get this cell's polygon
               // int nPoly = m_pRasterGrid->m_Cell[nXCliff][nYCliff].nGetPolygonID();
               //
               // LogStream << endl;
               // LogStream << "*****************************" << endl;
               // LogStream << m_ulIter << ": after cliff collapse on nPoly = " << nPoly << " total consolidated sand = " << dTmpSandCons * m_dCellArea << " total unconsolidated sand = " << dTmpSandUncons * m_dCellArea << endl;
               // LogStream << m_ulIter << ": total consolidated sand lost this iteration =  " << (m_dStartIterConsSandAllCells - dTmpSandCons) * m_dCellArea << endl;
               // LogStream << m_ulIter << ": total unconsolidated sand added this iteration =  " << (dTmpSandUncons - m_dStartIterUnconsSandAllCells) * m_dCellArea << endl;
               //
               // double dTmpAllPolySandErosion = 0;
               // double dTmpAllPolySandDeposition = 0;
               // for (unsigned int n = 0; n < m_pVCoastPolygon.size(); n++)
               // {
               // double dTmpSandErosion = m_pVCoastPolygon[n]->dGetCliffCollapseErosionSand() * m_dCellArea ;
               // double dTmpSandDeposition = m_pVCoastPolygon[n]->dGetCliffCollapseUnconsSandDeposition() * m_dCellArea ;
               //
               // LogStream << m_ulIter << ": polygon = " << m_pVCoastPolygon[n]->nGetPolygonCoastID() << " sand erosion = " << dTmpSandErosion << " sand deposition = " << dTmpSandDeposition << endl;
               //
               // dTmpAllPolySandErosion += dTmpSandErosion;
               // dTmpAllPolySandDeposition += dTmpSandDeposition;
               // }
               //
               // LogStream << "-------------------------------------------" << endl;
               // LogStream << m_ulIter << ": all polygons, sand erosion = " << dTmpAllPolySandErosion << " sand deposition = " << dTmpAllPolySandDeposition << endl;
               // LogStream << "*****************************" << endl;
               // }
               // // DEBUG CODE ============================================================================================================================================
            }
         }
      }
   }

   if (m_nLogFileDetail >= LOG_FILE_ALL)
      LogStream << m_ulIter << ": \ttotal cliff collapse (m^3) = " << (m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons + m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons + m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea << " (fine = " << (m_dThisIterCliffCollapseErosionFineUncons + m_dThisIterCliffCollapseErosionFineCons) * m_dCellArea << ", sand = " << (m_dThisIterCliffCollapseErosionSandUncons + m_dThisIterCliffCollapseErosionSandCons) * m_dCellArea << ", coarse = " << (m_dThisIterCliffCollapseErosionCoarseUncons + m_dThisIterCliffCollapseErosionCoarseCons) * m_dCellArea << "), talus deposition (m^3) = " << (m_dThisIterUnconsSandCliffDeposition + m_dThisIterUnconsCoarseCliffDeposition) * m_dCellArea << " (sand = " << m_dThisIterUnconsSandCliffDeposition * m_dCellArea << ", coarse = " << m_dThisIterUnconsSandCliffDeposition * m_dCellArea << ")" << endl << endl;

   return RTN_OK;
}

//===============================================================================================================================
//! Simulates cliff collapse on a single cell. Collapse happens when when a notch which is incised into the cell's consolidated sediment layer exceeds a critical horizontal incision. This routine updates the cliff object, the cell 'under' the cliff object, and the polygon which contains the cliff object
//===============================================================================================================================
int CSimulation::nDoCliffCollapse(int const nCoast, CRWCliff* pCliff, double& dFineCollapse, double& dSandCollapse, double& dCoarseCollapse, int& nNotchLayer, double& dPreCollapseCellElev, double& dPostCollapseCellElevNoTalus)
{
   // Get the cliff cell's grid coords
   int const nX = pCliff->pPtiGetCellMarkedAsCliff()->nGetX();
   int const nY = pCliff->pPtiGetCellMarkedAsCliff()->nGetY();

   // Get this cell's polygon
   int const nPoly = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
   if (nPoly == INT_NODATA)
   {
      // This cell isn't in a polygon
      LogStream << m_ulIter << ": " << WARN << "in nDoCliffCollapse(), [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} is not in a polygon" << endl;
      return RTN_ERR_CLIFF_NOT_IN_POLYGON;
   }

   // Get a pointer to the polygon
   CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

   // Get the elevation of the apex of the notch from the cliff object
   double const dNotchElev = pCliff->dGetNotchApexElev();

   // Get the index of the layer containing the notch (layer 0 being just above basement)
   nNotchLayer = m_pRasterGrid->m_Cell[nX][nY].nGetLayerAtElev(dNotchElev);

   // Safety check: is the notch elevation above the top of the sediment? If so, do no more
   if (nNotchLayer == ELEV_ABOVE_SEDIMENT_TOP)
   {
      LogStream << m_ulIter << ": in nDoCliffCollapse(), [" << nX << "][" << nY << "] has nNotchLayer above sediment top" << endl;
      return RTN_OK;
   }

   // More safety checks
   if (nNotchLayer == ELEV_IN_BASEMENT)
   {
      LogStream << m_ulIter << ": " << WARN << "in nDoCliffCollapse(), [" << nX << "][" << nY << "] nNotchLayer is in basement" << endl;
      return RTN_ERR_CLIFF_NOTCH;
   }

   if (nNotchLayer < 0)
   {
      LogStream << m_ulIter << ": " << WARN << "in nDoCliffCollapse(), [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "} notch layer = " << nNotchLayer << ", dNotchElev = " << dNotchElev << " m_dNotchApexAboveMHW = " << m_dNotchApexAboveMHW << " dPreCollapseCellElev = " << dPreCollapseCellElev << endl;
      return RTN_ERR_CLIFF_NOTCH;
   }

   // Notch layer is OK, so flag the coastline cliff object as having collapsed
   pCliff->SetCliffCollapsed();

   int const nTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();

   // Safety check
   if (nTopLayer == INT_NODATA)
   {
      LogStream << m_ulIter << ": " << WARN << "in nDoCliffCollapse(), [" << nX << "][" << nY << "] nTopLayer = " << nTopLayer << endl;
      return RTN_ERR_NO_TOP_LAYER;
   }

   // Set flags to say that the notch layer, and all layers above it, have changed
   for (int nLayer = nNotchLayer; nLayer <= nTopLayer; nLayer++)
   {
      m_bConsChangedThisIter[nLayer] = true;
      m_bUnconsChangedThisIter[nLayer] = true;
   }

   // Get the pre-collapse cliff elevation
   dPreCollapseCellElev = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevOmitTalus();

   // Now calculate the vertical depth of sediment lost in this cliff collapse, note that this includes the sediment which filled the notch before any incision took place
   double dAvailable = 0;
   double dFineConsLost = 0;
   double dFineUnconsLost = 0;
   double dSandConsLost = 0;
   double dSandUnconsLost = 0;
   double dCoarseConsLost = 0;
   double dCoarseUnconsLost = 0;

   // Update the cell's sediment. If there are sediment layers above the notched layer, we must remove sediment from the whole depth of each layer
   for (int n = nTopLayer; n > nNotchLayer; n--)
   {
      // Start with the unconsolidated sediment
      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetFineDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchFineLost();

      if (dAvailable > 0)
      {
         dFineCollapse += dAvailable;
         dFineUnconsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetFineDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetNotchFineLost(0);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetSandDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchSandLost();

      if (dAvailable > 0)
      {
         dSandCollapse += dAvailable;
         dSandUnconsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetSandDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetNotchSandLost(0);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetCoarseDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();

      if (dAvailable > 0)
      {
         dCoarseCollapse += dAvailable;
         dCoarseUnconsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetCoarseDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetUnconsolidatedSediment()->SetNotchCoarseLost(0);
      }

      // Now get the consolidated sediment
      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetFineDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchFineLost();

      if (dAvailable > 0)
      {
         dFineCollapse += dAvailable;
         dFineConsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetFineDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetNotchFineLost(0);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetSandDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchSandLost();

      if (dAvailable > 0)
      {
         dSandCollapse += dAvailable;
         dSandConsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetSandDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetNotchSandLost(0);
      }

      dAvailable = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetCoarseDepth() - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->dGetNotchCoarseLost();

      if (dAvailable > 0)
      {
         dCoarseCollapse += dAvailable;
         dCoarseConsLost += dAvailable;
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetCoarseDepth(0);
         m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(n)->pGetConsolidatedSediment()->SetNotchCoarseLost(0);
      }
   }

   // Now calculate the sediment lost from the consolidated layer into which the erosional notch was incised
   double const dNotchLayerTop = m_pRasterGrid->m_Cell[nX][nY].dCalcLayerElev(nNotchLayer);
   double const dNotchLayerThickness = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->dGetTotalThickness();
   double const dNotchLayerFracRemoved = (dNotchLayerTop - dNotchElev) / dNotchLayerThickness;

   // Sort out the notched layer's sediment, both consolidated and unconsolidated, for this cell. First the unconsolidated sediment
   double dFineDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetFineDepth();
   dAvailable = dFineDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchFineLost();

   if (dAvailable > 0)
   {
      // Some unconsolidated fine sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      dFineUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetFineDepth(dFineDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetNotchFineLost(0);
   }

   double dSandDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
   dAvailable = dSandDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchSandLost();

   if (dAvailable > 0)
   {
      // Some unconsolidated sand sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      dSandUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetSandDepth(dSandDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetNotchSandLost(0);
   }

   double dCoarseDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();
   dAvailable = dCoarseDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->dGetNotchCoarseLost();

   if (dAvailable > 0)
   {
      // Some unconsolidated coarse sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      dCoarseUnconsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetCoarseDepth(dCoarseDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetUnconsolidatedSediment()->SetNotchCoarseLost(0);
   }

   // Do the same for the consolidated sediment
   dFineDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetFineDepth();
   dAvailable = dFineDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchFineLost();

   if (dAvailable > 0)
   {
      // Some consolidated fine sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dFineCollapse += dLost;
      dFineConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetFineDepth(dFineDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetNotchFineLost(0);
   }

   dSandDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetSandDepth();
   dAvailable = dSandDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchSandLost();

   if (dAvailable > 0)
   {
      // Some consolidated sand sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dSandCollapse += dLost;
      dSandConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetSandDepth(dSandDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetNotchSandLost(0);
   }

   dCoarseDepth = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetCoarseDepth();
   dAvailable = dCoarseDepth - m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->dGetNotchCoarseLost();

   if (dAvailable > 0)
   {
      // Some consolidated coarse sediment is available for collapse
      double const dLost = dAvailable * dNotchLayerFracRemoved;
      dCoarseCollapse += dLost;
      dCoarseConsLost += dLost;
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetCoarseDepth(dCoarseDepth - dLost);
      m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer)->pGetConsolidatedSediment()->SetNotchCoarseLost(0);
   }

   // Update the cell's totals for cliff collapse erosion
   m_pRasterGrid->m_Cell[nX][nY].IncrCliffCollapseErosion(dFineCollapse, dSandCollapse, dCoarseCollapse);

   // Update the cell's layer elevations (pre talus deposition) and d50
   m_pRasterGrid->m_Cell[nX][nY].CalcAllLayerElevsAndD50();

   // Get the post-collapse cell elevation (talus will be deposited above this)
   dPostCollapseCellElevNoTalus = m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevOmitTalus();

   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": cliff collapse at [" << nX << "][" << nY << "], before talus deposition: original cell elevation = " << dPreCollapseCellElev << ", new cell elevation = " << dPostCollapseCellElevNoTalus << ", change in elevation = " << dPreCollapseCellElev - dPostCollapseCellElevNoTalus << endl;

   // And update the cell's sea depth
   // m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

   // LogStream << m_ulIter << ": cell [" << nX << "][" << nY << "] after removing sediment, dGetVolEquivSedTopElev() = " << m_pRasterGrid->m_Cell[nX][nY].dGetVolEquivSedTopElev() << ", dGetSedimentTopElevOmitTalus() = " << m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElevOmitTalus() << endl << endl;

   // Update this-polygon totals: add to the depths of cliff collapse erosion for this polygon
   pPolygon->AddCliffCollapseErosionFine(dFineCollapse);
   pPolygon->AddCliffCollapseToSuspensionFine(dFineCollapse);
   pPolygon->AddCliffCollapseErosionSand(dSandCollapse);
   pPolygon->AddCliffCollapseErosionCoarse(dCoarseCollapse);

   // And update the this-timestep totals and the grand totals for the number of cells with cliff collapse
   m_nNThisIterCliffCollapse++;
   m_nNTotCliffCollapse++;

   // Add to this-iteration totals of fine sediment (consolidated and unconsolidated) eroded via cliff collapse
   m_dThisIterCliffCollapseErosionFineUncons += dFineUnconsLost;
   m_dThisIterCliffCollapseErosionFineCons += dFineConsLost;

   // Also add to the total suspended load. Note that this addition to the suspended load has not yet been added to all cells, this happens in nUpdateGrid()
   m_dThisIterFineSedimentToSuspension += (dFineConsLost + dFineUnconsLost);

   // Add to this-iteration totals of sand and coarse sediment (consolidated and unconsolidated) eroded via cliff collapse
   m_dThisIterCliffCollapseErosionSandUncons += dSandUnconsLost;
   m_dThisIterCliffCollapseErosionSandCons += dSandConsLost;
   m_dThisIterCliffCollapseErosionCoarseUncons += dCoarseUnconsLost;
   m_dThisIterCliffCollapseErosionCoarseCons += dCoarseConsLost;

   // Save the timestep at which cliff collapse occurred
   m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffCollapseTimestep(m_ulIter);

   // Reset cell cliff info
   m_pRasterGrid->m_Cell[nX][nY].pGetLandform()->SetCliffNotchIncisionDepth(m_dCellSide);

   // Final safety check
   int const nNewTopLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopLayerAboveBasement();
   if (nNewTopLayer == INT_NODATA)
      return RTN_ERR_NO_TOP_LAYER;

   return RTN_OK;
}

//===============================================================================================================================
//! Calculate the incision (if any) of a cliff notch, assuming a linear decrease in incision with distance downwards from notch apex. Then do the incision
//===============================================================================================================================
void CSimulation::IncreaseCliffNotchIncision(CRWCliff* pCliff, double const dWaveEnergy)
{
   // Get the apex elevation of the cliff notch, and the cut-off elevation (if thi-iteration SWL is below this, there is no incision)
   double const dNotchElev = pCliff->dGetNotchApexElev();
   double const dCutoffElev = dNotchElev - CLIFF_NOTCH_CUTOFF_DISTANCE;

   // Get the coast number, and the coastline point, of the cliff
   int const nCoast = pCliff->nGetCoast();
   int const nCoastPoint = pCliff->nGetPointOnCoast();

   // And get the wave runup at this point
   double const dRunup = m_VCoast[nCoast].dGetRunUp(nCoastPoint);
   double const dWaveElev = m_dThisIterSWL + dRunup;

   if (dWaveElev < dCutoffElev )
      // SWL is below the cutoff elevation, so no notch incision
      return;

   // We have some notch incision
   double dWeight;
   if (dWaveElev > dNotchElev)
      // Should not happen: safety check
      dWeight = 1;
   else
      // Assume a linear decrease in incision with distance downwards from notch apex
      dWeight = 1 - ((dNotchElev - dWaveElev) / CLIFF_NOTCH_CUTOFF_DISTANCE);

   // Calculate this-timestep cliff notch erosion (is a length in external CRS units). Note that only consolidated sediment can have a cliff notch
   double const dNotchIncision = dWeight * dWaveEnergy / m_dCliffErosionResistance;

   // Deepen the cliff object's erosional notch as a result of wave energy during this timestep. Note that notch deepening may be constrained, since this-timestep notch extension cannot exceed the length (i.e. cellside minus notch depth) of sediment remaining on the cell
   pCliff->IncreaseNotchIncision(dNotchIncision);

   // LogStream << m_ulIter << ": dRunup = " << dRunup << " dWaveElev = " << dWaveElev << " dWeight = " << dWeight << " dNotchIncision = " << dNotchIncision << endl;
}

//===============================================================================================================================
//! Deposit the unconsolidated sediment from cliff collapse as talus on the cell on which collapse occurred.
//===============================================================================================================================
int CSimulation::nDoCliffCollapseTalusDeposition(int const nCoast, CRWCliff const* pCliff, double const dSandFromCollapse, double const dCoarseFromCollapse, int const nNotchLayer)
{
   // Check: is there some sand- or coarse-sized sediment to deposit?
   if ((dSandFromCollapse + dCoarseFromCollapse) < SEDIMENT_ELEV_TOLERANCE)
      return RTN_OK;

   LogStream << m_ulIter << "; coast " << nCoast << " in nDoCliffCollapseTalusDeposition() dSandFromCollapse = " << dSandFromCollapse << " dCoarseFromCollapse = " << dCoarseFromCollapse << endl;

   // Get the cliff cell's grid coords
   int const nX = pCliff->pPtiGetCellMarkedAsCliff()->nGetX();
   int const nY = pCliff->pPtiGetCellMarkedAsCliff()->nGetY();

   // Get a pointer to the layer in which the notch was incised
   CRWCellLayer* pLayer = m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nNotchLayer);

   // And get a pointer to the cell layer's talus object
   CRWCellTalus* pTalus = pLayer->pGetOrCreateTalus();

   // Add the sediment from the collapse to the talus object for this layer
   pTalus->AddSandDepth(dSandFromCollapse);
   pTalus->AddCoarseDepth(dCoarseFromCollapse);

   return RTN_OK;
}
