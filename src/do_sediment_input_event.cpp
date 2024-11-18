/*!
 *
 * \file do_sediment_input_event.cpp
 * \brief Deposits sediment onto the grid
 * \details TODO 001 A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo
 *
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
#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
using std::begin;
using std::end;
using std::find;

#include "cme.h"
#include "cell.h"
#include "coast.h"
#include "sediment_input_event.h"

//===============================================================================================================================
//! Check to see if we have any sediment input events this timestep, if so then do the event(s)
//===============================================================================================================================
int CSimulation::nCheckForSedimentInputEvent(void)
{
   // Go through all sediment input events, check for any this timestep
   int nEvents = static_cast<int>(m_pVSedInputEvent.size());

   for (int n = 0; n < nEvents; n++)
   {
      if (m_pVSedInputEvent[n]->ulGetEventTimeStep() == m_ulIter)
      {
         m_bSedimentInputThisIter = true;

         int nRet = nDoSedimentInputEvent(n);
         if (nRet != RTN_OK)
            return nRet;
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Do a sediment input event
//===============================================================================================================================
int CSimulation::nDoSedimentInputEvent(int const nEvent)
{
   // Get values for the sediment input event
   int nLocID = m_pVSedInputEvent[nEvent]->nGetLocationID();
   double
       dFineSedVol = m_pVSedInputEvent[nEvent]->dGetFineSedVol(),
       dSandSedVol = m_pVSedInputEvent[nEvent]->dGetSandSedVol(),
       dCoarseSedVol = m_pVSedInputEvent[nEvent]->dGetCoarseSedVol(),
       dLen = m_pVSedInputEvent[nEvent]->dGetLen(),
       dWidth = m_pVSedInputEvent[nEvent]->dGetWidth();
      //  dThick = m_pVSedInputEvent[nEvent]->dGetThick();

   if (m_bSedimentInputAtPoint || m_bSedimentInputAtCoast)
   {
      // The sediment input event is at a user-specified location, or in a block at the nearest point on a coast to a user-specified location. So get the location from values read from the shapefile
      int
          nEvents = static_cast<int>(m_VnSedimentInputLocationID.size()),
          nPointGridX = -1,
          nPointGridY = -1;

      for (int n = 0; n < nEvents; n++)
      {
         if (m_VnSedimentInputLocationID[n] == nLocID)
         {
            nPointGridX = nRound(m_VdSedimentInputLocationX[n]);
            nPointGridY = nRound(m_VdSedimentInputLocationY[n]);
         }
      }

      if (nPointGridX == -1)
         // Should never get here
         return RTN_ERR_SEDIMENT_INPUT_EVENT;

      // Is this sediment input event at a pre-specified point, or at a block on a coast, or along a line intersecting with a coast?
      if (m_bSedimentInputAtPoint)
      {
         // Sediment input is at a pre-specified point
         LogStream << "Sediment input event " << nEvent << " at pre-specified point [" << nPointGridX << "][" << nPointGridY << "] = {" << dGridXToExtCRSX(nPointGridX) << ", " << dGridYToExtCRSY(nPointGridY) << "] with Location ID " << nLocID;

         int nTopLayer = m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].nGetTopLayerAboveBasement();

         // Add to this cell's unconsolidated sediment
         double dFineDepth = dFineSedVol / m_dCellArea;
         m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddFineDepth(dFineDepth);
         m_dThisiterUnconsFineInput += dFineDepth;

         double dSandDepth = dSandSedVol / m_dCellArea;
         m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSandDepth(dSandDepth);
         m_dThisiterUnconsSandInput += dSandDepth;

         double dCoarseDepth = dCoarseSedVol / m_dCellArea;
         m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarseDepth(dCoarseDepth);
         m_dThisiterUnconsCoarseInput += dCoarseDepth;

         // And update the cell's total
         m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddToTotSedimentInputDepth(dFineDepth + dSandDepth + dCoarseDepth);

         LogStream << ", depth of fine sediment added = " << dFineDepth << " m, depth of sand sediment added = " << dSandDepth << " m, depth of coarse sediment added = " << dCoarseDepth << " m" << endl;
      }
      else if (m_bSedimentInputAtCoast)
      {
         // Is in a sediment block, seaward from a coast
         LogStream << "Sediment input event " << nEvent << " with Location ID " << nLocID << " at closest point on coast to [" << nPointGridX << "][" << nPointGridY << "] = {" << dGridXToExtCRSX(nPointGridX) << ", " << dGridYToExtCRSY(nPointGridY) << "]" << endl;

         // Find the closest point on the coastline
         CGeom2DIPoint PtiCoastPoint = PtiFindClosestCoastPoint(nPointGridX, nPointGridY);

         int
             nCoastX = PtiCoastPoint.nGetX(),
             nCoastY = PtiCoastPoint.nGetY();

         LogStream << "Closest coast point is at [" << nCoastX << "][" << nCoastY << "] = {" << dGridXToExtCRSX(nCoastX) << ", " << dGridYToExtCRSY(nCoastY) << "}, along-coast width of sediment block = " << dWidth << " m, coast-normal length of sediment block = " << dLen << " m" << endl;

         int
             nCoast = m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLandform()->nGetCoast(),
             nCoastPoint = m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLandform()->nGetPointOnCoast(),
             nHalfWidth = nRound(dWidth / m_dCellSide),
             nLength = nRound(dLen / m_dCellSide),
             nCoastLen = m_VCoast[nCoast].nGetCoastlineSize(),
             nCoastXBefore = nCoastX,
             nCoastYBefore = nCoastY,
             nCoastXAfter = nCoastX,
             nCoastYAfter = nCoastY;

         if (nCoastPoint > 0)
         {
            nCoastXBefore = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint - 1)->nGetX();
            nCoastYBefore = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint - 1)->nGetY();
         }

         if (nCoastPoint < nCoastLen - 1)
         {
            nCoastXAfter = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint + 1)->nGetX();
            nCoastYAfter = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPoint + 1)->nGetY();
         }

         int
             nCoastHand = m_VCoast[nCoast].nGetSeaHandedness(),
             nPerpHand = LEFT_HANDED;

         if (nCoastHand == LEFT_HANDED)
            nPerpHand = RIGHT_HANDED;
         
         vector<CGeom2DIPoint> VPoints;

         // If size of plume is smaller than a cell, choose the current coast point
         if (nHalfWidth == 0)
         {
            VPoints.push_back(PtiCoastPoint);
         }
         else
         {
            // It is larger than a cell
            vector<int> VnCentrePointsXOffset, VnCentrePointsYOffset;
            
            for (int m = 0; m < nHalfWidth; m++)
            {
               if (m == 0)
               {
                  VPoints.push_back(CGeom2DIPoint(nCoastX, nCoastY));
                  for (int n = 1; n < nLength; n++)
                  {
                     CGeom2DIPoint PtiTmp = PtiGetPerpendicular(nCoastXBefore, nCoastYBefore, nCoastXAfter, nCoastYAfter, n, nPerpHand);

                     if (bIsWithinValidGrid(&PtiTmp))
                     {
                        VPoints.push_back(PtiTmp);

                        VnCentrePointsXOffset.push_back(PtiTmp.nGetX() - nCoastX);
                        VnCentrePointsYOffset.push_back(PtiTmp.nGetY() - nCoastY);
                     }
                  }
               }
               else
               {
                  int
                      nCoastPointInBlockBefore = nCoastPoint - m,
                      nCoastPointInBlockAfter = nCoastPoint + m;

                  if (nCoastPointInBlockBefore >= 0)
                  {
                     int
                         nCoastXInBlockBefore = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPointInBlockBefore)->nGetX(),
                         nCoastYInBlockBefore = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPointInBlockBefore)->nGetY();

                     if (bIsWithinValidGrid(nCoastXInBlockBefore, nCoastYInBlockBefore))
                     {
                        VPoints.push_back(CGeom2DIPoint(nCoastXInBlockBefore, nCoastYInBlockBefore));

                        for (unsigned int n = 0; n < VnCentrePointsXOffset.size(); n++)
                        {
                           int
                               nXTmp = nCoastXInBlockBefore + VnCentrePointsXOffset[n],
                               nYTmp = nCoastYInBlockBefore + VnCentrePointsYOffset[n];

                           if (bIsWithinValidGrid(nXTmp, nYTmp))
                              VPoints.push_back(CGeom2DIPoint(nXTmp, nYTmp));
                        }
                     }
                  }

                  if (nCoastPointInBlockAfter < nCoastLen)
                  {
                     int
                         nCoastXInBlockAfter = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPointInBlockAfter)->nGetX(),
                         nCoastYInBlockAfter = m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nCoastPointInBlockAfter)->nGetY();

                     if (bIsWithinValidGrid(nCoastXInBlockAfter, nCoastYInBlockAfter))
                     {
                        VPoints.push_back(CGeom2DIPoint(nCoastXInBlockAfter, nCoastYInBlockAfter));

                        for (unsigned int n = 0; n < VnCentrePointsXOffset.size(); n++)
                        {
                           int
                               nXTmp = nCoastXInBlockAfter + VnCentrePointsXOffset[n],
                               nYTmp = nCoastYInBlockAfter + VnCentrePointsYOffset[n];

                           if (bIsWithinValidGrid(nXTmp, nYTmp))
                           {
                              CGeom2DIPoint PtiTmp(nXTmp, nYTmp);

                              // Is this cell already in the array?
                              if (find(begin(VPoints), end(VPoints), PtiTmp) == end(VPoints))
                                 // No it isn't
                                 VPoints.push_back(PtiTmp);
                           }
                        }
                     }
                  }
               }
            }

            // // DEBUG CODE ===============================================
            //       LogStream << endl;
            //       unsigned int m = 0;
            //       for (unsigned int n = 0; n < VPoints.size(); n++)
            //       {
            //          LogStream << "[" << VPoints[n].nGetX() << ", " << VPoints[n].nGetY() << "] ";
            //          m++;
            //          if (m == VnCentrePointsXOffset.size() + 1)
            //          {
            //             LogStream << endl;
            //             m = 0;
            //          }
            //       }
            //       LogStream << endl;
            // // DEBUG CODE ===============================================
         }

         // OK we now know which cells are part of the sediment block, and so will receive sediment input. Next calculate the volume per cell
         double
             dFineDepth = dFineSedVol / m_dCellArea,
             dSandDepth = dSandSedVol / m_dCellArea,
             dCoarseDepth = dCoarseSedVol / m_dCellArea;

         LogStream << "Total depth of fine sediment added = " << dFineDepth << " m, total depth of sand sediment added = " << dSandDepth << " m, total depth of coarse sediment added = " << dCoarseDepth << " m" << endl;

         size_t nArea = VPoints.size();
         double
             dArea = static_cast<double>(nArea),
             dFineDepthPerCell = dFineDepth / dArea,
             dSandDepthPerCell = dSandDepth / dArea,
             dCoarseDepthPerCell = dCoarseDepth / dArea;

         // OK, so finally: put some sediment onto each cell in the sediment block
         int nTopLayer = m_pRasterGrid->m_Cell[nPointGridX][nPointGridY].nGetTopLayerAboveBasement();
         for (unsigned int n = 0; n < nArea; n++)
         {
            // Add to this cell's unconsolidated sediment
            m_pRasterGrid->m_Cell[VPoints[n].nGetX()][VPoints[n].nGetY()].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddFineDepth(dFineDepthPerCell);
            m_dThisiterUnconsFineInput += dFineDepth;

            m_pRasterGrid->m_Cell[VPoints[n].nGetX()][VPoints[n].nGetY()].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSandDepth(dSandDepthPerCell);
            m_dThisiterUnconsSandInput += dSandDepth;

            m_pRasterGrid->m_Cell[VPoints[n].nGetX()][VPoints[n].nGetY()].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarseDepth(dCoarseDepthPerCell);
            m_dThisiterUnconsCoarseInput += dCoarseDepth;

            // And update the cell's total
            m_pRasterGrid->m_Cell[VPoints[n].nGetX()][VPoints[n].nGetY()].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddToTotSedimentInputDepth(dFineDepth + dSandDepth + dCoarseDepth);
         }
      }
   }
   else if (m_bSedimentInputAlongLine)
   {
      // The sediment input event is where a line intersects a coast. So get the line from values read from the shapefile
      int nPoints = static_cast<int>(m_VnSedimentInputLocationID.size());
      vector<int> VnLineGridX, VnLineGridY;

      for (int n = 0; n < nPoints; n++)
      {
         if (m_VnSedimentInputLocationID[n] == nLocID)
         {
            VnLineGridX.push_back(nRound(m_VdSedimentInputLocationX[n]));
            VnLineGridY.push_back(nRound(m_VdSedimentInputLocationY[n]));
         }
      }

      // Should never get here
      if (VnLineGridX.size() == 0)
         return RTN_ERR_SEDIMENT_INPUT_EVENT;

      LogStream << "Sediment input event " << nEvent << " at line/coast intersection for line with ID " << nLocID;

      // Now go along the line, joining each pair of points by a straight line
      int
          nCoastX = -1,
          nCoastY = -1;

      for (unsigned int n = 0; n < VnLineGridX.size() - 1; n++)
      {
         int
             nXStart = VnLineGridX[n],
             nXEnd = VnLineGridX[n + 1],
             nYStart = VnLineGridY[n],
             nYEnd = VnLineGridY[n + 1];

         // Interpolate between pairs of points using a simple DDA line algorithm, see http://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm) Note that Bresenham's algorithm gave occasional gaps
         double
             dXInc = nXEnd - nXStart,
             dYInc = nYEnd - nYStart,
             dLength = tMax(tAbs(dXInc), tAbs(dYInc));

         dXInc /= dLength;
         dYInc /= dLength;

         double
             dX = nXStart,
             dY = nYStart;

         // Process each interpolated point
         for (int m = 0; m <= nRound(dLength); m++)
         {
            int
                nX = static_cast<int>(dX),
                nY = static_cast<int>(dY);

            // Have we hit a coastline cell?
            if (bIsWithinValidGrid(nX, nY) && m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               nCoastX = nX;
               nCoastY = nY;
               break;
            }

            // Two diagonal(ish) raster lines can cross each other without any intersection, so must also test an adjacent cell for intersection (does not matter which adjacent cell)
            if (bIsWithinValidGrid(nX, nY + 1) && m_pRasterGrid->m_Cell[nX][nY + 1].bIsCoastline())
            {
               nCoastX = nX;
               nCoastY = nY + 1;
               break;
            }

            // Increment for next time
            dX += dXInc;
            dY += dYInc;
         }

         // Did we find an intersection?
         if (nCoastX != -1)
            break;
      }

      // Did we find an intersection?
      if (nCoastX == -1)
      {
         // Nope
         LogStream << endl
                   << "No intersection found" << endl;
         return RTN_ERR_SEDIMENT_INPUT_EVENT;
      }

      // OK we have an intersection of the line and coast
      LogStream << ", intersection is at [" << nCoastX << "][" << nCoastY << "] = {" << dGridXToExtCRSX(nCoastX) << ", " << dGridYToExtCRSY(nCoastY) << "}" << endl;

      int nTopLayer = m_pRasterGrid->m_Cell[nCoastX][nCoastY].nGetTopLayerAboveBasement();

      // Add to this cell's unconsolidated sediment
      double dFineDepth = dFineSedVol / m_dCellArea;
      m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddFineDepth(dFineDepth);
      m_dThisiterUnconsFineInput += dFineDepth;

      double dSandDepth = dSandSedVol / m_dCellArea;
      m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddSandDepth(dSandDepth);
      m_dThisiterUnconsSandInput += dSandDepth;

      double dCoarseDepth = dCoarseSedVol / m_dCellArea;
      m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddCoarseDepth(dCoarseDepth);
      m_dThisiterUnconsCoarseInput += dCoarseDepth;

      // And update the cell's total
      m_pRasterGrid->m_Cell[nCoastX][nCoastY].pGetLayerAboveBasement(nTopLayer)->pGetUnconsolidatedSediment()->AddToTotSedimentInputDepth(dFineDepth + dSandDepth + dCoarseDepth);

      LogStream << "Depth of fine sediment added = " << dFineDepth << " m, depth of sand sediment added = " << dSandDepth << " m, depth of coarse sediment added = " << dCoarseDepth << " m" << endl;
   }

   return RTN_OK;
}
