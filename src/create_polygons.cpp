/*!
 *
 * \file create_polygons.cpp
 * \brief Creates coast polygons for sediment transport calcs
 * \details TODO 001 A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2025
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
using std::endl;

#include <string>
using std::to_string;

#include <stack>
using std::stack;

#include "cme.h"
#include "simulation.h"
#include "coast.h"

//===============================================================================================================================
//! Create polygons, and marks the polygon boundaries on the raster grid
//===============================================================================================================================
int CSimulation::nCreateAllPolygons(void)
{
   // Global polygon count
   m_nGlobalPolygonID = -1;

   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // Do this for every point on the coastline
      int nPolygon = -1;          // This-coast-only polygon ID
      int nNode = -1;
      int nPrevProfileCoastPoint = -1;
      int nPrevProfile = -1;

      for (int nCoastPoint = 0; nCoastPoint < m_VCoast[nCoast].nGetCoastlineSize(); nCoastPoint++)
      {
         int nThisProfile = m_VCoast[nCoast].nGetProfileNumber(nCoastPoint);
         CGeomProfile* pThisProfile = m_VCoast[nCoast].pGetProfile(nThisProfile);

         if ((nThisProfile != INT_NODATA) && (pThisProfile->bOKIncStartAndEndOfCoast()))
         {
            // There is a valid profile at this coast point
            if (nPrevProfileCoastPoint >= 0)
            {
               // Calculate half the along-coast distance (in coast points) between this profile and the previous (i.e. up-coast) profile
               int nDist = (nCoastPoint - nPrevProfileCoastPoint) / 2;

               // OK, set the node point in the coast object. We do this now, instead of earlier on, since some profiles (i.e. polygon boundaries) may have been marked as invalid
               int nNodePoint = nCoastPoint - nDist;
               m_VCoast[nCoast].SetPolygonNode(nNodePoint, ++nNode);

               // Get the grid CRS co-ordinates of the coast node
               CGeom2DIPoint PtiNode = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(nNodePoint);

               // Get the previous profile, also set some defaults (assuming for now that this polygon is not approximately triangular i.e. both normals do not meet)
               CGeomProfile* pPrevProfile = m_VCoast[nCoast].pGetProfile(nPrevProfile);
               int nPrevProfileEnd = pPrevProfile->nGetProfileSize()-1;
               int nThisProfileEnd = pThisProfile->nGetProfileSize()-1;
               bool bMeetsAtAPoint = false;
               CGeom2DPoint PtCoastwardTip;

//                // DEBUG CODE =============================
//                CGeom2DPoint PtPrevEndTmp = *pPrevProfile->pPtGetPointInProfile(nPrevProfileEnd);
//                double
//                   dXTmp = dExtCRSXToGridX(PtPrevEndTmp.dGetX()),
//                   dYTmp = dExtCRSYToGridY(PtPrevEndTmp.dGetY());
//                assert(bIsWithinValidGrid(dXTmp, dYTmp));
//
//                CGeom2DPoint PtThisEndTmp = *pThisProfile->pPtGetPointInProfile(nThisProfileEnd);
//                dXTmp = dExtCRSXToGridX(PtThisEndTmp.dGetX());
//                dYTmp = dExtCRSYToGridY(PtThisEndTmp.dGetY());
//                assert(bIsWithinValidGrid(dXTmp, dYTmp));
//                // DEBUG CODE =============================

               // Now check to see if the two normals do meet i.e. if they are coincident
               if (pThisProfile->bFindProfileInCoincidentProfiles(nPrevProfile))
               {
                  // Yes they do meet
                  bMeetsAtAPoint = true;

                  // Find the most coastward point at which this normal and the previous normal touch. If they do not touch, the polygon requires a 'joining line'
                  pThisProfile->GetMostCoastwardSharedLineSegment(nPrevProfile, nThisProfileEnd, nPrevProfileEnd);
                  if (nThisProfileEnd == -1)
                  {
                     LogStream << m_ulIter << ": " << ERR << "profile " << nPrevProfile << " should be coincident with profile " << nThisProfile << " but was not found" << endl;
                     return RTN_ERR_BAD_MULTILINE;
                  }

                  PtCoastwardTip = *pThisProfile->pPtGetPointInProfile(nThisProfileEnd);
               }

               // Next, calculate the polygon's boundary (external CRS), do this in an anti-clockwise sequence
               vector<CGeom2DPoint> PtVBoundary;

               // Start appending points: begin at the node point, then move down-coast as far as the down-coast (this) normal
               for (int i = nNodePoint; i <= nCoastPoint; i++)
                  PtVBoundary.push_back(*m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(i));

               // Use the penultimate coastline point as the start point for the point-in-polygon search later, during flood fill
               int nPointInPolygonStartPoint = static_cast<int>(PtVBoundary.size()) - 2;

               // Append the points in the down-coast normal. Omit the last point of this normal if the the most seaward point of the this normal, and the most seaward point of the up-coast (previous) normal are the same
               int nFinishPoint = nThisProfileEnd;

               if (bMeetsAtAPoint)
                  nFinishPoint--;

               for (int i = 0; i <= nFinishPoint; i++)
               {
                  CGeom2DPoint PtThis = *pThisProfile->pPtGetPointInProfile(i);
                  if (! (PtThis == &PtVBoundary.back()))
                     PtVBoundary.push_back(PtThis);
               }

               // Append the points in the up-coast (previous) normal, in reverse order
               for (int i = nPrevProfileEnd; i >= 0; i--)
               {
                  CGeom2DPoint PtThis = *pPrevProfile->pPtGetPointInProfile(i);
                  if (! (PtThis == &PtVBoundary.back()))
                     PtVBoundary.push_back(PtThis);
               }

               // Append the points from the remaining bit of coast, moving down-coast to finish at the node point. Note that we must include the node point here, in order to obtain a closed polygon
               for (int i = nPrevProfileCoastPoint; i <= nNodePoint; i++)
               {
                  CGeom2DPoint PtThis = *m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(i);
                  if (! (PtThis == &PtVBoundary.back()))
                     PtVBoundary.push_back(PtThis);
               }

               // Now identify the 'anti-node', this is the seaward point 'opposite' the polygon's coastal node
               CGeom2DIPoint PtiAntiNode;
               if (bMeetsAtAPoint)
                  PtiAntiNode = PtiExtCRSToGridRound(&PtCoastwardTip);
               else
               {
                  CGeom2DPoint PtAvg = PtAverage(pThisProfile->pPtGetPointInProfile(nThisProfileEnd), pPrevProfile->pPtGetPointInProfile(nPrevProfileEnd));
                  PtiAntiNode = PtiExtCRSToGridRound(&PtAvg);
               }

               // Safety check

               // Create the coast's polygon object
               m_VCoast[nCoast].CreatePolygon(++m_nGlobalPolygonID, ++nPolygon, nNodePoint, &PtiNode, &PtiAntiNode, nPrevProfile, nThisProfile, &PtVBoundary, nPrevProfileEnd+1, nThisProfileEnd+1, nPointInPolygonStartPoint);

               // Get a pointer to this polygon object
               CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPolygon);

               // And store this pointer for simulation-wide access
               m_pVCoastPolygon.push_back(pPolygon);

               // Now rasterize the polygon boundaries: first, the coastline. This is necessary so that sand/coarse sediment derived from platform erosion of the coast cells is correctly added to the containing polygon's unconsolidated sediment
               for (int i = nPrevProfileCoastPoint; i <= nCoastPoint; i++)
               {
                  CGeom2DIPoint PtiToMark = *m_VCoast[nCoast].pPtiGetCellMarkedAsCoastline(i);
                  m_pRasterGrid->m_Cell[PtiToMark.nGetX()][PtiToMark.nGetY()].SetPolygonID(m_nGlobalPolygonID);
               }

               // Do the upcoast normal profile (does the whole length, including any shared line segments. So some cells are marked twice, however this is not a problem)
               int nCellsInProfile = pPrevProfile->nGetNumCellsInProfile();
               for (int i = 0; i < nCellsInProfile; i++)
               {
                  CGeom2DIPoint PtiToMark = *pPrevProfile->pPtiGetCellInProfile(i);
                  m_pRasterGrid->m_Cell[PtiToMark.nGetX()][PtiToMark.nGetY()].SetPolygonID(m_nGlobalPolygonID);
               }

               // Do the downcoast normal profile (again does the whole length, including any shared line segments)
               nCellsInProfile = pThisProfile->nGetNumCellsInProfile();
               for (int i = 0; i < nCellsInProfile; i++)
               {
                  CGeom2DIPoint PtiToMark = *pThisProfile->pPtiGetCellInProfile(i);
                  m_pRasterGrid->m_Cell[PtiToMark.nGetX()][PtiToMark.nGetY()].SetPolygonID(m_nGlobalPolygonID);
               }

               // If the polygon doesn't meet at a point at its seaward end, also need to rasterize the 'joining line'
               if (! bMeetsAtAPoint)
               {
                  CGeom2DPoint
                     PtUpCoastNormalEnd = *pPrevProfile->pPtGetPointInProfile(nPrevProfileEnd),
                     PtDownCoastNormalEnd = *pThisProfile->pPtGetPointInProfile(nThisProfileEnd);

                  RasterizePolygonJoiningLine(&PtUpCoastNormalEnd, &PtDownCoastNormalEnd);

//                   pPolygon->SetNotPointed();
               }
            }

            nPrevProfileCoastPoint = nCoastPoint;
            nPrevProfile = nThisProfile;
         }
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Puts a polygon 'joining line' (the line which is the seaward boundary of the polygon, if the polygon doesn't meet at a point) onto the raster grid
//===============================================================================================================================
void CSimulation::RasterizePolygonJoiningLine(CGeom2DPoint const* pPt1, CGeom2DPoint const* pPt2)
{
   // The start point of the line, must convert from the external CRS to grid CRS
   double dXStart = dExtCRSXToGridX(pPt1->dGetX());
   double dYStart = dExtCRSYToGridY(pPt1->dGetY());

   // The end point of the line, again convert from the external CRS to grid CRS
   double dXEnd = dExtCRSXToGridX(pPt2->dGetX());
   double dYEnd = dExtCRSYToGridY(pPt2->dGetY());

   // Safety check, in case the two points are identical (can happen due to rounding errors)
   if ((bFPIsEqual(dXStart, dXEnd, TOLERANCE)) && (bFPIsEqual(dYStart, dYEnd, TOLERANCE)))
      return;

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
      int nX = nRound(dX);
      int nY = nRound(dY);

      // Safety check
      if (! bIsWithinValidGrid(nX, nY))
         KeepWithinValidGrid(nRound(dXStart), nRound(dYStart), nX, nY);

      // Mark this point on the raster grid
      m_pRasterGrid->m_Cell[nX][nY].SetPolygonID(m_nGlobalPolygonID);

      // And increment for next time
      dX += dXInc;
      dY += dYInc;
   }
}

//===============================================================================================================================
//! Marks cells of the raster grid that are within each coastal polygon. The flood fill code used here is adapted from an example by Lode Vandevenne (http://lodev.org/cgtutor/floodfill.html#Scanline_Floodfill_Algorithm_With_Stack) 
//===============================================================================================================================
void CSimulation::MarkPolygonCells(void)
{
   // Do this for each coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      // Do this for every coastal polygon
      for (int nPoly = 0; nPoly < m_VCoast[nCoast].nGetNumPolygons(); nPoly++)
      {
         int nCellsInPolygon = 0;
         double dTotDepth = 0;
         double dStoredUnconsFine = 0;
         double dStoredUnconsSand = 0;
         double dStoredUnconsCoarse = 0;
         double dStoredConsFine = 0;
         double dStoredConsSand = 0;
         double dStoredConsCoarse = 0;
         double dSedimentInputFine = 0;
         double dSedimentInputSand = 0;
         double dSedimentInputCoarse = 0;

         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);
         int nPolyID = pPolygon->nGetGlobalID();

         // Create an empty stack
         stack<CGeom2DIPoint> PtiStack;

         // Since the polygon's vector boundary does not coincide exactly with the polygon's raster boundary, and the point-in-polygon check gives an indeterminate result if the point is exactly on the polygon's boundary, for safety we must construct a vector 'inner buffer' which is smaller than, and inside, the vector boundary
         int nHand = m_VCoast[nCoast].nGetSeaHandedness();
         int nSize = pPolygon->nGetBoundarySize();
         vector<CGeom2DPoint> PtVInnerBuffer;

         for (int i = 0; i < nSize-1; i++)
         {
            int j = i + 1;
            if (i == nSize-2)       // We must ignore the duplicated node point
               j = 0;
            
            CGeom2DPoint PtThis = *pPolygon->pPtGetBoundaryPoint(i);
            CGeom2DPoint PtNext = *pPolygon->pPtGetBoundaryPoint(j);
               
            // Safety check
            if (PtThis == PtNext)
               continue;
               
            CGeom2DPoint PtBuffer = PtGetPerpendicular(&PtThis, &PtNext, m_dCellSide, nHand);

            PtVInnerBuffer.push_back(PtBuffer);
         }

//          // DEBUG STUFF
//          LogStream << endl << m_ulIter << ":: coast " << nCoast << ", polygon " << nPoly << endl;
//          LogStream << "Boundary\t\t\tBuffer" << endl;
//          for (int i = 0; i < pPolygon->nGetBoundarySize()-1; i++)
//             LogStream << "{" << pPolygon->pPtGetBoundaryPoint(i)->dGetX() << ", " << pPolygon->pPtGetBoundaryPoint(i)->dGetY() << "}\t{" << PtVInnerBuffer[i].dGetX() << ", " << PtVInnerBuffer[i].dGetY() << "}" << endl;
//          LogStream << endl;

         // Now look for a point that is within the polygon as a start point for the flood fill. For a first attempt, calculate the polygon's centroid
         CGeom2DPoint PtStart = pPolygon->PtGetCentroid();

         // Is the centroid within the inner buffer?
         if (! bIsWithinPolygon(&PtStart, &PtVInnerBuffer))
         {
            // No, it is not: the polygon must be a concave polygon. So keep looking for a point which is definitely inside the polygon, using an alternative method
            PtStart = PtFindPointInPolygon(&PtVInnerBuffer, pPolygon->nGetPointInPolygonSearchStartPoint());
         }

         // Safety check (PtFindPointInPolygon() returns CGeom2DPoint(DBL_NODATA, DBL_NODATA) if it cannot find a valid start point)
         if (bFPIsEqual(PtStart.dGetX(), DBL_NODATA, TOLERANCE))
         {
            LogStream << m_ulIter << ": " << ERR << "could not find a flood fill start point for coast " << nCoast << ", polygon " << nPoly << endl;
            break;
         }

         // We have a flood fill start point which is definitely within the polygon so push this point onto the stack
         CGeom2DIPoint PtiStart = PtiExtCRSToGridRound(&PtStart);                // Grid CRS
         PtiStack.push(PtiStart);

//          LogStream << m_ulIter << ": filling polygon " << nPoly << " from [" << PtiStart.nGetX() << "][" << PtiStart.nGetY() << "] = {" << dGridCentroidXToExtCRSX(PtiStart.nGetX()) << ", " << dGridCentroidYToExtCRSY(PtiStart.nGetY()) << "}" << endl;

         // Then do the flood fill: loop until there are no more cell co-ordinates on the stack
         while (! PtiStack.empty())
         {
            CGeom2DIPoint Pti = PtiStack.top();
            PtiStack.pop();

            int nX = Pti.nGetX();
            int nY = Pti.nGetY();
               
            while ((nX >= 0) && (m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID() == INT_NODATA))
               nX--;

            nX++;
            
            bool bSpanAbove = false;
            bool bSpanBelow = false;

            while ((nX < m_nXGridSize) && (m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID() == INT_NODATA))
            {
               // Mark the cell as being in this polygon
               m_pRasterGrid->m_Cell[nX][nY].SetPolygonID(nPolyID);
//                LogStream << "[" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
               
               // Increment the running totals for this polygon
               // dCliffCollapseErosionFine += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionFine();
               // dCliffCollapseErosionSand += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionSand();
               // dCliffCollapseErosionCoarse += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseErosionCoarse();
               // dCliffCollapseTalusSand += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseSandTalusDeposition();
               // dCliffCollapseTalusCoarse += m_pRasterGrid->m_Cell[nX][nY].dGetThisIterCliffCollapseCoarseTalusDeposition();
               
               // Get the number of the highest layer with non-zero thickness
               int nThisLayer = m_pRasterGrid->m_Cell[nX][nY].nGetTopNonZeroLayerAboveBasement();
               
               // And increment some more running totals for this polygon TODO 066 should this be for ALL layers above the basement?
               dStoredUnconsFine += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetFineDepth();
               dStoredUnconsSand += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetSandDepth();
               dStoredUnconsCoarse += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetCoarseDepth();

               dStoredConsFine += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetFineDepth();
               dStoredConsSand += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetSandDepth();
               dStoredConsCoarse += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetConsolidatedSediment()->dGetCoarseDepth();

               // Add to the start-iteration total of suspended fine sediment within polygons
               m_dStartIterSuspFineInPolygons += m_pRasterGrid->m_Cell[nX][nY].dGetSuspendedSediment();

               // Add to the total of sediment derived from sediment input events
               dSedimentInputFine += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetFineSedimentInputDepth();
               dSedimentInputSand += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetSandSedimentInputDepth();
               dSedimentInputCoarse += m_pRasterGrid->m_Cell[nX][nY].pGetLayerAboveBasement(nThisLayer)->pGetUnconsolidatedSediment()->dGetCoarseSedimentInputDepth();

               nCellsInPolygon++;
               dTotDepth += m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth();

               if ((! bSpanAbove) && (nY > 0) && (m_pRasterGrid->m_Cell[nX][nY-1].nGetPolygonID() == INT_NODATA))
               {
                  PtiStack.push(CGeom2DIPoint(nX, nY-1));
                  bSpanAbove = true;
               }
               else if (bSpanAbove && (nY > 0) && (m_pRasterGrid->m_Cell[nX][nY-1].nGetPolygonID() != INT_NODATA))
               {
                  bSpanAbove = false;
               }

               if ((! bSpanBelow) && (nY < m_nYGridSize-1) && (m_pRasterGrid->m_Cell[nX][nY+1].nGetPolygonID() == INT_NODATA))
               {
                  PtiStack.push(CGeom2DIPoint(nX, nY+1));
                  bSpanBelow = true;
               }
               else if (bSpanBelow && (nY < m_nYGridSize-1) && (m_pRasterGrid->m_Cell[nX][nY+1].nGetPolygonID() != INT_NODATA))
               {
                  bSpanBelow = false;
               }

               nX++;
            }
         }
         
         // Store this polygon's stored unconsolidated sediment depths
         pPolygon->SetPreExistingUnconsFine(dStoredUnconsFine);
         pPolygon->SetPreExistingUnconsSand(dStoredUnconsSand);
         pPolygon->SetPreExistingUnconsCoarse(dStoredUnconsCoarse);

         // Store this polygon's values for unconsolidated sediment derived from sediment input event(s)
         pPolygon->SetSedimentInputUnconsFine(dSedimentInputFine);
         pPolygon->SetSedimentInputUnconsSand(dSedimentInputSand);
         pPolygon->SetSedimentInputUnconsCoarse(dSedimentInputCoarse);

         // Store this polygon's stored consolidated sediment depths
         pPolygon->SetPreExistingConsFine(dStoredConsFine);
         pPolygon->SetPreExistingConsSand(dStoredConsSand);
         pPolygon->SetPreExistingConsCoarse(dStoredConsCoarse);

         // Store the number of cells in the interior of the polygon
         pPolygon->SetNumCellsInPolygon(nCellsInPolygon);
         // LogStream << m_ulIter << ": N cells = " << nCellsInPolygon << " in polygon " << nPoly << endl;

         // Calculate the total volume of seawater on the polygon (m3) and store it
         double dSeaVolume = dTotDepth * m_dCellSide;
         pPolygon->SetSeawaterVolume(dSeaVolume);
      }
   }

//    // DEBUG CODE ===========================================
//    string strOutFile = m_strOutPath + "polygon_test_";
//    strOutFile += to_string(m_ulIter);
//    strOutFile += ".tif";
// 
//    GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
//    GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
//    pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
//    pDataSet->SetGeoTransform(m_dGeoTransform);
//    double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
//    int
//       n = 0,
//       nInPoly = 0,
//       nNotInPoly = 0;
//       
//    for (int nY = 0; nY < m_nYGridSize; nY++)
//    {
//       for (int nX = 0; nX < m_nXGridSize; nX++)
//       {
//          int nID = m_pRasterGrid->m_Cell[nX][nY].nGetPolygonID();
//          if (nID == INT_NODATA)
//             nNotInPoly++;
//          else
//             nInPoly++;
//          
//          // if (nX == 12)
//          // {
//          //    LogStream << m_ulIter << " " << "[" << nX << "][" << nY << "] is ";
//          //    if (nID == INT_NODATA)
//          //       LogStream << "NOT in poly" << endl;
//          //    else
//          //       LogStream << "in poly " << nID << endl;
//          // }       
//          // 
//          // if (nX == 72)
//          // {
//          //    LogStream << m_ulIter << " " << "[" << nX << "][" << nY << "] is ";
//          //    if (nID == INT_NODATA)
//          //       LogStream << "NOT in poly" << endl;
//          //    else
//          //       LogStream << "in poly " << nID << endl;
//          // }       
// 
//          pdRaster[n++] = nID;
//       }      
//    }
// 
//    GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
//    pBand->SetNoDataValue(m_dMissingValue);
//    int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
//    if (nRet == CE_Failure)
//       return;
// 
//    GDALClose(pDataSet);
//    delete[] pdRaster;
// 
//    // LogStream << m_ulIter << " Number of cells in a polygon = " << nInPoly << endl;
//    // LogStream << m_ulIter << " Number of cells not in any polygon = " << nNotInPoly << endl;
// 
//    // DEBUG CODE ===========================================
}


//===============================================================================================================================
//! For between-polygon potential sediment routing: find which are the adjacent polygons, and calc the length of the shared normal between this polygon and the adjacent polygons TODO 012 Will need to change this when length of coastline-normal profiles (and so polygon seaward length) is determined by depth of closure
//===============================================================================================================================
int CSimulation::nDoPolygonSharedBoundaries(void)
{
   // Do this for every coast
   for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
   {
      int nNumPolygons = m_VCoast[nCoast].nGetNumPolygons();

      // Do this for every coastal polygon
      for (int nPoly = 0; nPoly < nNumPolygons; nPoly++)
      {
         CGeomCoastPolygon* pPolygon = m_VCoast[nCoast].pGetPolygon(nPoly);

         vector<int>
            nVUpCoastAdjacentPolygon,
            nVDownCoastAdjacentPolygon;

         vector<double>
            dVUpCoastBoundaryShare,
            dVDownCoastBoundaryShare;

         // First deal with up-coast adjacent polygons
         if (nPoly == 0)
         {
            // We are at the start of the coastline, no other polygon is adjacent to the up-coast profile of the start-of-coast polygon
            nVUpCoastAdjacentPolygon.push_back(INT_NODATA);
            dVUpCoastBoundaryShare.push_back(1);

            // Store in the polygon
            pPolygon->SetUpCoastAdjacentPolygons(&nVUpCoastAdjacentPolygon);
            pPolygon->SetUpCoastAdjacentPolygonBoundaryShares(&dVUpCoastBoundaryShare);
         }
         else
         {
            // We are not at the start of the coastline, so there is at least one other polygon adjacent to the up-coast profile of this polygon
            int
               nProfile = pPolygon->nGetUpCoastProfile(),
               nPointsInProfile = pPolygon->nGetUpCoastProfileNumPointsUsed();

            double dUpCoastTotBoundaryLen = 0;

            CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

            for (int nPoint = 0; nPoint < nPointsInProfile-1; nPoint++)
            {
               CGeom2DPoint
                  PtStart = *pProfile->pPtGetPointInProfile(nPoint),
                  PtEnd = *pProfile->pPtGetPointInProfile(nPoint+1);

               // Calculate the length of this segment of the normal profile. Note that it should not be zero, since we checked for duplicate points when creating profiles
               double dDistBetween = dGetDistanceBetween(&PtStart, &PtEnd);

               // Find out which polygon is adjacent to each line segment of the polygon's up-coast profile boundary. The basic approach used is to count the number of coincident profiles in each line segment, and (because we are going up-coast) subtract this number from 'this' polygon's number. However, some of these coincident profiles may be invalid, so we must count only the valid co-incident profiles
               int
                  nNumCoinc = pProfile->nGetNumCoincidentProfilesInLineSegment(nPoint),
                  nNumValidCoinc = 0;

               for (int nCoinc = 0; nCoinc < nNumCoinc; nCoinc++)
               {
                  int nProf = pProfile->nGetCoincidentProfileForLineSegment(nPoint, nCoinc);

                  // Safety check
                  if (nProf == -1)
                     continue;

                  CGeomProfile const* pProf = m_VCoast[nCoast].pGetProfile(nProf);

                  if (pProf->bProfileOK())
                     nNumValidCoinc++;
               }

               // First stab at calculating the number of the adjacent polygon
               int nAdj = nPoly - nNumValidCoinc;

               // However, if 'this' polygon is close to the start of the coastline, we get polygon numbers below zero i.e. beyond the start of the coastline. If this happens, set the adjacent polygon to 'off-edge'
               if (nAdj < 0)
                  nAdj = INT_NODATA;

               nVUpCoastAdjacentPolygon.push_back(nAdj);

               dUpCoastTotBoundaryLen += dDistBetween;
               dVUpCoastBoundaryShare.push_back(dDistBetween);
            }

            // Calculate the up-coast boundary share
            for (unsigned int n = 0; n < dVUpCoastBoundaryShare.size(); n++)
            {
               dVUpCoastBoundaryShare[n] /= dUpCoastTotBoundaryLen;
               
               // // Safety check
               // dVUpCoastBoundaryShare[n] = tMin(dVUpCoastBoundaryShare[n], 1.0);
            }

            // Store in the polygon
            pPolygon->SetUpCoastAdjacentPolygons(&nVUpCoastAdjacentPolygon);
            pPolygon->SetUpCoastAdjacentPolygonBoundaryShares(&dVUpCoastBoundaryShare);
         }

         // Now deal with down-coast adjacent polygons
         if (nPoly == nNumPolygons-1)
         {
            // We are at the end of the coastline, no other polygon is adjacent to the down-coast profile of the end-of-coast polygon
            nVDownCoastAdjacentPolygon.push_back(INT_NODATA);
            dVDownCoastBoundaryShare.push_back(1);

            // Store in the polygon
            pPolygon->SetDownCoastAdjacentPolygons(&nVDownCoastAdjacentPolygon);
            pPolygon->SetDownCoastAdjacentPolygonBoundaryShares(&dVDownCoastBoundaryShare);
         }
         else
         {
            // We are not at the end of the coastline, so there is at least one other polygon adjacent to the down-coast profile of this polygon
            int nProfile = pPolygon->nGetDownCoastProfile();
            int nPointsInProfile = pPolygon->nGetDownCoastProfileNumPointsUsed();

            double dDownCoastTotBoundaryLen = 0;

            CGeomProfile* pProfile = m_VCoast[nCoast].pGetProfile(nProfile);

            for (int nPoint = 0; nPoint < nPointsInProfile-1; nPoint++)
            {
               CGeom2DPoint PtStart = *pProfile->pPtGetPointInProfile(nPoint);
               CGeom2DPoint PtEnd = *pProfile->pPtGetPointInProfile(nPoint+1);

               // Calculate the length of this segment of the normal profile. Note that it should not be zero, since we checked for duplicate points when creating profiles
               double dDistBetween = dGetDistanceBetween(&PtStart, &PtEnd);

               // Find out which polygon is adjacent to each line segment of the polygon's down-coast profile boundary. The basic approach used is to count the number of coincident profiles in each line segment, and (because we are going down-coast) add this number to 'this' polygon's number. However, some of these coincident profiles may be invalid, so we must count only the valid co-incident profiles
               int nNumCoinc = pProfile->nGetNumCoincidentProfilesInLineSegment(nPoint);
               int nNumValidCoinc = 0;

               for (int nCoinc = 0; nCoinc < nNumCoinc; nCoinc++)
               {
                  int nProf = pProfile->nGetCoincidentProfileForLineSegment(nPoint, nCoinc);

                  // Safety check
                  if (nProf == -1)
                     continue;

                  CGeomProfile const* pProf = m_VCoast[nCoast].pGetProfile(nProf);

                  if (pProf->bProfileOK())
                     nNumValidCoinc++;
               }

               // First stab at calculating the number of the adjacent polygon
               int nAdj = nPoly + nNumValidCoinc;

               // However, if 'this' polygon is close to the end of the coastline, we get polygon numbers greater than the number of polygons i.e. beyond the end of the coastline. If this happens, set the adjacent polygon to 'off-edge'
               if (nAdj >= nNumPolygons)
                  nAdj = INT_NODATA;

               nVDownCoastAdjacentPolygon.push_back(nAdj);

               dDownCoastTotBoundaryLen += dDistBetween;
               dVDownCoastBoundaryShare.push_back(dDistBetween);
            }

            // Calculate the down-coast boundary share
            for (unsigned int n = 0; n < dVDownCoastBoundaryShare.size(); n++)
            {
               dVDownCoastBoundaryShare[n] /= dDownCoastTotBoundaryLen;
               
               // // Safety check
               // dVDownCoastBoundaryShare[n] = tMin(dVDownCoastBoundaryShare[n], 1.0);
            }

            // Store in the polygon
            pPolygon->SetDownCoastAdjacentPolygons(&nVDownCoastAdjacentPolygon);
            pPolygon->SetDownCoastAdjacentPolygonBoundaryShares(&dVDownCoastBoundaryShare);
         }

         // Finally, calculate the distance between the coast node and the antinode of the polygon
         double dPolygonSeawardLen = dGetDistanceBetween(pPolygon->pPtiGetNode(), pPolygon->pPtiGetAntiNode());

         // And store it
         m_VCoast[nCoast].AppendPolygonLength(dPolygonSeawardLen);

//          // DEBUG CODE ======================================================
//          assert(dVUpCoastBoundaryShare.size() == nVUpCoastAdjacentPolygon.size());
//          assert(dVDownCoastBoundaryShare.size() == nVDownCoastAdjacentPolygon.size());
//
//          //             LogStream << m_ulIter << ": polygon = " << nPoly << (pPolygon->bIsPointed() ? " IS TRIANGULAR" : "") << endl;
//          LogStream << m_ulIter << ": coast " << nCoast << " polygon " << nPoly << endl;
//
//          LogStream << "\tThere are " << nVUpCoastAdjacentPolygon.size() << " UP-COAST adjacent polygon(s) = ";
//          for (unsigned int n = 0; n < nVUpCoastAdjacentPolygon.size(); n++)
//             LogStream << nVUpCoastAdjacentPolygon[n] << " ";
//          LogStream << endl;
//
//          LogStream << "\tThere are " << nVDownCoastAdjacentPolygon.size() << " DOWN-COAST adjacent polygon(s) = ";
//          for (unsigned int n = 0; n < nVDownCoastAdjacentPolygon.size(); n++)
//             LogStream << nVDownCoastAdjacentPolygon[n] << " ";
//          LogStream << endl;
//
//          LogStream << "\tUP-COAST boundary share(s) = ";
//          for (unsigned int n = 0; n < dVUpCoastBoundaryShare.size(); n++)
//             LogStream << dVUpCoastBoundaryShare[n] << " ";
//          LogStream << endl;
// //          LogStream << "\tTotal UP-COAST boundary length = " << dUpCoastTotBoundaryLen << endl;
//
//          LogStream << "\tDOWN-COAST boundary share(s) = ";
//          for (unsigned int n = 0; n < dVDownCoastBoundaryShare.size(); n++)
//             LogStream << dVDownCoastBoundaryShare[n] << " ";
//          LogStream << endl;
// //          LogStream << "\tTotal DOWN-COAST boundary length = " << dDownCoastTotBoundaryLen << endl;
//          // DEBUG CODE ======================================================
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Determines whether a point is within a polygon: however if the point is exactly on the edge of the polygon, then the result is indeterminate. Modified from code at http://alienryderflex.com/polygon/, our thanks to Darel Rex Finley (DarelRex@gmail.com)
//===============================================================================================================================
bool CSimulation::bIsWithinPolygon(CGeom2DPoint const* pPtStart, vector<CGeom2DPoint> const* pPtPoints)
{
   bool bOddNodes = false;

   int
      nPolyCorners = static_cast<int>(pPtPoints->size()),
      j = nPolyCorners-1;

   double
      dX = pPtStart->dGetX(),
      dY = pPtStart->dGetY();

   for (int i = 0; i < nPolyCorners; i++)
   {
      double
         dCorneriX = pPtPoints->at(i).dGetX(),
         dCorneriY = pPtPoints->at(i).dGetY(),
         dCornerjX = pPtPoints->at(j).dGetX(),
         dCornerjY = pPtPoints->at(j).dGetY();

      if ((dCorneriY < dY && dCornerjY >= dY) || (dCornerjY < dY && dCorneriY >= dY))
      {
         if (dCorneriX + (dY - dCorneriY) / (dCornerjY - dCorneriY) * (dCornerjX - dCorneriX) < dX)
         {
            bOddNodes = ! bOddNodes;
         }
      }

      j = i;
   }

  return bOddNodes;
}

//===============================================================================================================================
//! Finds a point in a polygon: is guaranteed to succeed, as every strictly closed polygon has at least one triangle that is completely contained within the polygon. Derived from an algorithm at http://stackoverflow.com/questions/9797448/get-a-point-inside-the-polygon
//===============================================================================================================================
CGeom2DPoint CSimulation::PtFindPointInPolygon(vector<CGeom2DPoint> const* pPtPoints, int const nStartPoint)
{
   int nPolySize = static_cast<int>(pPtPoints->size());
   int nOffSet = 0;
   CGeom2DPoint PtStart;

   do
   {
      // Choose three consecutive points from the polygon
      vector <CGeom2DPoint> nVTestPoints;
      for (int n = 0; n < 3; n++)
      {
         int nIndex = n + nStartPoint + nOffSet;
         if (nIndex > nPolySize-1)
            nIndex -= nPolySize;

         // Safety check
         if (nIndex < 0)
            return CGeom2DPoint(DBL_NODATA, DBL_NODATA);

         nVTestPoints.push_back(pPtPoints->at(nIndex));
      }

      // Increment ready for next time
      nOffSet++;

      // Safety check
      if (nOffSet >= (nPolySize + 3))
         return CGeom2DPoint(DBL_NODATA, DBL_NODATA);

      // Check if the halfway point between the first and the third point is inside the polygon
      PtStart = PtAverage(&nVTestPoints[0], &nVTestPoints[2]);
   }
   while (! bIsWithinPolygon(&PtStart, pPtPoints));

   return PtStart;
}
