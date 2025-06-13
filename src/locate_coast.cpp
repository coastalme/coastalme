/*!

   \file locate_coast.cpp
   \brief Finds the coastline on the raster grid
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
// #include <assert.h>
#include <cfloat>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <iomanip>
using std::resetiosflags;
using std::setiosflags;
using std::setprecision;
using std::setw;

#include <string>
using std::to_string;

#include <stack>
using std::stack;

#include "cme.h"
#include "i_line.h"
#include "line.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"

//===============================================================================================================================
//! First find all connected sea areas, then locate the vector coastline(s), then put these onto the raster grid
//===============================================================================================================================
int CSimulation::nLocateSeaAndCoasts(int& nValidCoast)
{
   // Find all connected sea cells
   FindAllSeaCells();

   // Find every coastline on the raster grid, mark raster cells, then create the vector coastline
   int nRet = nTraceAllCoasts(nValidCoast);

   if (nRet != RTN_OK)
      return nRet;

   // Have we created any coasts?
   if (m_VCoast.empty())
   {
      cerr << m_ulIter << ": " << ERR << "no coastline located: this iteration SWL = " << m_dThisIterSWL << ", maximum DEM top surface elevation = " << m_dThisIterTopElevMax << ", minimum DEM top surface elevation = " << m_dThisIterTopElevMin << endl;
      return RTN_ERR_NOCOAST;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Finds and flags all sea areas which have at least one cell at a grid edge (i.e. does not flag 'inland' seas)
//===============================================================================================================================
void CSimulation::FindAllSeaCells(void)
{
   // Go along the list of edge cells
   for (unsigned int n = 0; n < m_VEdgeCell.size(); n++)
   {
      if (m_bOmitSearchNorthEdge && m_VEdgeCellEdge[n] == NORTH)
         continue;

      if (m_bOmitSearchSouthEdge && m_VEdgeCellEdge[n] == SOUTH)
         continue;

      if (m_bOmitSearchWestEdge && m_VEdgeCellEdge[n] == WEST)
         continue;

      if (m_bOmitSearchEastEdge && m_VEdgeCellEdge[n] == EAST)
         continue;

      int nX = m_VEdgeCell[n].nGetX();
      int nY = m_VEdgeCell[n].nGetY();

      // if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth() == 0))
      if ((m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))

         // This edge cell is below SWL and sea depth remains set to zero
         FloodFillSea(nX, nY);
   }
}

//===============================================================================================================================
//! Flood-fills all sea cells starting from a given cell. The flood fill code used here is adapted from an example by Lode Vandevenne (http://lodev.org/cgtutor/floodfill.html#Scanline_Floodfill_Algorithm_With_Stack)
//===============================================================================================================================
void CSimulation::FloodFillSea(int const nXStart, int const nYStart)
{
   // For safety check
   int nRoundLoopMax = m_nXGridSize * m_nYGridSize;

   // Create an empty stack
   stack<CGeom2DIPoint> PtiStack;

   // Start at the given edge cell, push this onto the stack
   PtiStack.push(CGeom2DIPoint(nXStart, nYStart));

   // Then do the flood fill loop until there are no more cell coordinates on the stack
   int nRoundLoop = 0;

   while (! PtiStack.empty())
   {
      // Safety check
      if (nRoundLoop++ > nRoundLoopMax)
         break;

      CGeom2DIPoint Pti = PtiStack.top();
      PtiStack.pop();

      int nX = Pti.nGetX();
      int nY = Pti.nGetY();

      while ((nX >= 0) && (! m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY].bIsInundated()))
         nX--;

      nX++;

      bool bSpanAbove = false;
      bool bSpanBelow = false;

      while ((nX < m_nXGridSize) && (! m_pRasterGrid->m_Cell[nX][nY].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY].bIsInundated()) && (bFPIsEqual(m_pRasterGrid->m_Cell[nX][nY].dGetSeaDepth(), 0.0, TOLERANCE)))
      {
         // Set the sea depth for this cell
         m_pRasterGrid->m_Cell[nX][nY].SetSeaDepth();

         CRWCellLandform* pLandform = m_pRasterGrid->m_Cell[nX][nY].pGetLandform();
         int nCat = pLandform->nGetLFCategory();

         // Have we had sediment input here?
         if ((nCat == LF_CAT_SEDIMENT_INPUT) || (nCat == LF_CAT_SEDIMENT_INPUT_SUBMERGED) || (nCat == LF_CAT_SEDIMENT_INPUT_NOT_SUBMERGED))
         {
            if (m_pRasterGrid->m_Cell[nX][nY].bIsInundated())
            {
               pLandform->SetLFCategory(LF_CAT_SEDIMENT_INPUT_SUBMERGED);
               m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();

               // Set this sea cell to have deep water (off-shore) wave orientation and height, will change this later for cells closer to the shoreline if we have on-shore waves
               m_pRasterGrid->m_Cell[nX][nY].SetWaveValuesToDeepWaterWaveValues();
            }

            else
            {
               pLandform->SetLFCategory(LF_CAT_SEDIMENT_INPUT_NOT_SUBMERGED);
            }
         }

         else
         {
            // No sediment input here, just mark as sea
            m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();
            pLandform->SetLFCategory(LF_CAT_SEA);

            // Set this sea cell to have deep water (off-shore) wave orientation and height, will change this later for cells closer to the shoreline if we have on-shore waves
            m_pRasterGrid->m_Cell[nX][nY].SetWaveValuesToDeepWaterWaveValues();
         }

         // Now sort out the x-y extremities of the contiguous sea for the bounding box (used later in wave propagation)
         if (nX < m_nXMinBoundingBox)
            m_nXMinBoundingBox = nX;

         if (nX > m_nXMaxBoundingBox)
            m_nXMaxBoundingBox = nX;

         if (nY < m_nYMinBoundingBox)
            m_nYMinBoundingBox = nY;

         if (nY > m_nYMaxBoundingBox)
            m_nYMaxBoundingBox = nY;

         // Update count
         m_ulThisIterNumSeaCells++;

         if ((! bSpanAbove) && (nY > 0) && (! m_pRasterGrid->m_Cell[nX][nY-1].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY-1].bIsInundated()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY-1));
            bSpanAbove = true;
         }

         else if (bSpanAbove && (nY > 0) && (! m_pRasterGrid->m_Cell[nX][nY-1].bBasementElevIsMissingValue()) && (! m_pRasterGrid->m_Cell[nX][nY-1].bIsInundated()))
         {
            bSpanAbove = false;
         }

         if ((! bSpanBelow) && (nY < m_nYGridSize-1) && (! m_pRasterGrid->m_Cell[nX][nY+1].bBasementElevIsMissingValue()) && (m_pRasterGrid->m_Cell[nX][nY+1].bIsInundated()))
         {
            PtiStack.push(CGeom2DIPoint(nX, nY+1));
            bSpanBelow = true;
         }

         else if (bSpanBelow && (nY < m_nYGridSize-1) && (! m_pRasterGrid->m_Cell[nX][nY+1].bBasementElevIsMissingValue()) && (! m_pRasterGrid->m_Cell[nX][nY+1].bIsInundated()))
         {
            bSpanBelow = false;
         }

         nX++;
      }
   }

   // // DEBUG CODE ===========================================================================================================
   // string strOutFile = m_strOutPath + "is_contiguos_sea_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";

   // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
   // int n = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   //    for (int nX = 0; nX < m_nXGridSize; nX++)
   //    {
   //       pdRaster[n++] = m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea();
   //    }
   // }

   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_dMissingValue);
   // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   // if (nRet == CE_Failure)
   //    return;

   // GDALClose(pDataSet);
   // delete[] pdRaster;

   // string strOutFile = m_strOutPath + "is_inundated_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";

   // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];

   // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);

   // pdRaster = new double[m_nXGridSize * m_nYGridSize];
   // int n = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   //    for (int nX = 0; nX < m_nXGridSize; nX++)
   //    {
   //       pdRaster[n++] = m_pRasterGrid->m_Cell[nX][nY].bIsInundated();
   //    }
   // }

   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_dMissingValue);
   // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   // if (nRet == CE_Failure)
   //    return;

   // GDALClose(pDataSet);
   // delete[] pdRaster;
   // // DEBUG CODE ===========================================================================================================

   //    LogStream << m_ulIter << ": flood fill of sea from [" << nXStart << "][" << nYStart << "] = {" << dGridCentroidXToExtCRSX(nXStart) << ", " << dGridCentroidYToExtCRSY(nYStart) << "} with SWL = " << m_dThisIterSWL << ", " << m_ulThisIterNumSeaCells << " of " << m_ulNumCells << " cells now marked as sea (" <<  fixed << setprecision(3) << 100.0 * m_ulThisIterNumSeaCells / m_ulNumCells << " %)" << endl;

   //    LogStream << " m_nXMinBoundingBox = " << m_nXMinBoundingBox << " m_nXMaxBoundingBox = " << m_nXMaxBoundingBox << " m_nYMinBoundingBox = " << m_nYMinBoundingBox << " m_nYMaxBoundingBox = " << m_nYMaxBoundingBox << endl;
}

//===============================================================================================================================
//! Locates all the potential coastline start points on the edges of the raster grid, then traces vector coastline(s) from these start points
//===============================================================================================================================
int CSimulation::nTraceAllCoasts(int& nValidCoast)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Tracing coast" << endl;

   vector<bool> VbPossibleStartCellLHEdge;
   vector<bool> VbTraced;
   vector<int> VnSearchDirection;
   vector<CGeom2DIPoint> V2DIPossibleStartCell;

   // Go along the list of edge cells and look for possible coastline start cells
   for (unsigned int n = 0; n < m_VEdgeCell.size() - 1; n++)
   {
      if (m_bOmitSearchNorthEdge && (m_VEdgeCellEdge[n] == NORTH || m_VEdgeCellEdge[n + 1] == NORTH))
         continue;

      if (m_bOmitSearchSouthEdge && (m_VEdgeCellEdge[n] == SOUTH || m_VEdgeCellEdge[n + 1] == SOUTH))
         continue;

      if (m_bOmitSearchWestEdge && (m_VEdgeCellEdge[n] == WEST || m_VEdgeCellEdge[n + 1] == WEST))
         continue;

      if (m_bOmitSearchEastEdge && (m_VEdgeCellEdge[n] == EAST || m_VEdgeCellEdge[n + 1] == EAST))
         continue;

      int nXThis = m_VEdgeCell[n].nGetX();
      int nYThis = m_VEdgeCell[n].nGetY();
      int nXNext = m_VEdgeCell[n + 1].nGetX();
      int nYNext = m_VEdgeCell[n + 1].nGetY();

      // Get "Is it sea?" information for 'this' and 'next' cells
      bool bThisCellIsSea = m_pRasterGrid->m_Cell[nXThis][nYThis].bIsInContiguousSea();
      bool bNextCellIsSea = m_pRasterGrid->m_Cell[nXNext][nYNext].bIsInContiguousSea();

      // Are we at a coast?
      if ((! bThisCellIsSea) && bNextCellIsSea)
      {
         // 'This' cell is just inland, has it already been flagged as a possible start for a coastline (even if this subsequently 'failed' as a coastline)?
         if (! m_pRasterGrid->m_Cell[nXThis][nYThis].bIsPossibleCoastStartCell())
         {
            // It has not, so flag it
            m_pRasterGrid->m_Cell[nXThis][nYThis].SetPossibleCoastStartCell();

            if (m_nLogFileDetail >= LOG_FILE_ALL)
               LogStream << m_ulIter << ": Flagging [" << nXThis << "][" << nYThis << "] = {" << dGridCentroidXToExtCRSX(nXThis) << ", " << dGridCentroidYToExtCRSY(nYThis) << "} as possible coast start cell, with left_handed edge" << endl;

            // And save it
            V2DIPossibleStartCell.push_back(CGeom2DIPoint(nXThis, nYThis));
            VbPossibleStartCellLHEdge.push_back(true);
            VnSearchDirection.push_back(nGetOppositeDirection(m_VEdgeCellEdge[n]));
            VbTraced.push_back(false);
         }
      }

      else if (bThisCellIsSea && (! bNextCellIsSea))
      {
         // The 'next' cell is just inland, has it already been flagged as a possible start for a coastline (even if this subsequently 'failed' as a coastline)?
         if (! m_pRasterGrid->m_Cell[nXNext][nYNext].bIsPossibleCoastStartCell())
         {
            // It has not, so flag it
            m_pRasterGrid->m_Cell[nXNext][nYNext].SetPossibleCoastStartCell();

            if (m_nLogFileDetail >= LOG_FILE_ALL)
               LogStream << m_ulIter << ": Flagging [" << nXNext << "][" << nYNext << "] = {" << dGridCentroidXToExtCRSX(nXNext) << ", " << dGridCentroidYToExtCRSY(nYNext) << "} as possible coast start cell, with right_handed edge" << endl;

            // And save it
            V2DIPossibleStartCell.push_back(CGeom2DIPoint(nXNext, nYNext));
            VbPossibleStartCellLHEdge.push_back(false);
            VnSearchDirection.push_back(nGetOppositeDirection(m_VEdgeCellEdge[n + 1]));
            VbTraced.push_back(false);
         }
      }
   }

   for (unsigned int n = 0; n < V2DIPossibleStartCell.size(); n++)
   {
      if (! VbTraced[n])
      {
         int nRet = 0;

         if (VbPossibleStartCellLHEdge[n])
         {
            nRet = nTraceCoastLine(n, VnSearchDirection[n], LEFT_HANDED, &VbTraced, &V2DIPossibleStartCell);
         }

         else
         {
            nRet = nTraceCoastLine(n, VnSearchDirection[n], RIGHT_HANDED, &VbTraced, &V2DIPossibleStartCell);
         }

         if (nRet == RTN_OK)
         {
            // We have a valid coastline starting from this possible start cell
            VbTraced[n] = true;
            nValidCoast++;
         }
      }
   }

   if (nValidCoast > 0)
      return RTN_OK;

   else
   {
      cerr << m_ulIter << ": no valid coasts found" << endl;
      return RTN_ERR_TRACING_COAST;
   }
}

//===============================================================================================================================
//! Traces a coastline (which is defined to be just above still water level) on the grid using the 'wall follower' rule for maze traversal (http://en.wikipedia.org/wiki/Maze_solving_algorithm#Wall_follower). The vector coastlines are then smoothed
//===============================================================================================================================
int CSimulation::nTraceCoastLine(unsigned int const nTraceFromStartCellIndex, int const nStartSearchDirection, int const nHandedness, vector<bool>* pVbTraced, vector<CGeom2DIPoint> const* pV2DIPossibleStartCell)
{
   bool bHitStartCell = false;
   bool bAtCoast = false;
   bool bHasLeftStartEdge = false;
   bool bTooLong = false;
   bool bOffEdge = false;
   bool bRepeating = false;

   int nStartX = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetX();
   int nStartY = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetY();
   int nX = nStartX;
   int nY = nStartY;
   int nSearchDirection = nStartSearchDirection;
   int nRoundLoop = -1;
   //       nThisLen = 0;
   //       nLastLen = 0,
   //       nPreLastLen = 0;

   // Temporary coastline as integer points (grid CRS)
   CGeomILine ILTempGridCRS;

   // Mark the start cell as coast and add it to the vector object
   m_pRasterGrid->m_Cell[nStartX][nStartY].SetAsCoastline(true);
   CGeom2DIPoint PtiStart(nStartX, nStartY);
   ILTempGridCRS.Append(&PtiStart);

   // Start at this grid-edge point and trace the rest of the coastline using the 'wall follower' rule for maze traversal, trying to keep next to cells flagged as sea
   do
   {
      //       // DEBUG CODE ==============================================================================================================
      //       LogStream << "Now at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
      //       LogStream << "ILTempGridCRS is now:" << endl;
      //       for (int n = 0; n < ILTempGridCRS.nGetSize(); n++)
      //          LogStream << "[" << ILTempGridCRS[n].nGetX() << "][" << ILTempGridCRS[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[n].nGetY()) << "}" << endl;
      //       LogStream <<  "=================" << endl;
      //       // DEBUG CODE ==============================================================================================================

      // Safety check
      if (++nRoundLoop > m_nCoastMax)
      {
         bTooLong = true;

         //          LogStream << m_ulIter << ": abandoning coastline tracing from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, exceeded maximum search length (" << m_nCoastMax << ")" << endl;

         //          for (int n = 0; n < ILTempGridCRS.nGetSize(); n++)
         //             LogStream << "[" << ILTempGridCRS[n].nGetX() << "][" << ILTempGridCRS[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[n].nGetY()) << "}" << endl;
         //          LogStream << endl;

         break;
      }

      // Another safety check
      if ((nRoundLoop > 10) && (ILTempGridCRS.nGetSize() < 2))
      {
         // We've been 10 times round the loop but the coast is still less than 2 coastline points in length, so we must be repeating
         bRepeating = true;

         break;
      }

      // OK so far: so have we left the start edge?
      if (! bHasLeftStartEdge)
      {
         // We have not yet left the start edge
         if (((nStartSearchDirection == SOUTH) && (nY > nStartY)) || ((nStartSearchDirection == NORTH) && (nY < nStartY)) ||
               ((nStartSearchDirection == EAST) && (nX > nStartX)) || ((nStartSearchDirection == WEST) && (nX < nStartX)))
            bHasLeftStartEdge = true;

         // Flag this cell to ensure that it is not chosen as a coastline start cell later
         m_pRasterGrid->m_Cell[nX][nY].SetPossibleCoastStartCell();
         //          LogStream << "Flagging [" << nX << "][" << nY << "] as possible coast start cell NOT YET LEFT EDGE" << endl;
      }

      // Leave the loop if the vector coastline has left the start edge, then we find a coast cell which is a possible start cell from which a coastline has not yet been traced
      //       if (bHasLeftStartEdge && bAtCoast)
      {
         for (unsigned int nn = 0; nn < pVbTraced->size(); nn++)
         {
            if ((nn != nTraceFromStartCellIndex) && (!pVbTraced->at(nn)))
            {
               //                LogStream << "[" << pV2DIPossibleStartCell->at(nn).nGetX() << "][" << pV2DIPossibleStartCell->at(nn).nGetY() << "]" << endl;

               if (bAtCoast && (nX == pV2DIPossibleStartCell->at(nn).nGetX()) && (nY == pV2DIPossibleStartCell->at(nn).nGetY()))
               {
                  if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                     LogStream << m_ulIter << ": Valid coastline found, traced from [" << nStartX << "][" << nStartY << "] and hit another start cell at [" << nX << "][" << nY << "]" << endl;

                  pVbTraced->at(nn) = true;
                  bHitStartCell = true;
                  break;
               }
            }
         }

         //          LogStream << endl;
      }

      if (bHitStartCell)
         break;

      // OK now sort out the next iteration of the search
      int nXSeaward = 0;
      int nYSeaward = 0;
      int nSeawardNewDirection = 0;
      int nXStraightOn = 0;
      int nYStraightOn = 0;
      int nXAntiSeaward = 0;
      int nYAntiSeaward = 0;
      int nAntiSeawardNewDirection = 0;
      int nXGoBack = 0;
      int nYGoBack = 0;
      int nGoBackNewDirection = 0;

      CGeom2DIPoint Pti(nX, nY);

      // Set up the variables
      switch (nHandedness)
      {
         case RIGHT_HANDED:

            // The sea is to the right-hand side of the coast as we traverse it. We are just inland, so we need to keep heading right to find the sea
            switch (nSearchDirection)
            {
               case NORTH:
                  // The sea is towards the RHS (E) of the coast, so first try to go right (to the E)
                  nXSeaward = nX + 1;
                  nYSeaward = nY;
                  nSeawardNewDirection = EAST;

                  // If can't do this, try to go straight on (to the N)
                  nXStraightOn = nX;
                  nYStraightOn = nY - 1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (W)
                  nXAntiSeaward = nX - 1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = WEST;

                  // As a last resort, go back (to the S)
                  nXGoBack = nX;
                  nYGoBack = nY + 1;
                  nGoBackNewDirection = SOUTH;

                  break;

               case EAST:
                  // The sea is towards the RHS (S) of the coast, so first try to go right (to the S)
                  nXSeaward = nX;
                  nYSeaward = nY + 1;
                  nSeawardNewDirection = SOUTH;

                  // If can't do this, try to go straight on (to the E)
                  nXStraightOn = nX + 1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (N)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY - 1;
                  nAntiSeawardNewDirection = NORTH;

                  // As a last resort, go back (to the W)
                  nXGoBack = nX - 1;
                  nYGoBack = nY;
                  nGoBackNewDirection = WEST;

                  break;

               case SOUTH:
                  // The sea is towards the RHS (W) of the coast, so first try to go right (to the W)
                  nXSeaward = nX - 1;
                  nYSeaward = nY;
                  nSeawardNewDirection = WEST;

                  // If can't do this, try to go straight on (to the S)
                  nXStraightOn = nX;
                  nYStraightOn = nY + 1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (E)
                  nXAntiSeaward = nX + 1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = EAST;

                  // As a last resort, go back (to the N)
                  nXGoBack = nX;
                  nYGoBack = nY - 1;
                  nGoBackNewDirection = NORTH;

                  break;

               case WEST:
                  // The sea is towards the RHS (N) of the coast, so first try to go right (to the N)
                  nXSeaward = nX;
                  nYSeaward = nY - 1;
                  nSeawardNewDirection = NORTH;

                  // If can't do this, try to go straight on (to the W)
                  nXStraightOn = nX - 1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (S)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY + 1;
                  nAntiSeawardNewDirection = SOUTH;

                  // As a last resort, go back (to the E)
                  nXGoBack = nX + 1;
                  nYGoBack = nY;
                  nGoBackNewDirection = EAST;

                  break;
            }

            break;

         case LEFT_HANDED:

            // The sea is to the left-hand side of the coast as we traverse it. We are just inland, so we need to keep heading left to find the sea
            switch (nSearchDirection)
            {
               case NORTH:
                  // The sea is towards the LHS (W) of the coast, so first try to go left (to the W)
                  nXSeaward = nX - 1;
                  nYSeaward = nY;
                  nSeawardNewDirection = WEST;

                  // If can't do this, try to go straight on (to the N)
                  nXStraightOn = nX;
                  nYStraightOn = nY - 1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (E)
                  nXAntiSeaward = nX + 1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = EAST;

                  // As a last resort, go back (to the S)
                  nXGoBack = nX;
                  nYGoBack = nY + 1;
                  nGoBackNewDirection = SOUTH;

                  break;

               case EAST:
                  // The sea is towards the LHS (N) of the coast, so first try to go left (to the N)
                  nXSeaward = nX;
                  nYSeaward = nY - 1;
                  nSeawardNewDirection = NORTH;

                  // If can't do this, try to go straight on (to the E)
                  nXStraightOn = nX + 1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (S)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY + 1;
                  nAntiSeawardNewDirection = SOUTH;

                  // As a last resort, go back (to the W)
                  nXGoBack = nX - 1;
                  nYGoBack = nY;
                  nGoBackNewDirection = WEST;

                  break;

               case SOUTH:
                  // The sea is towards the LHS (E) of the coast, so first try to go left (to the E)
                  nXSeaward = nX + 1;
                  nYSeaward = nY;
                  nSeawardNewDirection = EAST;

                  // If can't do this, try to go straight on (to the S)
                  nXStraightOn = nX;
                  nYStraightOn = nY + 1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (W)
                  nXAntiSeaward = nX - 1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = WEST;

                  // As a last resort, go back (to the N)
                  nXGoBack = nX;
                  nYGoBack = nY - 1;
                  nGoBackNewDirection = NORTH;

                  break;

               case WEST:
                  // The sea is towards the LHS (S) of the coast, so first try to go left (to the S)
                  nXSeaward = nX;
                  nYSeaward = nY + 1;
                  nSeawardNewDirection = SOUTH;

                  // If can't do this, try to go straight on (to the W)
                  nXStraightOn = nX - 1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (N)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY - 1;
                  nAntiSeawardNewDirection = NORTH;

                  // As a last resort, go back (to the E)
                  nXGoBack = nX + 1;
                  nYGoBack = nY;
                  nGoBackNewDirection = EAST;

                  break;
            }

            break;
      }

      // Now do the actual search for this timestep: first try going in the direction of the sea. Is this seaward cell still within the grid?
      if (bIsWithinValidGrid(nXSeaward, nYSeaward))
      {
         // It is, so check if the cell in the seaward direction is a sea cell
         if (m_pRasterGrid->m_Cell[nXSeaward][nYSeaward].bIsInContiguousSea())
         {
            // There is sea in this seaward direction, so we are on the coast
            bAtCoast = true;

            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((bIsInterventionCell(nX, nY)) && (m_pRasterGrid->m_Cell[nX][nY].dGetInterventionTopElev() >= m_dThisIterSWL))
               {
                  // It is, so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  ILTempGridCRS.Append(&Pti);
               }

               else if (m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() >= m_dThisIterSWL)
               {
                  // The sediment top is above SWL so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  ILTempGridCRS.Append(&Pti);
               }
            }
         }

         else
         {
            // The seaward cell is not a sea cell, so we will move to it next time
            nX = nXSeaward;
            nY = nYSeaward;

            // And set a new search direction, to keep turning seaward
            nSearchDirection = nSeawardNewDirection;
            continue;
         }
      }

      // OK, we couldn't move seaward (but we may have marked the current cell as coast) so next try to move straight on. Is this straight-ahead cell still within the grid?
      if (bIsWithinValidGrid(nXStraightOn, nYStraightOn))
      {
         // It is, so check if there is sea immediately in front
         if (m_pRasterGrid->m_Cell[nXStraightOn][nYStraightOn].bIsInContiguousSea())
         {
            // Sea is in front, so we are on the coast
            bAtCoast = true;

            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((bIsInterventionCell(nX, nY)) && (m_pRasterGrid->m_Cell[nX][nY].dGetInterventionTopElev() >= m_dThisIterSWL))
               {
                  // It is, so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  ILTempGridCRS.Append(&Pti);
               }

               else if (m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() >= m_dThisIterSWL)
               {
                  // The sediment top is above SWL so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  ILTempGridCRS.Append(&Pti);
               }
            }
         }

         else
         {
            // The straight-ahead cell is not a sea cell, so we will move to it next time
            nX = nXStraightOn;
            nY = nYStraightOn;

            // The search direction remains unchanged
            continue;
         }
      }

      // Couldn't move either seaward or straight on (but we may have marked the current cell as coast) so next try to move in the anti-seaward direction. Is this anti-seaward cell still within the grid?
      if (bIsWithinValidGrid(nXAntiSeaward, nYAntiSeaward))
      {
         // It is, so check if there is sea in this anti-seaward cell
         if (m_pRasterGrid->m_Cell[nXAntiSeaward][nYAntiSeaward].bIsInContiguousSea())
         {
            // There is sea on the anti-seaward side, so we are on the coast
            bAtCoast = true;

            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsCoastline())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((bIsInterventionCell(nX, nY)) && (m_pRasterGrid->m_Cell[nX][nY].dGetInterventionTopElev() >= m_dThisIterSWL))
               {
                  // It is, so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  ILTempGridCRS.Append(&Pti);
               }

               else if (m_pRasterGrid->m_Cell[nX][nY].dGetSedimentTopElev() >= m_dThisIterSWL)
               {
                  // The sediment top is above SWL so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(true);
                  ILTempGridCRS.Append(&Pti);
               }
            }
         }

         else
         {
            // The anti-seaward cell is not a sea cell, so we will move to it next time
            nX = nXAntiSeaward;
            nY = nYAntiSeaward;

            // And set a new search direction, to keep turning seaward
            nSearchDirection = nAntiSeawardNewDirection;
            continue;
         }
      }

      // Could not move to the seaward side, move straight ahead, or move to the anti-seaward side, so we must be in a single-cell dead end! As a last resort, turn round and move back to where we just came from, but first check that this is a valid cell
      if (bIsWithinValidGrid(nXGoBack, nYGoBack))
      {
         nX = nXGoBack;
         nY = nYGoBack;

         // And change the search direction
         nSearchDirection = nGoBackNewDirection;
      }

      else
      {
         // Our final choice is not a valid cell, so give up
         bOffEdge = true;
         break;
      }
   } while (true);

   // OK, we have finished tracing this coastline on the grid. But is the coastline too long or too short?
   int nCoastSize = ILTempGridCRS.nGetSize();

   if (bOffEdge)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": ignoring temporary coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since hit off-edge cell at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, coastline size is " << nCoastSize << endl;

      // Unmark these cells as coast cells
      for (int n = 0; n < nCoastSize; n++)
         m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsCoastline(false);

      return RTN_ERR_TRACING_COAST;
   }

   if (bTooLong)
   {
      // Around loop too many times, so abandon this coastline
      if (m_nLogFileDetail >= LOG_FILE_ALL)
      {
         LogStream << m_ulIter << ": ignoring temporary coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since round loop " << nRoundLoop << " times (m_nCoastMax = " << m_nCoastMax << "), coastline size is " << nCoastSize;

         if (nCoastSize > 0)
            LogStream << ", it ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";

         LogStream << endl;
      }

      // Unmark these cells as coast cells
      for (int n = 0; n < nCoastSize; n++)
         m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsCoastline(false);

      return RTN_ERR_TRACING_COAST;
   }

   if (bRepeating)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
      {
         LogStream << m_ulIter << ": ignoring temporary coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since repeating, coastline size is " << nCoastSize;

         if (nCoastSize > 0)
            LogStream << ", it ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";

         LogStream << endl;
      }

      // Unmark these cells as coast cells
      for (int n = 0; n < nCoastSize; n++)
         m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsCoastline(false);

      return RTN_ERR_TRACING_COAST;
   }

   if (nCoastSize == 0)
   {
      // Zero-length coastline, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": abandoning zero-length coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}" << endl;

      return RTN_ERR_TRACING_COAST;
   }

   if (nCoastSize < m_nCoastMin)
   {
      // The vector coastline is too small, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ": ignoring temporary coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "} since size (" << nCoastSize << ") is less than minimum (" << m_nCoastMin << ")" << endl;

      // Unmark these cells as coast cells
      for (int n = 0; n < nCoastSize; n++)
         m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsCoastline(false);

      return RTN_ERR_TRACING_COAST;
   }

   // OK this new coastline is fine
   int nEndX = nX;
   int nEndY = nY;
   int nCoastEndX = ILTempGridCRS[nCoastSize - 1].nGetX();
   int nCoastEndY = ILTempGridCRS[nCoastSize - 1].nGetY();

   if ((nCoastEndX != nEndX) || (nCoastEndY != nEndY))
   {
      // The grid-edge cell at nEndX, nEndY is not already at end of ILTempGridCRS. But is the final cell in ILTempGridCRS already at the edge of the grid?
      if (! m_pRasterGrid->m_Cell[nCoastEndX][nCoastEndY].bIsBoundingBoxEdge())
      {
         // The final cell in ILTempGridCRS is not a grid-edge cell, so add the grid-edge cell and mark the cell as coastline
         ILTempGridCRS.Append(nEndX, nEndY);
         nCoastSize++;

         m_pRasterGrid->m_Cell[nEndX][nEndY].SetAsCoastline(true);
      }
   }

   // Need to specify start edge and end edge for smoothing routines
   int nStartEdge = m_pRasterGrid->m_Cell[nStartX][nStartY].nGetBoundingBoxEdge();
   int nEndEdge = m_pRasterGrid->m_Cell[nEndX][nEndY].nGetBoundingBoxEdge();

   // Next, convert the grid coordinates in ILTempGridCRS (integer values stored as doubles) to external CRS coordinates (which will probably be non-integer, again stored as doubles). This is done now, so that smoothing is more effective
   CGeomLine LTempExtCRS;

   for (int j = 0; j < nCoastSize; j++)
   {
      LTempExtCRS.Append(dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()), dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()));
   }

   // Now do some smoothing of the vector output, if desired
   if (m_nCoastSmooth == SMOOTH_RUNNING_MEAN)
      LTempExtCRS = LSmoothCoastRunningMean(&LTempExtCRS);

   else if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
      LTempExtCRS = LSmoothCoastSavitzkyGolay(&LTempExtCRS, nStartEdge, nEndEdge);

   //    // DEBUG CODE ==================================================================================================
   //    LogStream << "==================================" << endl;
   //    for (int j = 0; j < nCoastSize; j++)
   //    {
   //       LogStream << "{" << dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()) << "}" << "\t{" << LTempExtCRS.dGetXAt(j) << ", " << LTempExtCRS.dGetYAt(j) << "}" << endl;
   //    }
   //    LogStream << "==================================" << endl;
   //    // DEBUG CODE ==================================================================================================

   // Create a new coastline object and append to it the vector of coastline objects
   CRWCoast CoastTmp(this);
   m_VCoast.push_back(CoastTmp);
   int nCoast = static_cast<int>(m_VCoast.size()) - 1;

   // Set the coastline (Ext CRS)
   m_VCoast[nCoast].SetCoastlineExtCRS(&LTempExtCRS);

   // Set the coastline (Grid CRS)
   m_VCoast[nCoast].SetCoastlineGridCRS(&ILTempGridCRS);

   //    CGeom2DPoint PtLast(DBL_MIN, DBL_MIN);
   //    for (int j = 0; j < nCoastSize; j++)
   //    {
   //       // Store the smoothed points (in external CRS) in the coast's m_LCoastlineExtCRS object, also append dummy values to the other attribute vectors
   //       if (PtLast != &LTempExtCRS[j])        // Avoid duplicate points
   //       {
   //          m_VCoast[nCoast].AppendPointToCoastlineExtCRS(LTempExtCRS[j].dGetX(), LTempExtCRS[j].dGetY());
   //
   //          // Also store the locations of the corresponding unsmoothed points (in raster grid CRS) in the coast's m_ILCellsMarkedAsCoastline vector
   //          m_VCoast[nCoast].AppendCellMarkedAsCoastline(&ILTempGridCRS[j]);
   //       }
   //
   //       PtLast = LTempExtCRS[j];
   //    }

   // Set values for the coast's other attributes: set the coast's handedness, and start and end edges
   m_VCoast[nCoast].SetSeaHandedness(nHandedness);
   m_VCoast[nCoast].SetStartEdge(nStartEdge);
   m_VCoast[nCoast].SetEndEdge(nEndEdge);

   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
   {
      LogStream << m_ulIter << ": Coast " << nCoast << " created, from [" << nStartX << "][" << nStartY << "] to [" << nEndX << "][" << nEndY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to {" << dGridCentroidXToExtCRSX(nEndX) << ", " << dGridCentroidYToExtCRSY(nEndY) << "} with " << nCoastSize << " points, handedness = " << (nHandedness == LEFT_HANDED ? "left" : "right") << endl;

      LogStream << m_ulIter << ": Smoothed coastline " << nCoast << " runs from {" << LTempExtCRS[0].dGetX() << ", " << LTempExtCRS[0].dGetY() << "} to {" << LTempExtCRS[nCoastSize - 1].dGetX() << ", " << LTempExtCRS[nCoastSize - 1].dGetY() << "} i.e. from the ";

      if (nStartEdge == NORTH)
         LogStream << "north";

      else if (nStartEdge == SOUTH)
         LogStream << "south";

      else if (nStartEdge == WEST)
         LogStream << "west";

      else if (nStartEdge == EAST)
         LogStream << "east";

      LogStream << " edge to the ";

      if (nEndEdge == NORTH)
         LogStream << "north";

      else if (nEndEdge == SOUTH)
         LogStream << "south";

      else if (nEndEdge == WEST)
         LogStream << "west";

      else if (nEndEdge == EAST)
         LogStream << "east";

      LogStream << " edge" << endl;
   }

   //    LogStream << "-----------------" << endl;
   //    for (int kk = 0; kk < m_VCoast.back().nGetCoastlineSize(); kk++)
   //       LogStream << kk << " [" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX() << "][" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY()) << "}" << endl;
   //    LogStream << "-----------------" << endl;

   // Next calculate the curvature of the vector coastline
   DoCoastCurvature(nCoast, nHandedness);

   // Calculate values for the coast's flux orientation vector
   CalcCoastTangents(nCoast);

   // And create the vector of pointers to coastline-normal objects
   m_VCoast[nCoast].CreateProfilesAtCoastPoints();

   return RTN_OK;
}

//===============================================================================================================================
//! First find all connected sea areas, then locate the vector coastline(s), then put these onto the raster grid
//===============================================================================================================================
int CSimulation::nLocateFloodAndCoasts(void)
{
   // Find all connected sea cells
   FindAllInundatedCells();

   // Find every coastline on the raster grid, mark raster cells, then create the vector coastline
   int nRet = nTraceAllFloodCoasts();

   if (nRet != RTN_OK)
      return nRet;

   // Have we created any coasts?
   switch (m_nLevel)
   {
      case 0: // WAVESETUP + SURGE:
      {
         if (m_VFloodWaveSetupSurge.empty())
         {
            cerr << m_ulIter << ": " << ERR << "no flood coastline located: this iteration SWL = " << m_dThisIterSWL << ", maximum DEM top surface elevation = " << m_dThisIterTopElevMax << ", minimum DEM top surface elevation = " << m_dThisIterTopElevMin << endl;
            return RTN_ERR_NOCOAST;
         }

         break;
      }

      case 1: // WAVESETUP + SURGE + RUNUP:
      {
         if (m_VFloodWaveSetupSurgeRunup.empty())
         {
            cerr << m_ulIter << ": " << ERR << "no flood coastline located: this iteration SWL = " << m_dThisIterSWL << ", maximum DEM top surface elevation = " << m_dThisIterTopElevMax << ", minimum DEM top surface elevation = " << m_dThisIterTopElevMin << endl;
            return RTN_ERR_NOCOAST;
         }

         break;
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Finds and flags all sea areas which have at least one cell at a grid edge (i.e. does not flag 'inland' seas)
//===============================================================================================================================
int CSimulation::FindAllInundatedCells(void)
{
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         m_pRasterGrid->m_Cell[nX][nY].UnSetCheckFloodCell();
         m_pRasterGrid->m_Cell[nX][nY].UnSetInContiguousFlood();
         m_pRasterGrid->m_Cell[nX][nY].SetAsFloodLine(false);
      }
   }

   // Go along the list of edge cells
   for (unsigned int n = 0; n < m_VEdgeCell.size(); n++)
   {
      if (m_bOmitSearchNorthEdge && m_VEdgeCellEdge[n] == NORTH)
         continue;

      if (m_bOmitSearchSouthEdge && m_VEdgeCellEdge[n] == SOUTH)
         continue;

      if (m_bOmitSearchWestEdge && m_VEdgeCellEdge[n] == WEST)
         continue;

      if (m_bOmitSearchEastEdge && m_VEdgeCellEdge[n] == EAST)
         continue;

      int nX = m_VEdgeCell[n].nGetX();
      int nY = m_VEdgeCell[n].nGetY();

      if ((! m_pRasterGrid->m_Cell[nX][nY].bIsCellFloodCheck()) && (m_pRasterGrid->m_Cell[nX][nY].bIsInundated()))
      {
         // This edge cell is below SWL and sea depth remains set to zero
         FloodFillLand(nX, nY);
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Use the sealevel, wave set-up and run-up to evaluate flood hydraulically connected TODO 007 Not clear why we need this. We already have a flood fill sea routine: every cell that isn't sea is land
//===============================================================================================================================
void CSimulation::FloodFillLand(int const nXStart, int const nYStart)
{
   // The flood is at a user-specified location. So get the location from values read from the shapefile
   long unsigned int nLocIDs = m_VdFloodLocationX.size();
   double dDiffTotWaterLevel = 0;

   double dAuxWaterLevelDiff = 0;

   if (nLocIDs == 0)
   {
      int pointCounter = 0;

      for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
      {
         int nCoastSize = m_VCoast[nCoast].pLGetCoastlineExtCRS()->nGetSize();

         for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
         {
            dAuxWaterLevelDiff = m_VCoast[nCoast].dGetLevel(nCoastPoint, m_nLevel);

            if (! isnan(dAuxWaterLevelDiff))
            {
               if (tAbs(dAuxWaterLevelDiff) < 1)       // Limiting the maximum value that can be found (dAuxWaterLevelDiff != DBL_NODATA)
               {
                  pointCounter++;
                  dDiffTotWaterLevel += dAuxWaterLevelDiff;
               }
            }
         }
      }

      if (pointCounter > 0)
         dDiffTotWaterLevel /= pointCounter;

      else
         dDiffTotWaterLevel = 0;
   }

   else
   {
      for (long unsigned int n = 0; n < nLocIDs; n++)
      {
         double dPointGridXExtCRS = m_VdFloodLocationX[n];
         double dPointGridYExtCRS = m_VdFloodLocationY[n];
         double dMinDiffTotWaterLevelAtCoast = 1e10;

         for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
         {
            int nCoastSize = m_VCoast[nCoast].pLGetCoastlineExtCRS()->nGetSize();
            double dMinDistSquare = 1e10;

            for (int nCoastPoint = 0; nCoastPoint < nCoastSize; nCoastPoint++)
            {
               double dCoastPointXExtCRS = m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint)->dGetX();
               double dCoastPointYExtCRS = m_VCoast[nCoast].pPtGetCoastlinePointExtCRS(nCoastPoint)->dGetY();

               double dDistSquare = (dCoastPointXExtCRS - dPointGridXExtCRS) * (dCoastPointXExtCRS - dPointGridXExtCRS) + (dCoastPointYExtCRS - dPointGridYExtCRS) * (dCoastPointYExtCRS - dPointGridYExtCRS);

               if (dDistSquare < dMinDistSquare)
               {
                  dAuxWaterLevelDiff = m_VCoast[nCoast].dGetLevel(nCoastPoint, m_nLevel);

                  if (! isnan(dAuxWaterLevelDiff))
                  {
                     dMinDistSquare = dDistSquare;
                     dMinDiffTotWaterLevelAtCoast = dAuxWaterLevelDiff;
                  }
               }
            }
         }

         dDiffTotWaterLevel += dMinDiffTotWaterLevelAtCoast;
      }

      dDiffTotWaterLevel /= static_cast<double>(nLocIDs);
   }

   m_dThisIterDiffTotWaterLevel = dDiffTotWaterLevel;

   switch (m_nLevel)
   {
      case 0: // WAVESETUP + STORMSURGE:
         m_dThisIterDiffWaveSetupSurgeWaterLevel = m_dThisIterDiffTotWaterLevel;
         break;

      case 1: // WAVESETUP + STORMSURGE + RUNUP:
         m_dThisIterDiffWaveSetupSurgeRunupWaterLevel = m_dThisIterDiffTotWaterLevel;
         break;
   }

   // Create an empty stack
   stack<CGeom2DIPoint> PtiStackFlood;

   // Start at the given edge cell, push this onto the stack
   PtiStackFlood.push(CGeom2DIPoint(nXStart, nYStart));

   // Then do the flood fill loop until there are no more cell coordinates on the stack
   while (! PtiStackFlood.empty())
   {
      CGeom2DIPoint Pti = PtiStackFlood.top();
      PtiStackFlood.pop();

      int nX = Pti.nGetX();
      int nY = Pti.nGetY();

      while (nX >= 0)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsCellFloodCheck())
            break;

         if (! m_pRasterGrid->m_Cell[nX][nY].bIsElevLessThanWaterLevel())
            break;

         nX--;
      }

      nX++;

      bool bSpanAbove = false;
      bool bSpanBelow = false;

      while (nX < m_nXGridSize)
      {
         if (m_pRasterGrid->m_Cell[nX][nY].bIsCellFloodCheck())
            break;

         if (! m_pRasterGrid->m_Cell[nX][nY].bIsElevLessThanWaterLevel())
            break;

         // Flood this cell
         m_pRasterGrid->m_Cell[nX][nY].SetCheckFloodCell();
         m_pRasterGrid->m_Cell[nX][nY].SetInContiguousFlood();

         switch (m_nLevel)
         {
            case 0: // WAVESETUP + STORMSURGE:
               m_pRasterGrid->m_Cell[nX][nY].SetFloodBySetupSurge();
               break;

            case 1: // WAVESETUP + STORMSURGE + RUNUP:
               m_pRasterGrid->m_Cell[nX][nY].SetFloodBySetupSurgeRunup();
               break;
         }

         if ((! bSpanAbove) && (nY > 0) && (m_pRasterGrid->m_Cell[nX][nY - 1].bIsElevLessThanWaterLevel()) && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsCellFloodCheck()))
         {
            PtiStackFlood.push(CGeom2DIPoint(nX, nY - 1));
            bSpanAbove = true;
         }

         else if (bSpanAbove && (nY > 0) && (! m_pRasterGrid->m_Cell[nX][nY - 1].bIsElevLessThanWaterLevel()))
         {
            bSpanAbove = false;
         }

         if ((! bSpanBelow) && (nY < m_nYGridSize - 1) && (m_pRasterGrid->m_Cell[nX][nY + 1].bIsElevLessThanWaterLevel()) && (! m_pRasterGrid->m_Cell[nX][nY + 1].bIsCellFloodCheck()))
         {
            PtiStackFlood.push(CGeom2DIPoint(nX, nY + 1));
            bSpanBelow = true;
         }

         else if (bSpanBelow && (nY < m_nYGridSize - 1) && (! m_pRasterGrid->m_Cell[nX][nY + 1].bIsElevLessThanWaterLevel()))
         {
            bSpanBelow = false;
         }

         nX++;
      }
   }
}

//===============================================================================================================================
//! Locates all the potential coastline start points on the edges of the raster grid, then traces vector coastline(s) from these start points
//===============================================================================================================================
int CSimulation::nTraceAllFloodCoasts(void)
{
   vector<bool> VbPossibleStartCellLHEdge;
   vector<bool> VbTraced;
   vector<int> VnSearchDirection;
   vector<CGeom2DIPoint> V2DIPossibleStartCell;

   // Go along the list of edge cells and look for possible coastline start cells
   for (unsigned int n = 0; n < m_VEdgeCell.size() - 1; n++)
   {
      if (m_bOmitSearchNorthEdge && (m_VEdgeCellEdge[n] == NORTH || m_VEdgeCellEdge[n + 1] == NORTH))
         continue;

      if (m_bOmitSearchSouthEdge && (m_VEdgeCellEdge[n] == SOUTH || m_VEdgeCellEdge[n + 1] == SOUTH))
         continue;

      if (m_bOmitSearchWestEdge && (m_VEdgeCellEdge[n] == WEST || m_VEdgeCellEdge[n + 1] == WEST))
         continue;

      if (m_bOmitSearchEastEdge && (m_VEdgeCellEdge[n] == EAST || m_VEdgeCellEdge[n + 1] == EAST))
         continue;

      int nXThis = m_VEdgeCell[n].nGetX();
      int nYThis = m_VEdgeCell[n].nGetY();
      int nXNext = m_VEdgeCell[n + 1].nGetX();
      int nYNext = m_VEdgeCell[n + 1].nGetY();

      // Get "Is it sea?" information for 'this' and 'next' cells
      bool bThisCellIsSea = m_pRasterGrid->m_Cell[nXThis][nYThis].bIsInContiguousSeaArea();
      bool bNextCellIsSea = m_pRasterGrid->m_Cell[nXNext][nYNext].bIsInContiguousSeaArea();

      // Are we at a coast?
      if ((! bThisCellIsSea) && bNextCellIsSea)
      {
         // 'This' cell is just inland, has it already been flagged as a possible start for a coastline (even if this subsequently 'failed' as a coastline)?
         // if (! m_pRasterGrid->m_Cell[nXThis][nYThis].bIsPossibleCoastStartCell())
         {
            // It has not, so flag it
            m_pRasterGrid->m_Cell[nXThis][nYThis].SetPossibleFloodStartCell();

            // And save it
            V2DIPossibleStartCell.push_back(CGeom2DIPoint(nXThis, nYThis));
            VbPossibleStartCellLHEdge.push_back(true);
            VnSearchDirection.push_back(nGetOppositeDirection(m_VEdgeCellEdge[n]));
            VbTraced.push_back(false);
         }
      }

      else if (bThisCellIsSea && (! bNextCellIsSea))
      {
         // The 'next' cell is just inland, has it already been flagged as a possible start for a coastline (even if this subsequently 'failed' as a coastline)?
         // if (! m_pRasterGrid->m_Cell[nXNext][nYNext].bIsPossibleCoastStartCell())
         {
            // It has not, so flag it
            m_pRasterGrid->m_Cell[nXNext][nYNext].SetPossibleFloodStartCell();

            // And save it
            V2DIPossibleStartCell.push_back(CGeom2DIPoint(nXNext, nYNext));
            VbPossibleStartCellLHEdge.push_back(false);
            VnSearchDirection.push_back(nGetOppositeDirection(m_VEdgeCellEdge[n + 1]));
            VbTraced.push_back(false);
         }
      }
   }

   bool bAtLeastOneCoastTraced = false;

   for (unsigned int n = 0; n < V2DIPossibleStartCell.size(); n++)
   {
      if (! VbTraced[n])
      {
         int nRet = 0;

         if (VbPossibleStartCellLHEdge[n])
         {
            nRet = nTraceFloodCoastLine(n, VnSearchDirection[n], LEFT_HANDED, &VbTraced, &V2DIPossibleStartCell);
         }

         else
         {
            nRet = nTraceFloodCoastLine(n, VnSearchDirection[n], RIGHT_HANDED, &VbTraced, &V2DIPossibleStartCell);
         }

         if (nRet == RTN_OK)
         {
            // We have a valid coastline starting from this possible start cell
            VbTraced[n] = true;
            bAtLeastOneCoastTraced = true;
         }
      }
   }

   if (bAtLeastOneCoastTraced)
      return RTN_OK;

   else
      return RTN_ERR_TRACING_COAST;
}

//===============================================================================================================================
//! Traces a coastline (which is defined to be just above still water level) on the grid using the 'wall follower' rule for maze traversal (http://en.wikipedia.org/wiki/Maze_solving_algorithm#Wall_follower). The vector coastlines are then smoothed
//===============================================================================================================================
int CSimulation::nTraceFloodCoastLine(unsigned int const nTraceFromStartCellIndex, int const nStartSearchDirection, int const nHandedness, vector<bool>* pVbTraced, vector<CGeom2DIPoint> const* pV2DIPossibleStartCell)
{
   bool bHitStartCell = false;
   bool bAtCoast = false;
   bool bHasLeftStartEdge = false;
   bool bTooLong = false;
   bool bOffEdge = false;
   bool bRepeating = false;

   int nStartX = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetX();
   int nStartY = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetY();
   int nX = nStartX;
   int nY = nStartY;
   int nSearchDirection = nStartSearchDirection;
   int nRoundLoop = -1;
   //       nThisLen = 0;
   //       nLastLen = 0,
   //       nPreLastLen = 0;

   // Temporary coastline as integer points (grid CRS)
   CGeomILine ILTempGridCRS;

   // Mark the start cell as coast and add it to the vector object
   m_pRasterGrid->m_Cell[nStartX][nStartY].SetAsFloodLine(true);
   CGeom2DIPoint PtiStart(nStartX, nStartY);
   ILTempGridCRS.Append(&PtiStart);

   // Start at this grid-edge point and trace the rest of the coastline using the 'wall follower' rule for maze traversal, trying to keep next to cells flagged as sea
   do
   {
      //       // DEBUG CODE ==============================================================================================================
      //       LogStream << "Now at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
      //       LogStream << "ILTempGridCRS is now:" << endl;
      //       for (int n = 0; n < ILTempGridCRS.nGetSize(); n++)
      //          LogStream << "[" << ILTempGridCRS[n].nGetX() << "][" << ILTempGridCRS[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[n].nGetY()) << "}" << endl;
      //       LogStream <<  "=================" << endl;
      //       // DEBUG CODE ==============================================================================================================

      // Safety check
      if (++nRoundLoop > m_nCoastMax)
      {
         bTooLong = true;

         //          LogStream << m_ulIter << ": abandoning coastline tracing from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, exceeded maximum search length (" << m_nCoastMax << ")" << endl;

         //          for (int n = 0; n < ILTempGridCRS.nGetSize(); n++)
         //             LogStream << "[" << ILTempGridCRS[n].nGetX() << "][" << ILTempGridCRS[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[n].nGetY()) << "}" << endl;
         //          LogStream << endl;

         break;
      }

      // Another safety check
      if ((nRoundLoop > 10) && (ILTempGridCRS.nGetSize() < 2))
      {
         // We've been 10 times round the loop but the coast is still less than 2 coastline points in length, so we must be repeating
         bRepeating = true;

         break;
      }

      // OK so far: so have we left the start edge?
      if (! bHasLeftStartEdge)
      {
         // We have not yet left the start edge
         if (((nStartSearchDirection == SOUTH) && (nY > nStartY)) || ((nStartSearchDirection == NORTH) && (nY < nStartY)) ||
               ((nStartSearchDirection == EAST) && (nX > nStartX)) || ((nStartSearchDirection == WEST) && (nX < nStartX)))
            bHasLeftStartEdge = true;

         // Flag this cell to ensure that it is not chosen as a coastline start cell later
         m_pRasterGrid->m_Cell[nX][nY].SetPossibleFloodStartCell();
         //          LogStream << "Flagging [" << nX << "][" << nY << "] as possible coast start cell NOT YET LEFT EDGE" << endl;
      }

      // Leave the loop if the vector coastline has left the start edge, then we find a coast cell which is a possible start cell from which a coastline has not yet been traced
      //       if (bHasLeftStartEdge && bAtCoast)
      {
         for (unsigned int nn = 0; nn < pVbTraced->size(); nn++)
         {
            if ((nn != nTraceFromStartCellIndex) && (!pVbTraced->at(nn)))
            {
               //                LogStream << "[" << pV2DIPossibleStartCell->at(nn).nGetX() << "][" << pV2DIPossibleStartCell->at(nn).nGetY() << "]" << endl;

               if (bAtCoast && (nX == pV2DIPossibleStartCell->at(nn).nGetX()) && (nY == pV2DIPossibleStartCell->at(nn).nGetY()))
               {
                  if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                     LogStream << m_ulIter << ": valid flood coastline found, traced from [" << nStartX << "][" << nStartY << "] and hit another start cell at [" << nX << "][" << nY << "]" << endl;

                  pVbTraced->at(nn) = true;
                  bHitStartCell = true;
                  break;
               }
            }
         }

         //          LogStream << endl;
      }

      if (bHitStartCell)
         break;

      // OK now sort out the next iteration of the search
      int nXSeaward = 0;
      int nYSeaward = 0;
      int nSeawardNewDirection = 0;
      int nXStraightOn = 0;
      int nYStraightOn = 0;
      int nXAntiSeaward = 0;
      int nYAntiSeaward = 0;
      int nAntiSeawardNewDirection = 0;
      int nXGoBack = 0;
      int nYGoBack = 0;
      int nGoBackNewDirection = 0;

      CGeom2DIPoint Pti(nX, nY);

      // Set up the variables
      switch (nHandedness)
      {
         case RIGHT_HANDED:

            // The sea is to the right-hand side of the coast as we traverse it. We are just inland, so we need to keep heading right to find the sea
            switch (nSearchDirection)
            {
               case NORTH:
                  // The sea is towards the RHS (E) of the coast, so first try to go right (to the E)
                  nXSeaward = nX + 1;
                  nYSeaward = nY;
                  nSeawardNewDirection = EAST;

                  // If can't do this, try to go straight on (to the N)
                  nXStraightOn = nX;
                  nYStraightOn = nY - 1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (W)
                  nXAntiSeaward = nX - 1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = WEST;

                  // As a last resort, go back (to the S)
                  nXGoBack = nX;
                  nYGoBack = nY + 1;
                  nGoBackNewDirection = SOUTH;

                  break;

               case EAST:
                  // The sea is towards the RHS (S) of the coast, so first try to go right (to the S)
                  nXSeaward = nX;
                  nYSeaward = nY + 1;
                  nSeawardNewDirection = SOUTH;

                  // If can't do this, try to go straight on (to the E)
                  nXStraightOn = nX + 1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (N)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY - 1;
                  nAntiSeawardNewDirection = NORTH;

                  // As a last resort, go back (to the W)
                  nXGoBack = nX - 1;
                  nYGoBack = nY;
                  nGoBackNewDirection = WEST;

                  break;

               case SOUTH:
                  // The sea is towards the RHS (W) of the coast, so first try to go right (to the W)
                  nXSeaward = nX - 1;
                  nYSeaward = nY;
                  nSeawardNewDirection = WEST;

                  // If can't do this, try to go straight on (to the S)
                  nXStraightOn = nX;
                  nYStraightOn = nY + 1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (E)
                  nXAntiSeaward = nX + 1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = EAST;

                  // As a last resort, go back (to the N)
                  nXGoBack = nX;
                  nYGoBack = nY - 1;
                  nGoBackNewDirection = NORTH;

                  break;

               case WEST:
                  // The sea is towards the RHS (N) of the coast, so first try to go right (to the N)
                  nXSeaward = nX;
                  nYSeaward = nY - 1;
                  nSeawardNewDirection = NORTH;

                  // If can't do this, try to go straight on (to the W)
                  nXStraightOn = nX - 1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the LHS (S)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY + 1;
                  nAntiSeawardNewDirection = SOUTH;

                  // As a last resort, go back (to the E)
                  nXGoBack = nX + 1;
                  nYGoBack = nY;
                  nGoBackNewDirection = EAST;

                  break;
            }

            break;

         case LEFT_HANDED:

            // The sea is to the left-hand side of the coast as we traverse it. We are just inland, so we need to keep heading left to find the sea
            switch (nSearchDirection)
            {
               case NORTH:
                  // The sea is towards the LHS (W) of the coast, so first try to go left (to the W)
                  nXSeaward = nX - 1;
                  nYSeaward = nY;
                  nSeawardNewDirection = WEST;

                  // If can't do this, try to go straight on (to the N)
                  nXStraightOn = nX;
                  nYStraightOn = nY - 1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (E)
                  nXAntiSeaward = nX + 1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = EAST;

                  // As a last resort, go back (to the S)
                  nXGoBack = nX;
                  nYGoBack = nY + 1;
                  nGoBackNewDirection = SOUTH;

                  break;

               case EAST:
                  // The sea is towards the LHS (N) of the coast, so first try to go left (to the N)
                  nXSeaward = nX;
                  nYSeaward = nY - 1;
                  nSeawardNewDirection = NORTH;

                  // If can't do this, try to go straight on (to the E)
                  nXStraightOn = nX + 1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (S)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY + 1;
                  nAntiSeawardNewDirection = SOUTH;

                  // As a last resort, go back (to the W)
                  nXGoBack = nX - 1;
                  nYGoBack = nY;
                  nGoBackNewDirection = WEST;

                  break;

               case SOUTH:
                  // The sea is towards the LHS (E) of the coast, so first try to go left (to the E)
                  nXSeaward = nX + 1;
                  nYSeaward = nY;
                  nSeawardNewDirection = EAST;

                  // If can't do this, try to go straight on (to the S)
                  nXStraightOn = nX;
                  nYStraightOn = nY + 1;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (W)
                  nXAntiSeaward = nX - 1;
                  nYAntiSeaward = nY;
                  nAntiSeawardNewDirection = WEST;

                  // As a last resort, go back (to the N)
                  nXGoBack = nX;
                  nYGoBack = nY - 1;
                  nGoBackNewDirection = NORTH;

                  break;

               case WEST:
                  // The sea is towards the LHS (S) of the coast, so first try to go left (to the S)
                  nXSeaward = nX;
                  nYSeaward = nY + 1;
                  nSeawardNewDirection = SOUTH;

                  // If can't do this, try to go straight on (to the W)
                  nXStraightOn = nX - 1;
                  nYStraightOn = nY;

                  // If can't do either of these, try to go anti-seaward i.e. towards the RHS (N)
                  nXAntiSeaward = nX;
                  nYAntiSeaward = nY - 1;
                  nAntiSeawardNewDirection = NORTH;

                  // As a last resort, go back (to the E)
                  nXGoBack = nX + 1;
                  nYGoBack = nY;
                  nGoBackNewDirection = EAST;

                  break;
            }

            break;
      }

      // Now do the actual search for this timestep: first try going in the direction of the sea. Is this seaward cell still within the grid?
      if (bIsWithinValidGrid(nXSeaward, nYSeaward))
      {
         // It is, so check if the cell in the seaward direction is a sea cell
         if (m_pRasterGrid->m_Cell[nXSeaward][nYSeaward].bIsInContiguousSeaArea())
         {
            // There is sea in this seaward direction, so we are on the coast
            bAtCoast = true;

            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsFloodLine())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((bIsInterventionCell(nX, nY)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsElevLessThanWaterLevel()))
               {
                  // It is, so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsFloodLine(true);
                  ILTempGridCRS.Append(&Pti);
               }

               else if (! m_pRasterGrid->m_Cell[nX][nY].bIsElevLessThanWaterLevel())
               {
                  // The sediment top is above SWL so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsFloodLine(true);
                  ILTempGridCRS.Append(&Pti);
               }
            }
         }

         else
         {
            // The seaward cell is not a sea cell, so we will move to it next time
            nX = nXSeaward;
            nY = nYSeaward;

            // And set a new search direction, to keep turning seaward
            nSearchDirection = nSeawardNewDirection;
            continue;
         }
      }

      // OK, we couldn't move seaward (but we may have marked the current cell as coast) so next try to move straight on. Is this straight-ahead cell still within the grid?
      if (bIsWithinValidGrid(nXStraightOn, nYStraightOn))
      {
         // It is, so check if there is sea immediately in front
         if (m_pRasterGrid->m_Cell[nXStraightOn][nYStraightOn].bIsInContiguousSeaArea())
         {
            // Sea is in front, so we are on the coast
            bAtCoast = true;

            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsFloodLine())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((bIsInterventionCell(nX, nY)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsElevLessThanWaterLevel()))
               {
                  // It is, so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsFloodLine(true);
                  ILTempGridCRS.Append(&Pti);
               }

               else if (! m_pRasterGrid->m_Cell[nX][nY].bIsElevLessThanWaterLevel())
               {
                  // The sediment top is above SWL so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsFloodLine(true);
                  ILTempGridCRS.Append(&Pti);
               }
            }
         }

         else
         {
            // The straight-ahead cell is not a sea cell, so we will move to it next time
            nX = nXStraightOn;
            nY = nYStraightOn;

            // The search direction remains unchanged
            continue;
         }
      }

      // Couldn't move either seaward or straight on (but we may have marked the current cell as coast) so next try to move in the anti-seaward direction. Is this anti-seaward cell still within the grid?
      if (bIsWithinValidGrid(nXAntiSeaward, nYAntiSeaward))
      {
         // It is, so check if there is sea in this anti-seaward cell
         if (m_pRasterGrid->m_Cell[nXAntiSeaward][nYAntiSeaward].bIsInContiguousSeaArea())
         {
            // There is sea on the anti-seaward side, so we are on the coast
            bAtCoast = true;

            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->m_Cell[nX][nY].bIsFloodLine())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((bIsInterventionCell(nX, nY)) && (! m_pRasterGrid->m_Cell[nX][nY].bIsElevLessThanWaterLevel()))
               {
                  // It is, so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsFloodLine(true);
                  ILTempGridCRS.Append(&Pti);
               }

               else if (! m_pRasterGrid->m_Cell[nX][nY].bIsElevLessThanWaterLevel())
               {
                  // The sediment top is above SWL so mark as coast and add it to the vector object
                  m_pRasterGrid->m_Cell[nX][nY].SetAsFloodLine(true);
                  ILTempGridCRS.Append(&Pti);
               }
            }
         }

         else
         {
            // The anti-seaward cell is not a sea cell, so we will move to it next time
            nX = nXAntiSeaward;
            nY = nYAntiSeaward;

            // And set a new search direction, to keep turning seaward
            nSearchDirection = nAntiSeawardNewDirection;
            continue;
         }
      }

      // Could not move to the seaward side, move straight ahead, or move to the anti-seaward side, so we must be in a single-cell dead end! As a last resort, turn round and move back to where we just came from, but first check that this is a valid cell
      if (bIsWithinValidGrid(nXGoBack, nYGoBack))
      {
         nX = nXGoBack;
         nY = nYGoBack;

         // And change the search direction
         nSearchDirection = nGoBackNewDirection;
      }

      else
      {
         // Our final choice is not a valid cell, so give up
         bOffEdge = true;
         break;
      }
   } while (true);

   // OK, we have finished tracing this coastline on the grid. But is the coastline too long or too short?
   int nCoastSize = ILTempGridCRS.nGetSize();

   if (bOffEdge)
   {
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ": ignoring temporary flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since hit off-edge cell at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, coastline size is " << nCoastSize << endl;

      // Unmark these cells as coast cells
      for (int n = 0; n < nCoastSize; n++)
         m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsFloodLine(false);

      return RTN_ERR_TRACING_COAST;
   }

   if (bTooLong)
   {
      // Around loop too many times, so abandon this coastline
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      {
         LogStream << m_ulIter << ": ignoring temporary flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since round loop " << nRoundLoop << " times (m_nCoastMax = " << m_nCoastMax << "), coastline size is " << nCoastSize;

         if (nCoastSize > 0)
            LogStream << ", it ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";

         LogStream << endl;
      }

      // Unmark these cells as coast cells
      for (int n = 0; n < nCoastSize; n++)
         m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsFloodLine(false);

      return RTN_ERR_TRACING_COAST;
   }

   if (bRepeating)
   {
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
      {
         LogStream << m_ulIter << ": ignoring temporary flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since repeating, coastline size is " << nCoastSize;

         if (nCoastSize > 0)
            LogStream << ", it ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";

         LogStream << endl;
      }

      // Unmark these cells as coast cells
      for (int n = 0; n < nCoastSize; n++)
         m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsFloodLine(false);

      return RTN_ERR_TRACING_COAST;
   }

   if (nCoastSize == 0)
   {
      // Zero-length coastline, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ": abandoning zero-length flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}" << endl;

      return RTN_ERR_TRACING_COAST;
   }

   if (nCoastSize < m_nCoastMin)
   {
      // The vector coastline is too small, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ": ignoring temporary flood coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "} since size (" << nCoastSize << ") is less than minimum (" << m_nCoastMin << ")" << endl;

      // Unmark these cells as coast cells
      for (int n = 0; n < nCoastSize; n++)
         m_pRasterGrid->m_Cell[ILTempGridCRS[n].nGetX()][ILTempGridCRS[n].nGetY()].SetAsFloodLine(false);

      return RTN_ERR_TRACING_COAST;
   }

   // OK this new coastline is fine
   int nEndX = nX;
   int nEndY = nY;
   int nCoastEndX = ILTempGridCRS[nCoastSize - 1].nGetX();
   int nCoastEndY = ILTempGridCRS[nCoastSize - 1].nGetY();

   if ((nCoastEndX != nEndX) || (nCoastEndY != nEndY))
   {
      // The grid-edge cell at nEndX, nEndY is not already at end of ILTempGridCRS. But is the final cell in ILTempGridCRS already at the edge of the grid?
      if (! m_pRasterGrid->m_Cell[nCoastEndX][nCoastEndY].bIsBoundingBoxEdge())
      {
         // The final cell in ILTempGridCRS is not a grid-edge cell, so add the grid-edge cell and mark the cell as coastline
         ILTempGridCRS.Append(nEndX, nEndY);
         nCoastSize++;

         m_pRasterGrid->m_Cell[nEndX][nEndY].SetAsFloodLine(true);
      }
   }

   // Need to specify start edge and end edge for smoothing routines
   // int
   //     nStartEdge = m_pRasterGrid->m_Cell[nStartX][nStartY].nGetBoundingBoxEdge(),
   //     nEndEdge = m_pRasterGrid->m_Cell[nEndX][nEndY].nGetBoundingBoxEdge();

   // Next, convert the grid coordinates in ILTempGridCRS (integer values stored as doubles) to external CRS coordinates (which will probably be non-integer, again stored as doubles). This is done now, so that smoothing is more effective
   CGeomLine LTempExtCRS;

   for (int j = 0; j < nCoastSize; j++)
   {
      LTempExtCRS.Append(dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()), dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()));
   }

   // Now do some smoothing of the vector output, if desired
   // if (m_nCoastSmooth == SMOOTH_RUNNING_MEAN)
   //    LTempExtCRS = LSmoothCoastRunningMean(&LTempExtCRS);
   // else if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
   //    LTempExtCRS = LSmoothCoastSavitzkyGolay(&LTempExtCRS, nStartEdge, nEndEdge);

   //    // DEBUG CODE ==================================================================================================
   //    LogStream << "==================================" << endl;
   //    for (int j = 0; j < nCoastSize; j++)
   //    {
   //       LogStream << "{" << dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()) << "}" << "\t{" << LTempExtCRS.dGetXAt(j) << ", " << LTempExtCRS.dGetYAt(j) << "}" << endl;
   //    }
   //    LogStream << "==================================" << endl;
   //    // DEBUG CODE ==================================================================================================

   // Create a new coastline object and append to it the vector of coastline objects
   CRWCoast CoastTmp(this);
   int nCoast;

   switch (m_nLevel)
   {
      case 0:
         m_VFloodWaveSetupSurge.push_back(CoastTmp);
         nCoast = static_cast<int>(m_VFloodWaveSetupSurge.size()) - 1;
         m_VFloodWaveSetupSurge[nCoast].SetCoastlineExtCRS(&LTempExtCRS);
         break;

      case 1:
         m_VFloodWaveSetupSurgeRunup.push_back(CoastTmp);
         nCoast = static_cast<int>(m_VFloodWaveSetupSurgeRunup.size()) - 1;
         m_VFloodWaveSetupSurgeRunup[nCoast].SetCoastlineExtCRS(&LTempExtCRS);
   }

   return RTN_OK;
}
