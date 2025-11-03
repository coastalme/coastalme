/*!
   \file locate_coast.cpp
   \brief Finds the coastline on the raster grid using constrained flood-fill
   \details Implements a seed-based constrained flood-fill algorithm that identifies sea cells
            by expanding from known sea locations (grid edges) through cells below still water level.
            This prevents incorrect flooding of isolated depressions behind elevation barriers.
            See DEPRESSION_FLOODING_BUG_FIX.md for detailed algorithm documentation.
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

#include <iostream>
using std::cerr;
using std::endl;
using std::ios;

#include <ios>
using std::fixed;

#include <iomanip>
using std::setprecision;

#include <stack>
using std::stack;

#include "cme.h"
#include "2di_point.h"
#include "i_line.h"
#include "line.h"
#include "simulation.h"
#include "raster_grid.h"
#include "coast.h"

//===============================================================================================================================
//! First find all connected sea areas, then locate the vector coastline(s), then put these onto the raster grid
//===============================================================================================================================
int CSimulation::nLocateSeaAndCoasts(void)
{
   // Find all connected sea cells
   FindAllSeaCells();

   // Find every coastline on the raster grid, mark raster cells, then create the vector coastline
   int const nRet = nTraceAllCoasts();
   if (nRet != RTN_OK)
      return nRet;

   // Have we created any coasts?
   if (m_VCoast.empty())
   {
      cerr << m_ulIter << ": " << ERR << "no coastline located: this iteration SWL = " << m_dThisIterSWL << ", maximum DEM top surface elevation = " << m_dThisIterTopElevMax << ", minimum DEM top surface elevation = " << m_dThisIterTopElevMin << endl;

      return RTN_ERR_NO_COAST;
   }

   // Is this the highest SWL so far? If so, save this for all coasts
   if (m_bHighestSWLSoFar)
   {
      m_VHighestSWLCoastLine.clear();

      for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
      {
         CGeomLine LCoast;

         LCoast = *m_VCoast[nCoast].pLGetCoastlineExtCRS();
         m_VHighestSWLCoastLine.push_back(LCoast);
      }
   }

   // Is this the lowest SWL so far? If so, save this for all coasts
   if (m_bLowestSWLSoFar)
   {
      m_VLowestSWLCoastLine.clear();

      for (int nCoast = 0; nCoast < static_cast<int>(m_VCoast.size()); nCoast++)
      {
         CGeomLine LCoast;

         LCoast = *m_VCoast[nCoast].pLGetCoastlineExtCRS();
         m_VLowestSWLCoastLine.push_back(LCoast);
      }
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Finds and flags all sea areas which have at least one cell at a grid edge (i.e. does not flag 'inland' seas)
//===============================================================================================================================
void CSimulation::FindAllSeaCells(void)
{
   // Seed-based constrained flood-fill algorithm
   // Seeds can come from either: (1) shapefile points, or (2) grid edge cells
   // Phase 1: Identify valid seed points (parallelized if using edge cells)
   // Phase 2: Sequential flood-fill from each seed with elevation constraints

   vector<CGeom2DIPoint> vSeeds;  // Will hold all valid seed points

   // Check if we have seed points from shapefile
   if (! m_VSeaFloodSeedPoint.empty())
   {
      // Use seed points from shapefile
      // Validate each seed point: must be below SWL and not yet processed
      for (const auto& seedPoint : m_VSeaFloodSeedPoint)
      {
         int const nX = seedPoint.nGetX();
         int const nY = seedPoint.nGetY();

         // Verify point is inundated and hasn't been processed yet
         if ((m_pRasterGrid->Cell(nX, nY).bIsInundated()) &&
             (bFPIsEqual(m_pRasterGrid->Cell(nX, nY).dGetSeaDepth(), 0.0, TOLERANCE)))
         {
            vSeeds.push_back(seedPoint);
         }
         else
         {
            // Seed point is not valid - log warning
            if (! m_pRasterGrid->Cell(nX, nY).bIsInundated())
            {
               LogStream << WARN << "seed point at grid [" << nX << "][" << nY << "] is above SWL (elev = "
                         << m_pRasterGrid->Cell(nX, nY).dGetSedimentTopElev() << "m, SWL = "
                         << m_dThisIterSWL << "m), skipping" << endl;
            }
         }
      }

      LogStream << "Using " << vSeeds.size() << " seed point" << (vSeeds.size() != 1 ? "s" : "")
                << " from shapefile (out of " << m_VSeaFloodSeedPoint.size() << " total)" << endl;
   }
   else
   {
      // No shapefile seeds - use grid edge cells (original behavior)
         for (unsigned int n = 0; n < m_VEdgeCell.size(); n++)
         {
            // Check edge omission flags
            if (m_bOmitSearchNorthEdge && m_VEdgeCellEdge[n] == NORTH)
               continue;

            if (m_bOmitSearchSouthEdge && m_VEdgeCellEdge[n] == SOUTH)
               continue;

            if (m_bOmitSearchWestEdge && m_VEdgeCellEdge[n] == WEST)
               continue;

            if (m_bOmitSearchEastEdge && m_VEdgeCellEdge[n] == EAST)
               continue;

            int const nX = m_VEdgeCell[n].nGetX();
            int const nY = m_VEdgeCell[n].nGetY();

            // Check if this is a valid seed: below SWL, not yet processed, and has valid data
            if ((m_pRasterGrid->Cell(nX, nY).bIsInundated()) &&
                (bFPIsEqual(m_pRasterGrid->Cell(nX, nY).dGetSeaDepth(), 0.0, TOLERANCE)))
               vSeeds.push_back(CGeom2DIPoint(nX, nY));
         }

   }

   // Phase 2: Sequential flood-fill from each seed
   // This is inherently sequential due to shared state updates
   for (const auto& seed : vSeeds)
   {
      CellByCellFillSea(seed.nGetX(), seed.nGetY());
   }
}

//===============================================================================================================================
//! Cell-by-cell fills all sea cells starting from a given seed point, using constrained flood-fill that only expands through cells below SWL (still water level). This prevents incorrect flooding of isolated depressions behind elevation barriers.
//===============================================================================================================================
void CSimulation::CellByCellFillSea(int const nXStart, int const nYStart)
{
   // For safety check
   int const nRoundLoopMax = m_nXGridSize * m_nYGridSize;

   // Create visited tracking array to avoid revisiting cells
   vector<vector<bool>> VVbVisited(m_nXGridSize, vector<bool>(m_nYGridSize, false));

   // Create an empty stack
   stack<CGeom2DIPoint> PtiStack;

   // Start at the given seed cell, push this onto the stack
   PtiStack.push(CGeom2DIPoint(nXStart, nYStart));
   VVbVisited[nXStart][nYStart] = true;

   // Do the cell-by-cell fill loop until there are no more cell coordinates on the stack
   int nRoundLoop = 0;

   while (! PtiStack.empty())
   {
      // Safety check
      if (nRoundLoop++ > nRoundLoopMax)
         break;

      CGeom2DIPoint const Pti = PtiStack.top();
      PtiStack.pop();

      int const nX = Pti.nGetX();
      int const nY = Pti.nGetY();

      // Check if this cell should be flooded (below SWL and not yet processed)
      double const dElev = m_pRasterGrid->Cell(nX, nY).dGetSedimentTopElev();

      if (dElev < m_dThisIterSWL && bFPIsEqual(m_pRasterGrid->Cell(nX, nY).dGetSeaDepth(), 0.0, TOLERANCE))
      {
         // Set the sea depth for this cell
         m_pRasterGrid->Cell(nX, nY).SetSeaDepth();

         CRWCellLandform* pLandform = m_pRasterGrid->Cell(nX, nY).pGetLandform();
         int const nCat = pLandform->nGetLFCategory();

         // Have we had sediment input here?
         if ((nCat == LF_SEDIMENT_INPUT_CONSOLIDATED) || (nCat == LF_SEDIMENT_INPUT_UNCONSOLIDATED))
         {
            if (m_pRasterGrid->Cell(nX, nY).bIsInundated())
            {
               m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();

               // Set this sea cell to have deep water (off-shore) wave orientation and height
               m_pRasterGrid->Cell(nX, nY).SetWaveValuesToDeepWaterWaveValues();
            }
         }
         else
         {
            // No sediment input here, just mark as sea
            m_pRasterGrid->m_Cell[nX][nY].SetInContiguousSea();
            pLandform->SetLFCategory(LF_SEA);

            // Set this sea cell to have deep water (off-shore) wave orientation and height
            m_pRasterGrid->Cell(nX, nY).SetWaveValuesToDeepWaterWaveValues();
         }

         // Update bounding box for wave propagation
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

         // Check all 4 neighbors (N, S, E, W) and add to stack if they meet criteria
         int const VnDX[] = {0, 1, 0, -1};  // Neighbor offsets: N, E, S, W
         int const VnDY[] = {-1, 0, 1, 0};

         for (int dir = 0; dir < 4; dir++)
         {
            int const nXN = nX + VnDX[dir];
            int const nYN = nY + VnDY[dir];

            // Check bounds
            if (nXN < 0 || nXN >= m_nXGridSize || nYN < 0 || nYN >= m_nYGridSize)
               continue;

            // Check if already visited
            if (VVbVisited[nXN][nYN])
               continue;

            // Check for missing data
            if (m_pRasterGrid->Cell(nXN, nYN).bBasementElevIsMissingValue())
               continue;

            // Get neighbor elevation
            double const dElevN = m_pRasterGrid->Cell(nXN, nYN).dGetSedimentTopElev();

            // KEY CONSTRAINT: Only expand through cells below SWL
            // This prevents flooding isolated depressions behind barriers
            if (dElevN < m_dThisIterSWL)
            {
               PtiStack.push(CGeom2DIPoint(nXN, nYN));
               VVbVisited[nXN][nYN] = true;
            }
         }
      }
   }

   // // DEBUG CODE ===========================================================================================================
   // string strOutFile = m_strOutPath + "is_contiguous_sea_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   //
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
   //    pdRaster[n++] = m_pRasterGrid->Cell(nX, nY).bIsInContiguousSea();
   //    }
   // }
   //
   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_dMissingValue);
   // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   // if (nRet == CE_Failure)
   // return;
   //
   // GDALClose(pDataSet);
   // delete[] pdRaster;
   // // DEBUG CODE ===========================================================================================================

   // // DEBUG CODE ===========================================================================================================
   // string strOutFile = m_strOutPath + "is_inundated_";
   // strOutFile += to_string(m_ulIter);
   // strOutFile += ".tif";
   //
   // GDALDriver* pDriver = GetGDALDriverManager()->GetDriverByName("gtiff");
   // GDALDataset* pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   // double* pdRaster = new double[m_nXGridSize * m_nYGridSize];
   //
   // pDataSet = pDriver->Create(strOutFile.c_str(), m_nXGridSize, m_nYGridSize, 1, GDT_Float64, m_papszGDALRasterOptions);
   // pDataSet->SetProjection(m_strGDALBasementDEMProjection.c_str());
   // pDataSet->SetGeoTransform(m_dGeoTransform);
   //
   // pdRaster = new double[m_nXGridSize * m_nYGridSize];
   // int n = 0;
   // for (int nY = 0; nY < m_nYGridSize; nY++)
   // {
   //    for (int nX = 0; nX < m_nXGridSize; nX++)
   //    {
   //       pdRaster[n++] = m_pRasterGrid->Cell(nX, nY).bIsInundated();
   //    }
   // }
   //
   // GDALRasterBand* pBand = pDataSet->GetRasterBand(1);
   // pBand = pDataSet->GetRasterBand(1);
   // pBand->SetNoDataValue(m_dMissingValue);
   // int nRet = pBand->RasterIO(GF_Write, 0, 0, m_nXGridSize, m_nYGridSize, pdRaster, m_nXGridSize, m_nYGridSize, GDT_Float64, 0, 0, NULL);
   // if (nRet == CE_Failure)
   // return;
   //
   // GDALClose(pDataSet);
   // delete[] pdRaster;
   // // DEBUG CODE ===========================================================================================================

   // // DEBUG CODE ===========================================================================================================
   // LogStream << m_ulIter << ": cell-by-cell fill of sea from [" << nXStart << "][" << nYStart << "] = {" << dGridCentroidXToExtCRSX(nXStart) << ", " << dGridCentroidYToExtCRSY(nYStart) << "} with SWL = " << m_dThisIterSWL << ", " << m_ulThisIterNumSeaCells << " of " << m_ulNumCells << " cells now marked as sea (" <<  fixed << setprecision(3) << 100.0 * m_ulThisIterNumSeaCells / m_ulNumCells << " %)" << endl;
   //
   // LogStream << " m_nXMinBoundingBox = " << m_nXMinBoundingBox << " m_nXMaxBoundingBox = " << m_nXMaxBoundingBox << " m_nYMinBoundingBox = " << m_nYMinBoundingBox << " m_nYMaxBoundingBox = " << m_nYMaxBoundingBox << endl;
   // // DEBUG CODE ===========================================================================================================
}

//===============================================================================================================================
//! Locates all the potential coastline start/finish points on the edges of the raster grid, then traces vector coastline(s) from these start points
//===============================================================================================================================
int CSimulation::nTraceAllCoasts(void)
{
   if (m_nLogFileDetail >= LOG_FILE_MIDDLE_DETAIL)
      LogStream << m_ulIter << ": Tracing coasts" << endl;


//    // TEST ================================================================
//    int const BUFFER = 10;
//    int const DUMMY_COAST_NUMBER = 99;
//    int nValidCoast = -1;
//    int nXCoastMin = tMax(m_nXMinBoundingBox + BUFFER, 0);
//    int nXCoastMax = tMin(m_nXMaxBoundingBox + BUFFER, m_nXGridSize);
//    int nYCoastMin = tMax(m_nYMinBoundingBox + BUFFER, 0);
//    int nYCoastMax = tMin(m_nYMaxBoundingBox + BUFFER, m_nYGridSize);
//
//    for (int nX = nXCoastMin; nX < nXCoastMax; nX++)
//    {
//       for (int nY = nYCoastMin; nY < nYCoastMax; nY++)
//       {
//          for (int nSearchDirection = NORTH; nSearchDirection <= NORTH_WEST; nSearchDirection++)
//          {
//             int nXAdj;
//             int nYAdj;
//
//             switch (nSearchDirection)
//             {
//             case NORTH:
//                nXAdj = nX - 1;
//                nYAdj = nY;
//
//                if (bIsWithinValidGrid(nXAdj, nYAdj))
//                {
//                   if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
//                   {
//                      m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(DUMMY_COAST_NUMBER);
//                      break;
//                   }
//                }
//
//                break;
//
//             case NORTH_EAST:
//                nXAdj = nX;
//                nYAdj = nY - 1;
//
//                if (bIsWithinValidGrid(nXAdj, nYAdj))
//                {
//                   if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
//                   {
//                      m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(DUMMY_COAST_NUMBER);
//                      break;
//                   }
//                }
//
//                break;
//
//             case EAST:
//                nXAdj = nX;
//                nYAdj = nY - 1;
//
//                if (bIsWithinValidGrid(nXAdj, nYAdj))
//                {
//                   if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
//                   {
//                      m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(DUMMY_COAST_NUMBER);
//                      break;
//                   }
//                }
//
//                break;
//
//             case SOUTH_EAST:
//                nXAdj = nX + 1;
//                nYAdj = nY;
//
//                if (bIsWithinValidGrid(nXAdj, nYAdj))
//                {
//                   if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
//                   {
//                      m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(DUMMY_COAST_NUMBER);
//                      break;
//                   }
//                }
//
//                break;
//
//             case SOUTH:
//                nXAdj = nX + 1;
//                nYAdj = nY;
//
//                if (bIsWithinValidGrid(nXAdj, nYAdj))
//                {
//                   if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
//                   {
//                      m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(DUMMY_COAST_NUMBER);
//                      break;
//                   }
//                }
//
//                break;
//
//             case SOUTH_WEST:
//                nXAdj = nX + 1;
//                nYAdj = nY;
//
//                if (bIsWithinValidGrid(nXAdj, nYAdj))
//                {
//                   if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
//                   {
//                      m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(DUMMY_COAST_NUMBER);
//                      break;
//                   }
//                }
//
//                break;
//
//             case WEST:
//                nXAdj = nX;
//                nYAdj = nY + 1;
//
//                if (bIsWithinValidGrid(nXAdj, nYAdj))
//                {
//                   if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
//                   {
//                      m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(DUMMY_COAST_NUMBER);
//                      break;
//                   }
//                }
//
//                break;
//
//             case NORTH_WEST:
//                nXAdj = nX;
//                nYAdj = nY + 1;
//
//                if (bIsWithinValidGrid(nXAdj, nYAdj))
//                {
//                   if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].bIsInContiguousSea())
//                   {
//                      m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(DUMMY_COAST_NUMBER);
//                      break;
//                   }
//                }
//
//                break;
//             }
//          }
//       }
//    }
//
//    // Now go along the list of edge cells, look for DUMMY_COAST_NUMBER
//    bool bEdgeFound = false;
//    do
//    {
//       for (unsigned int n = 0; n < m_VEdgeCell.size(); n++)
//       {
//          if (m_bOmitSearchNorthEdge && m_VEdgeCellEdge[n] == NORTH)
//             continue;
//
//          if (m_bOmitSearchSouthEdge && m_VEdgeCellEdge[n] == SOUTH)
//             continue;
//
//          if (m_bOmitSearchWestEdge && m_VEdgeCellEdge[n] == WEST)
//             continue;
//
//          if (m_bOmitSearchEastEdge && m_VEdgeCellEdge[n] == EAST)
//             continue;
//
//          int const nX = m_VEdgeCell[n].nGetX();
//          int const nY = m_VEdgeCell[n].nGetY();
//
//          if (m_pRasterGrid->m_Cell[nX][nY].nGetCoastline() == DUMMY_COAST_NUMBER)
//          {
//             bEdgeFound = true;
//             nValidCoast++;
//
//             // Set this edge cell
//             m_pRasterGrid->m_Cell[nX][nY].SetAsCoastline(nValidCoast);
//
//             // Create the coast vector
//             CGeomILine ILTempGridCRS;
//             CGeom2DIPoint const PtiStart(nX, nY);
//             ILTempGridCRS.Append(&PtiStart);
//
//             bool bAdjFound = false;
//
//             // Now look for other cells
//             do
//             {
//                for (int nSearchDirection = NORTH; nSearchDirection <= NORTH_WEST; nSearchDirection++)
//                {
//                   int nXAdj;
//                   int nYAdj;
//
//                   switch (nSearchDirection)
//                   {
//                   case NORTH:
//                      nXAdj = nX - 1;
//                      nYAdj = nY;
//
//                      if (bIsWithinValidGrid(nXAdj, nYAdj))
//                      {
//                         if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline() == DUMMY_COAST_NUMBER)
//                         {
//                            m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(nValidCoast);
//
//                            CGeom2DIPoint const PtiPoint(nX, nY);
//                            ILTempGridCRS.Append(&PtiPoint);
//
//                            bAdjFound = true;
//
//                            break;
//                         }
//                      }
//
//                      break;
//
//                   case NORTH_EAST:
//                      nXAdj = nX;
//                      nYAdj = nY - 1;
//
//                      if (bIsWithinValidGrid(nXAdj, nYAdj))
//                      {
//                         if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline() == DUMMY_COAST_NUMBER)
//                         {
//                            m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(nValidCoast);
//
//                            CGeom2DIPoint const PtiPoint(nX, nY);
//                            ILTempGridCRS.Append(&PtiPoint);
//
//                            bAdjFound = true;
//
//                            break;
//                         }
//                      }
//
//                      break;
//
//                   case EAST:
//                      nXAdj = nX;
//                      nYAdj = nY - 1;
//
//                      if (bIsWithinValidGrid(nXAdj, nYAdj))
//                      {
//                         if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline() == DUMMY_COAST_NUMBER)
//                         {
//                            m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(nValidCoast);
//
//                            CGeom2DIPoint const PtiPoint(nX, nY);
//                            ILTempGridCRS.Append(&PtiPoint);
//
//                            bAdjFound = true;
//
//                            break;
//                         }
//                      }
//
//                      break;
//
//                   case SOUTH_EAST:
//                      nXAdj = nX + 1;
//                      nYAdj = nY;
//
//                      if (bIsWithinValidGrid(nXAdj, nYAdj))
//                      {
//                         if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline() == DUMMY_COAST_NUMBER)
//                         {
//                            m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(nValidCoast);
//
//                            CGeom2DIPoint const PtiPoint(nX, nY);
//                            ILTempGridCRS.Append(&PtiPoint);
//
//                            bAdjFound = true;
//
//                            break;
//                         }
//                      }
//
//                      break;
//
//                   case SOUTH:
//                      nXAdj = nX + 1;
//                      nYAdj = nY;
//
//                      if (bIsWithinValidGrid(nXAdj, nYAdj))
//                      {
//                         if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline() == DUMMY_COAST_NUMBER)
//                         {
//                            m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(nValidCoast);
//
//                            CGeom2DIPoint const PtiPoint(nX, nY);
//                            ILTempGridCRS.Append(&PtiPoint);
//
//                            bAdjFound = true;
//
//                            break;
//                         }
//                      }
//
//                      break;
//
//                   case SOUTH_WEST:
//                      nXAdj = nX + 1;
//                      nYAdj = nY;
//
//                      if (bIsWithinValidGrid(nXAdj, nYAdj))
//                      {
//                         if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline() == DUMMY_COAST_NUMBER)
//                         {
//                            m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(nValidCoast);
//
//                            CGeom2DIPoint const PtiPoint(nX, nY);
//                            ILTempGridCRS.Append(&PtiPoint);
//
//                            bAdjFound = true;
//
//                            break;
//                         }
//                      }
//
//                      break;
//
//                   case WEST:
//                      nXAdj = nX;
//                      nYAdj = nY + 1;
//
//                      if (bIsWithinValidGrid(nXAdj, nYAdj))
//                      {
//                         if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline() == DUMMY_COAST_NUMBER)
//                         {
//                            m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(nValidCoast);
//
//                            CGeom2DIPoint const PtiPoint(nX, nY);
//                            ILTempGridCRS.Append(&PtiPoint);
//
//                            bAdjFound = true;
//
//                            break;
//                         }
//                      }
//
//                      break;
//
//                   case NORTH_WEST:
//                      nXAdj = nX;
//                      nYAdj = nY + 1;
//
//                      if (bIsWithinValidGrid(nXAdj, nYAdj))
//                      {
//                         if (m_pRasterGrid->m_Cell[nXAdj][nYAdj].nGetCoastline() == DUMMY_COAST_NUMBER)
//                         {
//                            m_pRasterGrid->m_Cell[nXAdj][nYAdj].SetAsCoastline(nValidCoast);
//
//                            CGeom2DIPoint const PtiPoint(nX, nY);
//                            ILTempGridCRS.Append(&PtiPoint);
//
//                            bAdjFound = true;
//
//                            break;
//                         }
//                      }
//
//                      break;
//                   }
//                }
//
//             } while (bAdjFound);
//
//
//
//
//
//
//          }
//
//
//       }
//    }
//    while (bEdgeFound);
//
//
//
//
//
//
//
//
//    // ============================================================*/


   int const TOOCLOSE = 1;
   int nValidCoast = 0;
   vector<bool> VbPossibleStartCellLHEdge;
   vector<bool> VbTraced;
   vector<int> VnSearchDirection;
   vector<CGeom2DIPoint> V2DIPossibleStartCell;

   // Go along the list of edge cells and look for possible coastline start/finish cells
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

      int const nXThis = m_VEdgeCell[n].nGetX();
      int const nYThis = m_VEdgeCell[n].nGetY();
      int const nXNext = m_VEdgeCell[n + 1].nGetX();
      int const nYNext = m_VEdgeCell[n + 1].nGetY();

      // Get "Is it sea?" information for 'this' and 'next' cells
      bool const bThisCellIsSea = m_pRasterGrid->Cell(nXThis, nYThis).bIsInContiguousSea();
      bool const bNextCellIsSea = m_pRasterGrid->Cell(nXNext, nYNext).bIsInContiguousSea();

      // Is one cell land and the sea?
      if ((! bThisCellIsSea) && bNextCellIsSea)
      {
         // Yes, we are at a possible coastline start cell with 'this' cell just above the sea. Has 'this' cell already been flagged as a possible start for a coastline (even if this subsequently 'failed' as a coastline)?
         if (m_pRasterGrid->m_Cell[nXThis][nYThis].bIsPossibleCoastStartCell())
            continue;

         // Is 'this' cell too close to a cell that has previously been flagged as a possible coast start cell? Search first in one direction
         bool bTooClose = false;
         for (int nn = 1; nn <= TOOCLOSE; nn++)
         {
            int const nTmp = n + nn;
            if (nTmp == static_cast<int>(m_VEdgeCell.size()))
               break;

            int const nXTmp = m_VEdgeCell[nTmp].nGetX();
            int const nYTmp = m_VEdgeCell[nTmp].nGetY();

            if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsPossibleCoastStartCell())
            {
               bTooClose = true;
               break;
            }
         }
         if (bTooClose)
            continue;

         // Now search in the other direction
         for (int nn = 1; nn <= TOOCLOSE; nn++)
         {
            int const nTmp = n - nn;
            if (nTmp < 0)
               break;

            int const nXTmp = m_VEdgeCell[nTmp].nGetX();
            int const nYTmp = m_VEdgeCell[nTmp].nGetY();

            if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsPossibleCoastStartCell())
            {
               bTooClose = true;
               break;
            }
         }
         if (bTooClose)
            continue;

         // All OK, so flag 'this' cell
         m_pRasterGrid->m_Cell[nXThis][nYThis].SetPossibleCoastStartCell();

         if (m_nLogFileDetail >= LOG_FILE_ALL)
            LogStream << m_ulIter << ": \tflagging [" << nXThis << "][" << nYThis << "] = {" << dGridCentroidXToExtCRSX(nXThis) << ", " << dGridCentroidYToExtCRSY(nYThis) << "} as possible coast start cell (left_handed edge)" << endl;

         // And save it as a possible coastline start cell
         V2DIPossibleStartCell.push_back(CGeom2DIPoint(nXThis, nYThis));
         VbPossibleStartCellLHEdge.push_back(true);
         VnSearchDirection.push_back(nGetOppositeDirection(m_VEdgeCellEdge[n]));
         VbTraced.push_back(false);
      }
      else if (bThisCellIsSea && (! bNextCellIsSea))
      {
         // We are at a possible coastline start cell with the 'next' cell just above the sea. Has the 'next' cell already been flagged as a possible start for a coastline (even if this subsequently 'failed' as a coastline)?
         if (m_pRasterGrid->m_Cell[nXNext][nYNext].bIsPossibleCoastStartCell())
            continue;

         // Is the 'next' cell too close to a cell that has previously been flagged as a possible coast start cell? Search first in one direction
         bool bTooClose = false;
         for (int nn = 1; nn <= TOOCLOSE; nn++)
         {
            int const nTmp = n + 1 + nn;
            if (nTmp == static_cast<int>(m_VEdgeCell.size()))
               break;

            int const nXTmp = m_VEdgeCell[nTmp].nGetX();
            int const nYTmp = m_VEdgeCell[nTmp].nGetY();

            if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsPossibleCoastStartCell())
            {
               bTooClose = true;
               break;
            }
         }
         if (bTooClose)
            continue;

         // Now search in the other direction
         for (int nn = 1; nn <= TOOCLOSE; nn++)
         {
            int const nTmp = n + 1 - nn;
            if (nTmp < 0)
               break;

            int const nXTmp = m_VEdgeCell[nTmp].nGetX();
            int const nYTmp = m_VEdgeCell[nTmp].nGetY();

            if (m_pRasterGrid->m_Cell[nXTmp][nYTmp].bIsPossibleCoastStartCell())
            {
               bTooClose = true;
               break;
            }
         }
         if (bTooClose)
            continue;

         // All OK, so flag the 'next' cell
         m_pRasterGrid->m_Cell[nXNext][nYNext].SetPossibleCoastStartCell();

         if (m_nLogFileDetail >= LOG_FILE_ALL)
            LogStream << m_ulIter << ": \tflagging [" << nXNext << "][" << nYNext << "] = {" << dGridCentroidXToExtCRSX(nXNext) << ", " << dGridCentroidYToExtCRSY(nYNext) << "} as possible coast start cell (right_handed edge)" << endl;

         // And save it as a possible coastline start cell
         V2DIPossibleStartCell.push_back(CGeom2DIPoint(nXNext, nYNext));
         VbPossibleStartCellLHEdge.push_back(false);
         VnSearchDirection.push_back(nGetOppositeDirection(m_VEdgeCellEdge[n + 1]));
         VbTraced.push_back(false);
      }
   }

   // Any possible coastline start/finish cells found?
   if (V2DIPossibleStartCell.size() == 0)
   {
      LogStream << m_ulIter << ": no coastline start/finish points found after grid edges searched.";

      if (m_bOmitSearchNorthEdge || m_bOmitSearchSouthEdge || m_bOmitSearchWestEdge || m_bOmitSearchEastEdge)
      {
         LogStream << " Note that the following grid edges were not searched: " << (m_bOmitSearchNorthEdge ? "N " : "") << (m_bOmitSearchSouthEdge ? "S " : "") << (m_bOmitSearchWestEdge ? "W " : "") << (m_bOmitSearchEastEdge ? "E " : "");
      }

      LogStream << endl;

      return RTN_ERR_NO_START_FINISH_POINTS_TRACING_COAST;
   }

   // Some possible coastline start/finish points were found
   // LogStream << m_ulIter << ": " << V2DIPossibleStartCell.size() << " possible coastline start/finish points found" << endl;

   // // Are any of the possible start/finsh points adjacent?
   // vector<bool> VbToRemove(V2DIPossibleStartCell.size(), false);
   // for (int nThisPoint = 0; nThisPoint < static_cast<int>(V2DIPossibleStartCell.size()); nThisPoint++)
   // {
   //    for (int nOtherPoint = 0; nOtherPoint < static_cast<int>(V2DIPossibleStartCell.size()); nOtherPoint++)
   //    {
   //       if ((nThisPoint == nOtherPoint) || VbToRemove[nThisPoint] || VbToRemove[nOtherPoint])
   //          continue;
   //
   //       if (bIsAdjacentEdgeCell(&V2DIPossibleStartCell[nThisPoint], &V2DIPossibleStartCell[nOtherPoint]))
   //          VbToRemove[nOtherPoint] = true;
   //    }
   // }

   // // Remove each start/finish point which has been flagged as adjacent
   // for (int nPoint = 0; nPoint < static_cast<int>(VbToRemove.size()); nPoint++)
   // {
   //    if (VbToRemove[nPoint])
   //    {
   //       V2DIPossibleStartCell.erase(V2DIPossibleStartCell.begin() + nPoint);
   //       VbPossibleStartCellLHEdge.erase(VbPossibleStartCellLHEdge.begin() + nPoint);
   //       VnSearchDirection.erase(VnSearchDirection.begin() + nPoint);
   //       VbTraced.erase(VbTraced.begin() + nPoint);
   //    }
   // }

   // All OK, now trace from each of these possible start/finish points
   for (int n = static_cast<int>(V2DIPossibleStartCell.size())-1; n >= 0; n--)
   {
      if (! VbTraced[n])
      {
         int nRet = 0;

         if (VbPossibleStartCellLHEdge[n])
            nRet = nTraceCoastLine(n, VnSearchDirection[n], LEFT_HANDED, &VbTraced, &V2DIPossibleStartCell);
         else
            nRet = nTraceCoastLine(n, VnSearchDirection[n], RIGHT_HANDED, &VbTraced, &V2DIPossibleStartCell);

         if (nRet == RTN_OK)
         {
            // We have a valid coastline starting from this possible start cell
            VbTraced[n] = true;
            nValidCoast++;
         }
      }
   }

   if (nValidCoast == 0)
   {
      // No valid coasts found so try again, this time working through the possible start/finish points in reverse order
      for (int n = 0; n < static_cast<int>(VbTraced.size()); n++)
         VbTraced[n] = false;

      for (int n = 0; n < static_cast<int>(V2DIPossibleStartCell.size()); n++)
      {
         if (! VbTraced[n])
         {
            int nRet = 0;

            if (VbPossibleStartCellLHEdge[n])
               nRet = nTraceCoastLine(n, VnSearchDirection[n], LEFT_HANDED, &VbTraced, &V2DIPossibleStartCell);
            else
               nRet = nTraceCoastLine(n, VnSearchDirection[n], RIGHT_HANDED, &VbTraced, &V2DIPossibleStartCell);

            if (nRet == RTN_OK)
            {
               // We have a valid coastline starting from this possible start cell
               VbTraced[n] = true;
               nValidCoast++;
            }
         }
      }
   }

   if (nValidCoast == 0)
   {
      // Still no valid coasts found, so we have to give up
      cerr << m_ulIter << ": no valid coasts found, see " << m_strLogFile << " for more information" << endl;
      return RTN_ERR_NO_VALID_COAST;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Traces a coastline (which is defined to be just above still water level) on the grid using the 'wall follower' rule for maze traversal (http://en.wikipedia.org/wiki/Maze_solving_algorithm#Wall_follower). The resulting vector coastline is then smoothed
//===============================================================================================================================
int CSimulation::nTraceCoastLine(unsigned int const nTraceFromStartCellIndex, int const nStartSearchDirection, int const nHandedness, vector<bool>* pVbTraced, vector<CGeom2DIPoint> const* pV2DIPossibleStartCell)
{
   bool bHitStartCell = false;
   bool bAtCoast = false;
   bool bHasLeftStartEdge = false;
   bool bTooLong = false;
   bool bOffEdge = false;
   bool bRepeating = false;

   int const nStartX = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetX();
   int const nStartY = pV2DIPossibleStartCell->at(nTraceFromStartCellIndex).nGetY();
   int nX = nStartX;
   int nY = nStartY;
   int nSearchDirection = nStartSearchDirection;
   int nRoundLoop = -1;
   // nThisLen = 0;
   // nLastLen = 0,
   // nPreLastLen = 0;

   // Temporary coastline as integer points (grid CRS)
   CGeomILine ILTempGridCRS;

   // Add the start cell to the vector
   CGeom2DIPoint const PtiStart(nStartX, nStartY);
   ILTempGridCRS.Append(&PtiStart);

   // Start at this grid-edge point and trace the rest of the coastline using the 'wall follower' rule for maze traversal, trying to keep next to cells flagged as sea
   do
   {
      //       // DEBUG CODE ==============================================================================================================
      // LogStream << "Now at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;
      // LogStream << "ILTempGridCRS is now:" << endl;
      // for (int n = 0; n < ILTempGridCRS.nGetSize(); n++)
      // LogStream << "[" << ILTempGridCRS[n].nGetX() << "][" << ILTempGridCRS[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[n].nGetY()) << "}" << endl;
      // LogStream <<  "=================" << endl;
      //       // DEBUG CODE ==============================================================================================================

      // Safety check
      if (++nRoundLoop > m_nCoastMax)
      {
         bTooLong = true;

         LogStream << m_ulIter << ": \tabandoning possible coastline, traced from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, exceeded maximum search length (" << m_nCoastMax << ")" << endl;

         // for (int n = 0; n < ILTempGridCRS.nGetSize(); n++)
         // LogStream << "[" << ILTempGridCRS[n].nGetX() << "][" << ILTempGridCRS[n].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[n].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[n].nGetY()) << "}" << endl;
         // LogStream << endl;

         break;
      }

      // Another safety check
      if ((nRoundLoop > 10) && (ILTempGridCRS.nGetSize() < 2))
      {
         // We've been 10 times round the loop but the coast is still less than 2 coastline points in length, so we must be repeating
         bRepeating = true;

         LogStream << m_ulIter << ": \tabandoning possible coastline, traced from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, is looping" << endl;

         break;
      }

      // OK so far: so have we left the start edge?
      if (! bHasLeftStartEdge)
      {
         // We have not yet left the start edge
         if (((nStartSearchDirection == SOUTH) && (nY > nStartY)) || ((nStartSearchDirection == NORTH) && (nY < nStartY)) || ((nStartSearchDirection == EAST) && (nX > nStartX)) || ((nStartSearchDirection == WEST) && (nX < nStartX)))
            bHasLeftStartEdge = true;

         // Flag this cell to ensure that it is not chosen as a coastline start cell later
         m_pRasterGrid->Cell(nX, nY).SetPossibleCoastStartCell();
         // LogStream << "Flagging [" << nX << "][" << nY << "] as possible coast start cell NOT YET LEFT EDGE" << endl;
      }

      // If the vector coastline has left the start edge, and we hit a possible coast start point from which a coastline has not yet been traced, then leave the loop
      // LogStream << "bHasLeftStartEdge = " << bHasLeftStartEdge << " bAtCoast = " << bAtCoast << endl;
      if (bHasLeftStartEdge && bAtCoast)
      {
         for (unsigned int nn = 0; nn < pVbTraced->size(); nn++)
         {
            bool const bTraced = pVbTraced->at(nn);
            if ((nn != nTraceFromStartCellIndex) && (! bTraced))
            {
               int const nXPoss = pV2DIPossibleStartCell->at(nn).nGetX();
               int const nYPoss = pV2DIPossibleStartCell->at(nn).nGetY();

               // LogStream << "In 'Leave the edge' loop for [" << nX << "][" << nY << "] bTraced = " << bTraced << " nn = " << nn << " nTraceFromStartCellIndex = " << nTraceFromStartCellIndex << " possible start cell = [" << nXPoss << "][" << nYPoss << "]" << endl;

               if ((nX == nXPoss) && (nY == nYPoss))
               {
                  if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
                     LogStream << m_ulIter << ": \tpossible coastline found, traced from [" << nStartX << "][" << nStartY << "]  = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}, ended at possible coast start cell at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}" << endl;

                  pVbTraced->at(nn) = true;
                  bHitStartCell = true;
                  break;
               }
            }
         }
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

      CGeom2DIPoint const Pti(nX, nY);

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
         if (m_pRasterGrid->Cell(nXSeaward, nYSeaward).bIsInContiguousSea())
         {
            // There is sea in this seaward direction, so we are on the coast
            bAtCoast = true;

            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->Cell(nX, nY).bIsCoastline())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((bIsInterventionCell(nX, nY)) && (m_pRasterGrid->Cell(nX, nY).dGetInterventionTopElev() >= m_dThisIterSWL))
               {
                  // It is, so add it to the vector
                  ILTempGridCRS.AppendIfNotPrevious(&Pti);
               }
               else if (m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() >= m_dThisIterSWL)
               {
                  // The sediment top (inc any talus) is above SWL so add it to the vector object
                  ILTempGridCRS.AppendIfNotPrevious(&Pti);
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
         if (m_pRasterGrid->Cell(nXStraightOn, nYStraightOn).bIsInContiguousSea())
         {
            // Sea is in front, so we are on the coast
            bAtCoast = true;

            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->Cell(nX, nY).bIsCoastline())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((bIsInterventionCell(nX, nY)) && (m_pRasterGrid->Cell(nX, nY).dGetInterventionTopElev() >= m_dThisIterSWL))
               {
                  // It is, so add it to the vector object
                  ILTempGridCRS.AppendIfNotPrevious(&Pti);
               }
               else if (m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() >= m_dThisIterSWL)
               {
                  // The sediment top (inc any talus) is above SWL so add it to the vector object
                  ILTempGridCRS.AppendIfNotPrevious(&Pti);
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
         if (m_pRasterGrid->Cell(nXAntiSeaward, nYAntiSeaward).bIsInContiguousSea())
         {
            // There is sea on the anti-seaward side, so we are on the coast
            bAtCoast = true;

            // Has the current cell already marked been marked as a coast cell?
            if (! m_pRasterGrid->Cell(nX, nY).bIsCoastline())
            {
               // Not already marked, is this an intervention cell with the top above SWL?
               if ((bIsInterventionCell(nX, nY)) && (m_pRasterGrid->Cell(nX, nY).dGetInterventionTopElev() >= m_dThisIterSWL))
               {
                  // It is, so add it to the vector object
                  ILTempGridCRS.AppendIfNotPrevious(&Pti);
               }
               else if (m_pRasterGrid->m_Cell[nX][nY].dGetAllSedTopElevIncTalus() >= m_dThisIterSWL)
               {
                  // The sediment top (inc any talus) is above SWL so add it to the vector object
                  ILTempGridCRS.AppendIfNotPrevious(&Pti);
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

   // OK, we have a coastline. So is the coastline too long or too short?
   int nCoastSize = ILTempGridCRS.nGetSize();

   if (bOffEdge)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": \t**** TEST abandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since hit off-edge cell at [" << nX << "][" << nY << "] = {" << dGridCentroidXToExtCRSX(nX) << ", " << dGridCentroidYToExtCRSY(nY) << "}, coastline size is " << nCoastSize << endl;

      // return RTN_ERR_IGNORING_COAST;
   }

   if (bTooLong)
   {
      // Around loop too many times, so abandon this coastline
      if (m_nLogFileDetail >= LOG_FILE_ALL)
      {
         LogStream << m_ulIter << ": \tabandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since round loop " << nRoundLoop << " times, coastline size is " << nCoastSize;

         if (nCoastSize > 0)
            LogStream << ", ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";

         LogStream << endl;
      }

      return RTN_ERR_TOO_LONG_TRACING_COAST;
   }

   if (bRepeating)
   {
      if (m_nLogFileDetail >= LOG_FILE_ALL)
      {
         LogStream << m_ulIter << ": abandon possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} since repeating, coastline size is " << nCoastSize;

         if (nCoastSize > 0)
            LogStream << ", it ended at [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "}";

         LogStream << endl;
      }

      return RTN_ERR_REPEATING_WHEN_TRACING_COAST;
   }

   if (nCoastSize == 0)
   {
      // Zero-length coastline, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_ALL)
         LogStream << m_ulIter << ": abandoning zero-length coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "}" << endl;

      return RTN_ERR_ZERO_LENGTH_COAST;
   }

   if (nCoastSize < m_nCoastMin)
   {
      // The vector coastline is too small, so abandon it
      if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
         LogStream << m_ulIter << ": \tabandoning possible coastline from [" << nStartX << "][" << nStartY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to [" << ILTempGridCRS[nCoastSize - 1].nGetX() << "][" << ILTempGridCRS[nCoastSize - 1].nGetY() << "] = {" << dGridCentroidXToExtCRSX(ILTempGridCRS[nCoastSize - 1].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[nCoastSize - 1].nGetY()) << "} since size (" << nCoastSize << ") is less than minimum (" << m_nCoastMin << ")" << endl;

      return RTN_ERR_COAST_TOO_SMALL;
   }

   // OK this new coastline is fine
   int const nEndX = nX;
   int const nEndY = nY;
   int const nCoastEndX = ILTempGridCRS[nCoastSize - 1].nGetX();
   int const nCoastEndY = ILTempGridCRS[nCoastSize - 1].nGetY();

   if ((nCoastEndX != nEndX) || (nCoastEndY != nEndY))
   {
      // The grid-edge cell at nEndX, nEndY is not already at end of ILTempGridCRS. But is the final cell in ILTempGridCRS already at the edge of the grid?
      if (! m_pRasterGrid->Cell(nCoastEndX, nCoastEndY).bIsBoundingBoxEdge())
      {
         // The final cell in ILTempGridCRS is not a grid-edge cell, so add the grid-edge cell and mark the cell as coastline
         ILTempGridCRS.AppendIfNotPrevious(nEndX, nEndY);
         nCoastSize++;
      }
   }

   // Need to specify start edge and end edge for smoothing routines
   int const nStartEdge = m_pRasterGrid->Cell(nStartX, nStartY).nGetBoundingBoxEdge();
   int const nEndEdge = m_pRasterGrid->Cell(nEndX, nEndY).nGetBoundingBoxEdge();

   // Next, convert the grid coordinates in ILTempGridCRS (integer values stored as doubles) to external CRS coordinates (which will probably be non-integer, again stored as doubles). This is done now, so that smoothing is more effective
   CGeomLine LTempExtCRS;

   for (int j = 0; j < nCoastSize; j++)
      LTempExtCRS.Append(dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()), dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()));

   // Now do some smoothing of the vector output, if desired
   if (m_nCoastSmooth == SMOOTH_RUNNING_MEAN)
      LTempExtCRS = LSmoothCoastRunningMean(&LTempExtCRS);
   else if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
      LTempExtCRS = LSmoothCoastSavitzkyGolay(&LTempExtCRS, nStartEdge, nEndEdge);

   //    // DEBUG CODE ==================================================================================================
   // LogStream << "==================================" << endl;
   // for (int j = 0; j < nCoastSize; j++)
   // {
   // LogStream << "{" << dGridCentroidXToExtCRSX(ILTempGridCRS[j].nGetX()) << ", " << dGridCentroidYToExtCRSY(ILTempGridCRS[j].nGetY()) << "}" << "\t{" << LTempExtCRS.dGetXAt(j) << ", " << LTempExtCRS.dGetYAt(j) << "}" << endl;
   // }
   // LogStream << "==================================" << endl;
   //    // DEBUG CODE ==================================================================================================

   // Create a new coastline object and append to it the vector of coastline objects
   CRWCoast const CoastTmp(this);
   m_VCoast.push_back(CoastTmp);
   int const nCoast = static_cast<int>(m_VCoast.size()) - 1;

   // Now mark the coastline on the grid
   for (int n = 0; n < nCoastSize; n++)
      m_pRasterGrid->Cell(ILTempGridCRS[n].nGetX(), ILTempGridCRS[n].nGetY()).SetAsCoastline(nCoast);

   // Set the coastline (Ext CRS)
   m_VCoast[nCoast].SetCoastlineExtCRS(&LTempExtCRS);

   // Set the coastline (Grid CRS)
   m_VCoast[nCoast].SetCoastlineGridCRS(&ILTempGridCRS);

   // CGeom2DPoint PtLast(DBL_MIN, DBL_MIN);
   // for (int j = 0; j < nCoastSize; j++)
   // {
   //       // Store the smoothed points (in external CRS) in the coast's m_LCoastlineExtCRS object, also append dummy values to the other attribute vectors
   // if (PtLast != &LTempExtCRS[j])        // Avoid duplicate points
   // {
   // m_VCoast[nCoast].AppendPointToCoastlineExtCRS(LTempExtCRS[j].dGetX(), LTempExtCRS[j].dGetY());
   //
   //          // Also store the locations of the corresponding unsmoothed points (in raster grid CRS) in the coast's m_ILCellsMarkedAsCoastline vector
   // m_VCoast[nCoast].AppendCellMarkedAsCoastline(&ILTempGridCRS[j]);
   // }
   //
   // PtLast = LTempExtCRS[j];
   // }

   // Set values for the coast's other attributes: set the coast's handedness, and start and end edges
   m_VCoast[nCoast].SetSeaHandedness(nHandedness);
   m_VCoast[nCoast].SetStartEdge(nStartEdge);
   m_VCoast[nCoast].SetEndEdge(nEndEdge);

   if (m_nLogFileDetail >= LOG_FILE_HIGH_DETAIL)
   {
      LogStream << m_ulIter << ": \tvalid coast " << nCoast << " created, from [" << nStartX << "][" << nStartY << "] to [" << nEndX << "][" << nEndY << "] = {" << dGridCentroidXToExtCRSX(nStartX) << ", " << dGridCentroidYToExtCRSY(nStartY) << "} to {" << dGridCentroidXToExtCRSX(nEndX) << ", " << dGridCentroidYToExtCRSY(nEndY) << "} with " << nCoastSize << " points, handedness = " << (nHandedness == LEFT_HANDED ? "left" : "right") << endl;

      LogStream << m_ulIter << ": \tsmoothed coastline " << nCoast << " runs from {" << LTempExtCRS[0].dGetX() << ", " << LTempExtCRS[0].dGetY() << "} to {" << LTempExtCRS[nCoastSize - 1].dGetX() << ", " << LTempExtCRS[nCoastSize - 1].dGetY() << "} i.e. from the ";

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

   // LogStream << "-----------------" << endl;
   // for (int kk = 0; kk < m_VCoast.back().nGetCoastlineSize(); kk++)
   // LogStream << kk << " [" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX() << "][" << m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY() << "] = {" << dGridCentroidXToExtCRSX(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetX()) << ", " << dGridCentroidYToExtCRSY(m_VCoast.back().pPtiGetCellMarkedAsCoastline(kk)->nGetY()) << "}" << endl;
   // LogStream << "-----------------" << endl;

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
   int const nRet = nTraceAllFloodCoasts();

   if (nRet != RTN_OK)
      return nRet;

   // Have we created any coasts?
   switch (m_nLevel)
   {
      case 0:      // WAVESETUP + SURGE:
      {
         if (m_VFloodWaveSetupSurge.empty())
         {
            cerr << m_ulIter << ": " << ERR << "no flood coastline located: this iteration SWL = " << m_dThisIterSWL << ", maximum DEM top surface elevation = " << m_dThisIterTopElevMax << ", minimum DEM top surface elevation = " << m_dThisIterTopElevMin << endl;
            return RTN_ERR_NO_COAST;
         }

         break;
      }

      case 1:      // WAVESETUP + SURGE + RUNUP:
      {
         if (m_VFloodWaveSetupSurgeRunup.empty())
         {
            cerr << m_ulIter << ": " << ERR << "no flood coastline located: this iteration SWL = " << m_dThisIterSWL << ", maximum DEM top surface elevation = " << m_dThisIterTopElevMax << ", minimum DEM top surface elevation = " << m_dThisIterTopElevMin << endl;
            return RTN_ERR_NO_COAST;
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
         m_pRasterGrid->Cell(nX, nY).UnSetCheckFloodCell();
         m_pRasterGrid->Cell(nX, nY).UnSetInContiguousFlood();
         m_pRasterGrid->Cell(nX, nY).SetAsFloodline(false);
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

      int const nX = m_VEdgeCell[n].nGetX();
      int const nY = m_VEdgeCell[n].nGetY();

      if ((! m_pRasterGrid->Cell(nX, nY).bIsCellFloodCheck()) && (m_pRasterGrid->Cell(nX, nY).bIsInundated()))
      {
         // This edge cell is below SWL and sea depth remains set to zero
         FloodFillLand(nX, nY);
      }
   }

   return RTN_OK;
}
