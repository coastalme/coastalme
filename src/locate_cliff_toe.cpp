/*!

  \file locate_cliff_toe.cpp
  \brief Locates and traces the cliff toe on the raster grid
  \details TODO 001 A more detailed description of these routines.
  \author David Favis-Mortlock
  \author Andres Payo
  \date 2025
  \copyright GNU General Public License

*/

/* ==============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

   This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave,
Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include <assert.h>
#include <cfloat>

#include <cpl_port.h>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <iomanip>
using std::resetiosflags;
using std::setiosflags;
using std::setprecision;
using std::setw;

#include <sstream>
using std::stringstream;

#include "cme.h"
#include "raster_grid.h"
#include "simulation.h"

/*===============================================================================================================================

 Calculates slope at every cell throughout the grid using finite difference
method

===============================================================================================================================*/
void CSimulation::nCalcSlopeAtAllCells(void) {
  for (int nX = 0; nX < m_nXGridSize; nX++) {
    for (int nY = 0; nY < m_nYGridSize; nY++) {

      // Now lets look at surounding cells
      if (nX > 0 && nX < m_nXGridSize - 1 && nY > 0 && nY < m_nYGridSize - 1) {
        double dElevLeft =
            m_pRasterGrid->m_Cell[nX - 1][nY].dGetOverallTopElev();
        double dElevRight =
            m_pRasterGrid->m_Cell[nX + 1][nY].dGetOverallTopElev();
        double dElevUp = m_pRasterGrid->m_Cell[nX][nY - 1].dGetOverallTopElev();
        double dElevDown =
            m_pRasterGrid->m_Cell[nX][nY + 1].dGetOverallTopElev();
        // Calculate slope using finite difference method
        double dSlopeX = (dElevRight - dElevLeft) / (2.0 * m_dCellSide);
        double dSlopeY = (dElevDown - dElevUp) / (2.0 * m_dCellSide);
        double dSlope = sqrt(dSlopeX * dSlopeX + dSlopeY * dSlopeY);
        m_pRasterGrid->m_Cell[nX][nY].SetSlope(dSlope);
      }
    }
  }
}

/*===============================================================================================================================

Highligts cells with slope greater than threshold

===============================================================================================================================*/
void CSimulation::nLocateCliffCell(void) {
  // TODO: Allow user input of cliff threshold
  double cliff_slope_limit{0.3};
  for (int nX = 0; nX < m_nXGridSize; nX++) {
    for (int nY = 0; nY < m_nYGridSize; nY++) {
      double dSlope = m_pRasterGrid->m_Cell[nX][nY].dGetSlope();
      if (dSlope >= cliff_slope_limit) {
        m_pRasterGrid->m_Cell[nX][nY].SetAsCliff(true);
      }
    }
  }
}

/*===============================================================================================================================

Traces the seaward extent of cliff cells using wall following algorithm

===============================================================================================================================*/
void CSimulation::nTraceSeawardCliffEdge(void) {
  // Clear previous cliff edges
  m_VCliffEdge.clear();

  // Find all possible cliff edge start points by scanning for cliff/non-cliff
  // transitions
  vector<CGeom2DIPoint> V2DIPossibleStartCell;
  vector<bool> VbPossibleStartCellHandedness;
  vector<int> VnSearchDirection;

  // Scan all cells to find cliff toe edges by checking elevation patterns
  for (int nX = 2; nX < m_nXGridSize - 2; nX++) {
    for (int nY = 2; nY < m_nYGridSize - 2; nY++) {
      if (m_pRasterGrid->m_Cell[nX][nY].bIsCliff()) {

        // East direction (check if this is a seaward-facing cliff toe)
        if (!m_pRasterGrid->m_Cell[nX][nY + 1].bIsCliff()) {
          V2DIPossibleStartCell.push_back(CGeom2DIPoint(nX, nY));
          VbPossibleStartCellHandedness.push_back(true);
          VnSearchDirection.push_back(EAST);
        }

        // South direction
        if (!m_pRasterGrid->m_Cell[nX + 1][nY].bIsCliff()) {
          V2DIPossibleStartCell.push_back(CGeom2DIPoint(nX, nY));
          VbPossibleStartCellHandedness.push_back(true);
          VnSearchDirection.push_back(SOUTH);
        }

        // West direction
        if (!m_pRasterGrid->m_Cell[nX][nY - 1].bIsCliff()) {
          V2DIPossibleStartCell.push_back(CGeom2DIPoint(nX, nY));
          VbPossibleStartCellHandedness.push_back(true);
          VnSearchDirection.push_back(WEST);
        }

        // North direction
        if (!m_pRasterGrid->m_Cell[nX - 1][nY].bIsCliff()) {
          V2DIPossibleStartCell.push_back(CGeom2DIPoint(nX, nY));
          VbPossibleStartCellHandedness.push_back(true);
          VnSearchDirection.push_back(NORTH);
        }
      }
    }
  }

  // Create a 2D array to track which cells have been used in cliff edge tracing
  vector<vector<bool>> bUsedInCliffTrace(m_nXGridSize,
                                         vector<bool>(m_nYGridSize, false));

  int nCliffEdgesTraced = 0;

  // For each possible start point, attempt to trace a cliff edge
  for (size_t nStartPoint = 0; nStartPoint < V2DIPossibleStartCell.size();
       nStartPoint++) {
    int nXStart = V2DIPossibleStartCell[nStartPoint].nGetX();
    int nYStart = V2DIPossibleStartCell[nStartPoint].nGetY();

    // Skip if this cell has already been used in another cliff trace
    if (bUsedInCliffTrace[nXStart][nYStart]) {
      continue;
    }

    // Begin cliff edge tracing using wall following algorithm
    vector<CGeom2DIPoint> VCliffEdge;
    bool bHandedness = VbPossibleStartCellHandedness[nStartPoint];
    int nSearchDirection = VnSearchDirection[nStartPoint];

    int nX = nXStart;
    int nY = nYStart;
    int nLength = 0;
    const int nMaxLength = m_nXGridSize * m_nYGridSize; // Safety limit

    do {
      // Add current point to cliff edge
      VCliffEdge.push_back(CGeom2DIPoint(nX, nY));
      bUsedInCliffTrace[nX][nY] = true;
      nLength++;

      // Find next cell using wall following algorithm (right-hand rule)
      int nXSeaward, nYSeaward, nXStraightOn, nYStraightOn;
      int nXAntiSeaward, nYAntiSeaward, nXGoBack, nYGoBack;

      // Calculate candidate positions based on search direction
      switch (nSearchDirection) {
      case NORTH:
        nXSeaward = nX + 1;
        nYSeaward = nY; // Right (East)
        nXStraightOn = nX;
        nYStraightOn = nY - 1; // Straight (North)
        nXAntiSeaward = nX - 1;
        nYAntiSeaward = nY; // Left (West)
        nXGoBack = nX;
        nYGoBack = nY + 1; // Back (South)
        break;
      case EAST:
        nXSeaward = nX;
        nYSeaward = nY + 1; // Right (South)
        nXStraightOn = nX + 1;
        nYStraightOn = nY; // Straight (East)
        nXAntiSeaward = nX;
        nYAntiSeaward = nY - 1; // Left (North)
        nXGoBack = nX - 1;
        nYGoBack = nY; // Back (West)
        break;
      case SOUTH:
        nXSeaward = nX - 1;
        nYSeaward = nY; // Right (West)
        nXStraightOn = nX;
        nYStraightOn = nY + 1; // Straight (South)
        nXAntiSeaward = nX + 1;
        nYAntiSeaward = nY; // Left (East)
        nXGoBack = nX;
        nYGoBack = nY - 1; // Back (North)
        break;
      case WEST:
        nXSeaward = nX;
        nYSeaward = nY - 1; // Right (North)
        nXStraightOn = nX - 1;
        nYStraightOn = nY; // Straight (West)
        nXAntiSeaward = nX;
        nYAntiSeaward = nY + 1; // Left (South)
        nXGoBack = nX + 1;
        nYGoBack = nY; // Back (East)
        break;
      default:
        nXSeaward = nXStraightOn = nXAntiSeaward = nXGoBack = nX;
        nYSeaward = nYStraightOn = nYAntiSeaward = nYGoBack = nY;
        break;
      }

      // Try to move using wall following priority: seaward, straight,
      // anti-seaward, back
      bool bFoundNextCell = false;

      // 1. Try seaward (right turn)
      if (bIsWithinValidGrid(nXSeaward, nYSeaward) &&
          m_pRasterGrid->m_Cell[nXSeaward][nYSeaward].bIsCliff()) {
        nX = nXSeaward;
        nY = nYSeaward;
        // Update search direction (turn right)
        switch (nSearchDirection) {
        case NORTH:
          nSearchDirection = EAST;
          break;
        case EAST:
          nSearchDirection = SOUTH;
          break;
        case SOUTH:
          nSearchDirection = WEST;
          break;
        case WEST:
          nSearchDirection = NORTH;
          break;
        }
        bFoundNextCell = true;
      }
      // 2. Try straight ahead
      else if (bIsWithinValidGrid(nXStraightOn, nYStraightOn) &&
               m_pRasterGrid->m_Cell[nXStraightOn][nYStraightOn].bIsCliff()) {
        nX = nXStraightOn;
        nY = nYStraightOn;
        // Direction stays the same
        bFoundNextCell = true;
      }
      // 3. Try anti-seaward (left turn)
      else if (bIsWithinValidGrid(nXAntiSeaward, nYAntiSeaward) &&
               m_pRasterGrid->m_Cell[nXAntiSeaward][nYAntiSeaward].bIsCliff()) {
        nX = nXAntiSeaward;
        nY = nYAntiSeaward;
        // Update search direction (turn left)
        switch (nSearchDirection) {
        case NORTH:
          nSearchDirection = WEST;
          break;
        case EAST:
          nSearchDirection = NORTH;
          break;
        case SOUTH:
          nSearchDirection = EAST;
          break;
        case WEST:
          nSearchDirection = SOUTH;
          break;
        }
        bFoundNextCell = true;
      }
      // 4. Try going back (U-turn)
      else if (bIsWithinValidGrid(nXGoBack, nYGoBack) &&
               m_pRasterGrid->m_Cell[nXGoBack][nYGoBack].bIsCliff()) {
        nX = nXGoBack;
        nY = nYGoBack;
        // Update search direction (U-turn)
        nSearchDirection = nGetOppositeDirection(nSearchDirection);
        bFoundNextCell = true;
      }

      if (!bFoundNextCell) {
        break; // No valid next cell found, end this cliff edge trace
      }

    } while ((nX != nXStart || nY != nYStart) && nLength < nMaxLength);

    // Only keep cliff edges that have a reasonable length
    if (VCliffEdge.size() > 2) {
      nCliffEdgesTraced++;

      // Convert grid coordinates to external CRS for smoothing
      CGeomLine CliffEdgeExtCRS;
      for (const auto &point : VCliffEdge) {
        CliffEdgeExtCRS.Append(dGridCentroidXToExtCRSX(point.nGetX()), dGridCentroidYToExtCRSY(point.nGetY()));
        bUsedInCliffTrace[point.nGetX()][point.nGetY()] = true;
      }

      // Apply the same smoothing as used for coastlines
      if (m_nCoastSmooth == SMOOTH_RUNNING_MEAN)
        CliffEdgeExtCRS = LSmoothCoastRunningMean(&CliffEdgeExtCRS);
      else if (m_nCoastSmooth == SMOOTH_SAVITZKY_GOLAY)
        CliffEdgeExtCRS = LSmoothCoastSavitzkyGolay(&CliffEdgeExtCRS, 0, 0);

      // Store the smoothed cliff edge in external CRS
      m_VCliffEdge.push_back(CliffEdgeExtCRS);
    }
  }
}

/*===============================================================================================================================

Removes small cliff islands below a given number of cellsusing flood fill
algorithm

===============================================================================================================================*/
void CSimulation::nRemoveSmallCliffIslands(int const dMinCliffCellThreshold) {
  // Create a 2D array to track which cells have been visited during flood fill
  vector<vector<bool>> bVisited(m_nXGridSize,
                                vector<bool>(m_nYGridSize, false));

  // Vector to store cells that belong to small cliff islands to be removed
  vector<pair<int, int>> VSmallIslandCells;

  // Loop through all cells to find unvisited cliff cells
  for (int nX = 0; nX < m_nXGridSize; nX++) {
    for (int nY = 0; nY < m_nYGridSize; nY++) {
      // Check if this is an unvisited cliff cell
      if (!bVisited[nX][nY] && m_pRasterGrid->m_Cell[nX][nY].bIsCliff()) {
        // Found the start of a new cliff region - use flood fill to find all
        // connected cliff cells
        vector<pair<int, int>> VCurrentCliffRegion;

        // Stack for iterative flood fill algorithm
        vector<pair<int, int>> VStack;
        VStack.push_back(std::make_pair(nX, nY));

        // Flood fill to find all connected cliff cells
        while (!VStack.empty()) {
          pair<int, int> currentCell = VStack.back();
          VStack.pop_back();

          int nCurX = currentCell.first;
          int nCurY = currentCell.second;

          // Skip if already visited or out of bounds
          if (nCurX < 0 || nCurX >= m_nXGridSize || nCurY < 0 ||
              nCurY >= m_nYGridSize || bVisited[nCurX][nCurY] ||
              !m_pRasterGrid->m_Cell[nCurX][nCurY].bIsCliff()) {
            continue;
          }

          // Mark as visited and add to current cliff region
          bVisited[nCurX][nCurY] = true;
          VCurrentCliffRegion.push_back(std::make_pair(nCurX, nCurY));

          // Add neighboring cells to stack (8-connectivity: N, NE, E, SE, S,
          // SW, W, NW)
          VStack.push_back(std::make_pair(nCurX - 1, nCurY));     // North
          VStack.push_back(std::make_pair(nCurX - 1, nCurY + 1)); // Northeast
          VStack.push_back(std::make_pair(nCurX, nCurY + 1));     // East
          VStack.push_back(std::make_pair(nCurX + 1, nCurY + 1)); // Southeast
          VStack.push_back(std::make_pair(nCurX + 1, nCurY));     // South
          VStack.push_back(std::make_pair(nCurX + 1, nCurY - 1)); // Southwest
          VStack.push_back(std::make_pair(nCurX, nCurY - 1));     // West
          VStack.push_back(std::make_pair(nCurX - 1, nCurY - 1)); // Northwest
        }

        // Calculate area of this cliff region (number of cells * cell area)
        int dCliffRegionArea = VCurrentCliffRegion.size();

        // If area is below threshold, mark all cells in this region for removal
        if (dCliffRegionArea < dMinCliffCellThreshold) {
          VSmallIslandCells.insert(VSmallIslandCells.end(),
                                   VCurrentCliffRegion.begin(),
                                   VCurrentCliffRegion.end());
        }
      }
    }
  }

  // Remove cliff designation from all small island cells
  for (const auto &cell : VSmallIslandCells) {
    m_pRasterGrid->m_Cell[cell.first][cell.second].SetAsCliff(false);
  }
}

/*===============================================================================================================================

 Locates and traces the cliff toe

===============================================================================================================================*/
int CSimulation::nLocateCliffToe(void) {
  // First step: calculate slope at every cell throughout the grid
  nCalcSlopeAtAllCells();
  nLocateCliffCell();
  nRemoveSmallCliffIslands(10);
  nTraceSeawardCliffEdge();

  // TODO: Additional implementation will be added later

  return RTN_OK;
}
