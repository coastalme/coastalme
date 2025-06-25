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
  double cliff_slope_limit{0.4};
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

Highligts cells with slope greater than threshold

===============================================================================================================================*/
void CSimulation::nTraceSeawardCliffEdge(void) {
  // TODO: Additional implementation will be added later
  ;
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
  nRemoveSmallCliffIslands(5);
  nTraceSeawardCliffEdge();

  // TODO: Additional implementation will be added later

  return RTN_OK;
}
