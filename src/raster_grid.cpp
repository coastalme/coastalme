/*!

   \file raster_grid.cpp
   \brief CGeomRasterGrid routines
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
#include <cstdio>

#include <new>
using std::bad_alloc;

#include "cme.h"
#include "cell.h"
#include "simulation.h"
#include "raster_grid.h"

CGeomRasterGrid* CGeomCell::m_pGrid = NULL; // Initialise m_pGrid, the static member of CGeomCell

//! Constructor
CGeomRasterGrid::CGeomRasterGrid(CSimulation* pSimIn)
    : m_dD50Fine(0),
      m_dD50Sand(0),
      m_dD50Coarse(0),
      m_pSim(pSimIn),
      m_nXSize(0),
      m_nYSize(0),
      m_CellData(NULL)
{
}

//! Destructor
CGeomRasterGrid::~CGeomRasterGrid(void)
{
   // Free the single contiguous cell array
   delete[] m_CellData;
}

//! Returns a pointer to the simulation object
CSimulation* CGeomRasterGrid::pGetSim(void)
{
   return m_pSim;
}

//! Creates a single contiguous CGeomCell array for optimal cache performance
int CGeomRasterGrid::nCreateGrid(void)
{
   // Store grid dimensions
   m_nXSize = m_pSim->nGetGridXMax();
   m_nYSize = m_pSim->nGetGridYMax();

   // Calculate total cells needed
   int const nTotalCells = m_nXSize * m_nYSize;

   // TODO 038 Do better error handling if insufficient memory
   try
   {
      // Single contiguous allocation - much faster than array-of-arrays
      // Memory layout: all cells in a single block, indexed as [nX * m_nYSize + nY]
      m_CellData = new CGeomCell[nTotalCells];
   }

   catch (bad_alloc&)
   {
      // Uh-oh, not enough memory
      return RTN_ERR_MEMALLOC;
   }

   // Initialize the CGeomCell shared pointer to the CGeomRasterGrid object
   CGeomCell::m_pGrid = this;

   return RTN_OK;
}
