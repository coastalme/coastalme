/*!
 *
 * \file raster_grid.cpp
 * \brief CGeomRasterGrid routines
 * \details TODO A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2017
 * \copyright GNU General Public License
 *
 */

/*===============================================================================================================================

 This file is part of CoastalME, the Coastal Modelling Environment.

 CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include "cme.h"
#include "raster_grid.h"


CGeomRasterGrid* CGeomCell::m_pGrid = NULL;          // Initialise m_pGrid, the static member of CGeomCell


CGeomRasterGrid::CGeomRasterGrid(CSimulation* pSimIn)
: m_dD50Fine(0),
  m_dD50Sand(0),
  m_dD50Coarse(0),
  m_pSim(pSimIn)
{
}


CGeomRasterGrid::~CGeomRasterGrid(void)
{
}


CSimulation* CGeomRasterGrid::pGetSim(void)
{
   return m_pSim;
}


// CGeomCell* CGeomRasterGrid::pGetCell(int const nX, int const nY)
// {
//    return &m_Cell[nX][nY];
// }


int CGeomRasterGrid::nCreateGrid(void)
{
   // Create the 2D vector CGeomCell array
   int
      nXMax = m_pSim->nGetGridXMax(),
      nYMax = m_pSim->nGetGridYMax();

   // TODO Check if we don't have enough memory, if so return RTN_ERR_MEMALLOC
   m_Cell.resize(nXMax);
   for (int nX = 0; nX < nXMax; nX++)
      m_Cell[nX].resize(nYMax);

   // Initialize the CGeomCell shared pointer to the CGeomRasterGrid object
   CGeomCell::m_pGrid = this;

   return RTN_OK;
}

