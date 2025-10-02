/*!

   \class CGeomRasterGrid
   \brief Geometry cass used to represent the raster grid of cell objects
   \details TODO 001 This is a more detailed description of the CGeomRasterGrid class.
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License
   \file raster_grid.h
   \brief Contains CGeomRasterGrid definitions

*/

#ifndef RASTERGRID_H
#define RASTERGRID_H
/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
#include <vector>

class CGeomCell;        // Forward declaration
class CSimulation;      // Ditto

class CGeomRasterGrid
{
   //! The CSimulation class is a friend of the CGeomRasterGrid class
   friend class CSimulation;

   //! The CGeomProfile class is a friend of the CGeomRasterGrid class
   friend class CGeomProfile;

   private:
   //! The d50 of fine-sized  sediment
   double m_dD50Fine;

   //! The d50 of sand-sized sediment
   double m_dD50Sand;

   //! The d50 of coarse-sized sediment
   double m_dD50Coarse;

   //! A pointer to the CSimulation object
   CSimulation* m_pSim;

   //! The 1D array of m_Cell objects stored as a flat vector for optimal cache performance
   std::vector<CGeomCell> m_Cell;

   //! Grid dimensions
   int m_nXSize;
   int m_nYSize;

   protected:
   public:
   explicit CGeomRasterGrid(CSimulation*);
   ~CGeomRasterGrid(void);

   CSimulation* pGetSim(void);
   // CGeomCell* pGetCell(int const, int const);
   int nCreateGrid(void);

   //! Inline accessor for backward compatibility - compiler will optimize to zero overhead
   inline CGeomCell& Cell(int const nX, int const nY)
   {
      return m_Cell[nY * m_nXSize + nX];
   }

   //! Const accessor for read-only access
   inline CGeomCell const& Cell(int const nX, int const nY) const
   {
      return m_Cell[nY * m_nXSize + nX];
   }

   //! Direct pointer access for hot loops
   inline CGeomCell* CellData(void)
   {
      return m_Cell.data();
   }

   //! Calculate 1D index from 2D coordinates
   inline int nGetIndex(int const nX, int const nY) const
   {
      return nY * m_nXSize + nX;
   }
};
#endif      // RASTERGRID_H
