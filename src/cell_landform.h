/*!

   \class CRWCellLandform
   \brief Real-world class used to represent the landform of a cell
   \details TODO 001 This is a more detailed description of the CCRWCellLandform class.
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License
   \file cell_landform.h
   \brief Contains CRWCellLandform definitions

*/

#ifndef CELL_LANDFORM_H
#define CELL_LANDFORM_H
/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   ===============================================================================================================================*/
class CRWCellLandform
{
private:
   //! Landform category for this cell
   int m_nCategory;

   //! Landform subcategory for this cell
   int m_nSubCategory;

   //! Coast on which this landform sits (if any)
   int m_nCoast;

   //! Point on coast on which this landform sits (if any)
   int m_nPointOnCoast;

   //! Accumulate wave energy for this landform on this cell
   double m_dAccumWaveEnergy;

   //! The m_uLFData will hold landform data: currently, only cliffs are considered
   union
   {
      struct
      {
         //! Currently unused for all landforms except cliffs
         int m_nDummy;
      }
      m_sBeachData;

      struct
      {
         //! Cliff notch base elevation
         double m_dNotchBaseElev;

         //! Cliff notch incised depth
         double m_dNotchDepth;

         //! Cliff notch depth remaining
         double m_dRemaining;
      }
      m_sCliffData;

   }
   m_uLFData;

protected:

public:
   CRWCellLandform();
   ~CRWCellLandform(void);

   void SetLFCategory(int const);
   int nGetLFCategory(void) const;
   void SetLFSubCategory(int const);
   int nGetLFSubCategory(void) const;
   void SetCoast(int const);
   int nGetCoast(void) const;
   void SetPointOnCoast(int const);
   int nGetPointOnCoast(void) const;
   void SetAccumWaveEnergy(double const);
   double dGetAccumWaveEnergy(void) const;
   void SetCliffNotchBaseElev(double const);
   double dGetCliffNotchBaseElev(void) const;
   void SetCliffNotchDepth(double const);
   double dGetCliffNotchDepth(void) const;
   void SetCliffRemaining(double const);
   // double dGetCliffRemaining(void) const;
};
#endif // CELL_LANDFORM_H
