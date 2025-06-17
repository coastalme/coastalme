/*!

   \file cell_landform.cpp
   \brief CRWCellLandform routines
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
#include "cme.h"
#include "cell_landform.h"

//! Constructor
CRWCellLandform::CRWCellLandform()
   :  m_nCategory(LF_NONE),
      m_nSubCategory(LF_NONE),
      m_nCoast(-1),
      m_nPointOnCoast(-1),
      m_dAccumWaveEnergy(0)
{
}

//! Destructor
CRWCellLandform::~CRWCellLandform(void)
{
}

//! Set the landform category
void CRWCellLandform::SetLFCategory(int const nClassIn)
{
   m_nCategory = nClassIn;
}

//! Get the landform category
int CRWCellLandform::nGetLFCategory(void) const
{
   return m_nCategory;
}

//! Set the both the landform sub-category, and the landform category
void CRWCellLandform::SetLFSubCategory(int const nClassIn)
{
   m_nSubCategory = nClassIn;

   if ((nClassIn == LF_SUBCAT_CLIFF_ON_COASTLINE) || (nClassIn == LF_SUBCAT_CLIFF_INLAND))
      m_nCategory = LF_CAT_CLIFF;

   else if ((nClassIn == LF_SUBCAT_DRIFT_TALUS) || (nClassIn == LF_SUBCAT_DRIFT_BEACH) || (nClassIn == LF_SUBCAT_DRIFT_MIXED))
      m_nCategory = LF_CAT_DRIFT;
}

//! Get the landform sub-category
int CRWCellLandform::nGetLFSubCategory(void) const
{
   return m_nSubCategory;
}

//! Set the coast number
void CRWCellLandform::SetCoast(int const nCoastIn)
{
   m_nCoast = nCoastIn;
}

//! Get the coast number
int CRWCellLandform::nGetCoast(void) const
{
   return m_nCoast;
}

//! Set the number of the point on the coastline
void CRWCellLandform::SetPointOnCoast(int const nPointOnCoastIn)
{
   m_nPointOnCoast = nPointOnCoastIn;
}

//! Set the number of the point on the coastline
int CRWCellLandform::nGetPointOnCoast(void) const
{
   return m_nPointOnCoast;
}

//! Set accumulated wave energy
void CRWCellLandform::SetAccumWaveEnergy(double const dEnergyIn)
{
   m_dAccumWaveEnergy = dEnergyIn;
}

//! Get accumulated wave energy
double CRWCellLandform::dGetAccumWaveEnergy(void) const
{
   return m_dAccumWaveEnergy;
}

//! Set cliff notch base elevation
void CRWCellLandform::SetCliffNotchBaseElev(double const dElevIn)
{
   m_uLFData.m_sCliffData.m_dNotchBaseElev = dElevIn;
}

//! Get cliff notch base elevation
double CRWCellLandform::dGetCliffNotchBaseElev(void) const
{
   return m_uLFData.m_sCliffData.m_dNotchBaseElev;
}

//! Set the cliff notch overhang which remains on this cell
void CRWCellLandform::SetCliffNotchDepth(double const dLenIn)
{
   m_uLFData.m_sCliffData.m_dNotchDepth = dLenIn;
}

//! Get the cliff notch overhang which remains on this cell
double CRWCellLandform::dGetCliffNotchDepth(void) const
{
   return m_uLFData.m_sCliffData.m_dNotchDepth;
}

//! Set the cliff depth remaining on this cell
void CRWCellLandform::SetCliffRemaining(double const dLenIn)
{
   m_uLFData.m_sCliffData.m_dRemaining = dLenIn;
}

// //! Get the cliff depth remaining on this cell
// double CRWCellLandform::dGetCliffRemaining(void) const
// {
//    return m_uLFData.m_sCliffData.m_dRemaining;
// }


