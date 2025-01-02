/*!
 *
 * \file coast_landform.cpp
 * \brief CACoastLandform routines
 * \details TODO 001 A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2025
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
#include "coast_landform.h"

//! Constructor with initialization list
CACoastLandform::CACoastLandform(void)
:  m_nCoast(0),
   m_nPointOnCoast(0),
   m_nCategory(LF_NONE),
   m_dTotAccumWaveEnergy(0),
   pCoast(NULL)
{
}

//! Destructor
CACoastLandform::~CACoastLandform(void)
{
}

//! Get the number of the coast on which this coast landform sits
int CACoastLandform::nGetCoast(void) const
{
   return m_nCoast;
}

//! Get the point on the coast on which this coast landform sits
int CACoastLandform::nGetPointOnCoast(void) const
{
   return m_nPointOnCoast;
}

// void CACoastLandform::SetLandFormCategory(int const nCategoryIn)
// {
//    m_nCategory = nCategoryIn;
// }

//! Get the landform category
int CACoastLandform::nGetLandFormCategory(void) const
{
   return m_nCategory;
}

//! Get the grid co-ordinates of the cell on which this coast landform sits
CGeom2DIPoint* CACoastLandform::pPtiGetCellMarkedAsLF(void) const
{
   return pCoast->pPtiGetCellMarkedAsCoastline(m_nPointOnCoast);
}

// void CACoastLandform::SetTotAccumWaveEnergy(double const dWaveEnergy)
// {
//    m_dTotAccumWaveEnergy = dWaveEnergy;
// }

//! Increment total accumulated wave energy
void CACoastLandform::IncTotAccumWaveEnergy(double const dWaveEnergy)
{
   m_dTotAccumWaveEnergy += dWaveEnergy;
}

//! Get total accumulated wave energy
double CACoastLandform::dGetTotAccumWaveEnergy(void) const
{
   return m_dTotAccumWaveEnergy;
}
