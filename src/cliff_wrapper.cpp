/*!

   \file cliff.cpp
   \brief CRWCliff routines
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
#include <assert.h>

#include <iostream>
// using std::cout;
// using std::cerr;
// using std::endl;
using std::ios;

#include "cme.h"
#include "cliff_wrapper.h"
#include "coast.h"

//! Constructor with seven parameters and an intialization list
CRWCliffWrap::CRWCliffWrap(CRWCoast* pCoastIn, int const nCoast, int const nPointOnCoast, double const dCellSide, double const dNotchDepthIn, double const dNotchElevIn, double const dAccumWaveEnergyIn)
{
   m_bCliffHasCollapsed = false;

   pCoast = pCoastIn;

   m_nCoast = nCoast;
   m_nPointOnCoast = nPointOnCoast;
   m_nCategory = LF_CAT_CLIFF;

   m_dMaxDepth = dCellSide;
   m_dNotchDepth = dNotchDepthIn;
   m_dNotchBaseElev = dNotchElevIn;
   m_dTotAccumWaveEnergy = dAccumWaveEnergyIn;
   // assert(m_dRemaining >=0);
}

//! Destructor
CRWCliffWrap::~CRWCliffWrap(void)
{
}

//! Returns the value of the cliff collapse switch
bool CRWCliffWrap::bHasCollapsed(void) const
{
   return m_bCliffHasCollapsed;
}

//! Flags the cliff as having collapsed
void CRWCliffWrap::SetCliffCollapsed(void)
{
   m_bCliffHasCollapsed = true;
}

//! Returns the elevation of the base of the erosional notch
double CRWCliffWrap::dGetNotchBaseElev(void) const
{
   return m_dNotchBaseElev;
}

//! Sets the elevation of the base of the erosional notch
void CRWCliffWrap::SetNotchBaseElev(double const dNewElev)
{
   m_dNotchBaseElev = dNewElev;
}

//! Returns the length (in external CRS units) of the cliff's remaining sediment 'behind' the erosional notch
double CRWCliffWrap::dGetRemaining(void) const
{
   return (m_dMaxDepth - m_dNotchDepth);
}

// //! Sets the horizontal depth of the cliff's erosional notch
// void CRWCliff::SetNotchDepth(double const dLenIn)
// {
// m_dNotchDepth = dLenIn;
// }

//! Returns the horizontal depth of the cliff's erosional notch (the 'overhang')
double CRWCliffWrap::dGetNotchDepth(void) const
{
   return m_dNotchDepth;
}

//! Returns true if the horizontal depth of the erosional notch exceeds the critical notch overhang
bool CRWCliffWrap::bReadyToCollapse(double const dThresholdNotchDepth) const
{
   if (m_dNotchDepth >= dThresholdNotchDepth)
      return true;

   else
      return false;
}

//! Increases the XY-plane length (in external CRS units) of the erosional notch, measured inland from the side of the cell that touches the sea
void CRWCliffWrap::DeepenErosionalNotch(double const dLenIn)
{
   m_dNotchDepth += dLenIn;

   // Constrain the notch depth, it cannot be greater than the max notch depth
   m_dNotchDepth = tMin(m_dNotchDepth, m_dMaxDepth);

   // assert((m_dMaxDepth - m_dNotchDepth) >=0);
}

//! Instantiates the pure virtual function in the abstract parent class, so that CRWCliff is not an abstract class
void CRWCliffWrap::Display(void)
{
}
