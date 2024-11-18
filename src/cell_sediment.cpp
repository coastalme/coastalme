/*!
 *
 * \file cell_sediment.cpp
 * \brief CRWCellSediment routines
 * \details TODO 001 A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2024
 * \copyright GNU General Public License
 *
 */

/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
//#include <assert.h>

#include "cme.h"
#include "cell_sediment.h"


//! CRWCellSediment constructor, initialization list sets all internal values to zero
CRWCellSediment::CRWCellSediment(void)
:  m_dFine(0),
   m_dNotchFineLost(0),
   m_dSand(0),
   m_dNotchSandLost(0),
   m_dCoarse(0),
   m_dNotchCoarseLost(0),
   m_dTotSedimentInput(0)
{
}

//! CRWCellSediment copy constructor
CRWCellSediment::CRWCellSediment(CRWCellSediment const& OtherSediment)
{
   m_dFine                = OtherSediment.m_dFine;
   m_dNotchFineLost       = OtherSediment.m_dNotchFineLost;
   m_dSand                = OtherSediment.m_dSand;
   m_dNotchSandLost       = OtherSediment.m_dNotchSandLost;
   m_dCoarse              = OtherSediment.m_dCoarse;
   m_dNotchCoarseLost     = OtherSediment.m_dNotchCoarseLost;
   m_dTotSedimentInput    = OtherSediment.m_dTotSedimentInput;
}

//! Sets one CRWCellSediment object equal to another
CRWCellSediment& CRWCellSediment::operator= (const CRWCellSediment& OtherSediment)
{
   // This copies all fields from one object to another
   m_dFine                = OtherSediment.m_dFine;
   m_dNotchFineLost       = OtherSediment.m_dNotchFineLost;
   m_dSand                = OtherSediment.m_dSand;
   m_dNotchSandLost       = OtherSediment.m_dNotchSandLost;
   m_dCoarse              = OtherSediment.m_dCoarse;
   m_dNotchCoarseLost     = OtherSediment.m_dNotchCoarseLost;
   m_dTotSedimentInput    = OtherSediment.m_dTotSedimentInput;

   return (*this);
}

//! Sets this sediment layer object's fine sediment depth equivalent. Note no checks here to see if new equiv depth is sensible (e.g. non-negative)
void CRWCellSediment::SetFineDepth(double const dNewSedDepth)
{
   m_dFine = dNewSedDepth;
}

//! Returns the fine sediment depth equivalent for this sediment layer object
double CRWCellSediment::dGetFineDepth(void) const
{
   return m_dFine;
}

//! Adds fine sediment (depth equivalent) to this sediment layer object's fine sediment
void CRWCellSediment::AddFineDepth(double const dSedDepthToAdd)
{
   m_dFine += dSedDepthToAdd;
}

//! Sets this sediment layer object's sand sediment depth equivalent. Note no checks here to see if new equiv depth is sensible (e.g. non-negative)
void CRWCellSediment::SetSandDepth(double const dNewSedDepth)
{
   m_dSand = dNewSedDepth;
//    assert(m_dSand >= 0);
}

//! Returns the sand sediment depth equivalent for this sediment layer
double CRWCellSediment::dGetSandDepth(void) const
{
   return m_dSand;
}

//! Adds sand sediment (depth equivalent) to this sediment layer object's sand sediment
void CRWCellSediment::AddSandDepth(double const dSedDepthToAdd)
{
   m_dSand += dSedDepthToAdd;
}

//! Sets this sediment layer object's coarse sediment depth equivalent. Note no checks here to see if new equiv depth is sensible (e.g. non-negative)
void CRWCellSediment::SetCoarseDepth(double const dNewSedDepth)
{
   m_dCoarse = dNewSedDepth;
}

//! Returns the coarse sediment depth equivalent for this sediment layer object
double CRWCellSediment::dGetCoarseDepth(void) const
{
   return m_dCoarse;
}

//! Adds coarse sediment (depth equivalent) to this sediment layer object's coarse sediment
void CRWCellSediment::AddCoarseDepth(double const dSedDepthToAdd)
{
   m_dCoarse += dSedDepthToAdd;
}

//! Sets the depth equivalent of fine sediment lost by notch incision
void CRWCellSediment::SetNotchFineLost(double const dDepthIn)
{
   m_dNotchFineLost = dDepthIn;
}

// //! Increments the depth equivalent of fine sediment lost by notch incision
// void CRWCellSediment::IncrNotchFineLost(double const dDepthIn)
// {
//    m_dNotchFineLost += dDepthIn;
// //    assert(m_dNotchFineLost <= m_dFine);
// }

//! Gets the depth equivalent of fine sediment lost by notch incision
double CRWCellSediment::dGetNotchFineLost(void) const
{
   return m_dNotchFineLost;
}

//! Sets the depth equivalent of sand sediment lost by notch incision
void CRWCellSediment::SetNotchSandLost(double const dDepthIn)
{
   m_dNotchSandLost = dDepthIn;
}

// //! Increments the depth equivalent of sand sediment lost by notch incision
// void CRWCellSediment::IncrNotchSandLost(double const dDepthIn)
// {
//    m_dNotchSandLost += dDepthIn;
// //    assert(m_dNotchSandLost <= m_dSand);
// }

//! Gets the depth equivalent of sand sediment lost by notch incision
double CRWCellSediment::dGetNotchSandLost(void) const
{
   return m_dNotchSandLost;
}

//! Sets the depth equivalent of coarse sediment lost by notch incision
void CRWCellSediment::SetNotchCoarseLost(double const dDepthIn)
{
   m_dNotchCoarseLost = dDepthIn;
}

// //! Increments the depth equivalent of coarse sediment lost by notch incision
// void CRWCellSediment::IncrNotchCoarseLost(double const dDepthIn)
// {
//    m_dNotchCoarseLost += dDepthIn;
// //    assert(m_dNotchCoarseLost <= m_dCoarse);
// }

//! Gets the depth equivalent of coarse sediment lost by notch incision
double CRWCellSediment::dGetNotchCoarseLost(void) const
{
   return m_dNotchCoarseLost;
}

//! Adds to the total (since the start of simulation) depth of sediment input, all size classes
void CRWCellSediment::AddToTotSedimentInputDepth(double const dDepth)
{
   m_dTotSedimentInput += dDepth;
}

//! Returns the total (since the start of simulation) depth of sediment input, all size classes
double CRWCellSediment::dGetTotSedimentInputDepth(void) const
{
   return m_dTotSedimentInput;
}
