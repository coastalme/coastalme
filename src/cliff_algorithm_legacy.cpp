/*!

   \file cliff_algorithm_legacy.cpp
   \brief CCliffAlgorithmLegacy implementation
   \details Legacy cliff collapse algorithm wrapper implementation
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License

*/

/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/

#include "cliff_algorithm_legacy.h"

//! Constructor
CCliffAlgorithmLegacy::CCliffAlgorithmLegacy()
{
   m_strAlgorithmName = "Legacy Original Algorithm";
}

//! Initialize the algorithm
bool CCliffAlgorithmLegacy::Initialize(const string& strAlgorithmName)
{
   m_strAlgorithmName = strAlgorithmName;
   return true;
}

//! Process a single timestep
CCliffResults CCliffAlgorithmLegacy::ProcessTimestep(const CCliffData& cliffData)
{
   CCliffResults results;
   
   // TODO: This will be implemented to call the original CoastalME cliff collapse logic
   // For now, return empty results
   results.NewCliffToe = cliffData.CliffToe;
   results.NewCliffTop = cliffData.CliffTop;
   
   return results;
}

//! Reset algorithm state
void CCliffAlgorithmLegacy::Reset()
{
   // Nothing to reset for legacy algorithm
}

//! Get algorithm name
string CCliffAlgorithmLegacy::GetAlgorithmName() const
{
   return m_strAlgorithmName;
}

//! Get algorithm description
string CCliffAlgorithmLegacy::GetAlgorithmDescription() const
{
   return "Legacy wrapper for the original CoastalME cliff collapse algorithm for comparison purposes.";
}