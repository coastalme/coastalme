/*!

   \file cliff_algorithm_simple_notch.cpp
   \brief CCliffAlgorithmSimpleNotch implementation
   \details Simple notch-based cliff collapse algorithm implementation
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

#include "cliff_algorithm_simple_notch.h"
#include <sstream>
#include <cmath>
#include <algorithm>

using std::stringstream;
using std::min;
using std::max;

//! Constructor
CCliffAlgorithmSimpleNotch::CCliffAlgorithmSimpleNotch()
{
   // Initialize default parameters (these would typically come from config)
   m_dCliffErosionResistance = 1000.0;     // J/m (default value)
   m_dNotchDepthAtCollapse = 0.5;          // m (default value) 
   m_dMaxNotchDepth = 1.0;                 // m (default value)
   m_strAlgorithmName = "Simple Notch Algorithm";
}

//! Initialize the algorithm
bool CCliffAlgorithmSimpleNotch::Initialize(const string& strAlgorithmName)
{
   m_strAlgorithmName = strAlgorithmName;
   
   // Clear any existing state
   m_mapNotchDepth.clear();
   m_mapAccumulatedWaveEnergy.clear();
   
   return true;
}

//! Process a single timestep
CCliffResults CCliffAlgorithmSimpleNotch::ProcessTimestep(const CCliffData& cliffData)
{
   CCliffResults results;
   
   // Generate unique key for this cliff point
   string strKey = GenerateCliffKey(cliffData.nCoast, cliffData.nPointOnCoast);
   
   // Get current notch depth for this cliff point
   double dCurrentNotchDepth = 0.0;
   if (m_mapNotchDepth.find(strKey) != m_mapNotchDepth.end())
   {
      dCurrentNotchDepth = m_mapNotchDepth[strKey];
   }
   
   // Calculate notch erosion for this timestep
   double dNotchErosion = CalculateNotchErosion(cliffData.dWaveEnergy, cliffData.dTimeStep);
   
   // Update notch depth, constrained by maximum depth
   double dNewNotchDepth = min(dCurrentNotchDepth + dNotchErosion, m_dMaxNotchDepth);
   m_mapNotchDepth[strKey] = dNewNotchDepth;
   
   // Update accumulated wave energy
   double dCurrentEnergy = 0.0;
   if (m_mapAccumulatedWaveEnergy.find(strKey) != m_mapAccumulatedWaveEnergy.end())
   {
      dCurrentEnergy = m_mapAccumulatedWaveEnergy[strKey];
   }
   m_mapAccumulatedWaveEnergy[strKey] = dCurrentEnergy + cliffData.dWaveEnergy;
   
   // Set basic erosion results
   results.dErosionVolume = dNotchErosion * (cliffData.CliffTop.dGetY() - cliffData.CliffToe.dGetY()); // Approximate volume
   
   // Check if collapse should occur
   if (dNewNotchDepth >= m_dNotchDepthAtCollapse)
   {
      results.bCollapseOccurred = true;
      
      // Calculate collapse volumes
      double dCollapseDepth = dNewNotchDepth;
      CalculateCollapseVolumes(cliffData, dCollapseDepth, 
                              results.dCollapseVolume, results.dCollapseVolumeFine,
                              results.dCollapseVolumeSand, results.dCollapseVolumeCoarse);
      
      // Reset notch depth after collapse
      m_mapNotchDepth[strKey] = 0.0;
      
      // Update cliff positions after collapse (simplified - move toe inland)
      double dCliffHeight = cliffData.CliffTop.dGetY() - cliffData.CliffToe.dGetY();
      double dCliffWidth = cliffData.CliffTop.dGetX() - cliffData.CliffToe.dGetX();
      
      // Calculate retreat based on collapse depth
      double dRetreatsRatio = dCollapseDepth / dCliffHeight;
      double dRetreatDistance = dRetreatsRatio * dCliffWidth;
      
      // Move cliff toe inland
      results.NewCliffToe.SetX(cliffData.CliffToe.dGetX() + dRetreatDistance);
      results.NewCliffToe.SetY(cliffData.CliffToe.dGetY());
      
      // Move cliff top inland proportionally
      results.NewCliffTop.SetX(cliffData.CliffTop.dGetX() + dRetreatDistance * 0.5);
      results.NewCliffTop.SetY(cliffData.CliffTop.dGetY());
   }
   else
   {
      // No collapse - maintain current positions
      results.NewCliffToe = cliffData.CliffToe;
      results.NewCliffTop = cliffData.CliffTop;
   }
   
   return results;
}

//! Reset algorithm state
void CCliffAlgorithmSimpleNotch::Reset()
{
   m_mapNotchDepth.clear();
   m_mapAccumulatedWaveEnergy.clear();
}

//! Get algorithm name
string CCliffAlgorithmSimpleNotch::GetAlgorithmName() const
{
   return m_strAlgorithmName;
}

//! Get algorithm description
string CCliffAlgorithmSimpleNotch::GetAlgorithmDescription() const
{
   return "Simple notch-based cliff collapse algorithm that accumulates wave energy to erode notches and triggers collapse when critical depth is reached.";
}

//! Generate unique key for cliff point
string CCliffAlgorithmSimpleNotch::GenerateCliffKey(int nCoast, int nPointOnCoast) const
{
   stringstream ss;
   ss << nCoast << "_" << nPointOnCoast;
   return ss.str();
}

//! Calculate notch erosion for this timestep
double CCliffAlgorithmSimpleNotch::CalculateNotchErosion(double dWaveEnergy, double dTimeStep) const
{
   // Simple erosion model: notch erosion = wave energy / resistance
   // Convert timestep from hours to seconds for energy calculation
   double dTimeStepSeconds = dTimeStep * 3600.0;
   
   // Calculate erosion depth based on wave energy and resistance
   double dErosionDepth = (dWaveEnergy * dTimeStepSeconds) / m_dCliffErosionResistance;
   
   return max(0.0, dErosionDepth);
}

//! Calculate collapsed volumes based on cliff geometry
void CCliffAlgorithmSimpleNotch::CalculateCollapseVolumes(const CCliffData& cliffData, double dCollapseDepth,
                                                          double& dTotalVolume, double& dFineVolume,
                                                          double& dSandVolume, double& dCoarseVolume) const
{
   // Calculate cliff dimensions
   double dCliffHeight = cliffData.CliffTop.dGetY() - cliffData.CliffToe.dGetY();
   double dCliffWidth = cliffData.CliffTop.dGetX() - cliffData.CliffToe.dGetX();
   
   // Calculate total collapse volume (simplified triangular approximation)
   dTotalVolume = 0.5 * dCollapseDepth * dCliffHeight * dCliffWidth;
   
   // Simplified sediment size distribution (would normally come from cliff geology)
   dFineVolume = dTotalVolume * 0.3;    // 30% fine sediment
   dSandVolume = dTotalVolume * 0.5;    // 50% sand
   dCoarseVolume = dTotalVolume * 0.2;  // 20% coarse sediment
}