/*!

   \class CCliffAlgorithmSimpleNotch
   \brief Simple notch-based cliff collapse algorithm
   \details This algorithm implements a simple notch erosion and collapse mechanism based on wave energy
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file cliff_algorithm_simple_notch.h
   \brief Contains CCliffAlgorithmSimpleNotch definitions

*/

#ifndef CLIFF_ALGORITHM_SIMPLE_NOTCH_H
#define CLIFF_ALGORITHM_SIMPLE_NOTCH_H
/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/

#include "cliff_algorithm.h"
#include <map>

using std::map;

//! Simple notch-based cliff collapse algorithm
class CCliffAlgorithmSimpleNotch : public CCliffAlgorithm
{
   private:
      //! Algorithm name
      string m_strAlgorithmName;
      
      //! Map to store cliff state for each cliff point (keyed by coast_point string)
      map<string, double> m_mapNotchDepth;
      
      //! Map to store accumulated wave energy for each cliff point
      map<string, double> m_mapAccumulatedWaveEnergy;
      
      //! Cliff erosion resistance parameter (J/m)
      double m_dCliffErosionResistance;
      
      //! Critical notch depth for collapse (m)
      double m_dNotchDepthAtCollapse;
      
      //! Maximum possible notch depth (m)
      double m_dMaxNotchDepth;
      
      //! Generate unique key for cliff point
      string GenerateCliffKey(int nCoast, int nPointOnCoast) const;
      
      //! Calculate notch erosion for this timestep
      double CalculateNotchErosion(double dWaveEnergy, double dTimeStep) const;
      
      //! Calculate collapsed volumes based on cliff geometry
      void CalculateCollapseVolumes(const CCliffData& cliffData, double dCollapseDepth, 
                                    double& dTotalVolume, double& dFineVolume, 
                                    double& dSandVolume, double& dCoarseVolume) const;
      
   public:
      //! Constructor
      CCliffAlgorithmSimpleNotch();
      
      //! Destructor
      ~CCliffAlgorithmSimpleNotch() override = default;
      
      //! Initialize the algorithm
      bool Initialize(const string& strAlgorithmName) override;
      
      //! Process a single timestep
      CCliffResults ProcessTimestep(const CCliffData& cliffData) override;
      
      //! Reset algorithm state
      void Reset() override;
      
      //! Get algorithm name
      string GetAlgorithmName() const override;
      
      //! Get algorithm description
      string GetAlgorithmDescription() const override;
};

#endif // CLIFF_ALGORITHM_SIMPLE_NOTCH_H