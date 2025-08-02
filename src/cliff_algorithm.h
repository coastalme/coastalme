/*!

   \class CCliffAlgorithm
   \brief Abstract base class for cliff collapse algorithms
   \details This class provides the interface for pluggable cliff collapse algorithms in CoastalME.
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file cliff_algorithm.h
   \brief Contains CCliffAlgorithm definitions and data structures

*/

#ifndef CLIFF_ALGORITHM_H
#define CLIFF_ALGORITHM_H
/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/

#include <vector>
#include <string>
#include "2d_point.h"

using std::vector;
using std::string;

//! Data structure containing cliff geometry and environmental conditions for algorithm input
struct CCliffData
{
   //! Cliff toe position (x, y coordinates)
   CGeom2DPoint CliffToe;
   
   //! Cliff top position (x, y coordinates)
   CGeom2DPoint CliffTop;
   
   //! Still water level for current timestep (m)
   double dStillWaterLevel;
   
   //! Significant wave height at breaking (m)
   double dWaveHeight;
   
   //! Wave period (s)
   double dWavePeriod;
   
   //! Wave direction (degrees, 0 = North, clockwise)
   double dWaveDirection;
   
   //! Wave energy flux at breaking (J/m/s)
   double dWaveEnergy;
   
   //! Timestep duration (hours)
   double dTimeStep;
   
   //! Coast index for this cliff section
   int nCoast;
   
   //! Point index on coast for this cliff section
   int nPointOnCoast;
   
   //! Default constructor
   CCliffData() : dStillWaterLevel(0.0), dWaveHeight(0.0), dWavePeriod(0.0), 
                  dWaveDirection(0.0), dWaveEnergy(0.0), dTimeStep(0.0), 
                  nCoast(0), nPointOnCoast(0) {}
};

//! Data structure containing cliff algorithm results
struct CCliffResults
{
   //! Total eroded volume from cliff face (m³)
   double dErosionVolume;
   
   //! Volume of fine sediment eroded (m³)
   double dErosionVolumeFine;
   
   //! Volume of sand eroded (m³)
   double dErosionVolumeSand;
   
   //! Volume of coarse sediment eroded (m³)
   double dErosionVolumeCoarse;
   
   //! Flag indicating if cliff collapse occurred this timestep
   bool bCollapseOccurred;
   
   //! Total collapse volume if collapse occurred (m³)
   double dCollapseVolume;
   
   //! Volume of fine sediment collapsed (m³)
   double dCollapseVolumeFine;
   
   //! Volume of sand collapsed (m³)
   double dCollapseVolumeSand;
   
   //! Volume of coarse sediment collapsed (m³)
   double dCollapseVolumeCoarse;
   
   //! New cliff toe position after erosion/collapse
   CGeom2DPoint NewCliffToe;
   
   //! New cliff top position after erosion/collapse
   CGeom2DPoint NewCliffTop;
   
   //! Default constructor
   CCliffResults() : dErosionVolume(0.0), dErosionVolumeFine(0.0), dErosionVolumeSand(0.0), 
                     dErosionVolumeCoarse(0.0), bCollapseOccurred(false), dCollapseVolume(0.0),
                     dCollapseVolumeFine(0.0), dCollapseVolumeSand(0.0), dCollapseVolumeCoarse(0.0) {}
};

//! Abstract base class for cliff collapse algorithms
class CCliffAlgorithm
{
   public:
      //! Virtual destructor
      virtual ~CCliffAlgorithm() = default;
      
      //! Initialize the algorithm with configuration parameters
      virtual bool Initialize(const string& strAlgorithmName) = 0;
      
      //! Process a single timestep for cliff evolution
      virtual CCliffResults ProcessTimestep(const CCliffData& cliffData) = 0;
      
      //! Reset algorithm state (e.g., for new simulation run)
      virtual void Reset() = 0;
      
      //! Get algorithm name
      virtual string GetAlgorithmName() const = 0;
      
      //! Get algorithm description
      virtual string GetAlgorithmDescription() const = 0;
};

#endif // CLIFF_ALGORITHM_H