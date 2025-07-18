/*!

   \class CCliffCollapseManager
   \brief Manager class for cliff collapse algorithms
   \details This class coordinates between CoastalME simulation and cliff collapse algorithms
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file cliff_collapse_manager.h
   \brief Contains CCliffCollapseManager definitions

*/

#ifndef CLIFF_COLLAPSE_MANAGER_H
#define CLIFF_COLLAPSE_MANAGER_H
/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/

#include <memory>
#include <vector>
#include <string>
#include "cliff_algorithm.h"
#include "cliff_algorithm_factory.h"

using std::unique_ptr;
using std::vector;
using std::string;

// Forward declarations
class CSimulation;
class CRWCoast;

//! Manager class for cliff collapse algorithms
class CCliffCollapseManager
{
   private:
      //! Pointer to the main simulation object
      CSimulation* m_pSimulation;
      
      //! The cliff algorithm instance
      unique_ptr<CCliffAlgorithm> m_pCliffAlgorithm;
      
      //! Algorithm name
      string m_strAlgorithmName;
      
      //! Flag to indicate if cliff collapse is enabled
      bool m_bCliffCollapseEnabled;
      
      //! Extract cliff data from simulation for a specific coast point
      CCliffData ExtractCliffData(int nCoast, int nPointOnCoast) const;
      
      //! Apply cliff algorithm results back to simulation
      void ApplyCliffResults(const CCliffResults& results, int nCoast, int nPointOnCoast);
      
      //! Apply sediment deposition from cliff collapse
      void ApplyCollapseDeposition(const CCliffResults& results, int nCoast, int nPointOnCoast);
      
      //! Update cliff geometry after erosion/collapse
      void UpdateCliffGeometry(const CCliffResults& results, int nCoast, int nPointOnCoast);
      
   public:
      //! Constructor
      CCliffCollapseManager(CSimulation* pSimulation);
      
      //! Destructor
      ~CCliffCollapseManager() = default;
      
      //! Initialize the cliff collapse manager
      bool Initialize(const string& strAlgorithmName);
      
      //! Process cliff collapse for all coasts and points
      int ProcessAllCliffCollapse();
      
      //! Process cliff collapse for a single coast
      int ProcessCoastCliffCollapse(int nCoast);
      
      //! Process cliff collapse for a single point on a coast
      int ProcessPointCliffCollapse(int nCoast, int nPointOnCoast);
      
      //! Reset algorithm state
      void Reset();
      
      //! Enable/disable cliff collapse
      void SetCliffCollapseEnabled(bool bEnabled);
      
      //! Check if cliff collapse is enabled
      bool IsCliffCollapseEnabled() const;
      
      //! Get current algorithm name
      string GetAlgorithmName() const;
      
      //! Get algorithm description
      string GetAlgorithmDescription() const;
      
      //! Get list of available algorithms
      vector<string> GetAvailableAlgorithms() const;
      
      //! Switch to a different algorithm
      bool SwitchAlgorithm(const string& strNewAlgorithmName);
};

#endif // CLIFF_COLLAPSE_MANAGER_H