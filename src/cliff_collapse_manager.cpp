/*!

   \file cliff_collapse_manager.cpp
   \brief CCliffCollapseManager implementation
   \details Manager class implementation for cliff collapse algorithms
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

#include "cliff_collapse_manager.h"
#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "cliff.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//! Constructor
CCliffCollapseManager::CCliffCollapseManager(CSimulation* pSimulation)
   : m_pSimulation(pSimulation)
   , m_bCliffCollapseEnabled(false)
{
}

//! Initialize the cliff collapse manager
bool CCliffCollapseManager::Initialize(const string& strAlgorithmName)
{
   if (! m_pSimulation)
      return false;
   
   // Create the cliff algorithm
   m_pCliffAlgorithm = CCliffAlgorithmFactory::CreateAlgorithm(strAlgorithmName);
   if (! m_pCliffAlgorithm)
      return false;
   
   // Initialize the algorithm
   if (! m_pCliffAlgorithm->Initialize(strAlgorithmName))
      return false;
   
   m_strAlgorithmName = strAlgorithmName;
   m_bCliffCollapseEnabled = true;
   
   return true;
}

//! Process cliff collapse for all coasts and points
int CCliffCollapseManager::ProcessAllCliffCollapse()
{
   if (! m_bCliffCollapseEnabled || !m_pCliffAlgorithm)
      return RTN_OK;
   
   // Process cliff collapse for all coasts
   int nCoasts = static_cast<int>(m_pSimulation->m_VCoast.size());
   
   // Use OpenMP to parallelize coast processing
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (int nCoast = 0; nCoast < nCoasts; nCoast++)
   {
      ProcessCoastCliffCollapse(nCoast);
   }
   
   return RTN_OK;
}

//! Process cliff collapse for a single coast
int CCliffCollapseManager::ProcessCoastCliffCollapse(int nCoast)
{
   if (! m_bCliffCollapseEnabled || !m_pCliffAlgorithm)
      return RTN_OK;
   
   // Get the coast object
   if (nCoast >= static_cast<int>(m_pSimulation->m_VCoast.size()))
      return RTN_ERR_NOCOAST;
   
   CRWCoast* pCoast = &(m_pSimulation->m_VCoast[nCoast]);
   
   // Process cliff collapse for all points on this coast
   int nCoastPoints = pCoast->nGetCoastlineSize();
   
   // Use OpenMP to parallelize point processing within each coast
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (int nPoint = 0; nPoint < nCoastPoints; nPoint++)
   {
      ProcessPointCliffCollapse(nCoast, nPoint);
   }
   
   return RTN_OK;
}

//! Process cliff collapse for a single point on a coast
int CCliffCollapseManager::ProcessPointCliffCollapse(int nCoast, int nPointOnCoast)
{
   if (! m_bCliffCollapseEnabled || !m_pCliffAlgorithm)
      return RTN_OK;
   
   // Extract cliff data from simulation
   CCliffData cliffData = ExtractCliffData(nCoast, nPointOnCoast);
   
   // Skip if no cliff data available
   if (bFPIsEqual(cliffData.CliffToe.dGetX(), 0.0, TOLERANCE) && bFPIsEqual(cliffData.CliffToe.dGetY(), 0.0, TOLERANCE))
      return RTN_OK;
   
   // Process the cliff collapse algorithm
   CCliffResults results = m_pCliffAlgorithm->ProcessTimestep(cliffData);
   
   // Apply results back to simulation
   ApplyCliffResults(results, nCoast, nPointOnCoast);
   
   return RTN_OK;
}

//! Reset algorithm state
void CCliffCollapseManager::Reset()
{
   if (m_pCliffAlgorithm)
      m_pCliffAlgorithm->Reset();
}

//! Enable/disable cliff collapse
void CCliffCollapseManager::SetCliffCollapseEnabled(bool bEnabled)
{
   m_bCliffCollapseEnabled = bEnabled;
}

//! Check if cliff collapse is enabled
bool CCliffCollapseManager::IsCliffCollapseEnabled() const
{
   return m_bCliffCollapseEnabled;
}

//! Get current algorithm name
string CCliffCollapseManager::GetAlgorithmName() const
{
   if (m_pCliffAlgorithm)
      return m_pCliffAlgorithm->GetAlgorithmName();

   return m_strAlgorithmName;
}

//! Get algorithm description
string CCliffCollapseManager::GetAlgorithmDescription() const
{
   if (m_pCliffAlgorithm)
      return m_pCliffAlgorithm->GetAlgorithmDescription();
   return "";
}

//! Get list of available algorithms
vector<string> CCliffCollapseManager::GetAvailableAlgorithms() const
{
   return CCliffAlgorithmFactory::GetAvailableAlgorithms();
}

//! Switch to a different algorithm
bool CCliffCollapseManager::SwitchAlgorithm(const string& strNewAlgorithmName)
{
   // Create new algorithm
   auto pNewAlgorithm = CCliffAlgorithmFactory::CreateAlgorithm(strNewAlgorithmName);
   if (! pNewAlgorithm)
      return false;
   
   // Initialize new algorithm
   if (! pNewAlgorithm->Initialize(strNewAlgorithmName))
      return false;
   
   // Replace old algorithm
   m_pCliffAlgorithm = std::move(pNewAlgorithm);
   m_strAlgorithmName = strNewAlgorithmName;
   
   return true;
}

//! Extract cliff data from simulation for a specific coast point
CCliffData CCliffCollapseManager::ExtractCliffData(int nCoast, int nPointOnCoast) const
{
   CCliffData cliffData;
   
   // TODO: This will be implemented to extract actual cliff data from the simulation
   // For now, return empty data
   cliffData.nCoast = nCoast;
   cliffData.nPointOnCoast = nPointOnCoast;
   cliffData.dTimeStep = 1.0; // Default 1 hour timestep
   
   return cliffData;
}

//! Apply cliff algorithm results back to simulation
void CCliffCollapseManager::ApplyCliffResults(const CCliffResults& results, int nCoast, int nPointOnCoast)
{
   // Update cliff geometry
   UpdateCliffGeometry(results, nCoast, nPointOnCoast);
   
   // Apply collapse deposition if collapse occurred
   if (results.bCollapseOccurred)
   {
      ApplyCollapseDeposition(results, nCoast, nPointOnCoast);
   }
}

//! Apply sediment deposition from cliff collapse
void CCliffCollapseManager::ApplyCollapseDeposition(const CCliffResults& results, int nCoast, int nPointOnCoast)
{
   // TODO: Implement sediment deposition logic
   // This will distribute collapsed sediment volumes using the existing talus deposition system
}

//! Update cliff geometry after erosion/collapse
void CCliffCollapseManager::UpdateCliffGeometry(const CCliffResults& results, int nCoast, int nPointOnCoast)
{
   // TODO: Implement cliff geometry updates
   // This will update the cliff toe and top positions in the simulation
}
