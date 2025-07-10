/*!

   \file do_intervention.cpp
   \brief Checks for new interventions
   \details TODO 001 A more detailed description of these routines.
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License

*/

/* ==============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include "cme.h"
#include "simulation.h"
#include "raster_grid.h"
#include "cell.h"

using std::endl;

//===============================================================================================================================
//! Check for intervention failures and update intervention status
//===============================================================================================================================
int CSimulation::nUpdateIntervention(void)
{
   // Check for intervention failures based on trigger elevations
   int nRet = nCheckInterventionFailures();
   if (nRet != RTN_OK)
      return nRet;

   return RTN_OK;
}

//===============================================================================================================================
//! Check all intervention cells for failure conditions and remove failed interventions
//===============================================================================================================================
int CSimulation::nCheckInterventionFailures(void)
{
   int nInterventionsRemoved = 0;
   vector<CGeom2DIPoint> vFailedInterventions;

   // Scan all cells for interventions that should fail based on influence zone undermining
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         // Check if this cell is an intervention and if its influence zone has failed
         if (bIsInterventionCell(nX, nY) && bCheckInfluenceZoneFailure(nX, nY))
         {
            vFailedInterventions.push_back(CGeom2DIPoint(nX, nY));
         }
      }
   }

   // Group contiguous failed interventions and remove them as units
   vector<bool> vProcessed(vFailedInterventions.size(), false);
   
   for (size_t i = 0; i < vFailedInterventions.size(); i++)
   {
      if (vProcessed[i])
         continue;
         
      // Start a new intervention group
      vector<CGeom2DIPoint> vInterventionGroup;
      vector<size_t> vToProcess;
      vToProcess.push_back(i);
      
      // Find all contiguous intervention cells in this group
      while (!vToProcess.empty())
      {
         size_t nCurrent = vToProcess.back();
         vToProcess.pop_back();
         
         if (vProcessed[nCurrent])
            continue;
            
         vProcessed[nCurrent] = true;
         vInterventionGroup.push_back(vFailedInterventions[nCurrent]);
         
         int nX = vFailedInterventions[nCurrent].nGetX();
         int nY = vFailedInterventions[nCurrent].nGetY();
         
         // Check neighboring cells for more intervention failures (8-direction)
         int nXInc[8] = {-1, -1, -1,  0,  0,  1,  1,  1};
         int nYInc[8] = {-1,  0,  1, -1,  1, -1,  0,  1};
         
         for (int nDirection = 0; nDirection < 8; nDirection++)
         {
            int nXAdj = nX + nXInc[nDirection];
            int nYAdj = nY + nYInc[nDirection];
            
            if (!bIsWithinValidGrid(nXAdj, nYAdj))
               continue;
               
            // Look for this adjacent cell in our failed interventions list
            for (size_t j = 0; j < vFailedInterventions.size(); j++)
            {
               if (!vProcessed[j] && 
                   vFailedInterventions[j].nGetX() == nXAdj && 
                   vFailedInterventions[j].nGetY() == nYAdj)
               {
                  vToProcess.push_back(j);
                  break;
               }
            }
         }
      }
      
      // Remove this entire intervention group
      if (!vInterventionGroup.empty())
      {
         for (const auto& point : vInterventionGroup)
         {
            m_pRasterGrid->m_Cell[point.nGetX()][point.nGetY()].RemoveIntervention();
         }
         nInterventionsRemoved++;
         
         LogStream << m_ulIter << ": intervention group of " << vInterventionGroup.size() 
                   << " cells failed due to influence zone undermining and was removed" << endl;
      }
   }

   if (nInterventionsRemoved > 0)
   {
      LogStream << m_ulIter << ": total of " << nInterventionsRemoved 
                << " intervention groups failed due to influence zone undermining" << endl;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Set trigger depths for all intervention cells
//===============================================================================================================================
void CSimulation::SetInterventionTriggerDepths(double const dTriggerDepth)
{
   int nInterventionsWithTriggers = 0;
   
   for (int nX = 0; nX < m_nXGridSize; nX++)
   {
      for (int nY = 0; nY < m_nYGridSize; nY++)
      {
         if (bIsInterventionCell(nX, nY))
         {
            m_pRasterGrid->m_Cell[nX][nY].SetInterventionTriggerDepth(dTriggerDepth);
            nInterventionsWithTriggers++;
         }
      }
   }
   
   LogStream << "Set trigger depth of " << dTriggerDepth << " m for " 
             << nInterventionsWithTriggers << " intervention cells" << endl;
}

//===============================================================================================================================
//! Check if any cell within the influence zone of an intervention has ground lowering below trigger level
//===============================================================================================================================
bool CSimulation::bCheckInfluenceZoneFailure(int const nIntX, int const nIntY) const
{
   // Get the trigger elevation for this intervention cell
   double dTriggerElev = m_pRasterGrid->m_Cell[nIntX][nIntY].dGetInterventionTriggerElev();
   
   // If no valid trigger elevation is set, no failure
   if (bFPIsEqual(dTriggerElev, DBL_NODATA, TOLERANCE))
      return false;
   
   // Calculate search radius in grid cells (convert from meters to cells)
   int nSearchRadius = static_cast<int>(ceil(m_dInterventionInfluenceDistance / m_dCellSide));
   
   // Search within the influence zone around this intervention cell
   for (int nXOffset = -nSearchRadius; nXOffset <= nSearchRadius; nXOffset++)
   {
      for (int nYOffset = -nSearchRadius; nYOffset <= nSearchRadius; nYOffset++)
      {
         int nCheckX = nIntX + nXOffset;
         int nCheckY = nIntY + nYOffset;
         
         // Skip if outside grid bounds
         if (!bIsWithinValidGrid(nCheckX, nCheckY))
            continue;
         
         // Calculate actual distance from intervention cell to check cell
         double dDistance = sqrt((nXOffset * m_dCellSide) * (nXOffset * m_dCellSide) + 
                               (nYOffset * m_dCellSide) * (nYOffset * m_dCellSide));
         
         // Skip if outside influence distance
         if (dDistance > m_dInterventionInfluenceDistance)
            continue;
         
         // Check if this cell's ground elevation is below the trigger level
         double dCurrentElev = m_pRasterGrid->m_Cell[nCheckX][nCheckY].dGetSedimentTopElev();
         
         if (dCurrentElev < dTriggerElev)
         {
            return true; // Failure detected in influence zone
         }
      }
   }
   
   return false; // No failure detected in influence zone
}
