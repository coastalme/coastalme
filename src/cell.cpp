/*!
 *
 * \file cell.cpp
 * \brief CGeomCell routines
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
#include "cell.h"

//! Constructor with initialization list
CGeomCell::CGeomCell()
    : m_bInContiguousSea(false),
      m_bInContiguousFlood(false),        // TODO 007 What is this?
      m_bIsInActiveZone(false),
      m_bCoastline(false),
      m_bFloodLine(false),
      m_bWaveFlood(false),
      // m_bCheckCell(false),             // TODO 007 What is this?
      m_bCheckFloodCell(false),           // TODO 007 What is this?
      m_bShadowBoundary(false),
      m_bPossibleCoastStartCell(false),
      m_bPossibleFloodStartCell(false),   // TODO 007 What is this?
      m_nBoundingBoxEdge(NO_DIRECTION),
      m_nPolygonID(INT_NODATA),
      m_nCoastlineNormal(INT_NODATA),
      m_nShadowZoneNumber(0),
      m_nDownDriftZoneNumber(0),
      m_dLocalConsSlope(0),
      m_dBasementElevation(0),
      m_dSeaDepth(0),
      m_dTotSeaDepth(0),
      m_dWaveHeight(0),
      m_dTotWaveHeight(0),
      m_dWaveAngle(DBL_NODATA),
      m_dTotWaveAngle(DBL_NODATA),
      m_dDeepWaterWaveHeight(DBL_NODATA),
      m_dDeepWaterWaveAngle(DBL_NODATA),
      m_dBeachProtectionFactor(DBL_NODATA),
      m_dSuspendedSediment(0),
      m_dTotSuspendedSediment(0),
      m_dPotentialPlatformErosionThisIter(0),
      m_dTotPotentialPlatformErosion(0),
      m_dActualPlatformErosionThisIter(0),
      m_dTotActualPlatformErosion(0),
      m_dCliffCollapseFineThisIter(0),
      m_dCliffCollapseSandThisIter(0),
      m_dCliffCollapseCoarseThisIter(0),
      m_dTotFineCliffCollapse(0),
      m_dTotSandCliffCollapse(0),
      m_dTotCoarseCliffCollapse(0),
      m_dTalusSandDepositionThisIter(0),
      m_dTotTalusSandDeposition(0),
      m_dTalusCoarseDepositionThisIter(0),
      m_dTotTalusCoarseDeposition(0),
      m_dPotentialBeachErosionThisIter(0),
      m_dTotPotentialBeachErosion(0),
      m_dActualBeachErosionThisIter(0),
      m_dTotActualBeachErosion(0),
      m_dBeachDepositionThisIter(0),
      m_dTotBeachDeposition(0),
      m_dUnconsD50(0),
      m_dInterventionHeight(0)
{
   m_Landform.SetLFCategory(LF_NONE);
}

//! Destructor
CGeomCell::~CGeomCell(void)
{
}

//! Set the edge number if this cell is an edge bounding-box cell
void CGeomCell::SetBoundingBoxEdge(int const nDirection)
{
   m_nBoundingBoxEdge = nDirection;
}

//! Returns the number of the bounding-box edge, or NO_DIRECTION if it is not
int CGeomCell::nGetBoundingBoxEdge(void) const
{
   return m_nBoundingBoxEdge;
}

//! Is this an edge bounding-box cell?
bool CGeomCell::bIsBoundingBoxEdge(void) const
{
   return (m_nBoundingBoxEdge != NO_DIRECTION);
}

//! Set this cell as a sea cell
void CGeomCell::SetInContiguousSea(void)
{
   m_bInContiguousSea = true;
}

//! Is this a sea cell?
bool CGeomCell::bIsInContiguousSea(void) const
{
   return m_bInContiguousSea;
}

//! TODO 007 What do this do? Does it duplicate SetInContiguousSea()?
void CGeomCell::SetInContiguousFlood(void)
{
   m_bInContiguousFlood = true;
}

//! TODO 007 What does this do? Is it just the inverse of SetInContiguousSea()?
void CGeomCell::UnSetInContiguousFlood(void)
{
   m_bInContiguousFlood = false;
}

//! TODO 007 What does this do? Set this cell as flood by setup surger
void CGeomCell::SetFloodBySetupSurge(void)
{
   m_bFloodBySetupSurge = true;
}

//! TODO 007 What does this do? Is this cell flood by setup surge?
bool CGeomCell::bIsFloodBySetupSurge(void) const
{
   return m_bFloodBySetupSurge;
}

//! TODO 007 What does this do? Set this cell as flood by setup surge runup
void CGeomCell::SetFloodBySetupSurgeRunup(void)
{
   m_bFloodBySetupSurgeRunup = true;
}

//! TODO 007 What does this do? Is this cell flood by setup surge runup?
bool CGeomCell::bIsFloodBySetupSurgeRunup(void) const
{
   return m_bFloodBySetupSurgeRunup;
}

//! TODO 007 What does this do? Does it just duplicate bIsInContiguousSea()?
bool CGeomCell::bIsInContiguousFlood(void) const
{
   return m_bInContiguousFlood;
}

//! Sets a flag to show whether this cell is in the active zone
void CGeomCell::SetInActiveZone(bool const bFlag)
{
   m_bIsInActiveZone = bFlag;
}

//! Returns a flag which shows whether this cell is in the active zone
bool CGeomCell::bIsInActiveZone(void) const
{
   return m_bIsInActiveZone;
}

//! Sets a flag to show that this cell is a shadow zone boundary
void CGeomCell::SetShadowZoneBoundary(void)
{
   m_bShadowBoundary = true;
}

//! Returns a flag which shows whether this cell is a shadow zone boundary
bool CGeomCell::bIsShadowZoneBoundary(void) const
{
   return m_bShadowBoundary;
}

//! Sets a flag to show that this cell has been flagged as a possible start- or end-point for a coastline
void CGeomCell::SetPossibleCoastStartCell(void)
{
   m_bPossibleCoastStartCell = true;
}

//! Returns a flag which shows whether this cell has been flagged as a possible start- or end-point for a coastline
bool CGeomCell::bIsPossibleCoastStartCell(void) const
{
   return m_bPossibleCoastStartCell;
}

//! TODO 007 What is this for? Sets a flag to show that this cell has been flagged as a possible start-point for a coastline
void CGeomCell::SetPossibleFloodStartCell(void)
{
   m_bPossibleFloodStartCell = true;
}

// //! TODO 007 What is this for? Returns a flag which shows whether this cell has been flagged as a possible start- or end-point for a coastline
// bool CGeomCell::bIsPossibleFloodStartCell(void) const
// {
//    return m_bPossibleFloodStartCell;
// }

//! Returns true if this cell has had potential erosion this timestep
bool CGeomCell::bPotentialPlatformErosion(void) const
{
   return (m_dPotentialPlatformErosionThisIter > 0);
}

// bool CGeomCell::bActualPlatformErosion(void) const
// {
//    return (m_dActualPlatformErosionThisIter > 0);
// }

//! Marks this cell as 'under' a coastline
void CGeomCell::SetAsCoastline(bool const bNewFlag)
{
   m_bCoastline = bNewFlag;
}

//! Returns true if the cell is 'under' a coastline
bool CGeomCell::bIsCoastline(void) const
{
   return m_bCoastline;
}

//! Marks this cell is flood line
void CGeomCell::SetAsFloodLine(bool const bNewFlag)
{
   m_bFloodLine = bNewFlag;
}

//! Returns true if the cell is flood line
bool CGeomCell::bIsFloodLine(void) const
{
   return m_bFloodLine;
}

//! Marks this cell as 'under' a coastline-normal profile
void CGeomCell::SetProfile(int const nNormal)
{
   m_nCoastlineNormal = nNormal;
}

//! If this cell is 'under' a coastline-normal profile, returns the number of the profile. Otherwise it returns INT_NODATA
int CGeomCell::nGetProfile(void) const
{
   return m_nCoastlineNormal;
}

//! Returns true if this cell is 'under' a coastline normal
bool CGeomCell::bIsProfile(void) const
{
   if (m_nCoastlineNormal == INT_NODATA)
      return false;

   return true;
}

//! Sets the global ID number of the polygon which 'contains' this cell
void CGeomCell::SetPolygonID(int const nPolyID)
{
   m_nPolygonID = nPolyID;
}

//! Returns the global ID number of the polygon which 'contains' this cell (returns INT_NODATA if the cell is not 'in' a polygon)
int CGeomCell::nGetPolygonID(void) const
{
   return m_nPolygonID;
}

//! Set the number of the shadow zone that this cell is in
void CGeomCell::SetShadowZoneNumber(int const nCode)
{
   m_nShadowZoneNumber = nCode;
}

//! Gets the number of the shadow zone that this cell is in
int CGeomCell::nGetShadowZoneNumber(void) const
{
   return m_nShadowZoneNumber;
}

//! Returns true if this cell is in the shadow zone with number given by the parameter, false otherwise
bool CGeomCell::bIsinThisShadowZone(int const nZone) const
{
   if (m_nShadowZoneNumber == nZone)
      return true;

   return false;
}

//! Returns true if this cell is in any shadow zone, false otherwise
bool CGeomCell::bIsinAnyShadowZone(void) const
{
   if (m_nShadowZoneNumber != 0)
      return true;

   return false;
}

// Set this cell as flooded by swl + surge + setup + runup
// void CGeomCell::SetWaveFlood(void)
// {
//    m_bWaveFlood = true;
// }

// void CGeomCell::SetWaveSetup(int const dWaveSetup)
// {
//    m_dWaveSetup = dWaveSetup;
// }

// void CGeomCell::SetStormSurge(int const dStormSurge)
// {
//    m_dStormSurge = dStormSurge;
// }

// void CGeomCell::SetRunUp(int const dRunUp)
// {
//    m_dRunUp = dRunUp;
// }

// void CGeomCell::SetTotLevel(void) const
// {
//    m_dTotLevel = m_dWaveSetup + m_dStormSurge + m_dRunUp + m_pGrid->pGetSim()->CSimulation::dGetThisIterSWL();
// }

// int CGeomCell::nGetTotLevel(void) const
// {
//    return m_dTotLevel;
// }

//! Returns true if the top elevation of this cell (sediment plus any intervention) is less than this iteration's total water level
bool CGeomCell::bIsElevLessThanWaterLevel(void) const
{
   return ((m_VdAllHorizonTopElev.back() + m_dInterventionHeight) < (m_pGrid->pGetSim()->dGetThisIterTotWaterLevel() + m_pGrid->pGetSim()->dGetThisIterSWL()));
}

// //! Set this cell as checked TODO What is this used for?
// void CGeomCell::SetCheckCell(void)
// {
//    m_bCheckCell = true;
// }

// //! Returns true if this cell is checked, false otherwise TODO 007 What is this used for?
// bool CGeomCell::bIsCellCheck(void) const
// {
//    return m_bCheckCell;
// }

//! Set this cell as checked (flood switch)
void CGeomCell::SetCheckFloodCell(void)
{
   m_bCheckFloodCell = true;
}

//! Set the cell as not checked (flood switch)
void CGeomCell::UnSetCheckFloodCell(void)
{
   m_bCheckFloodCell = false;
}

//! Returns true if this cell is checked, false otherwise (flood switch)
bool CGeomCell::bIsCellFloodCheck(void) const
{
   return m_bCheckFloodCell;
}

//! Sets the down drift zone number
void CGeomCell::SetDownDriftZoneNumber(int const nCode)
{
   m_nDownDriftZoneNumber = nCode;
}

//! Gets the down drift zone number
int CGeomCell::nGetDownDriftZoneNumber(void) const
{
   return m_nDownDriftZoneNumber;
}

//! Returns a pointer to this cell's CRWCellLandform object
CRWCellLandform *CGeomCell::pGetLandform(void)
{
   return &m_Landform;
}

//! Sets the local slope of the consolidated sediment only
void CGeomCell::SetLocalConsSlope(double const dNewSlope)
{
   m_dLocalConsSlope = dNewSlope;
}

//! Returns the local slope of the consolidated sediment only
double CGeomCell::dGetLocalConsSlope(void) const
{
   return m_dLocalConsSlope;
}

//! Sets this cell's basement elevation
void CGeomCell::SetBasementElev(double const dNewElev)
{
   m_dBasementElevation = dNewElev;
}

//! Returns this cell's basement elevation
double CGeomCell::dGetBasementElev(void) const
{
   return (m_dBasementElevation);
}

//! Returns true if this cells's basement data is NODATA, is needed for irregularly-shaped DEMs
bool CGeomCell::bBasementElevIsMissingValue(void) const
{
   if (bFPIsEqual(m_dBasementElevation, m_pGrid->pGetSim()->CSimulation::dGetMissingValue(), TOLERANCE))
      return true;

   return false;
}

//! Returns the depth of seawater on this cell
double CGeomCell::dGetSeaDepth(void) const
{
   return (m_dSeaDepth);
}

//! Returns the total depth of seawater on this cell
double CGeomCell::dGetTotSeaDepth(void) const
{
   return (m_dTotSeaDepth);
}

//! Sets this cell's suspended sediment depth equivalent, it also increments the running total of suspended sediment depth equivalent
void CGeomCell::SetSuspendedSediment(double const dNewSedDepth)
{
   // Note no checks here to see if new equiv depth is sensible (e.g. non-negative)
   m_dSuspendedSediment = dNewSedDepth;
   m_dTotSuspendedSediment += dNewSedDepth;
}

//! Adds to this cell's suspended sediment depth equivalent, it also increments the running total of suspended sediment depth equivalent
void CGeomCell::AddSuspendedSediment(double const dIncSedDepth)
{
   // Note no checks here to see if increment equiv depth is sensible (e.g. non-negative)
   m_dSuspendedSediment += dIncSedDepth;
   m_dTotSuspendedSediment += dIncSedDepth;
}

//! Returns the suspended sediment depth equivalent on this cell
double CGeomCell::dGetSuspendedSediment(void) const
{
   return (m_dSuspendedSediment);
}

//! Returns the total suspended sediment depth equivalent on this cell
double CGeomCell::dGetTotSuspendedSediment(void) const
{
   return (m_dTotSuspendedSediment);
}

//! Returns the index of the topmost sediment layer (layer 0 being the one just above basement) with non-zero thickness. If there is no such layer, it returns NO_NONZERO_THICKNESS_LAYERS
int CGeomCell::nGetTopNonZeroLayerAboveBasement(void) const
{
   if (m_VLayerAboveBasement.empty())
      return INT_NODATA;

   int nTop = static_cast<int>(m_VLayerAboveBasement.size()) - 1;
   while (m_VLayerAboveBasement[nTop].dGetTotalThickness() <= 0)
   {
      if (--nTop < 0)
         return NO_NONZERO_THICKNESS_LAYERS;
   }

   return nTop;
}

//! Returns the index of the topmost sediment layer (layer 0 being the one just above basement), which could have zero thickness
int CGeomCell::nGetTopLayerAboveBasement(void) const
{
   if (m_VLayerAboveBasement.empty())
      return INT_NODATA;

   return static_cast<int>(m_VLayerAboveBasement.size()) - 1;
}

//! Returns the elevation of the top of the consolidated sediment only, for a given layer (layer 0 being the one just above basement)
double CGeomCell::dGetConsSedTopForLayerAboveBasement(int const nLayer) const
{
   // Note no check to see if nLayer < m_VLayerAboveBasement.size()
   double dTopElev = m_dBasementElevation;

   for (int n = 0; n < nLayer; n++)
   {
      dTopElev += m_VLayerAboveBasement[n].dGetUnconsolidatedThickness();
      dTopElev += m_VLayerAboveBasement[n].dGetConsolidatedThickness();
   }

   dTopElev += m_VLayerAboveBasement[nLayer].dGetConsolidatedThickness();

   return dTopElev;
}

//! Return a reference to the Nth sediment layer (layer 0 being just above basement)
CRWCellLayer *CGeomCell::pGetLayerAboveBasement(int const nLayer)
{
   // TODO 055 No check that nLayer < size()
   return &m_VLayerAboveBasement[nLayer];
}

// //! Returns the volume-equivalent elevation of the sediment's top surface for this cell (i.e. if there is a cliff notch, then lower the elevation by the notch's volume)
// double CGeomCell::dGetVolEquivSedTopElev(void) const
// {
//    double dTopElev = m_dBasementElevation;
//    for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
//    {
//       dTopElev += (m_VLayerAboveBasement[n].dGetUnconsolidatedThickness() - m_VLayerAboveBasement[n].dGetNotchUnconsolidatedLost());
//       dTopElev += (m_VLayerAboveBasement[n].dGetConsolidatedThickness() - m_VLayerAboveBasement[n].dGetNotchConsolidatedLost());
//    }
//
//    return dTopElev;
// }

//! Returns the true elevation of the sediment's top surface for this cell (if there is a cliff notch, ignore the missing volume)
double CGeomCell::dGetSedimentTopElev(void) const
{
   return m_VdAllHorizonTopElev.back();
}

//! Returns the true elevation of the sediment's top surface for this cell (if there is a cliff notch, ignore the missing volume) plus the height of any intervention
double CGeomCell::dGetSedimentPlusInterventionTopElev(void) const
{
   return m_VdAllHorizonTopElev.back() + m_dInterventionHeight;
}

//! Returns the highest elevation of the cell, which is either the sediment top elevation plus intervention height, or the sea surface elevation
double CGeomCell::dGetOverallTopElev(void) const
{
   return m_VdAllHorizonTopElev.back() + m_dInterventionHeight + m_dSeaDepth;
}

//! Returns true if the elevation of the sediment top surface for this cell (plus any intervention) is less than the grid's this-timestep still water elevation
bool CGeomCell::bIsInundated(void) const
{
   return ((m_VdAllHorizonTopElev.back() + m_dInterventionHeight) < m_pGrid->pGetSim()->CSimulation::dGetThisIterSWL());
}

//! Returns the sea surface elevation at current iteration
double CGeomCell::dGetThisIterSWL(void) const
{
   return m_pGrid->pGetSim()->CSimulation::dGetThisIterSWL();
}

//! Returns the total water level at current iteration
double CGeomCell::dGetThisIterTotWaterLevel(void) const
{
   return m_pGrid->pGetSim()->CSimulation::dGetThisIterTotWaterLevel();
}

// //! Returns true if the elevation of the sediment top surface for this cell is greater than or equal to the grid's this-timestep still water elevation. Also returns true if the cell has unconsolidated sediment on it and the elevation of the sediment top surface, minus a tolerance value, is less than the grid's this-timestep still water elevation
// bool CGeomCell::bIsSeaIncBeach(void) const
// {
//    if (m_bInContiguousSea)
//       // Sea
//       return true;
//
//    double
//        dWaterLevel = m_pGrid->pGetSim()->CSimulation::dGetThisIterSWL(),
//        dSedTop = m_VdAllHorizonTopElev.back();
//
//    // Beach
//    if ((m_VLayerAboveBasement.back().dGetUnconsolidatedThickness() > 0) && ((dSedTop - m_pGrid->pGetSim()->CSimulation::dGetMaxBeachElevAboveSWL()) < dWaterLevel))
//       return true;
//
//    return false;
// }

//! Returns the total thickness of fine consolidated sediment on this cell, minus the depth-equivalent of any cliff notch
double CGeomCell::dGetTotConsFineThickConsiderNotch(void) const
{
   double dTotThick = 0;
   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
   {
      CRWCellLayer m_Layer = m_VLayerAboveBasement[n];
      double dLayerThick = m_Layer.dGetFineConsolidatedThickness();
      double dNotchEquiv = m_Layer.pGetConsolidatedSediment()->dGetNotchFineLost();
      
      dTotThick += (dLayerThick - dNotchEquiv);
   }

   return dTotThick;   
}

//! Returns the total thickness of fine unconsolidated sediment on this cell
double CGeomCell::dGetTotUnconsFine(void) const
{
   double dTotThick = 0;
   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotThick += m_VLayerAboveBasement[n].dGetFineUnconsolidatedThickness();

   return dTotThick;   
}

//! Returns the total thickness of sand-sized consolidated sediment on this cell, minus the depth-equivalent of any cliff notch
double CGeomCell::dGetTotConsSandThickConsiderNotch(void) const
{
   double dTotThick = 0;
   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
   {
      CRWCellLayer m_Layer = m_VLayerAboveBasement[n];
      double dLayerThick = m_Layer.dGetSandConsolidatedThickness();
      double dNotchEquiv = m_Layer.pGetConsolidatedSediment()->dGetNotchSandLost();
      
      dTotThick += (dLayerThick - dNotchEquiv);
   }

   return dTotThick;   
}

//! Returns the total thickness of sand-sized unconsolidated sediment on this cell
double CGeomCell::dGetTotUnconsSand(void) const
{
   double dTotThick = 0;
   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotThick += m_VLayerAboveBasement[n].dGetSandUnconsolidatedThickness();

   return dTotThick;   
}

//! Returns the total thickness of coarse consolidated sediment on this cell, minus the depth-equivalent of any cliff notch
double CGeomCell::dGetTotConsCoarseThickConsiderNotch(void) const
{
   double dTotThick = 0;
   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
   {
      CRWCellLayer m_Layer = m_VLayerAboveBasement[n];
      double dLayerThick = m_Layer.dGetCoarseConsolidatedThickness();
      double dNotchEquiv = m_Layer.pGetConsolidatedSediment()->dGetNotchCoarseLost();
      
      dTotThick += (dLayerThick - dNotchEquiv);
   }

   return dTotThick;   
}

//! Returns the total thickness of coarse unconsolidated sediment on this cell
double CGeomCell::dGetTotUnconsCoarse(void) const
{
double dTotThick = 0;
   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotThick += m_VLayerAboveBasement[n].dGetCoarseUnconsolidatedThickness();

   return dTotThick;   
}

//! Returns the total thickness of consolidated sediment (all size classes) on this cell
double CGeomCell::dGetTotConsThickness(void) const
{
   double dTotThick = 0;
   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotThick += m_VLayerAboveBasement[n].dGetConsolidatedThickness();

   return dTotThick;
}

//! Returns the total thickness of unconsolidated sediment (all size classes) on this cell
double CGeomCell::dGetTotUnconsThickness(void) const
{
   double dTotThick = 0;
   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      dTotThick += m_VLayerAboveBasement[n].dGetUnconsolidatedThickness();

   return dTotThick;
}

//! Returns the total thickness of all sediment (all size classes) on this cell
double CGeomCell::dGetTotAllSedThickness(void) const
{
   return (this->dGetTotUnconsThickness() + this->dGetTotConsThickness());
}

//! Appends sediment layers
void CGeomCell::AppendLayers(int const nLayer)
{
   for (int i = 0; i < nLayer; i++)
      m_VLayerAboveBasement.push_back(CRWCellLayer());
}

//! For this cell, calculates the elevation of the top of every layer, and the d50 for the topmost unconsolidated sediment layer
void CGeomCell::CalcAllLayerElevsAndD50(void)
{
   m_VdAllHorizonTopElev.clear();
   m_VdAllHorizonTopElev.push_back(m_dBasementElevation);         // Elevation of top of the basement

   // Calculate the elevation of the top of all other layers
   int m = 0;
   for (unsigned int n = 0; n < m_VLayerAboveBasement.size(); n++)
      m_VdAllHorizonTopElev.push_back(m_VLayerAboveBasement[n].dGetTotalThickness() + m_VdAllHorizonTopElev[m++]); // Elevation of top of layer n

   // Now calculate the d50 of the topmost unconsolidated sediment layer with non-zero thickness. If there is no unconsolidated sediment, m_dUnconsD50 is set to DBL_NODATA
   m_dUnconsD50 = DBL_NODATA;
   for (int n = static_cast<int>(m_VLayerAboveBasement.size()) - 1; n >= 0; n--)
   {
      double dUnconsThick = m_VLayerAboveBasement[n].dGetUnconsolidatedThickness();
      if (dUnconsThick > 0)
      {
         // This is a layer with non-zero thickness of unconsolidated sediment
         CRWCellSediment const* pUnconsSedLayer = m_VLayerAboveBasement[n].pGetUnconsolidatedSediment();
         double dFineProp = pUnconsSedLayer->dGetFineDepth() / dUnconsThick;
         double dSandProp = pUnconsSedLayer->dGetSandDepth() / dUnconsThick;
         double dCoarseProp = pUnconsSedLayer->dGetCoarseDepth() / dUnconsThick;

         // Calculate d50 for the unconsolidated sediment
         m_dUnconsD50 = (dFineProp * m_pGrid->pGetSim()->dGetD50Fine()) + (dSandProp * m_pGrid->pGetSim()->dGetD50Sand()) + (dCoarseProp * m_pGrid->pGetSim()->dGetD50Coarse());

         break;
      }
   }
}

//! Given an elevation, this finds the index of the layer that contains that elevation (layer 0 being the first above the basement). Note that the elevation cannot exactly equal the elevation of the layer's top surface (this causes problems with e.g. cliff notches, which extend above this elevation)
int CGeomCell::nGetLayerAtElev(double const dElev) const
{
   /*! Returns ELEV_IN_BASEMENT if in basement, ELEV_ABOVE_SEDIMENT_TOP if higher than or equal to sediment top, or layer number (0 to n),  */
   if (dElev < m_VdAllHorizonTopElev[0])
      return ELEV_IN_BASEMENT;

   for (unsigned int nLayer = 1; nLayer < m_VdAllHorizonTopElev.size(); nLayer++)
   {
      if ((m_VLayerAboveBasement[nLayer - 1].dGetTotalThickness() > 0) && (dElev >= m_VdAllHorizonTopElev[nLayer - 1]) && (dElev <= m_VdAllHorizonTopElev[nLayer]))
         return (nLayer - 1);
   }

   return ELEV_ABOVE_SEDIMENT_TOP;
}

//! For this cell, calculates the elevation of the top of a given layer
double CGeomCell::dCalcLayerElev(const int nLayer)
{
   // Note no check to see if nLayer < m_VLayerAboveBasement.size()
   double dTopElev = m_dBasementElevation;

   for (int n = 0; n <= nLayer; n++)
      dTopElev += m_VLayerAboveBasement[n].dGetTotalThickness();

   return dTopElev;
}

//! Set potential (unconstrained) shore platform erosion and increment total shore platform potential erosion
void CGeomCell::SetPotentialPlatformErosion(double const dPotentialIn)
{
   m_dPotentialPlatformErosionThisIter = dPotentialIn;
   m_dTotPotentialPlatformErosion += dPotentialIn;
}

//! Get potential (unconstrained) shore platform erosion
double CGeomCell::dGetPotentialPlatformErosion(void) const
{
   return m_dPotentialPlatformErosionThisIter;
}

//! Get total potential (unconstrained) shore platform erosion
double CGeomCell::dGetTotPotentialPlatformErosion(void) const
{
   return m_dTotPotentialPlatformErosion;
}

//! Set this-timestep actual (constrained) shore platform erosion and increment total actual shore platform erosion
void CGeomCell::SetActualPlatformErosion(double const dThisActualErosion)
{
   m_dActualPlatformErosionThisIter = dThisActualErosion;
   m_dTotActualPlatformErosion += dThisActualErosion;
}

//! Get actual (constrained) shore platform erosion
double CGeomCell::dGetActualPlatformErosion(void) const
{
   return m_dActualPlatformErosionThisIter;
}

//! Get total actual (constrained) shore platform erosion
double CGeomCell::dGetTotActualPlatformErosion(void) const
{
   return m_dTotActualPlatformErosion;
}

//! Returns the depth of seawater on this cell if the sediment top is < SWL, or zero
void CGeomCell::SetSeaDepth(void)
{
   m_dSeaDepth = tMax(m_pGrid->pGetSim()->CSimulation::dGetThisIterSWL() - (m_VdAllHorizonTopElev.back() + m_dInterventionHeight), 0.0);
}

//! Initialise values for this cell
void CGeomCell::InitCell(void)
{
   m_bInContiguousSea =
   m_bInContiguousFlood =           // TODO 007 What is this?
   m_bCoastline =
   m_bFloodLine =
   m_bIsInActiveZone =
   // m_bEstimated =
   m_bShadowBoundary =
   m_bPossibleCoastStartCell =      // TODO 007 What is this?
   m_bPossibleFloodStartCell =
   m_bWaveFlood =
   // m_bCheckCell =                // TODO 007 What is this?
   m_bCheckFloodCell = false;

   m_nPolygonID =
   m_nCoastlineNormal = INT_NODATA;

   m_nShadowZoneNumber =
   m_nDownDriftZoneNumber = 0;

   m_dLocalConsSlope =
   m_dPotentialPlatformErosionThisIter =
   m_dActualPlatformErosionThisIter =
   m_dCliffCollapseFineThisIter =
   m_dCliffCollapseSandThisIter =
   m_dCliffCollapseCoarseThisIter =
   m_dTalusSandDepositionThisIter =
   m_dTotTalusSandDeposition = 
   m_dTalusCoarseDepositionThisIter = 
   m_dTotTalusCoarseDeposition = 
   m_dPotentialBeachErosionThisIter =
   m_dActualBeachErosionThisIter =
   m_dBeachDepositionThisIter =
   m_dSeaDepth =
   m_dWaveHeight =
   m_dWaveAngle = 0;

   m_dBeachProtectionFactor = DBL_NODATA;
}

//! Sets the wave height on this cell, also increments the total wave height
void CGeomCell::SetWaveHeight(double const dWaveHeight)
{
   m_dWaveHeight = dWaveHeight;
   m_dTotWaveHeight += dWaveHeight;

   //    if (m_dWaveHeight != DBL_NODATA)
   //       assert(m_dWaveHeight >= 0);
}

//! Returns the wave height on this cell
double CGeomCell::dGetWaveHeight(void) const
{
   return m_dWaveHeight;
}

//! Returns the total wave height on this cell
double CGeomCell::dGetTotWaveHeight(void) const
{
   return m_dTotWaveHeight;
}

//! Sets the wave orientation on this cell, also increments the total wave orientation
void CGeomCell::SetWaveAngle(double const dWaveAngle)
{
   m_dWaveAngle = dWaveAngle;
   m_dTotWaveAngle += dWaveAngle;
}

//! Returns the wave orientation on this cell
double CGeomCell::dGetWaveAngle(void) const
{
   return m_dWaveAngle;
}

//! Returns the total wave orientation on this cell
double CGeomCell::dGetTotWaveAngle(void) const
{
   return m_dTotWaveAngle;
}

//! Sets the deep water wave height on this cell
void CGeomCell::SetCellDeepWaterWaveHeight(double const dWaveHeight)
{
   m_dDeepWaterWaveHeight = dWaveHeight;
}

//! Returns the deep water wave height on this cell
double CGeomCell::dGetCellDeepWaterWaveHeight(void) const
{
   return m_dDeepWaterWaveHeight;
}

//! Sets the deep water wave orientation on this cell
void CGeomCell::SetCellDeepWaterWaveAngle(double const dWaveAngle)
{
   m_dDeepWaterWaveAngle = dWaveAngle;
}

//! Returns the deep water wave orientation on this cell
double CGeomCell::dGetCellDeepWaterWaveAngle(void) const
{
   return m_dDeepWaterWaveAngle;
}

//! Sets the deep water wave Period on this cell
void CGeomCell::SetCellDeepWaterWavePeriod(double const dWavePeriod)
{
   m_dDeepWaterWavePeriod = dWavePeriod;
}

//! Returns the deep water wave period on this cell
double CGeomCell::dGetCellDeepWaterWavePeriod(void) const
{
   return m_dDeepWaterWavePeriod;
}

//! Sets wave height to the deep water wave height value, and sets wave orientation to the deep water wave orientation value
void CGeomCell::SetWaveValuesToDeepWaterWaveValues(void)
{
   m_dWaveHeight = m_dDeepWaterWaveHeight;
   m_dWaveAngle = m_dDeepWaterWaveAngle;
   m_dWavePeriod = m_dDeepWaterWavePeriod;
}

// Sets this cell's beach protection factor
void CGeomCell::SetBeachProtectionFactor(double const dFactor)
{
   m_dBeachProtectionFactor = dFactor;
}

//! Returns this cell's beach protection factor
double CGeomCell::dGetBeachProtectionFactor(void) const
{
   return m_dBeachProtectionFactor;
}

//! Increments the fine, sand, and coarse depths of this-timestep cliff collapse on this cell, also increments the totals
void CGeomCell::IncrCliffCollapseErosion(double const dFineDepth, double const dSandDepth, double const dCoarseDepth)
{
   m_dCliffCollapseFineThisIter += dFineDepth;
   m_dCliffCollapseSandThisIter += dSandDepth;
   m_dCliffCollapseCoarseThisIter += dCoarseDepth;

   m_dTotFineCliffCollapse += dFineDepth;
   m_dTotSandCliffCollapse += dSandDepth;
   m_dTotCoarseCliffCollapse += dCoarseDepth;
}

//! Returns the depth of this-timestep fine-sized sediment cliff collapse on this cell
double CGeomCell::dGetThisIterCliffCollapseErosionFine(void) const
{
   return m_dCliffCollapseFineThisIter;
}

//! Returns the depth of this-timestep sand-sized sediment cliff collapse on this cell
double CGeomCell::dGetThisIterCliffCollapseErosionSand(void) const
{
   return m_dCliffCollapseSandThisIter;
}

//! Returns the depth of this-timestep coarse-sized sediment cliff collapse on this cell
double CGeomCell::dGetThisIterCliffCollapseErosionCoarse(void) const
{
   return m_dCliffCollapseCoarseThisIter;
}

//! Returns the running total depth of fine-sized sediment eroded by cliff collapse on this cell
double CGeomCell::dGetTotCliffCollapseFine(void) const
{
   return m_dTotFineCliffCollapse;
}

//! Returns the running total depth of sand-sized sediment eroded by cliff collapse on this cell
double CGeomCell::dGetTotCliffCollapseSand(void) const
{
   return m_dTotSandCliffCollapse;
}

//! Returns the running total depth of coarse-sized sediment eroded by cliff collapse on this cell
double CGeomCell::dGetTotCliffCollapseCoarse(void) const
{
   return m_dTotCoarseCliffCollapse;
}

//! Increments the depth of this-timestep sand-sized talus from cliff collapse on this cell, also increments the total
void CGeomCell::AddSandTalusDeposition(double const dDepth)
{
   m_dTalusSandDepositionThisIter += dDepth;
   m_dTotTalusSandDeposition += dDepth;
}

//! Increments the depth of this-timestep coarse-sized talus from cliff collapse on this cell, also increments the total
void CGeomCell::AddCoarseTalusDeposition(double const dDepth)
{
   m_dTalusCoarseDepositionThisIter += dDepth;
   m_dTotTalusCoarseDeposition += dDepth;
}

//! Returns the depth of this-timestep sand talus deposition from cliff collapse on this cell
double CGeomCell::dGetThisIterCliffCollapseSandTalusDeposition(void) const
{
   return m_dTalusSandDepositionThisIter;
}

//! Retuns the depth of this-timestep coarse talus deposition from cliff collapse on this cell
double CGeomCell::dGetThisIterCliffCollapseCoarseTalusDeposition(void) const
{
   return m_dTalusCoarseDepositionThisIter;
}

//! Returns the total depth of sand talus deposition from cliff collapse on this cell
double CGeomCell::dGetTotSandTalusDeposition(void) const
{
   return m_dTotTalusSandDeposition;
}

//! Returns the total depth of coarse talus deposition from cliff collapse on this cell
double CGeomCell::dGetTotCoarseTalusDeposition(void) const
{
   return m_dTotTalusCoarseDeposition;
}

//! Set potential (unconstrained) beach erosion and increment total beach potential erosion
void CGeomCell::SetPotentialBeachErosion(double const dPotentialIn)
{
   m_dPotentialBeachErosionThisIter = dPotentialIn;
   m_dTotPotentialBeachErosion += dPotentialIn;
}

//! Get potential (unconstrained) beach erosion
double CGeomCell::dGetPotentialBeachErosion(void) const
{
   return m_dPotentialBeachErosionThisIter;
}

//! Get total potential (supply-unconstrained) beach erosion
double CGeomCell::dGetTotPotentialBeachErosion(void) const
{
   return m_dTotPotentialBeachErosion;
}

//! Set this-timestep actual (supply-constrained) beach erosion and increment total actual beach erosion
void CGeomCell::SetActualBeachErosion(double const dThisActualErosion)
{
   m_dActualBeachErosionThisIter = dThisActualErosion;
   m_dTotActualBeachErosion += dThisActualErosion;
}

//! Get actual (supply-constrained) beach erosion
double CGeomCell::dGetActualBeachErosion(void) const
{
   return m_dActualBeachErosionThisIter;
}

//! Get total actual (supply-constrained) beach erosion
double CGeomCell::dGetTotActualBeachErosion(void) const
{
   return m_dTotActualBeachErosion;
}

// //! Returns true if there has been actual beach erosion this timestep
// bool CGeomCell::bActualBeachErosionThisIter(void) const
// {
//    return (m_dActualBeachErosionThisIter > 0 ? true : false);
// }

//! Increment this-timestep beach deposition, also increment total beach deposition
void CGeomCell::IncrBeachDeposition(double const dThisDeposition)
{
   m_dBeachDepositionThisIter += dThisDeposition;
   m_dTotBeachDeposition += dThisDeposition;
}

//! Get beach deposition
double CGeomCell::dGetBeachDeposition(void) const
{
   return m_dBeachDepositionThisIter;
}

//! Get beach erosion
double CGeomCell::dGetTotBeachDeposition(void) const
{
   return m_dTotBeachDeposition;
}

// //! Returns true if there has been beach deposition this timestep
// bool CGeomCell::bBeachDepositionThisIter(void) const
// {
//    return (m_dBeachDepositionThisIter > 0 ? true : false);
// }

//! Returns true only if this cell has had no deposition or erosion this timestep
bool CGeomCell::bBeachErosionOrDepositionThisIter(void) const
{
   if ((m_dActualBeachErosionThisIter > 0) || (m_dBeachDepositionThisIter > 0))
      return true;

   return false;
}

//! Returns the D50 of unconsolidated sediment on this cell
double CGeomCell::dGetUnconsD50(void) const
{
   return m_dUnconsD50;
}

//! Sets the landform category and subcategory for an intervention
void CGeomCell::SetInterventionClass(int const nSubCatCode)
{
   if (nSubCatCode != LF_NONE)
   {
      this->m_Landform.SetLFCategory(LF_CAT_INTERVENTION);

      if (nSubCatCode == IO_INTERVENTION_STRUCT)
         this->m_Landform.SetLFSubCategory(LF_SUBCAT_INTERVENTION_STRUCT);
      else if (nSubCatCode == IO_INTERVENTION_NON_STRUCT)
         this->m_Landform.SetLFSubCategory(LF_SUBCAT_INTERVENTION_NON_STRUCT);
   }
}

//! Gets the intervention class
int CGeomCell::nGetInterventionClass(void) const
{
   int nTmp = INT_NODATA;

   if (this->m_Landform.nGetLFCategory() == LF_CAT_INTERVENTION)
   {
      if (this->m_Landform.nGetLFSubCategory() == LF_SUBCAT_INTERVENTION_STRUCT)
         nTmp = IO_INTERVENTION_STRUCT;
      else if (this->m_Landform.nGetLFSubCategory() == LF_SUBCAT_INTERVENTION_NON_STRUCT)
         nTmp = IO_INTERVENTION_NON_STRUCT;
   }

   return nTmp;
}

//! Sets the intervention height
void CGeomCell::SetInterventionHeight(double const dHeight)
{
   m_dInterventionHeight = dHeight;
}

//! Returns the intervention height
double CGeomCell::dGetInterventionHeight(void) const
{
   return m_dInterventionHeight;
}

//! Returns the elevation of the top of the intervention, assuming it rests on the sediment-top surface
double CGeomCell::dGetInterventionTopElev(void) const
{
   return m_VdAllHorizonTopElev.back() + m_dInterventionHeight;
}
