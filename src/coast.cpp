/*!
 *
 * \file coast.cpp
 * \brief CRWCoast routines
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
#include <assert.h>

#include <vector>
#include <algorithm>

#include "cme.h"
#include "coast.h"
#include "line.h"
#include "i_line.h"

//! Constructor with initialization list
CRWCoast::CRWCoast(void)
    : m_nSeaHandedness(NULL_HANDED),
      m_nStartEdge(INT_NODATA),
      m_nEndEdge(INT_NODATA),
      m_dCurvatureDetailedMean(0),
      m_dCurvatureDetailedSTD(0),
      m_dCurvatureSmoothMean(0),
      m_dCurvatureSmoothSTD(0)
{
}

//! Destructor
CRWCoast::~CRWCoast(void)
{
   for (unsigned int i = 0; i < m_pVLandforms.size(); i++)
      delete m_pVLandforms[i];

   for (unsigned int i = 0; i < m_pVPolygon.size(); i++)
      delete m_pVPolygon[i];
}

//! Sets the handedness of the coast
void CRWCoast::SetSeaHandedness(int const nNewHandedness)
{
   m_nSeaHandedness = nNewHandedness;
}

//! Gers the handedness of the coast
int CRWCoast::nGetSeaHandedness(void) const
{
   return m_nSeaHandedness;
}

//! Sets the coast's start edge
void CRWCoast::SetStartEdge(int const nEdge)
{
   m_nStartEdge = nEdge;
}

//! Gets the coast's start edge
int CRWCoast::nGetStartEdge(void) const
{
   return m_nStartEdge;
}

//! Sets the coast's end edge
void CRWCoast::SetEndEdge(int const nEdge)
{
   m_nEndEdge = nEdge;
}

//! Gets the coast's end edge
int CRWCoast::nGetEndEdge(void) const
{
   return m_nEndEdge;
}

//! Given the vector line of a coast initializes, initializes coastline values (curvature, breaking wave height, wave angle, and flux orientation etc.)
void CRWCoast::SetCoastlineExtCRS(CGeomLine const* pLCoast)
{
   m_LCoastlineExtCRS = *pLCoast;

   int nLen = m_LCoastlineExtCRS.nGetSize();

   m_VnProfileNumber = vector<int>(nLen, INT_NODATA);
   m_VnPolygonNode = vector<int>(nLen, INT_NODATA);
   m_VnBreakingDistance = vector<int>(nLen, INT_NODATA);

   m_VdCurvatureDetailed = vector<double>(nLen, DBL_NODATA);
   m_VdCurvatureSmooth = vector<double>(nLen, DBL_NODATA);
   m_VdDeepWaterWaveHeight = vector<double>(nLen, DBL_NODATA);
   m_VdDeepWaterWaveAngle = vector<double>(nLen, DBL_NODATA);
   m_VdDeepWaterWavePeriod = vector<double>(nLen, DBL_NODATA);
   m_VdBreakingWaveHeight = vector<double>(nLen, DBL_NODATA);
   m_VdWaveSetupSurge = vector<double>(nLen, 0); // it is better to initiate with DBL_NODATA but some values are outside of range in the interpolation
   // m_VdStormSurge = vector<double>(nLen, DBL_NODATA);
   m_VdRunUp = vector<double>(nLen, 0);
   m_VdCoastWaveHeight = vector<double>(nLen, DBL_NODATA);
   m_VdBreakingWaveAngle = vector<double>(nLen, DBL_NODATA);
   m_VdDepthOfBreaking = vector<double>(nLen, DBL_NODATA);
   m_VdFluxOrientation = vector<double>(nLen, DBL_NODATA);
   m_VdWaveEnergyAtBreaking = vector<double>(nLen, 0);
}

// //! Appends a coastline point (in external CRS), also appends dummy values for curvature, breaking wave height, wave angle, and flux orientation
// void CRWCoast::AppendPointToCoastlineExtCRS(double const dX, double const dY)
// {
//    m_LCoastlineExtCRS.Append(dX, dY);
//
//    m_VnProfileNumber.push_back(INT_NODATA);
//    m_VnPolygonNode.push_back(INT_NODATA);
//    m_VnBreakingDistance.push_back(INT_NODATA);
//
//    m_VdCurvatureDetailed.push_back(DBL_NODATA);
//    m_VdCurvatureSmooth.push_back(DBL_NODATA);
//    m_VdDeepWaterWaveHeight.push_back(DBL_NODATA);
//    m_VdDeepWaterWaveAngle.push_back(DBL_NODATA);
//    m_VdBreakingWaveHeight.push_back(DBL_NODATA);
//    m_VdWaveSetupSurge.push_back(DBL_NODATA);
//    // m_VdStormSurge.push_back(DBL_NODATA);
//    m_VdRunUp.push_back(DBL_NODATA);
//    m_VdBreakingWaveAngle.push_back(DBL_NODATA);
//    m_VdDepthOfBreaking.push_back(DBL_NODATA);
//    m_VdFluxOrientation.push_back(DBL_NODATA);
//    m_VdWaveEnergyAtBreaking.push_back(0);
// }

// void CRWCoast::SetFloodWaveSetupPointExtCRS(CGeomLine const* pLCoast)
// {
//    m_LFloodWaveSetupExtCRS = *pLCoast;

//    // int nLen = m_LFloodWaveSetupLineExtCRS.nGetSize();
// }

// void CRWCoast::SetFloodWaveSetupSurgePointExtCRS(CGeomLine const* pLCoast)
// {
//    m_LFloodWaveSetupSurgeExtCRS = *pLCoast;

//    // int nLen = m_LFloodWaveSetupSurgeLineExtCRS.nGetSize();
// }

// void CRWCoast::SetFloodWaveSetupSurgeRunupPointExtCRS(CGeomLine const* pLCoast)
// {
//    m_LFloodWaveSetupSurgeRunupExtCRS = *pLCoast;

//    // int nLen = m_LFloodWaveSetupSurgeRunupLineExtCRS.nGetSize();
// }

//! Returns the coastline (external CRS)
CGeomLine *CRWCoast::pLGetCoastlineExtCRS(void)
{
   return &m_LCoastlineExtCRS;
}

//! Returns a given coast point in external CRS
CGeom2DPoint *CRWCoast::pPtGetCoastlinePointExtCRS(int const n)
{
   // Point is in external CRS TODO 055 No check to see that n is < m_LCoastlineExtCRS.Size()
   return &m_LCoastlineExtCRS[n];
}

// CGeomLine *CRWCoast::pLGetFloodWaveSetupExtCRS(void)
// {
//    return &m_LFloodWaveSetupExtCRS;
// }

// CGeom2DPoint *CRWCoast::pPtGetFloodWaveSetupPointExtCRS(int const n)
// {
//    // Point is in external CRS TODO 055 No check to see that n is < m_LCoastlineExtCRS.Size()
//    return &m_LFloodWaveSetupExtCRS[n];
// }

// CGeom2DPoint *CRWCoast::pPtGetFloodWaveSetupSurgePointExtCRS(int const n)
// {
//    // Point is in external CRS TODO 055 No check to see that n is < m_LCoastlineExtCRS.Size()
//    return &m_LFloodWaveSetupSurgeExtCRS[n];
// }

// CGeom2DPoint *CRWCoast::pPtGetFloodWaveSetupSurgeRunupPointExtCRS(int const n)
// {
//    // Point is in external CRS TODO 055 No check to see that n is < m_LCoastlineExtCRS.Size()
//    return &m_LFloodWaveSetupSurgeRunupExtCRS[n];
// }

//! Gets the size of the coastline
int CRWCoast::nGetCoastlineSize(void) const
{
   return m_LCoastlineExtCRS.nGetSize();
}

// void CRWCoast::DisplayCoastline(void)
// {
//    m_LCoastlineExtCRS.Display();
// }

//! Sets the co-ordinates (grid CRS) of the cells marked as coastline
void CRWCoast::SetCoastlineGridCRS(CGeomILine const* pILCoastCells)
{
   m_ILCellsMarkedAsCoastline = *pILCoastCells;
}

// void CRWCoast::AppendCellMarkedAsCoastline(CGeom2DIPoint const* pPti)
// {
//    m_ILCellsMarkedAsCoastline.Append(*pPti);
// }
//
// void CRWCoast::AppendCellMarkedAsCoastline(int const nX, int const nY)
// {
//    m_ILCellsMarkedAsCoastline.Append(CGeom2DIPoint(nX, nY));
// }

//! Returns the co-ordinates (grid CRS) of the cells marked as coastline
CGeom2DIPoint *CRWCoast::pPtiGetCellMarkedAsCoastline(int const n)
{
   // TODO 055 No check to see if n < size()
   return &m_ILCellsMarkedAsCoastline[n];
}

// int CRWCoast::nGetNCellsMarkedAsCoastline(void) const
// {
//    return m_ILCellsMarkedAsCoastline.size();
// }

// double CRWCoast::dGetCoastlineSegmentLength(int const m, int const n)
// {
//    // TODO 055 No check to see that m is < m_LCoastlineExtCRS.Size(), same for n
//    if (m == n)
//       return 0;
//
//    return hypot(m_LCoastlineExtCRS[n].dGetX() - m_LCoastlineExtCRS[m].dGetX(), m_LCoastlineExtCRS[n].dGetY() - m_LCoastlineExtCRS[m].dGetY());
// }

// double CRWCoast::dGetCoastlineLengthSoFar(int const n)
// {
//    // TODO 055 No check to see that n is < m_LCoastlineExtCRS.Size()
//    double dLen = 0;
//    for (int m = 0; m < n; m++)
//       dLen += dGetCoastlineSegmentLength(m, m+1);
//    return dLen;
// }

//! Returns the coastline number given a cell, or INT_NODATA if neither this cell or any of its neighbouring cells are 'under' a coastline. If it is a neighbouring cell that is under the coastline, then it also changes the cell that is supplied as an input parameter
int CRWCoast::nGetCoastPointGivenCell(CGeom2DIPoint* pPtiCell)
{
   for (int nCoastPoint = 0; nCoastPoint < m_ILCellsMarkedAsCoastline.nGetSize(); nCoastPoint++)
   {
      if (m_ILCellsMarkedAsCoastline[nCoastPoint] == pPtiCell)
      {
         return nCoastPoint;
      }
   }

   // This cell is not under a coastline, so try the adjacent cells
   int
       n = -1,
       nX = pPtiCell->nGetX(),
       nY = pPtiCell->nGetY(),
       nXAdj = 0,
       nYAdj = 0;

   while (n <= 7)
   {
      switch (++n)
      {
         case 0:
            nXAdj = nX;
            nYAdj = nY - 1;
            break;
         case 1:
            nXAdj = nX + 1;
            nYAdj = nY - 1;
            break;
         case 2:
            nXAdj = nX + 1;
            nYAdj = nY;
            break;
         case 3:
            nXAdj = nX + 1;
            nYAdj = nY + 1;
            break;
         case 4:
            nXAdj = nX;
            nYAdj = nY + 1;
            break;
         case 5:
            nXAdj = nX - 1;
            nYAdj = nY + 1;
            break;
         case 6:
            nXAdj = nX - 1;
            nYAdj = nY;
            break;
         case 7:
            nXAdj = nX - 1;
            nYAdj = nY - 1;
            break;
      }

      CGeom2DIPoint PtiTmp(nXAdj, nYAdj);
      for (int nCoastPoint = 0; nCoastPoint < m_ILCellsMarkedAsCoastline.nGetSize(); nCoastPoint++)
      {
         if (m_ILCellsMarkedAsCoastline[nCoastPoint] == &PtiTmp)
         {
            *pPtiCell = PtiTmp;
            return nCoastPoint;
         }
      }
   }

   return INT_NODATA;
}

//! Returns the detailed curvature for a coast point
double CRWCoast::dGetDetailedCurvature(int const nCoastPoint) const
{
   // TODO 055 No sanity check for nCoastPoint < m_VdCurvatureDetailed.Size()
   return m_VdCurvatureDetailed[nCoastPoint];
}

//! Sets the detailed curvature for a coast point
void CRWCoast::SetDetailedCurvature(int const nCoastPoint, double const dCurvature)
{
   // TODO 055 No check to see if nCoastPoint < m_VdCurvatureDetailed.size()
   m_VdCurvatureDetailed[nCoastPoint] = dCurvature;
}

//! Returns a vector of detailed curvature for all coast points
vector<double> *CRWCoast::pVGetDetailedCurvature(void)
{
   return &m_VdCurvatureDetailed;
}

//! Returns the smoothed curvature for a coast point
double CRWCoast::dGetSmoothCurvature(int const nCoastPoint) const
{
   // TODO 055 No sanity check for nCoastPoint < m_VdCurvatureSmooth.Size()
   return m_VdCurvatureSmooth[nCoastPoint];
}

//! Sets the smoothed curvature for a coast point
void CRWCoast::SetSmoothCurvature(int const nCoastPoint, double const dCurvature)
{
   // TODO 055 No check to see if nCoastPoint < m_VdCurvatureSmooth.size()
   m_VdCurvatureSmooth[nCoastPoint] = dCurvature;
}

//! Returns a vector of smoothed curvature for all coast points
vector<double> *CRWCoast::pVGetSmoothCurvature(void)
{
   return &m_VdCurvatureSmooth;
}

//! Sets the mean of the coast's detailed curvature
void CRWCoast::SetDetailedCurvatureMean(double const dMean)
{
   m_dCurvatureDetailedMean = dMean;
}

// //! Gets the mean of the coast's detailed curvature
// double CRWCoast::dGetDetailedCurvatureMean(void) const
// {
//    return m_dCurvatureDetailedMean;
// }

//! Sets the standard deviation of the coast's detailed curvature
void CRWCoast::SetDetailedCurvatureSTD(double const dSTD)
{
   m_dCurvatureDetailedSTD = dSTD;
}

// //! Gets the standard deviation of the coast's detailed curvature
// double CRWCoast::dGetDetailedCurvatureSTD(void) const
// {
//    return m_dCurvatureDetailedSTD;
// }

//! Sets the mean of the coast's smoothed curvature
void CRWCoast::SetSmoothCurvatureMean(double const dMean)
{
   m_dCurvatureSmoothMean = dMean;
}

//! Gets the mean of the coast's smoothed curvature
double CRWCoast::dGetSmoothCurvatureMean(void) const
{
   return m_dCurvatureSmoothMean;
}

//! Sets the standard deviation of the coast's smoothed curvature
void CRWCoast::SetSmoothCurvatureSTD(double const dSTD)
{
   m_dCurvatureSmoothSTD = dSTD;
}

//! Gets the standard deviation of the coast's smoothed curvature
double CRWCoast::dGetSmoothCurvatureSTD(void) const
{
   return m_dCurvatureSmoothSTD;
}

//! Returns a pointer to a profile
CGeomProfile* CRWCoast::pGetProfile(int const nProfile)
{
   // TODO 055 Maybe add a safety check? that nProfile < m_VProfile.size()
   return &m_VProfile[nProfile];
}

//! Appends a profile
void CRWCoast::AppendProfile(int const nCoastPoint, int const nProfile)
{
   CGeomProfile Profile(nCoastPoint);
   m_VProfile.push_back(Profile);

   m_VnProfileNumber[nCoastPoint] = nProfile;
}

// void CRWCoast::ReplaceProfile(int const nProfile, vector<CGeom2DPoint> const* pPtVProfileNew)
// {
//    // TODO 055 Maybe add a safety check? that nProfile < m_VProfile.size()
//    m_VProfile[nProfile].SetAllPointsInProfile(pPtVProfileNew);
// }

//! Returns the number of profiles on this coast
int CRWCoast::nGetNumProfiles(void) const
{
   return static_cast<int>(m_VProfile.size());
}

//! Returns true if the given coast point is the start point of a profile
bool CRWCoast::bIsProfileStartPoint(int const nCoastPoint) const
{
   // TODO 055 No sanity check for nCoastPoint < m_VnProfileNumber.Size()
   if (m_VnProfileNumber[nCoastPoint] != INT_NODATA)
      return true;

   return false;
}

//! Returns the number of the profile at the given coast point, returns INT_NODATA if there is no profile at this point
int CRWCoast::nGetProfileNumber(int const nCoastPoint) const
{
   return m_VnProfileNumber[nCoastPoint];
}

//! Creates an index containing the numbers of the coastline-normal profiles in along-coast sequence
void CRWCoast::CreateAlongCoastProfileIndex(void)
{
   for (int nCoastPoint = 0; nCoastPoint < m_LCoastlineExtCRS.nGetSize(); nCoastPoint++)
   {
      if (m_VnProfileNumber[nCoastPoint] != INT_NODATA)
         m_VnProfileCoastIndex.push_back(m_VnProfileNumber[nCoastPoint]);
   }
}

//! Returns the number of the coastline-normal profile which is at the given coast point in the along-coast sequence
int CRWCoast::nGetProfileFromAlongCoastProfileIndex(int const n) const
{

   return m_VnProfileCoastIndex[n];
}

//! Returns the number of the profile which is adjacent to and down-coast from the specified profile. It returns INT_NODATA if there is no valid up-coast profile
int CRWCoast::nGetDownCoastProfileNumber(int const nProfile) const
{
   for (unsigned int n = 0; n < m_VnProfileCoastIndex.size() - 1; n++)
   {
      if (nProfile == m_VnProfileCoastIndex[n])
         return m_VnProfileCoastIndex[n + 1];
   }

   // At end of m_VnProfileCoastIndex, so no down-coast profile
   return INT_NODATA;
}

// int CRWCoast::nGetAlongCoastlineIndexOfProfile(int const nProfile)
// {
//    // Returns the along-coastline index of a coastline-normal profile
//    for (int n = 0; n < m_VnProfileCoastIndex.size(); n++)
//       if (m_VnProfileCoastIndex[n] == nProfile)
//          return n;
//    return -1;
// }

//! Sets the deep water wave height for this coast point
void CRWCoast::SetCoastDeepWaterWaveHeight(int const nCoastPoint, double const dHeight)
{
   // TODO 055 No check to see if nCoastPoint < m_VdDeepWaterWaveHeight.size()
   m_VdDeepWaterWaveHeight[nCoastPoint] = dHeight;
}

// //! Gets the deep water wave height for this coast point
// double CRWCoast::dGetCoastDeepWaterWaveHeight(int const nCoastPoint) const
// {
//    // TODO 055 No check to see if nCoastPoint < m_VdDeepWaterWaveHeight.size()
//    return m_VdDeepWaterWaveHeight[nCoastPoint];
// }

//! Sets the deep water wave angle for this coast point
void CRWCoast::SetCoastDeepWaterWaveAngle(int const nCoastPoint, double const dOrientation)
{
   // TODO 055 No check to see if nCoastPoint < m_VdDeepWaterWaveAngle.size()
   m_VdDeepWaterWaveAngle[nCoastPoint] = dOrientation;
}

//! Gets the deep water wave angle for this coast point
double CRWCoast::dGetCoastDeepWaterWaveAngle(int const nCoastPoint) const
{
   // TODO 055 No check to see if nCoastPoint < m_VdDeepWaterWaveAngle.size()
   return m_VdDeepWaterWaveAngle[nCoastPoint];
}

//! Sets the deep water wave period for this coast point
void CRWCoast::SetCoastDeepWaterWavePeriod(int const nCoastPoint, double const dPeriod)
{
   m_VdDeepWaterWavePeriod[nCoastPoint] = dPeriod;
}

//! Gets the deep water wave period for this coast point
double CRWCoast::dGetCoastDeepWaterWavePeriod(int const nCoastPoint) const
{
   return m_VdDeepWaterWavePeriod[nCoastPoint];
}

//! Sets the breaking wave height for this coast point
void CRWCoast::SetBreakingWaveHeight(int const nCoastPoint, double const dHeight)
{
   // TODO 055 No check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   m_VdBreakingWaveHeight[nCoastPoint] = dHeight;
}

//! Gets the breaking wave height for this coast point
double CRWCoast::dGetBreakingWaveHeight(int const nCoastPoint) const
{
   // TODO 055 No check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   return m_VdBreakingWaveHeight[nCoastPoint];
}

//! Sets the wave setup surge for this coast point
void CRWCoast::SetWaveSetupSurge(int const nCoastPoint, double const dWaveSetup)
{
   m_VdWaveSetupSurge[nCoastPoint] = dWaveSetup;
}

//! Gets the wave setup surge for this coast point
double CRWCoast::dGetWaveSetupSurge(int const nCoastPoint) const
{
   return m_VdWaveSetupSurge[nCoastPoint];
}

// void CRWCoast::SetStormSurge(int const nCoastPoint, double const dStormSurge)
// {
//    m_VdStormSurge[nCoastPoint] = dStormSurge;
// }

// double CRWCoast::dGetStormSurge(int const nCoastPoint) const
// {
//    return m_VdStormSurge[nCoastPoint];
// }

//! Sets the wave runup for this coast point
void CRWCoast::SetRunUp(int const nCoastPoint, double const dRunUp)
{
   m_VdRunUp[nCoastPoint] = dRunUp;
}

//! Gets the wave runup for this coast point
double CRWCoast::dGetRunUp(int const nCoastPoint) const
{
   return m_VdRunUp[nCoastPoint];
}

//! Sets the wave level for this coast point
double CRWCoast::dGetLevel(int const nCoastPoint, int const level) const
{
   switch (level)
   {
      case 0: // WAVESETUPSURGE:
         return m_VdWaveSetupSurge[nCoastPoint];
         break;
      case 1: // WAVESETUPSURGE + RUNUP:
         return m_VdWaveSetupSurge[nCoastPoint] + m_VdRunUp[nCoastPoint];
         break;
      default:
         return 0;
   }
}

//! Sets the coast wave height for this coast point
void CRWCoast::SetCoastWaveHeight(int const nCoastPoint, double const dHeight)
{
   // TODO 055 No check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   m_VdCoastWaveHeight[nCoastPoint] = dHeight;
}

//! Gets the coast wave height for this coast point
double CRWCoast::dGetCoastWaveHeight(int const nCoastPoint) const
{
   // TODO 055 No check to see if nCoastPoint < m_VdBreakingWaveHeight.size()
   return m_VdCoastWaveHeight[nCoastPoint];
}

//! Sets the breaking wave angle for this coast point
void CRWCoast::SetBreakingWaveAngle(int const nCoastPoint, double const dOrientation)
{
   // TODO 055 No check to see if nCoastPoint < m_VdBreakingWaveAngle.size()
   m_VdBreakingWaveAngle[nCoastPoint] = dOrientation;
}

//! Gets the breaking wave angle for this coast point
double CRWCoast::dGetBreakingWaveAngle(int const nCoastPoint) const
{
   // TODO 055 No check to see if nCoastPoint < m_VdBreakingWaveAngle.size()
   return m_VdBreakingWaveAngle[nCoastPoint];
}

//! Sets the depth of breaking for this coast point
void CRWCoast::SetDepthOfBreaking(int const nCoastPoint, double const dDepth)
{
   // TODO 055 No check to see if nCoastPoint < m_VdDepthOfBreaking.size()
   m_VdDepthOfBreaking[nCoastPoint] = dDepth;
}

//! Gets the depth of breaking for this coast point
double CRWCoast::dGetDepthOfBreaking(int const nCoastPoint) const
{
   // TODO 055 No check to see if nCoastPoint < m_VdDepthOfBreaking.size()
   return m_VdDepthOfBreaking[nCoastPoint];
}

//! Sets the breaking distance for this coast point
void CRWCoast::SetBreakingDistance(int const nCoastPoint, int const nDist)
{
   // TODO 055 No check to see if nCoastPoint < m_VnBreakingDistance.size()
   m_VnBreakingDistance[nCoastPoint] = nDist;
}

//! Gets the breaking distance for this coast point
int CRWCoast::nGetBreakingDistance(int const nCoastPoint) const
{
   // TODO 055 No check to see if nCoastPoint < m_VnBreakingDistance.size()
   return m_VnBreakingDistance[nCoastPoint];
}

//! Sets the flux orientation for this coast point
void CRWCoast::SetFluxOrientation(int const nCoastPoint, double const dOrientation)
{
   // TODO 055 No check to see if nCoastPoint < m_VdFluxOrientation.size()
   m_VdFluxOrientation[nCoastPoint] = dOrientation;
}

//! Gets the flux orientation for this coast point
double CRWCoast::dGetFluxOrientation(int const nCoastPoint) const
{
   // TODO 055 No check to see if nCoastPoint < m_VdFluxOrientation.size()
   return m_VdFluxOrientation[nCoastPoint];
}

//! Sets the wave energy at breaking for this coast point
void CRWCoast::SetWaveEnergyAtBreaking(int const nCoastPoint, double const dEnergy)
{
   // TODO 055 No check to see if nCoastPoint < m_VdWaveEnergyAtBreaking.size()
   //    assert(isfinite(dEnergy));
   m_VdWaveEnergyAtBreaking[nCoastPoint] = dEnergy;
}

//! Gets the wave energy at breaking for this coast point
double CRWCoast::dGetWaveEnergyAtBreaking(int const nCoastPoint) const
{
   // TODO 055 No check to see if nCoastPoint < m_VdWaveEnergyAtBreaking.size()
   //    assert(isfinite(m_VdWaveEnergyAtBreaking[nCoastPoint]));
   return m_VdWaveEnergyAtBreaking[nCoastPoint];
}

//! Appends a coastal landform to this coast
void CRWCoast::AppendCoastLandform(CACoastLandform* pCoastLandform)
{
   m_pVLandforms.push_back(pCoastLandform);
}

//! Returns the coastal landform for a given coast point, or NULL if there is no coast landform here
CACoastLandform *CRWCoast::pGetCoastLandform(int const nCoastPoint)
{
   if (nCoastPoint < static_cast<int>(m_pVLandforms.size()))
      return m_pVLandforms[nCoastPoint];

   return NULL;
}

//! Sets a coast polygon node
void CRWCoast::SetPolygonNode(int const nPoint, int const nNode)
{
   // TODO 055 No check to see if nPoint < m_VnPolygonNode.size()
   m_VnPolygonNode[nPoint] = nNode;
}

//! Gets a coast polygon node
int CRWCoast::nGetPolygonNode(int const nPoint) const
{
   // TODO 055 No check to see if nPoint < m_VnPolygonNode.size()
   return m_VnPolygonNode[nPoint];
}

//! Creates a coast polygon
void CRWCoast::CreatePolygon(int const nGlobalID, int const nCoastID, int const nCoastPoint, CGeom2DIPoint const *PtiNode, CGeom2DIPoint const* PtiAntiNode, int const nProfileUpCoast, int const nProfileDownCoast, vector<CGeom2DPoint> const* pVIn, int const nPointsUpCoastProfile, int const nPointsDownCoastProfile, int const nPointInPolygonStartPoint)
{
   CGeomCoastPolygon* pPolygon = new CGeomCoastPolygon(nGlobalID, nCoastID, nCoastPoint, nProfileUpCoast, nProfileDownCoast, pVIn, nPointsUpCoastProfile, nPointsDownCoastProfile, PtiNode, PtiAntiNode, nPointInPolygonStartPoint);

   m_pVPolygon.push_back(pPolygon);
}

//! Returns the number of coast polygons
int CRWCoast::nGetNumPolygons(void) const
{
   return static_cast<int>(m_pVPolygon.size());
}

//! Returns a pointer to the coast polygon at a given coast point
CGeomCoastPolygon* CRWCoast::pGetPolygon(int const nPoly) const
{
   // TODO 055 No check to see if nPoint < m_VnPolygonNode.size()
   return m_pVPolygon[nPoly];
}

//! Appends to coast polygon length
void CRWCoast::AppendPolygonLength(const double dLength)
{
   m_VdPolygonLength.push_back(dLength);
}

//! Gets coast polygon length
double CRWCoast::dGetPolygonLength(int const nIndex) const
{
   // TODO 055 No check to see if nIndex < m_VdPolygonLength.size()
   return m_VdPolygonLength[nIndex];
}

//! Returns the number of shadow boundaries on this coast
int CRWCoast::nGetNumShadowBoundaries(void) const
{
   return static_cast<int>(m_LShadowBoundary.size());
}

//! Appends a shadow boundary to this coast
void CRWCoast::AppendShadowBoundary(CGeomLine const* pLBoundary)
{
   m_LShadowBoundary.push_back(*pLBoundary);
}

//! Returns a pointer to a shadow boundary
CGeomLine *CRWCoast::pGetShadowBoundary(int const n)
{
   // TODO 055 No check to see if n < m_LShadowBoundary.size()
   return &m_LShadowBoundary[n];
}

//! Returns the number of shadow zone downdrift boundaries on this coast
int CRWCoast::nGetNumShadowDowndriftBoundaries(void) const
{
   return static_cast<int>(m_LShadowDowndriftBoundary.size());
}

//! Appends a shadow zone downdrift boundary
void CRWCoast::AppendShadowDowndriftBoundary(CGeomLine const* pLBoundary)
{
   m_LShadowDowndriftBoundary.push_back(*pLBoundary);
}

//! Returns a pointer to a shadow zone downdrift boundary
CGeomLine* CRWCoast::pGetShadowDowndriftBoundary(int const n)
{
   // TODO 055 No check to see if n < m_LShadowDowndriftBoundary.size()
   return &m_LShadowDowndriftBoundary[n];
}
