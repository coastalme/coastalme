/*!
 *
 * \file profile.cpp
 * \brief CGeomProfile routines
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
// #include <assert.h>
#include <cmath>

#include <vector>

#include <algorithm>
using std::find;

#include "cme.h"
#include "profile.h"

//! Constructor with initialization list, requires one parameter (the coast point at which the profile starts)
CGeomProfile::CGeomProfile(int const nCoastPoint)
    : m_bStartOfCoast(false),
      m_bEndOfCoast(false),
      m_bHitLand(false),
      m_bHitCoast(false),
      m_bTooShort(false),
      m_bTruncated(false),
      m_bHitAnotherProfile(false),
      m_nNumCoastPoint(nCoastPoint),
      m_dDeepWaterWaveHeight(0),
      m_dDeepWaterWaveAngle(0),
      m_dDeepWaterWavePeriod(0)
{
}

//! Destructor
CGeomProfile::~CGeomProfile(void)
{
}

//! Returns the coast point at which the profile starts
int CGeomProfile::nGetNumCoastPoint(void) const
{
   return m_nNumCoastPoint;
}

//! Sets a switch to indicate whether this is a start-of-coast profile
void CGeomProfile::SetStartOfCoast(bool const bFlag)
{
   m_bStartOfCoast = bFlag;
}

//! Returns the switch to indicate whether this is a start-of-coast profile
bool CGeomProfile::bStartOfCoast(void) const
{
   return m_bStartOfCoast;
}

//! Sets a switch to indicate whether this is an end-of-coast profile
void CGeomProfile::SetEndOfCoast(bool const bFlag)
{
   m_bEndOfCoast = bFlag;
}

//! Returns the switch to indicate whether this is an end-of-coast profile
bool CGeomProfile::bEndOfCoast(void) const
{
   return m_bEndOfCoast;
}

// //! Sets a switch which indicates whether this profile has hit land
// void CGeomProfile::SetHitLand(bool const bFlag)
// {
//    m_bHitLand = bFlag;
// }

//! Returns the switch which indicates whether this profile has hit land
bool CGeomProfile::bHitLand(void) const
{
   return m_bHitLand;
}

//! Sets a switch which indicates whether this profile has hit a coast
void CGeomProfile::SetHitCoast(bool const bFlag)
{
   m_bHitCoast = bFlag;
}

//! Returns the switch which indicates whether this profile has hit a coast
bool CGeomProfile::bHitCoast(void) const
{
   return m_bHitCoast;
}

//! Sets a switch which indicates whether this profile is too short to be useful
void CGeomProfile::SetTooShort(bool const bFlag)
{
   m_bTooShort = bFlag;
}

// //! Returns the switch which indicates whether this profile is too short to be useful
// bool CGeomProfile::bTooShort(void) const
// {
//    return m_bTooShort;
// }

//! Sets a switch which indicates whether this profile is truncated
void CGeomProfile::SetTruncated(bool const bFlag)
{
   m_bTruncated = bFlag;
}

// //! Returns the switch which indicates whether this profile is truncated
// bool CGeomProfile::bTruncated(void) const
// {
//    return m_bTruncated;
// }

//! Sets a switch which indicates whether this profile hits another profile
void CGeomProfile::SetHitAnotherProfile(bool const bFlag)
{
   m_bHitAnotherProfile = bFlag;
}

//! Returns the switch which indicates whether this profile hits another profile
bool CGeomProfile::bHitAnotherProfile(void) const
{
   return m_bHitAnotherProfile;
}

//! Returns true if this is a problem-free profile, and is not a start-o-coast or an end-of-coast profile
bool CGeomProfile::bProfileOK(void) const
{
   // All profiles without problems, but not start- or end-of-coast profiles
   if ((! m_bStartOfCoast) &&
       (! m_bEndOfCoast) &&
       (! m_bHitLand) &&
       (! m_bHitCoast) &&
       (! m_bTooShort) &&
       (! m_bTruncated) &&
       (! m_bHitAnotherProfile))
      return true;

   return false;
}

//! Returns true if this is a problem-free profile (however it could still be a start-of-coast or an end-of-coast profile)
bool CGeomProfile::bOKIncStartAndEndOfCoast(void) const
{
   // All profiles without problems, including start- and end-of-coast profiles
   if ((! m_bHitLand) &&
       (! m_bHitCoast) &&
       (! m_bTooShort) &&
       (! m_bTruncated) &&
       (! m_bHitAnotherProfile))
      return true;

   return false;
}

// //! Returns true if this is a problem-free profile (however it could still be a start-of-coast profile)
// bool CGeomProfile::bOKIncStartOfCoast(void) const
// {
//    // All profiles without problems, including start-of-coast profile (but not end-of-coast profile)
//    if ((! m_bEndOfCoast) &&
//        (! m_bHitLand) &&
//        (! m_bHitCoast) &&
//        (! m_bTooShort) &&
//        (! m_bTruncated) &&
//        (! m_bHitAnotherProfile))
//       return true;
//
//    return false;
// }

//! Sets all points in the profile
void CGeomProfile::SetAllPointsInProfile(vector<CGeom2DPoint> const *VNewPoints)
{
   m_VPoints = *VNewPoints;
}

//! Sets a single point in the profile
void CGeomProfile::SetPointInProfile(int const nPoint, double const dNewX, double const dNewY)
{
   // TODO 055 No check to see if nPoint < m_VPoints,size()
   m_VPoints[nPoint] = CGeom2DPoint(dNewX, dNewY);
}

//! Appends a point to the profile
void CGeomProfile::AppendPointInProfile(double const dNewX, double const dNewY)
{
   m_VPoints.push_back(CGeom2DPoint(dNewX, dNewY));
}

//! Appends a point to the profile (overloaded version)
void CGeomProfile::AppendPointInProfile(CGeom2DPoint const* pPt)
{
   m_VPoints.push_back(*pPt);
}

//! Inserts an intersection into the profile
bool CGeomProfile::bInsertIntersection(double const dX, double const dY, int const nSeg)
{
   // Safety check
   if (nSeg >= nGetNumLineSegments())
      return false;

   vector<CGeom2DPoint>::iterator it;
   it = m_VPoints.begin();

   // Do the insertion
   m_VPoints.insert(it + nSeg + 1, CGeom2DPoint(dX, dY));

   // Now insert a line segment in the associated multi-line, this will inherit the profile/line seg details from the preceding line segment
   CGeomMultiLine::InsertLineSegment(nSeg);

   return true;
}

//! Truncates the profile
void CGeomProfile::TruncateProfile(int const nSize)
{
   m_VPoints.resize(nSize);
}

// void CGeomProfile::TruncateAndSetPointInProfile(int const nPoint, double const dNewX, double const dNewY)
// {
//    m_VPoints.resize(nPoint+1);
//    m_VPoints[nPoint] = CGeom2DPoint(dNewX, dNewY);
// }

// void CGeomProfile::ShowProfile(void) const
// {
//    for (int n = 0; n < m_VPoints.size(); n++)
//    {
//       cout << n << " [" << m_VPoints[n].dGetX() << "][" << m_VPoints[n].dGetY() << "]" << endl;
//    }
// }

//! Returns the number of points in the profile
int CGeomProfile::nGetProfileSize(void) const
{
   return static_cast<int>(m_VPoints.size());
}

//! Returns a single point in the profile
CGeom2DPoint *CGeomProfile::pPtGetPointInProfile(int const n)
{
   return &m_VPoints[n];
}

//! Returns a given point from the profile, and all points after this
vector<CGeom2DPoint> CGeomProfile::PtVGetThisPointAndAllAfter(int const nStart)
{
   return vector<CGeom2DPoint>(m_VPoints.begin() + nStart, m_VPoints.end());
}

//! Removes a line segment from the profile
void CGeomProfile::RemoveLineSegment(int const nPoint)
{
   m_VPoints.erase(m_VPoints.begin() + nPoint);
   CGeomMultiLine::RemoveLineSegment(nPoint);
}

//! Queries the profile: is the given point a profile point?
bool CGeomProfile::bIsPointInProfile(double const dX, double const dY)
{
   CGeom2DPoint Pt(dX, dY);
   auto it = find(m_VPoints.begin(), m_VPoints.end(), &Pt);
   if (it != m_VPoints.end())
      return true;
   else
      return false;
}

//! Queries the profile: is the given point a profile point? If so, then it also returns the number of the point in the profile
bool CGeomProfile::bIsPointInProfile(double const dX, double const dY, int &nPoint)
{
   CGeom2DPoint Pt(dX, dY);
   auto it = find(m_VPoints.begin(), m_VPoints.end(), &Pt);
   if (it != m_VPoints.end())
   {
      // Found, so return true and set nPoint to be the index of the point which was found
      nPoint = static_cast<int>(it - m_VPoints.begin());
      return true;
   }
   else
      return false;
}

// int CGeomProfile::nFindInsertionLineSeg(double const dInsertX, double const dInsertY)
// {
//    for (int n = 0; n < m_VPoints.back(); n++)
//    {
//       double
//          dThisX = m_VPoints[n].dGetX(),
//          dThisY = m_VPoints[n].dGetY(),
//          dNextX = m_VPoints[n+1].dGetX(),
//          dNextY = m_VPoints[n+1].dGetY();
//
//       bool
//          bBetweenX = false,
//          bBetweenY = false;
//
//       if (dNextX >= dThisX)
//       {
//          // Ascending
//          if ((dInsertX >= dThisX) && (dInsertX <= dNextX))
//             bBetweenX = true;
//       }
//       else
//       {
//          // Descending
//          if ((dInsertX >= dNextX) && (dInsertX <= dThisX))
//             bBetweenX = true;
//       }
//
//       if (dNextY >= dThisY)
//       {
//          // Ascending
//          if ((dInsertY >= dThisY) && (dInsertY <= dNextY))
//             bBetweenY = true;
//       }
//       else
//       {
//          // Descending
//          if ((dInsertY >= dNextY) && (dInsertY <= dThisY))
//             bBetweenY = true;
//       }
//
//       if (bBetweenX && bBetweenY)
//          return n;
//    }
//
//    return -1;
// }

// void CGeomProfile::AppendPointShared(bool const bShared)
// {
//    m_bVShared.push_back(bShared);
// }

// bool CGeomProfile::bPointShared(int const n) const
// {
//    // TODO 055 No check to see if n < size()
//    return m_bVShared[n];
// }

// void CGeomProfile::SetCoastPolyToLeft(int const n, int const nPoly)
// {
//    // TODO 055 No check to see if n < size()
//    m_VnCoastPolyToLeft[n] = nPoly;
// }

// int CGeomProfile::nGetCoastPolyToleft(int const n)
// {
//    // TODO 055 No check to see if n < size()
//    return m_VnCoastPolyToLeft[n];
// }

// void CGeomProfile::SetCoastPolyToRight(int const n, int const nPoly)
// {
//    // TODO 055 No check to see if n < size()
//    m_VnCoastPolyToRight[n] = nPoly;
// }

// int CGeomProfile::nGetCoastPolyToRight(int const n)
// {
//    // TODO 055 No check to see if n < size()
//    return m_VnCoastPolyToRight[n];
// }

//! Appends a cell to the profile
void CGeomProfile::AppendCellInProfile(CGeom2DIPoint const* pPti)
{
   // In grid CRS
   m_VCellInProfile.push_back(*pPti);
}

//! Appends a cell to the profile (overloaded version)
void CGeomProfile::AppendCellInProfile(int const nX, int const nY)
{
   // In grid CRS
   m_VCellInProfile.push_back(CGeom2DIPoint(nX, nY));
}

// void CGeomProfile::SetCellsInProfile(vector<CGeom2DIPoint>* VNewPoints)
// {
//    // In grid CRS
//    m_VCellInProfile = *VNewPoints;
// }

//! Returns all cells in the profile
vector<CGeom2DIPoint> *CGeomProfile::pPtiVGetCellsInProfile(void)
{
   // In grid CRS
   return &m_VCellInProfile;
}

//! Returns a single cell in the profile
CGeom2DIPoint *CGeomProfile::pPtiGetCellInProfile(int const n)
{
   // In grid CRS TODO 055 No check to see if n < size()
   return &m_VCellInProfile[n];
}

//! Returns the number of cells in the profile
int CGeomProfile::nGetNumCellsInProfile(void) const
{
   // In grid CRS
   return static_cast<int>(m_VCellInProfile.size());
}

// vector<CGeom2DPoint>* CGeomProfile::PtVGetCellsInProfileExtCRS(void)
// {
//    // In external CRS
//    return &m_VCellInProfileExtCRS;
// }

//! Appends a cell (specified in the external co-ordinate system) to the profile
void CGeomProfile::AppendCellInProfileExtCRS(double const dX, double const dY)
{
   // In external CRS
   m_VCellInProfileExtCRS.push_back(CGeom2DPoint(dX, dY));
}

//! Appends a cell (specified in the external co-ordinate system) to the profile (overloaded version)
void CGeomProfile::AppendCellInProfileExtCRS(CGeom2DPoint const* pPt)
{
   // In external CRS
   m_VCellInProfileExtCRS.push_back(*pPt);
}

//! Returns the index of the cell on this profile which has a sea depth which is just less than a given depth. If every cell on the profile has a sea depth which is less than the given depth it returns INT_NODATA
int CGeomProfile::nGetCellGivenDepth(CGeomRasterGrid const* pGrid, double const dDepthIn)
{
   int nIndex = INT_NODATA; // If not found, i.e. if every profile cell has sea depth less than dDepthIn

   for (unsigned int n = 0; n < m_VCellInProfile.size(); n++)
   {
      int
          nX = m_VCellInProfile[n].nGetX(),
          nY = m_VCellInProfile[n].nGetY();

      double dCellDepth = pGrid->m_Cell[nX][nY].dGetSeaDepth();
      if (dCellDepth >= dDepthIn)
      {
         nIndex = n;

         if (n > 0)
            nIndex = n - 1; // Grid CRS units

         break;
      }
   }

   if (nIndex == INT_NODATA)
   {
      nIndex = static_cast<int>(m_VCellInProfile.size()) - 1;
   }

   return nIndex;
}

//! Sets the deep-water wave height for this profile
void CGeomProfile::SetProfileDeepWaterWaveHeight(double const dWaveHeight)
{
   m_dDeepWaterWaveHeight = dWaveHeight;
}

//! Returns the deep-water wave height for this profile
double CGeomProfile::dGetProfileDeepWaterWaveHeight(void) const
{
   return m_dDeepWaterWaveHeight;
}

//! Sets the deep-water wave orientation for this profile
void CGeomProfile::SetProfileDeepWaterWaveAngle(double const dWaveAngle)
{
   m_dDeepWaterWaveAngle = dWaveAngle;
}

//! Returns the deep-water wave orientation for this profile
double CGeomProfile::dGetProfileDeepWaterWaveAngle(void) const
{
   return m_dDeepWaterWaveAngle;
}

//! Sets the deep-water wave Period for this profile
void CGeomProfile::SetProfileDeepWaterWavePeriod(double const dWavePeriod)
{
   m_dDeepWaterWavePeriod = dWavePeriod;
}

//! Returns the deep-water wave Period for this profile
double CGeomProfile::dGetProfileDeepWaterWavePeriod(void) const
{
   return m_dDeepWaterWavePeriod;
}
