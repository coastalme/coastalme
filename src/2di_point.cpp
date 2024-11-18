/*!
 *
 * \file 2di_point.cpp
 * \brief Geometry class used to represent 2D point objects with integer co-ordinates
 * \details The CGeom2DIPoint geometry class is used to represent 2D points where the x and y co-ordinates can only be integer values, e.g. points for which the x and y co-ordinates are in the raster-grid CRS (co-ordinate reference system)
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
#include "2di_point.h"
#include "cme.h"

//! Constructor with no parameters (the X and Y co-ordinates of the CGeom2DIPoint object are set to zero in an initialization list)
CGeom2DIPoint::CGeom2DIPoint(void)
:  nX(0),
   nY(0)
{
}

//! Constructor with two integer parameters, for the X and Y co-ordinates of the CGeom2DIPoint object
CGeom2DIPoint::CGeom2DIPoint(int const nNewX, int const nNewY)
:  nX(nNewX),
   nY(nNewY)
{
}

//! Returns the CGeom2DIPoint object's integer X co-ordinate
int CGeom2DIPoint::nGetX(void) const
{
   return nX;
}

//! Returns the CGeom2DIPoint object's integer Y co-ordinate
int CGeom2DIPoint::nGetY(void) const
{
   return nY;
}

//! Returns a reference to the CGeom2DIPoint object's integer X co-ordinate
int* CGeom2DIPoint::pnGetX(void)
{
   return &nX;
}

//! Returns a reference to the CGeom2DIPoint object's integer Y co-ordinate
int* CGeom2DIPoint::pnGetY(void)
{
   return &nY;
}

//! The integer parameter sets a value for the CGeom2DIPoint object's X co-ordinate
void CGeom2DIPoint::SetX(int const nNewX)
{
   nX = nNewX;
}

//! The integer parameter sets a value for the CGeom2DIPoint object's Y co-ordinate
void CGeom2DIPoint::SetY(int const nNewY)
{
   nY = nNewY;
}

//! The two integer parameters set values for the CGeom2DIPoint object's X and Y co-ordinates
void CGeom2DIPoint::SetXY(int const nNewX, int const nNewY)
{
   nX = nNewX;
   nY = nNewY;
}

//! The parameter is a pointer to a CGeom2DIPoint object, this is used to set values for the CGeom2DIPoint object's X and Y co-ordinates
// void CGeom2DIPoint::SetXY(CGeom2DIPoint const* Pti)
// {
//    nX = Pti->nGetX();
//    nY = Pti->nGetY();
// }

//! Adds the first integer parameter to the CGeom2DIPoint object's X co-ordinate, adds the second integer parameter to the CGeom2DIPoint object's Y co-ordinate
void CGeom2DIPoint::AddXAddY(int const nXToAdd, int const nYToAdd)
{
   nX += nXToAdd;
   nY += nYToAdd;
}

//! Adds the first double parameter (rounded) to the CGeom2DIPoint object's X co-ordinate, adds the second double parameter (rounded) to the CGeom2DIPoint object's Y co-ordinate
void CGeom2DIPoint::AddXAddY(double const dXToAdd, double const dYToAdd)
{
   nX += nRound(dXToAdd);
   nY += nRound(dYToAdd);
}

//! Divides the CGeom2DIPoint object's X co-ordinate by the first double parameter (rounded), divides the CGeom2DIPoint object's Y co-ordinate by the second double parameter (rounded)
void CGeom2DIPoint::DivXDivY(double const dXDiv, double const dYDiv)
{
   int
      nXDiv = nRound(dXDiv),
      nYDiv = nRound(dYDiv);

   // Check for zero division
   if (nXDiv != 0)
      nX /= nXDiv;

   if (nYDiv != 0)
      nY /= nYDiv;
}

//! Sets one CGeom2DIPoint object to be the same as another
void CGeom2DIPoint::operator= (CGeom2DIPoint const* pPti)
{
   nX = pPti->nGetX();
   nY = pPti->nGetY();
}

//! Compares two CGeom2DIPoint objects for equality
bool CGeom2DIPoint::operator== (CGeom2DIPoint const* pPti) const
{
   if ((pPti->nGetX() == nX) && (pPti->nGetY() == nY))
      return true;

   return false;
}

//! Compares two CGeom2DIPoint objects for equality
bool CGeom2DIPoint::operator== (CGeom2DIPoint Pti) const
{
   if ((Pti.nGetX() == nX) && (Pti.nGetY() == nY))
      return true;

   return false;
}

//! Compares two CGeom2DIPoint objects for inequality
bool CGeom2DIPoint::operator!= (CGeom2DIPoint const* pPti) const
{
   if ((pPti->nGetX() != nX) || (pPti->nGetY() != nY))
      return true;

   return false;
}

//! Compares two CGeom2DIPoint objects for inequality
bool CGeom2DIPoint::operator!= (CGeom2DIPoint Pti) const
{
   if ((Pti.nGetX() != nX) || (Pti.nGetY() != nY))
      return true;

   return false;
}
