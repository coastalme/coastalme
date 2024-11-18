/*!
\class CGeom2DPoint
\brief Geometry class used to represent 2D point objects with floating-point co-ordinates
\details The CGeom2DPoint geometry class is used to represent 2D points where the x and y co-ordinates are floating-point values, e.g. points for which the x and y co-ordinates are in the external CRS (co-ordinate reference system)
\author David Favis-Mortlock
\author Andres Payo
\date 2024
\copyright GNU General Public License

\file 2d_point.h
\brief Contains CGeom2DPoint definitions
*/

#ifndef C2DPOINT_H
#define C2DPOINT_H
/*===============================================================================================================================
This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
===============================================================================================================================*/
class CGeom2DPoint
{
private:
   //! The x co-ordinate
   double dX;
   
   //! The y co-ordinate
   double dY;

public:
   CGeom2DPoint(void);
   CGeom2DPoint(double const, double const);

   double dGetX(void) const;
   double dGetY(void) const;
   void SetX(double const);
   void SetY(double const);
//    void SetXY(double const, double const);
//    void SetXY(CGeom2DPoint const*);

   void operator= (CGeom2DPoint const*);
   bool operator== (CGeom2DPoint const*) const;
   bool operator== (CGeom2DPoint) const;
   bool operator!= (CGeom2DPoint const*) const;
};
#endif // C2DPOINT_H
