/*!
 *
 * \file interpolate.h
 * \brief Definitions of routines which return interpolated value at x from parallel arrays
 * \details TODO 001 A more detailed description of these routines.
 * \author Modified by David Favis-Mortlock and Andres Payo
 * \date 2025
 * \copyright GNU Lesser General Public License
 *
 */

#ifndef INTERPOLATE_H
#define INTERPOLATE_H
/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/
double dGetInterpolatedValue(vector<double> const*, vector<double> const*, double, bool);
double dGetInterpolatedValue(vector<int> const*, vector<double> const*, int, bool);
int nFindIndex(vector<double> const*, double const);
vector<double> VdInterpolateCShoreProfileOutput(vector<double> const*, vector<double> const*, vector<double> const*);
#endif // INTERPOLATE_H
