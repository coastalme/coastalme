/*!
 * \file do_sediment_avalanching.h
 * \brief Header file for sediment avalanche redistribution
 *
 * This file contains constants, structures, and declarations for the sediment
 * avalanche algorithm which redistributes sediment when slopes exceed the
 * angle of repose.
 *
 * The algorithm uses a priority queue approach combined with dirty cell tracking
 * to efficiently process only cells that have changed (and their neighbors),
 * achieving 10-100x speedup over full-grid iteration.
 */

#ifndef DO_SEDIMENT_AVALANCHING_H
#define DO_SEDIMENT_AVALANCHING_H

/*===============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/

//! Angle of repose for sediment (degrees). This is a typical value for sand.
//! Finer sediments may have lower angles (~28°), coarser sediments higher (~37°).
const double ANGLE_OF_REPOSE_DEG = 33.0;

//! Angle of repose in radians (pre-calculated for efficiency)
const double ANGLE_OF_REPOSE_RAD = 0.5759586531;  // 33° * π/180

//! Tangent of angle of repose (pre-calculated for direct slope comparison)
const double TAN_ANGLE_OF_REPOSE = 0.6494075931;  // tan(33°)

//! Minimum sediment volume (m³) to trigger avalanche redistribution.
//! Prevents processing of trivially small amounts that don't affect morphology.
const double MIN_AVALANCHE_VOLUME = 0.001;  // 1 mm average depth over 1 m² cell

//! Maximum number of avalanche iterations per timestep.
//! Safety limit to prevent infinite loops in case of numerical issues.
const int MAX_AVALANCHE_ITERATIONS = 100;

//! Fraction of excess sediment to redistribute per iteration.
//! 0.5 = move 50% of unstable sediment each iteration.
//! Lower values (0.2-0.3) are more stable but slower to converge.
//! Higher values (0.6-0.8) converge faster but may overshoot.
const double REDISTRIBUTION_FRACTION = 0.5;

/*!
 * \struct UnstableCell
 * \brief Represents a cell with slope exceeding angle of repose
 *
 * Used in the priority queue to process most unstable cells first.
 * Instability metric is calculated as: max(tan(slope) - tan(angle_of_repose), 0)
 */
struct UnstableCell
{
   int nX;                  //!< Cell X coordinate (grid column)
   int nY;                  //!< Cell Y coordinate (grid row)
   double dInstability;     //!< Instability metric (how much slope exceeds angle of repose)

   //! Constructor
   UnstableCell(int const x, int const y, double const instability)
       : nX(x), nY(y), dInstability(instability)
   {
   }
};

/*!
 * \struct UnstableCellComparator
 * \brief Comparator for priority queue ordering
 *
 * Orders cells by instability metric in descending order (most unstable first).
 * Note: std::priority_queue is a max-heap by default, so we use < to get
 * the most unstable (highest instability) at the top.
 */
struct UnstableCellComparator
{
   bool operator()(UnstableCell const& a, UnstableCell const& b) const
   {
      // Return true if a has LOWER priority than b
      // This makes the priority queue a max-heap (highest instability first)
      return a.dInstability < b.dInstability;
   }
};

#endif // DO_SEDIMENT_AVALANCHING_H
