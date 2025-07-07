/*!

   \class CCliffAlgorithmLegacy
   \brief Legacy cliff collapse algorithm wrapper
   \details This algorithm wraps the original CoastalME cliff collapse implementation for comparison
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file cliff_algorithm_legacy.h
   \brief Contains CCliffAlgorithmLegacy definitions

*/

#ifndef CLIFF_ALGORITHM_LEGACY_H
#define CLIFF_ALGORITHM_LEGACY_H
/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/

#include "cliff_algorithm.h"

//! Legacy cliff collapse algorithm wrapper
class CCliffAlgorithmLegacy : public CCliffAlgorithm
{
   private:
      //! Algorithm name
      string m_strAlgorithmName;
      
   public:
      //! Constructor
      CCliffAlgorithmLegacy();
      
      //! Destructor
      ~CCliffAlgorithmLegacy() override = default;
      
      //! Initialize the algorithm
      bool Initialize(const string& strAlgorithmName) override;
      
      //! Process a single timestep
      CCliffResults ProcessTimestep(const CCliffData& cliffData) override;
      
      //! Reset algorithm state
      void Reset() override;
      
      //! Get algorithm name
      string GetAlgorithmName() const override;
      
      //! Get algorithm description
      string GetAlgorithmDescription() const override;
};

#endif // CLIFF_ALGORITHM_LEGACY_H