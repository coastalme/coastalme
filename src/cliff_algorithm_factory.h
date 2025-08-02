/*!

   \class CCliffAlgorithmFactory
   \brief Factory class for creating cliff collapse algorithms
   \details This class provides a factory pattern for creating different cliff collapse algorithms
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License
   \file cliff_algorithm_factory.h
   \brief Contains CCliffAlgorithmFactory definitions

*/

#ifndef CLIFF_ALGORITHM_FACTORY_H
#define CLIFF_ALGORITHM_FACTORY_H
/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/

#include <memory>
#include <string>
#include <vector>
#include "cliff_algorithm.h"

using std::unique_ptr;
using std::string;
using std::vector;

//! Enumeration of available cliff collapse algorithms
enum class ECliffAlgorithm
{
   SIMPLE_NOTCH,        //! Simple notch-based algorithm (default)
   LEGACY_ORIGINAL      //! Original CoastalME algorithm for comparison
};

//! Factory class for creating cliff collapse algorithms
class CCliffAlgorithmFactory
{
   public:
      //! Create a cliff algorithm instance by name
      static unique_ptr<CCliffAlgorithm> CreateAlgorithm(const string& strAlgorithmName);
      
      //! Create a cliff algorithm instance by enum
      static unique_ptr<CCliffAlgorithm> CreateAlgorithm(ECliffAlgorithm algorithmType);
      
      //! Get list of available algorithm names
      static vector<string> GetAvailableAlgorithms();
      
      //! Check if algorithm name is valid
      static bool IsValidAlgorithm(const string& strAlgorithmName);
      
      //! Convert algorithm name to enum
      static ECliffAlgorithm GetAlgorithmEnum(const string& strAlgorithmName);
      
      //! Convert enum to algorithm name
      static string GetAlgorithmName(ECliffAlgorithm algorithmType);
};

#endif // CLIFF_ALGORITHM_FACTORY_H