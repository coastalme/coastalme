/*!

   \file cliff_algorithm_factory.cpp
   \brief CCliffAlgorithmFactory implementation
   \details Factory implementation for creating cliff collapse algorithms
   \author Wilf Chun
   \date 2025
   \copyright GNU General Public License

*/

/* ===============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

===============================================================================================================================*/

#include "cliff_algorithm_factory.h"
#include "cliff_algorithm_simple_notch.h"
#include "cliff_algorithm_legacy.h"

#include <algorithm>
#include <cctype>
#include <map>

using std::map;
using std::make_unique;
using std::transform;
using std::tolower;

//! Create a cliff algorithm instance by name
unique_ptr<CCliffAlgorithm> CCliffAlgorithmFactory::CreateAlgorithm(const string& strAlgorithmName)
{
   return CreateAlgorithm(GetAlgorithmEnum(strAlgorithmName));
}

//! Create a cliff algorithm instance by enum
unique_ptr<CCliffAlgorithm> CCliffAlgorithmFactory::CreateAlgorithm(ECliffAlgorithm algorithmType)
{
   switch (algorithmType)
   {
      case ECliffAlgorithm::SIMPLE_NOTCH:
         return make_unique<CCliffAlgorithmSimpleNotch>();
      
      case ECliffAlgorithm::LEGACY_ORIGINAL:
         return make_unique<CCliffAlgorithmLegacy>();
      
      default:
         // Default to simple notch algorithm
         return make_unique<CCliffAlgorithmSimpleNotch>();
   }
}

//! Get list of available algorithm names
vector<string> CCliffAlgorithmFactory::GetAvailableAlgorithms()
{
   return {"simple_notch", "legacy_original"};
}

//! Check if algorithm name is valid
bool CCliffAlgorithmFactory::IsValidAlgorithm(const string& strAlgorithmName)
{
   vector<string> availableAlgorithms = GetAvailableAlgorithms();
   string lowerName = strAlgorithmName;
   transform(lowerName.begin(), lowerName.end(), lowerName.begin(), [](unsigned char c){ return std::tolower(c); });
   
   for (const string& algorithm : availableAlgorithms)
   {
      if (algorithm == lowerName)
         return true;
   }
   return false;
}

//! Convert algorithm name to enum
ECliffAlgorithm CCliffAlgorithmFactory::GetAlgorithmEnum(const string& strAlgorithmName)
{
   string lowerName = strAlgorithmName;
   transform(lowerName.begin(), lowerName.end(), lowerName.begin(), [](unsigned char c){ return std::tolower(c); });
   
   if (lowerName == "simple_notch")
      return ECliffAlgorithm::SIMPLE_NOTCH;
   else if (lowerName == "legacy_original")
      return ECliffAlgorithm::LEGACY_ORIGINAL;
   else
      return ECliffAlgorithm::SIMPLE_NOTCH; // Default
}

//! Convert enum to algorithm name
string CCliffAlgorithmFactory::GetAlgorithmName(ECliffAlgorithm algorithmType)
{
   switch (algorithmType)
   {
      case ECliffAlgorithm::SIMPLE_NOTCH:
         return "simple_notch";
      case ECliffAlgorithm::LEGACY_ORIGINAL:
         return "legacy_original";
      default:
         return "simple_notch";
   }
}