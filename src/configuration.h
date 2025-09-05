/*!

   \file configuration.h
   \brief Unified configuration class for CoastalME simulation parameters
   \details Provides a single interface for accessing simulation parameters regardless of input format (.dat or YAML)
   \author David Favis-Mortlock
   \author Andres Payo
   \date 2025
   \copyright GNU General Public License

*/

/* ==============================================================================================================================

   This file is part of CoastalME, the Coastal Modelling Environment.

   CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <vector>

using std::string;
using std::vector;

//! Unified configuration class for CoastalME simulation parameters
class CConfiguration
{
   private:
   // Run Information
   string m_strRunName;
   int m_nLogFileDetail;
   bool m_bCSVPerTimestepResults;

   // Simulation timing
   string m_strStartDateTime;
   string m_strDuration;
   string m_strTimestep;
   vector<string> m_vecSaveTimes;
   int m_nRandomSeed;
   bool m_bUseSystemTimeForSeed;

   // GIS Output
   int m_nMaxSaveDigits;
   string m_strSaveDigitsMode;
   vector<string> m_vecRasterFiles;
   string m_strRasterFormat;
   bool m_bWorldFile;
   bool m_bScaleValues;
   vector<double> m_vecSliceElevations;
   vector<string> m_vecVectorFiles;
   string m_strVectorFormat;
   vector<string> m_vecTimeSeriesFiles;

   // Grid and Coastline
   int m_nCoastlineSmoothing;
   int m_nCoastlineSmoothingWindow;
   int m_nPolynomialOrder;
   string m_strOmitGridEdges;
   int m_nProfileSmoothingWindow;
   double m_dMaxLocalSlope;
   double m_dMaxBeachElevation;

   // Layers and Files
   int m_nNumLayers;
   string m_strBasementDEMFile;
   vector<string> m_vecUnconsFineFiles;
   vector<string> m_vecUnconsSandFiles;
   vector<string> m_vecUnconsCoarseFiles;
   vector<string> m_vecConsFineFiles;
   vector<string> m_vecConsSandFiles;
   vector<string> m_vecConsCoarseFiles;
   string m_strSuspendedSedFile;
   string m_strLandformFile;
   string m_strInterventionClassFile;
   string m_strInterventionHeightFile;

   // Hydrology
   int m_nWavePropagationModel;
   double m_dSeawaterDensity;
   double m_dInitialWaterLevel;
   double m_dFinalWaterLevel;
   bool m_bHasFinalWaterLevel;
   double m_dDeepWaterWaveHeight;
   string m_strWaveHeightTimeSeries;
   double m_dDeepWaterWaveOrientation;
   double m_dWavePeriod;
   string m_strTideDataFile;
   double m_dBreakingWaveRatio;

   // Sediment and Erosion
   bool m_bCoastPlatformErosion;
   double m_dPlatformErosionResistance;
   bool m_bBeachSedimentTransport;
   int m_nBeachTransportAtEdges;
   int m_nBeachErosionEquation;
   double m_dFineMedianSize;
   double m_dSandMedianSize;
   double m_dCoarseMedianSize;
   double m_dSedimentDensity;
   double m_dBeachSedimentPorosity;
   double m_dFineErosivity;
   double m_dSandErosivity;
   double m_dCoarseErosivity;
   double m_dTransportKLS;
   double m_dKamphuis;
   double m_dBermHeight;

   // Cliff parameters
   bool m_bCliffCollapse;
   double m_dCliffErosionResistance;
   double m_dNotchOverhang;
   double m_dNotchBase;
   double m_dCliffDepositionA;
   double m_dTalusWidth;
   double m_dMinTalusLength;
   double m_dMinTalusHeight;

   // Flood parameters
   bool m_bFloodInput;
   string m_strFloodCoastline;
   string m_strRunupEquation;
   string m_strFloodLocations;
   string m_strFloodInputLocation;

   // Sediment input parameters
   bool m_bSedimentInput;
   string m_strSedimentInputLocation;
   string m_strSedimentInputType;
   string m_strSedimentInputDetails;

   // Physics and Geometry
   double m_dGravitationalAcceleration;
   double m_dNormalSpacing;
   double m_dRandomFactor;
   double m_dNormalLength;
   double m_dStartDepthRatio;

   // Profile and Output Options
   bool m_bSaveProfileData;
   vector<int> m_vecProfileNumbers;
   vector<int> m_vecProfileTimesteps;
   bool m_bSaveParallelProfiles;
   bool m_bOutputErosionPotential;
   int m_nCurvatureWindow;

   // Cliff Edge Processing
   int m_nCliffEdgeSmoothing;
   int m_nCliffEdgeSmoothingWindow;
   int m_nCliffEdgePolynomialOrder;
   double m_dCliffSlopeLimit;

   public:
   CConfiguration();
   ~CConfiguration();

   // Setters for all parameters
   void SetRunName(string const& str)
   {
      m_strRunName = str;
   }
   void SetLogFileDetail(int n)
   {
      m_nLogFileDetail = n;
   }
   void SetCSVPerTimestepResults(bool b)
   {
      m_bCSVPerTimestepResults = b;
   }
   void SetStartDateTime(string const& str)
   {
      m_strStartDateTime = str;
   }
   void SetDuration(string const& str)
   {
      m_strDuration = str;
   }
   void SetTimestep(string const& str)
   {
      m_strTimestep = str;
   }
   void SetSaveTimes(vector<string> const& vec)
   {
      m_vecSaveTimes = vec;
   }
   void SetRandomSeed(int n)
   {
      m_nRandomSeed = n;
      m_bUseSystemTimeForSeed = false;
   }
   void UseSystemTimeForSeed()
   {
      m_bUseSystemTimeForSeed = true;
   }

   void SetMaxSaveDigits(int n)
   {
      m_nMaxSaveDigits = n;
   }
   void SetSaveDigitsMode(string const& str)
   {
      m_strSaveDigitsMode = str;
   }
   void SetRasterFiles(vector<string> const& vec)
   {
      m_vecRasterFiles = vec;
   }
   void SetRasterFormat(string const& str)
   {
      m_strRasterFormat = str;
   }
   void SetWorldFile(bool b)
   {
      m_bWorldFile = b;
   }
   void SetScaleValues(bool b)
   {
      m_bScaleValues = b;
   }
   void SetSliceElevations(vector<double> const& vec)
   {
      m_vecSliceElevations = vec;
   }
   void SetVectorFiles(vector<string> const& vec)
   {
      m_vecVectorFiles = vec;
   }
   void SetVectorFormat(string const& str)
   {
      m_strVectorFormat = str;
   }
   void SetTimeSeriesFiles(vector<string> const& vec)
   {
      m_vecTimeSeriesFiles = vec;
   }

   void SetCoastlineSmoothing(int n)
   {
      m_nCoastlineSmoothing = n;
   }
   void SetCoastlineSmoothingWindow(int n)
   {
      m_nCoastlineSmoothingWindow = n;
   }
   void SetPolynomialOrder(int n)
   {
      m_nPolynomialOrder = n;
   }
   void SetOmitGridEdges(string const& str)
   {
      m_strOmitGridEdges = str;
   }
   void SetProfileSmoothingWindow(int n)
   {
      m_nProfileSmoothingWindow = n;
   }
   void SetMaxLocalSlope(double d)
   {
      m_dMaxLocalSlope = d;
   }
   void SetMaxBeachElevation(double d)
   {
      m_dMaxBeachElevation = d;
   }

   void SetNumLayers(int n)
   {
      m_nNumLayers = n;
   }
   void SetBasementDEMFile(string const& str)
   {
      m_strBasementDEMFile = str;
   }
   void SetUnconsFineFiles(vector<string> const& vec)
   {
      m_vecUnconsFineFiles = vec;
   }
   void SetUnconsSandFiles(vector<string> const& vec)
   {
      m_vecUnconsSandFiles = vec;
   }
   void SetUnconsCoarseFiles(vector<string> const& vec)
   {
      m_vecUnconsCoarseFiles = vec;
   }
   void SetConsFineFiles(vector<string> const& vec)
   {
      m_vecConsFineFiles = vec;
   }
   void SetConsSandFiles(vector<string> const& vec)
   {
      m_vecConsSandFiles = vec;
   }
   void SetConsCoarseFiles(vector<string> const& vec)
   {
      m_vecConsCoarseFiles = vec;
   }
   void SetSuspendedSedFile(string const& str)
   {
      m_strSuspendedSedFile = str;
   }
   void SetLandformFile(string const& str)
   {
      m_strLandformFile = str;
   }
   void SetInterventionClassFile(string const& str)
   {
      m_strInterventionClassFile = str;
   }
   void SetInterventionHeightFile(string const& str)
   {
      m_strInterventionHeightFile = str;
   }

   void SetWavePropagationModel(int n)
   {
      m_nWavePropagationModel = n;
   }
   void SetSeawaterDensity(double d)
   {
      m_dSeawaterDensity = d;
   }
   void SetInitialWaterLevel(double d)
   {
      m_dInitialWaterLevel = d;
   }
   void SetFinalWaterLevel(double d)
   {
      m_dFinalWaterLevel = d;
      m_bHasFinalWaterLevel = true;
   }
   void SetDeepWaterWaveHeight(double d)
   {
      m_dDeepWaterWaveHeight = d;
   }
   void SetWaveHeightTimeSeries(string const& str)
   {
      m_strWaveHeightTimeSeries = str;
   }
   void SetDeepWaterWaveOrientation(double d)
   {
      m_dDeepWaterWaveOrientation = d;
   }
   void SetWavePeriod(double d)
   {
      m_dWavePeriod = d;
   }
   void SetTideDataFile(string const& str)
   {
      m_strTideDataFile = str;
   }
   void SetBreakingWaveRatio(double d)
   {
      m_dBreakingWaveRatio = d;
   }

   // Additional setters for comprehensive YAML support
   void SetCoastPlatformErosion(bool b)
   {
      m_bCoastPlatformErosion = b;
   }
   void SetPlatformErosionResistance(double d)
   {
      m_dPlatformErosionResistance = d;
   }
   void SetBeachSedimentTransport(bool b)
   {
      m_bBeachSedimentTransport = b;
   }
   void SetBeachTransportAtEdges(int n)
   {
      m_nBeachTransportAtEdges = n;
   }
   void SetBeachErosionEquation(int n)
   {
      m_nBeachErosionEquation = n;
   }
   void SetFineMedianSize(double d)
   {
      m_dFineMedianSize = d;
   }
   void SetSandMedianSize(double d)
   {
      m_dSandMedianSize = d;
   }
   void SetCoarseMedianSize(double d)
   {
      m_dCoarseMedianSize = d;
   }
   void SetSedimentDensity(double d)
   {
      m_dSedimentDensity = d;
   }
   void SetBeachSedimentPorosity(double d)
   {
      m_dBeachSedimentPorosity = d;
   }
   void SetFineErosivity(double d)
   {
      m_dFineErosivity = d;
   }
   void SetSandErosivity(double d)
   {
      m_dSandErosivity = d;
   }
   void SetCoarseErosivity(double d)
   {
      m_dCoarseErosivity = d;
   }
   void SetTransportKLS(double d)
   {
      m_dTransportKLS = d;
   }
   void SetKamphuis(double d)
   {
      m_dKamphuis = d;
   }
   void SetBermHeight(double d)
   {
      m_dBermHeight = d;
   }

   void SetCliffCollapse(bool b)
   {
      m_bCliffCollapse = b;
   }
   void SetCliffErosionResistance(double d)
   {
      m_dCliffErosionResistance = d;
   }
   void SetNotchOverhang(double d)
   {
      m_dNotchOverhang = d;
   }
   void SetNotchBase(double d)
   {
      m_dNotchBase = d;
   }
   void SetCliffDepositionA(double d)
   {
      m_dCliffDepositionA = d;
   }
   void SetTalusWidth(double d)
   {
      m_dTalusWidth = d;
   }
   void SetMinTalusLength(double d)
   {
      m_dMinTalusLength = d;
   }
   void SetMinTalusHeight(double d)
   {
      m_dMinTalusHeight = d;
   }

   void SetFloodInput(bool b)
   {
      m_bFloodInput = b;
   }
   void SetFloodCoastline(string const& str)
   {
      m_strFloodCoastline = str;
   }
   void SetRunupEquation(string const& str)
   {
      m_strRunupEquation = str;
   }
   void SetFloodLocations(string const& str)
   {
      m_strFloodLocations = str;
   }
   void SetFloodInputLocation(string const& str)
   {
      m_strFloodInputLocation = str;
   }

   void SetSedimentInput(bool b)
   {
      m_bSedimentInput = b;
   }
   void SetSedimentInputLocation(string const& str)
   {
      m_strSedimentInputLocation = str;
   }
   void SetSedimentInputType(string const& str)
   {
      m_strSedimentInputType = str;
   }
   void SetSedimentInputDetails(string const& str)
   {
      m_strSedimentInputDetails = str;
   }

   void SetGravitationalAcceleration(double d)
   {
      m_dGravitationalAcceleration = d;
   }
   void SetNormalSpacing(double d)
   {
      m_dNormalSpacing = d;
   }
   void SetRandomFactor(double d)
   {
      m_dRandomFactor = d;
   }
   void SetNormalLength(double d)
   {
      m_dNormalLength = d;
   }
   void SetStartDepthRatio(double d)
   {
      m_dStartDepthRatio = d;
   }

   void SetSaveProfileData(bool b)
   {
      m_bSaveProfileData = b;
   }
   void SetProfileNumbers(vector<int> const& vec)
   {
      m_vecProfileNumbers = vec;
   }
   void SetProfileTimesteps(vector<int> const& vec)
   {
      m_vecProfileTimesteps = vec;
   }
   void SetSaveParallelProfiles(bool b)
   {
      m_bSaveParallelProfiles = b;
   }
   void SetOutputErosionPotential(bool b)
   {
      m_bOutputErosionPotential = b;
   }
   void SetCurvatureWindow(int n)
   {
      m_nCurvatureWindow = n;
   }

   void SetCliffEdgeSmoothing(int n)
   {
      m_nCliffEdgeSmoothing = n;
   }
   void SetCliffEdgeSmoothingWindow(int n)
   {
      m_nCliffEdgeSmoothingWindow = n;
   }
   void SetCliffEdgePolynomialOrder(int n)
   {
      m_nCliffEdgePolynomialOrder = n;
   }
   void SetCliffSlopeLimit(double d)
   {
      m_dCliffSlopeLimit = d;
   }

   // Getters for all parameters
   string GetRunName() const
   {
      return m_strRunName;
   }
   int GetLogFileDetail() const
   {
      return m_nLogFileDetail;
   }
   bool GetCSVPerTimestepResults() const
   {
      return m_bCSVPerTimestepResults;
   }
   string GetStartDateTime() const
   {
      return m_strStartDateTime;
   }
   string GetDuration() const
   {
      return m_strDuration;
   }
   string GetTimestep() const
   {
      return m_strTimestep;
   }
   vector<string> GetSaveTimes() const
   {
      return m_vecSaveTimes;
   }
   int GetRandomSeed() const
   {
      return m_nRandomSeed;
   }
   bool UseSystemTimeForRandomSeed() const
   {
      return m_bUseSystemTimeForSeed;
   }

   int GetMaxSaveDigits() const
   {
      return m_nMaxSaveDigits;
   }
   string GetSaveDigitsMode() const
   {
      return m_strSaveDigitsMode;
   }
   vector<string> GetRasterFiles() const
   {
      return m_vecRasterFiles;
   }
   string GetRasterFormat() const
   {
      return m_strRasterFormat;
   }
   bool GetWorldFile() const
   {
      return m_bWorldFile;
   }
   bool GetScaleValues() const
   {
      return m_bScaleValues;
   }
   vector<double> GetSliceElevations() const
   {
      return m_vecSliceElevations;
   }
   vector<string> GetVectorFiles() const
   {
      return m_vecVectorFiles;
   }
   string GetVectorFormat() const
   {
      return m_strVectorFormat;
   }
   vector<string> GetTimeSeriesFiles() const
   {
      return m_vecTimeSeriesFiles;
   }

   int GetCoastlineSmoothing() const
   {
      return m_nCoastlineSmoothing;
   }
   int GetCoastlineSmoothingWindow() const
   {
      return m_nCoastlineSmoothingWindow;
   }
   int GetPolynomialOrder() const
   {
      return m_nPolynomialOrder;
   }
   string GetOmitGridEdges() const
   {
      return m_strOmitGridEdges;
   }
   int GetProfileSmoothingWindow() const
   {
      return m_nProfileSmoothingWindow;
   }
   double GetMaxLocalSlope() const
   {
      return m_dMaxLocalSlope;
   }
   double GetMaxBeachElevation() const
   {
      return m_dMaxBeachElevation;
   }

   int GetNumLayers() const
   {
      return m_nNumLayers;
   }
   string GetBasementDEMFile() const
   {
      return m_strBasementDEMFile;
   }
   vector<string> GetUnconsFineFiles() const
   {
      return m_vecUnconsFineFiles;
   }
   vector<string> GetUnconsSandFiles() const
   {
      return m_vecUnconsSandFiles;
   }
   vector<string> GetUnconsCoarseFiles() const
   {
      return m_vecUnconsCoarseFiles;
   }
   vector<string> GetConsFineFiles() const
   {
      return m_vecConsFineFiles;
   }
   vector<string> GetConsSandFiles() const
   {
      return m_vecConsSandFiles;
   }
   vector<string> GetConsCoarseFiles() const
   {
      return m_vecConsCoarseFiles;
   }
   string GetSuspendedSedFile() const
   {
      return m_strSuspendedSedFile;
   }
   string GetLandformFile() const
   {
      return m_strLandformFile;
   }
   string GetInterventionClassFile() const
   {
      return m_strInterventionClassFile;
   }
   string GetInterventionHeightFile() const
   {
      return m_strInterventionHeightFile;
   }

   int GetWavePropagationModel() const
   {
      return m_nWavePropagationModel;
   }
   double GetSeawaterDensity() const
   {
      return m_dSeawaterDensity;
   }
   double GetInitialWaterLevel() const
   {
      return m_dInitialWaterLevel;
   }
   double GetFinalWaterLevel() const
   {
      return m_dFinalWaterLevel;
   }
   bool HasFinalWaterLevel() const
   {
      return m_bHasFinalWaterLevel;
   }
   double GetDeepWaterWaveHeight() const
   {
      return m_dDeepWaterWaveHeight;
   }
   string GetWaveHeightTimeSeries() const
   {
      return m_strWaveHeightTimeSeries;
   }
   double GetDeepWaterWaveOrientation() const
   {
      return m_dDeepWaterWaveOrientation;
   }
   double GetWavePeriod() const
   {
      return m_dWavePeriod;
   }
   string GetTideDataFile() const
   {
      return m_strTideDataFile;
   }
   double GetBreakingWaveRatio() const
   {
      return m_dBreakingWaveRatio;
   }

   // Initialize with default values
   void InitializeDefaults();
};

#endif      // CONFIGURATION_H
