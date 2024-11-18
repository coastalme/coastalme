/*!
 *
 * \file read_input.cpp
 * \brief Reads non-GIS input files
 * \details TODO 001 A more detailed description of these routines.
 * \author David Favis-Mortlock
 * \author Andres Payo

 * \date 2024
 * \copyright GNU General Public License
 *
 */

/*==============================================================================================================================

This file is part of CoastalME, the Coastal Modelling Environment.

CoastalME is free software; you can redistribute it and/or modify it under the terms of the GNU General Public  License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

==============================================================================================================================*/
#include <stdlib.h> // for strtod()
#include <fstream>
using std::ifstream;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <string>
using std::to_string;

#include <algorithm>
using std::find;

#include "cme.h"
#include "simulation.h"
#include "sediment_input_event.h"

//===============================================================================================================================
//! The bReadIniFile member function reads the initialization file
//===============================================================================================================================
bool CSimulation::bReadIniFile(void)
{
   m_strCMEIni = m_strCMEDir;
   m_strCMEIni.append(CME_INI);

   // The .ini file is assumed to be in the CoastalME executable's directory
   string strFilePathName(m_strCMEIni);

   // Tell the user what is happening
   cout << READING_FILE_LOCATIONS << strFilePathName << endl;

   // Create an ifstream object
   ifstream InStream;

   // Try to open .ini file for input
   InStream.open(strFilePathName.c_str(), ios::in);

   // Did it open OK?
   if (! InStream.is_open())
   {
      // Error: cannot open .ini file for input
      cerr << ERR << "cannot open " << strFilePathName << " for input" << endl;
      return false;
   }

   int
      nLine = 0,
      i = 0;
   string strRec, strErr;

   while (getline(InStream, strRec))
   {
      nLine++;

      // Trim off leading and trailing whitespace
      strRec = strTrim(&strRec);

      // If it is a blank line or a comment then ignore it
      if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2))
      {
         // It isn't so increment counter
         i++;

         // Find the colon: note that lines MUST have a colon separating data from leading description portion
         size_t nPos = strRec.find(COLON);
         if (nPos == string::npos)
         {
            // Error: badly formatted (no colon)
            cerr << ERR << "on line " << nLine << ": badly formatted (no ':') in " << strFilePathName << endl << "'" << strRec << "'" << endl;
            return false;
         }

         if (nPos == strRec.size() - 1)
         {
            // Error: badly formatted (colon with nothing following)
            cerr << ERR << "on line " << nLine << ": badly formatted (nothing following ':') in " << strFilePathName << endl << "'" << strRec << "'" << endl;
            return false;
         }

         // Strip off leading portion (the bit up to and including the colon)
         string strRH = strRec.erase(0, nPos+1);

         // Remove leading whitespace
         strRH = strTrimLeft(&strRH);

         // Look for a trailing comment, if found then terminate string at that point and trim off any trailing whitespace
         nPos = strRH.rfind(QUOTE1);
         if (nPos != string::npos)
            strRH.resize(nPos);

         nPos = strRH.rfind(QUOTE2);
         if (nPos != string::npos)
            strRH.resize(nPos);

         // Remove trailing whitespace
         strRH = strTrimRight(&strRH);

         switch (i)
         {
         case 1:
            // The main input run-data filename
            if (strRH.empty())
               strErr = "line " + to_string(nLine) + ": path and name of main datafile";
            else
            {
               // First check that we don't already have an input run-data filename, e.g. one entered on the command-line
               if (m_strDataPathName.empty())
               {
                  // We don't: so first check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                  if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                     // It has an absolute path, so use it 'as is'
                     m_strDataPathName = strRH;
                  else
                  {
                     // It has a relative path, so prepend the CoastalME dir
                     m_strDataPathName = m_strCMEDir;
                     m_strDataPathName.append(strRH);
                  }
               }
            }
            break;

         case 2:
            // Path for CoastalME output
            if (strRH.empty())
               strErr = "line " + to_string(nLine) + ": path for CoastalME output";
            else
            {
               // Check for trailing slash on CoastalME output directory name (is vital)
               if (strRH[strRH.size() - 1] != PATH_SEPARATOR)
                  strRH.push_back(PATH_SEPARATOR);

               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                  // It is an absolute path, so use it 'as is'
                  m_strOutPath = strRH;
               else
               {
                  // It is a relative path, so prepend the CoastalME dir
                  m_strOutPath = m_strCMEDir;
                  m_strOutPath.append(strRH);
               }
            }
            break;

         case 3:
            // Email address, only useful if running under Linux/Unix
            if (! strRH.empty())
            {
               // Something was entered, do rudimentary check for valid email address
               if (strRH.find('@') == string::npos)
                  strErr = "line " + to_string(nLine) + ": email address for messages";
               else
                  m_strMailAddress = strRH;
            }
            break;
         }

         // Did an error occur?
         if (! strErr.empty())
         {
            // Error in input to initialisation file
            cerr << ERR << "reading " << strErr << " in " << strFilePathName << endl
                 << "'" << strRec << "'" << endl;
            InStream.close();

            return false;
         }
      }
   }

   InStream.close();
   return true;
}

//===============================================================================================================================
//! Reads the run details input file and does some initialization
//===============================================================================================================================
bool CSimulation::bReadRunDataFile(void)
{
   // Create an ifstream object
   ifstream InStream;

   // Try to open run details file for input
   InStream.open(m_strDataPathName.c_str(), ios::in);

   // Did it open OK?
   if (!InStream.is_open())
   {
      // Error: cannot open run details file for input
      cerr << ERR << "cannot open " << m_strDataPathName << " for input" << endl;
      return false;
   }

   int nLine = 0;
   int i = 0;
   size_t nPos;
   string strRec, strErr;

   while (getline(InStream, strRec))
   {
      nLine++;

      // Trim off leading and trailing whitespace
      strRec = strTrim(&strRec);

      // If it is a blank line or a comment then ignore it
      if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2))
      {
         // It isn't so increment counter
         i++;

         // Find the colon: note that lines MUST have a colon separating data from leading description portion
         nPos = strRec.find(COLON);
         if (nPos == string::npos)
         {
            // Error: badly formatted (no colon)
            cerr << ERR << "on line " << to_string(nLine) << "badly formatted (no ':') in " << m_strDataPathName << endl
                 << strRec << endl;
            return false;
         }

         // Strip off leading portion (the bit up to and including the colon)
         string strRH = strRec.erase(0, nPos+1);

         // Remove leading whitespace after the colon
         strRH = strTrimLeft(&strRH);

         // Look for trailing comments, if found then terminate string at that point and trim off any trailing whitespace
         bool bFound = true;
         while (bFound)
         {
            bFound = false;

            nPos = strRH.rfind(QUOTE1);
            if (nPos != string::npos)
            {
               strRH.resize(nPos);
               bFound = true;
            }

            nPos = strRH.rfind(QUOTE2);
            if (nPos != string::npos)
            {
               strRH.resize(nPos);
               bFound = true;
            }

            // Trim trailing spaces
            strRH = strTrimRight(&strRH);
         }

#ifdef _WIN32
         // For Windows, make sure has backslashes, not Unix-style slashes
         strRH = pstrChangeToBackslash(&strRH);
#endif
         bool
             bFirst = true;
         int
             nRet = 0,
             nHour = 0,
             nMin = 0,
             nSec = 0,
             nDay = 0,
             nMonth = 0,
             nYear = 0;
         double dMult = 0;
         string strTmp;
         vector<string> VstrTmp;

         switch (i)
         {
         // ---------------------------------------------- Run Information -----------------------------------------------------
         case 1:
            // Text output file names, don't change case
            if (strRH.empty())
               strErr = "line " + to_string(nLine) + ": output file names";
            else
            {
               m_strRunName = strRH;

               m_strOutFile = m_strOutPath;
               m_strOutFile.append(strRH);
               m_strOutFile.append(OUTEXT);

               m_strLogFile = m_strOutPath;
               m_strLogFile.append(strRH);
               m_strLogFile.append(LOGEXT);
            }
            break;

         case 2:
            // Content of log file, 0 = no log file, 1 = least detail, 3 = most detail
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for log file detail level '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nLogFileDetail = stoi(strRH);

            if ((m_nLogFileDetail < NO_LOG_FILE) || (m_nLogFileDetail > LOG_FILE_ALL))
               strErr = "line " + to_string(nLine) + ": log file detail level";
            break;

         case 3:
            // Get the start date/time of the simulation, format is [hh-mm-ss dd/mm/yyyy]
            VstrTmp = VstrSplit(&strRH, SPACE);

            // Both date and time here?
            if (VstrTmp.size() < 2)
            {
               strErr = "line " + to_string(nLine) + ": must have both date and time for simulation start in '" + m_strDataPathName + "'";
               break;
            }

            // OK, first sort out the time
            if (! bParseTime(&VstrTmp[0], nHour, nMin, nSec))
            {
               strErr = "line " + to_string(nLine) + ": could not understand simulation start time in '" + m_strDataPathName + "'";
               break;
            }

            // Next sort out the date
            if (! bParseDate(&VstrTmp[1], nDay, nMonth, nYear))
            {
               strErr = "line " + to_string(nLine) + ": could not understand simulation start date in '" + m_strDataPathName + "'";
               break;
            }

            // Store simulation start time and date
            m_nSimStartSec = nSec;
            m_nSimStartMin = nMin;
            m_nSimStartHour = nHour;
            m_nSimStartDay = nDay;
            m_nSimStartMonth = nMonth;
            m_nSimStartYear = nYear;

            break;

         case 4:
            // Duration of simulation (in hours, days, months, or years): sort out multiplier and user units, as used in the per-timestep output
            strRH = strToLower(&strRH);

            nRet = nDoSimulationTimeMultiplier(&strRH);
            if (nRet != RTN_OK)
            {
               strErr = "line " + to_string(nLine) + ": units for duration of simulation";
               break;
            }

            // And now calculate the duration of the simulation in hours: first find whitespace between the number and the unit
            nPos = strRH.rfind(SPACE);
            if (nPos == string::npos)
            {
               strErr = "line " + to_string(nLine) + ": format of duration simulation line";
               break;
            }

            // Cut off rh bit of string
            strRH.resize(nPos);

            // Remove trailing spaces
            strRH = strTrimRight(&strRH);

            // Calculate the duration of the simulation in hours
            m_dSimDuration = strtod(strRH.c_str(), NULL) * m_dDurationUnitsMult;

            if (m_dSimDuration <= 0)
               strErr = "line " + to_string(nLine) + ": duration of simulation must be > 0";
            break;

         case 5:
            // Timestep of simulation (in hours or days)
            strRH = strToLower(&strRH);

            dMult = dGetTimeMultiplier(&strRH);
            if (static_cast<int>(dMult) == TIME_UNKNOWN)
            {
               strErr = "line " + to_string(nLine) + ": units for simulation timestep";
               break;
            }

            // we have the multiplier, now calculate the timestep in hours: look for the whitespace between the number and unit
            nPos = strRH.rfind(SPACE);
            if (nPos == string::npos)
            {
               strErr = "line " + to_string(nLine) + ": format of simulation timestep";
               break;
            }

            // cut off rh bit of string
            strRH.resize(nPos);

            // remove trailing spaces
            strRH = strTrimRight(&strRH);

            // Check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for timestep '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dTimeStep = strtod(strRH.c_str(), NULL) * dMult; // in hours

            if (m_dTimeStep <= 0)
               strErr = "line " + to_string(nLine) + ": timestep of simulation must be > 0";

            if (m_dTimeStep >= m_dSimDuration)
               strErr = "line " + to_string(nLine) + ": timestep of simulation must be < the duration of the simulation";

            break;

         case 6:
            // Save interval(s)
            strRH = strToLower(&strRH);

            // First get the multiplier
            dMult = dGetTimeMultiplier(&strRH);
            if (static_cast<int>(dMult) == TIME_UNKNOWN)
            {
               strErr = "line " + to_string(nLine) + ": units for save intervals";
               break;
            }

            // Search for a space from the right-hand end of the string
            nPos = strRH.rfind(SPACE);
            if (nPos == string::npos)
            {
               strErr = "line " + to_string(nLine) + ": format of save times/intervals";
               break;
            }

            // Cut off rh bit of string
            strRH.resize(nPos);

            // Remove trailing spaces
            strRH = strTrimRight(&strRH);

            // Now find out whether we're dealing with a single regular save interval or a vector of save times: again check for a space
            nPos = strRH.find(SPACE);
            if (nPos != string::npos)
            {
               // There's another space, so we must have more than one number
               m_bSaveRegular = false;

               strRH += SPACE;

               do
               {
                  // Put the number into the array
                  if (m_nUSave > static_cast<int>(SAVEMAX) - 1)
                  {
                     strErr = "line " + to_string(nLine) + ": too many save intervals";
                     break;
                  }
                  m_dUSaveTime[m_nUSave++] = strtod(strRH.substr(0, nPos).c_str(), NULL) * dMult; // convert to hours

                  // Trim off the number and remove leading whitespace
                  strRH.erase(0, nPos);
                  strRH = strTrimLeft(&strRH);

                  // Now look for another space
                  nPos = strRH.find(SPACE);
               }
               //               while (nPos != string::npos);
               while (strRH.size() > 1);

               // Check that we have at least 2 numbers
               if (m_nUSave < 1)
               {
                  strErr = "line " + to_string(nLine) + ": must have at least two save times";
                  break;
               }

               if (m_dUSaveTime[0] < m_dTimeStep)
               {
                  strErr = "line " + to_string(nLine) + ": first save time cannot be less than timestep";
                  break;
               }

               // Put a dummy save interval as the last entry in the array: this is needed to stop problems at end of run
               m_dUSaveTime[m_nUSave] = m_dSimDuration + 1;
            }
            else
            {
               // There isn't a space, so we have only one number
               m_bSaveRegular = true;
               m_dRegularSaveInterval = strtod(strRH.c_str(), NULL) * dMult;        // Convert to hours
               if (m_dRegularSaveInterval < m_dTimeStep)
                  strErr = "line " + to_string(nLine) + ": save interval cannot be less than timestep";
               else
                  // Set up for first save
                  m_dRegularSaveTime = m_dRegularSaveInterval;
            }
            break;

         case 7:
            // Random number seed(s)
            m_ulRandSeed[0] = atol(strRH.c_str());
            if (0 == m_ulRandSeed[0])
            {
               strErr = "line " + to_string(nLine) + ": random number seed";
               break;
            }

            // TODO 039 Rewrite this, similar to reading raster slice elevations
            // OK, so find out whether we're dealing with a single seed or more than one: check for a space
            nPos = strRH.find(SPACE);
            if (nPos != string::npos)
            {
               // There's a space, so we must have more than one number
               int n = 0;
               do
               {
                  // Trim off the part before the first space then remove leading whitespace
                  strRH.erase(0, nPos);
                  strRH = strTrimLeft(&strRH);

                  // Put the number into the array
                  m_ulRandSeed[++n] = atol(strRH.c_str());

                  // Now look for another space
                  nPos = strRH.find(SPACE);
               } while ((n < NRNG) && (nPos != string::npos));
            }
            else
            {
               // Only one seed specified, so make all seeds the same
               for (int n = 1; n < NRNG; n++)
                  m_ulRandSeed[n] = m_ulRandSeed[n - 1];
            }
            break;

         case 8:
            // Max save digits for GIS output file names
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for max save digits for GIS output file names '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nGISMaxSaveDigits = stoi(strRH);

            if (m_nGISMaxSaveDigits < 2)
               strErr = "line " + to_string(nLine) + ": max save digits for GIS output file names must be > 1";
            break;

         case 9:
            // Save digits for GIS output sequential or iteration number?    [s = sequential, i = iteration]: s
            if (strRH.empty())
               strErr = "line " + to_string(nLine) + ": must specify save digits for GIS output as sequential or as iteration number";
            else
            {
               // Convert to lower case
               strRH = strToLower(&strRH);

               if (strRH.find('s') != string::npos)
               {
                  // First look for 's'
                  m_bGISSaveDigitsSequential = true;
               }
               else if (strRH.find('i') != string::npos)
               {
                  // Now look for 'i'
                  m_bGISSaveDigitsSequential = false;
               }
               else
               {
                  strErr = "line " + to_string(nLine) + ": invalid code for save digits for GIS output save number (must be s or i)";
                  break;
               }
            }
            break;


         case 10:
            // Raster GIS files to output
            if (strRH.empty())
               strErr = "line " + to_string(nLine) + ": must contain '" + RASTER_ALL_OUTPUT_CODE +"', or '" + RASTER_USUAL_OUTPUT_CODE + "', or at least one raster GIS output code";
            else
            {
               // Convert to lower case
               strRH = strToLower(&strRH);

               if (strRH.find(RASTER_ALL_OUTPUT_CODE) != string::npos)
               {
                  // Set switches for all GIS raster output. Some of these (e.g. all relating to fine sediment) are ignored if e.g. no fine sediment layers are read in
                  m_bSuspSedSave =
                  m_bAvgSuspSedSave =
                  m_bFineUnconsSedSave =
                  m_bFineConsSedSave =
                  m_bSandUnconsSedSave =
                  m_bSandConsSedSave =
                  m_bCoarseUnconsSedSave =
                  m_bCoarseConsSedSave =
                  m_bSedimentTopSurfSave =
                  m_bTopSurfSave =
                  m_bSeaDepthSave =
                  m_bWaveHeightSave =
                  m_bWaveAngleSave =
                  m_bPotentialPlatformErosionSave =
                  m_bActualPlatformErosionSave =
                  m_bTotalPotentialPlatformErosionSave =
                  m_bTotalActualPlatformErosionSave =
                  m_bPotentialBeachErosionSave =
                  m_bActualBeachErosionSave =
                  m_bTotalPotentialBeachErosionSave =
                  m_bTotalActualBeachErosionSave =
                  m_bBeachDepositionSave =
                  m_bTotalBeachDepositionSave =
                  m_bLandformSave =
                  m_bLocalSlopeSave =
                  m_bAvgSeaDepthSave =
                  m_bAvgWaveHeightSave =
                  m_bAvgWaveAngleSave =
                  m_bBeachProtectionSave =
                  m_bBasementElevSave =
                  m_bRasterCoastlineSave =
                  m_bRasterNormalSave =
                  m_bActiveZoneSave =
                  m_bCliffCollapseSave =
                  m_bTotCliffCollapseSave =
                  m_bCliffCollapseDepositionSave =
                  m_bTotCliffCollapseDepositionSave =
                  m_bRasterPolygonSave =
                  m_bPotentialPlatformErosionMaskSave =
                  m_bSeaMaskSave =
                  m_bBeachMaskSave =
                  m_bShadowZoneCodesSave =
                  m_bDeepWaterWaveAngleSave =
                  m_bDeepWaterWaveHeightSave =
                  m_bDeepWaterWavePeriodSave =
                  m_bPolygonUnconsSedUpOrDownDriftSave =
                  m_bPolygonUnconsSedGainOrLossSave = true;
               }
               else if (strRH.find(RASTER_USUAL_OUTPUT_CODE) != string::npos)
               {
                  // Set switches for usual GIS raster output. Again, some of these (e.g. all relating to fine sediment) are ignored if they are irrelevant
                  m_bSuspSedSave =
                  m_bAvgSuspSedSave =
                  m_bFineUnconsSedSave =
                  m_bFineConsSedSave =
                  m_bSandUnconsSedSave =
                  m_bSandConsSedSave =
                  m_bCoarseUnconsSedSave =
                  m_bCoarseConsSedSave =
                  m_bSedimentTopSurfSave =
                  m_bTopSurfSave =
                  m_bSeaDepthSave =
                  m_bWaveHeightSave =
                  m_bWaveAngleSave =
                  m_bPotentialPlatformErosionSave =
                  m_bActualPlatformErosionSave =
                  m_bTotalPotentialPlatformErosionSave =
                  m_bTotalActualPlatformErosionSave =
                  m_bPotentialBeachErosionSave =
                  m_bActualBeachErosionSave =
                  m_bTotalPotentialBeachErosionSave =
                  m_bTotalActualBeachErosionSave =
                  m_bBeachDepositionSave =
                  m_bTotalBeachDepositionSave =
                  m_bLandformSave =
                  m_bLocalSlopeSave =
                  m_bAvgWaveHeightSave =
                  m_bAvgWaveAngleSave =
                  m_bBeachProtectionSave =
                  m_bBasementElevSave =
                  m_bActiveZoneSave =
                  m_bCliffCollapseSave =
                  m_bTotCliffCollapseSave =
                  m_bCliffCollapseDepositionSave =
                  m_bTotCliffCollapseDepositionSave =
                  m_bRasterPolygonSave =
                  m_bShadowZoneCodesSave =
                  m_bPolygonUnconsSedUpOrDownDriftSave =
                  m_bPolygonUnconsSedGainOrLossSave = true;
               }
               else
               {
                  // We are not outputting either the "usual" or the "all" collections of raster GIS files, so set switches (and remove strings) for those optional files for which the user specified the code
                  if (strRH.find(RASTER_SEDIMENT_TOP_CODE) != string::npos)
                  {
                     m_bSedimentTopSurfSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_SEDIMENT_TOP_CODE);
                  }

                  if (strRH.find(RASTER_TOP_CODE) != string::npos)
                  {
                     m_bTopSurfSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOP_CODE);
                  }

                  if (strRH.find(RASTER_SEA_DEPTH_CODE) != string::npos)
                  {
                     m_bSeaDepthSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_SEA_DEPTH_CODE);
                  }

                  if (strRH.find(RASTER_WAVE_HEIGHT_CODE) != string::npos)
                  {
                     m_bWaveHeightSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_WAVE_HEIGHT_CODE);
                  }

                  if (strRH.find(RASTER_WAVE_ORIENTATION_CODE) != string::npos)
                  {
                     m_bWaveAngleSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_WAVE_ORIENTATION_CODE);
                  }

                  if (strRH.find(RASTER_WAVE_PERIOD_CODE) != string::npos)
                  {
                     m_bDeepWaterWavePeriodSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_WAVE_PERIOD_CODE);
                  }

                  if (strRH.find(RASTER_POTENTIAL_PLATFORM_EROSION_CODE) != string::npos)
                  {
                     m_bPotentialPlatformErosionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_POTENTIAL_PLATFORM_EROSION_CODE);
                  }

                  if (strRH.find(RASTER_ACTUAL_PLATFORM_EROSION_CODE) != string::npos)
                  {
                     m_bActualPlatformErosionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_ACTUAL_PLATFORM_EROSION_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_POTENTIAL_PLATFORM_EROSION_CODE) != string::npos)
                  {
                     m_bTotalPotentialPlatformErosionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_POTENTIAL_PLATFORM_EROSION_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_ACTUAL_PLATFORM_EROSION_CODE) != string::npos)
                  {
                     m_bTotalActualPlatformErosionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_ACTUAL_PLATFORM_EROSION_CODE);
                  }

                  if (strRH.find(RASTER_POTENTIAL_BEACH_EROSION_CODE) != string::npos)
                  {
                     m_bPotentialBeachErosionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_POTENTIAL_BEACH_EROSION_CODE);
                  }

                  if (strRH.find(RASTER_ACTUAL_BEACH_EROSION_CODE) != string::npos)
                  {
                     m_bActualBeachErosionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_ACTUAL_BEACH_EROSION_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_POTENTIAL_BEACH_EROSION_CODE) != string::npos)
                  {
                     m_bTotalPotentialBeachErosionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_POTENTIAL_BEACH_EROSION_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_ACTUAL_BEACH_EROSION_CODE) != string::npos)
                  {
                     m_bTotalActualBeachErosionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_ACTUAL_BEACH_EROSION_CODE);
                  }

                  if (strRH.find(RASTER_LANDFORM_CODE) != string::npos)
                  {
                     m_bLandformSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_LANDFORM_CODE);
                  }

                  if (strRH.find(RASTER_LOCAL_SLOPE_CODE) != string::npos)
                  {
                     m_bLocalSlopeSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_LOCAL_SLOPE_CODE);
                  }

                  if (strRH.find(RASTER_AVG_SEA_DEPTH_CODE) != string::npos)
                  {
                     m_bAvgSeaDepthSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_AVG_SEA_DEPTH_CODE);
                  }

                  if (strRH.find(RASTER_AVG_WAVE_HEIGHT_CODE) != string::npos)
                  {
                     m_bAvgWaveHeightSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_AVG_WAVE_HEIGHT_CODE);
                  }

                  if (strRH.find(RASTER_AVG_WAVE_ORIENTATION_CODE) != string::npos)
                  {
                     m_bAvgWaveAngleSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_AVG_WAVE_ORIENTATION_CODE);
                  }

                  if (strRH.find(RASTER_BEACH_PROTECTION_CODE) != string::npos)
                  {
                     m_bBeachProtectionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_BEACH_PROTECTION_CODE);
                  }

                  if (strRH.find(RASTER_BASEMENT_ELEVATION_CODE) != string::npos)
                  {
                     m_bBasementElevSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_BASEMENT_ELEVATION_CODE);
                  }

                  if (strRH.find(RASTER_SUSP_SED_CODE) != string::npos)
                  {
                     m_bSuspSedSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_SUSP_SED_CODE);
                  }

                  if (strRH.find(RASTER_AVG_SUSP_SED_CODE) != string::npos)
                  {
                     m_bAvgSuspSedSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_AVG_SUSP_SED_CODE);
                  }

                  if (strRH.find(RASTER_FINE_UNCONS_CODE) != string::npos)
                  {
                     m_bFineUnconsSedSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_FINE_UNCONS_CODE);
                  }

                  if (strRH.find(RASTER_SAND_UNCONS_CODE) != string::npos)
                  {
                     m_bSandUnconsSedSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_SAND_UNCONS_CODE);
                  }

                  if (strRH.find(RASTER_COARSE_UNCONS_CODE) != string::npos)
                  {
                     m_bCoarseUnconsSedSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_COARSE_UNCONS_CODE);
                  }

                  if (strRH.find(RASTER_FINE_CONS_CODE) != string::npos)
                  {
                     m_bFineConsSedSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_FINE_CONS_CODE);
                  }

                  if (strRH.find(RASTER_SAND_CONS_CODE) != string::npos)
                  {
                     m_bSandConsSedSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_SAND_CONS_CODE);
                  }

                  if (strRH.find(RASTER_COARSE_CONS_CODE) != string::npos)
                  {
                     m_bCoarseConsSedSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_COARSE_CONS_CODE);
                  }

                  if (strRH.find(RASTER_COAST_CODE) != string::npos)
                  {
                     m_bRasterCoastlineSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_COAST_CODE);
                  }

                  if (strRH.find(RASTER_COAST_NORMAL_CODE) != string::npos)
                  {
                     m_bRasterNormalSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_COAST_NORMAL_CODE);
                  }

                  if (strRH.find(RASTER_ACTIVE_ZONE_CODE) != string::npos)
                  {
                     m_bActiveZoneSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_ACTIVE_ZONE_CODE);
                  }

                  if (strRH.find(RASTER_CLIFF_COLLAPSE_EROSION_FINE_CODE) != string::npos)
                  {
                     m_bCliffCollapseSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_CLIFF_COLLAPSE_EROSION_FINE_CODE);
                  }

                  if (strRH.find(RASTER_CLIFF_COLLAPSE_EROSION_SAND_CODE) != string::npos)
                  {
                     m_bCliffCollapseSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_CLIFF_COLLAPSE_EROSION_SAND_CODE);
                  }
                  
                  if (strRH.find(RASTER_CLIFF_COLLAPSE_EROSION_COARSE_CODE) != string::npos)
                  {
                     m_bCliffCollapseSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_CLIFF_COLLAPSE_EROSION_COARSE_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_FINE_CODE) != string::npos)
                  {
                     m_bTotCliffCollapseSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_FINE_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_SAND_CODE) != string::npos)
                  {
                     m_bTotCliffCollapseSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_SAND_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE_CODE) != string::npos)
                  {
                     m_bTotCliffCollapseSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_CLIFF_COLLAPSE_EROSION_COARSE_CODE);
                  }

                  if (strRH.find(RASTER_CLIFF_COLLAPSE_DEPOSITION_SAND_CODE) != string::npos)
                  {
                     m_bCliffCollapseDepositionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_CLIFF_COLLAPSE_DEPOSITION_SAND_CODE);
                  }

                  if (strRH.find(RASTER_CLIFF_COLLAPSE_DEPOSITION_COARSE_CODE) != string::npos)
                  {
                     m_bCliffCollapseDepositionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_CLIFF_COLLAPSE_DEPOSITION_COARSE_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND_CODE) != string::npos)
                  {
                     m_bTotCliffCollapseDepositionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_SAND_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE_CODE) != string::npos)
                  {
                     m_bTotCliffCollapseDepositionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_CLIFF_COLLAPSE_DEPOSITION_COARSE_CODE);
                  }

                  if (strRH.find(RASTER_POLYGON_CODE) != string::npos)
                  {
                     m_bRasterPolygonSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_POLYGON_CODE);
                  }

                  if (strRH.find(RASTER_POTENTIAL_PLATFORM_EROSION_MASK_CODE) != string::npos)
                  {
                     m_bPotentialPlatformErosionMaskSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_POTENTIAL_PLATFORM_EROSION_MASK_CODE);
                  }

                  if (strRH.find(RASTER_INUNDATION_MASK_CODE) != string::npos)
                  {
                     m_bSeaMaskSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_INUNDATION_MASK_CODE);
                  }

                  if (strRH.find(RASTER_BEACH_MASK_CODE) != string::npos)
                  {
                     m_bBeachMaskSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_BEACH_MASK_CODE);
                  }

                  if (strRH.find(RASTER_INTERVENTION_CLASS_CODE) != string::npos)
                  {
                     m_bInterventionClassSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_INTERVENTION_CLASS_CODE);
                  }

                  if (strRH.find(RASTER_INTERVENTION_HEIGHT_CODE) != string::npos)
                  {
                     m_bInterventionHeightSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_INTERVENTION_HEIGHT_CODE);
                  }

                  if (strRH.find(RASTER_SHADOW_ZONE_CODE) != string::npos)
                  {
                     m_bShadowZoneCodesSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_SHADOW_ZONE_CODE);
                  }

                  if (strRH.find(RASTER_DEEP_WATER_WAVE_ORIENTATION_CODE) != string::npos)
                  {
                     m_bDeepWaterWaveAngleSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_DEEP_WATER_WAVE_ORIENTATION_CODE);
                  }

                  if (strRH.find(RASTER_DEEP_WATER_WAVE_HEIGHT_CODE) != string::npos)
                  {
                     m_bDeepWaterWaveHeightSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_DEEP_WATER_WAVE_HEIGHT_CODE);
                  }

                  if (strRH.find(RASTER_DEEP_WATER_WAVE_PERIOD_CODE) != string::npos)
                  {
                     m_bDeepWaterWavePeriodSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_DEEP_WATER_WAVE_PERIOD_CODE);
                  }

                  if (strRH.find(RASTER_POLYGON_UPDRIFT_OR_DOWNDRIFT_CODE) != string::npos)
                  {
                     m_bPolygonUnconsSedUpOrDownDriftSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_POLYGON_UPDRIFT_OR_DOWNDRIFT_CODE);
                  }

                  if (strRH.find(RASTER_POLYGON_GAIN_OR_LOSS_CODE) != string::npos)
                  {
                     m_bPolygonUnconsSedGainOrLossSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_POLYGON_GAIN_OR_LOSS_CODE);
                  }

                  if (strRH.find(RASTER_BEACH_DEPOSITION_CODE) != string::npos)
                  {
                     m_bBeachDepositionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_BEACH_DEPOSITION_CODE);
                  }

                  if (strRH.find(RASTER_TOTAL_BEACH_DEPOSITION_CODE) != string::npos)
                  {
                     m_bTotalBeachDepositionSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_TOTAL_BEACH_DEPOSITION_CODE);
                  }

                  if (strRH.find(RASTER_SEDIMENT_INPUT_EVENT_CODE) != string::npos)
                  {
                     m_bSedimentInputEventSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_SEDIMENT_INPUT_EVENT_CODE);
                  }

                  if (strRH.find(RASTER_SETUP_SURGE_FLOOD_MASK_CODE) != string::npos)
                  {
                     m_bSetupSurgeFloodMaskSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_SETUP_SURGE_FLOOD_MASK_CODE);
                  }
                  if (strRH.find(RASTER_SETUP_SURGE_RUNUP_FLOOD_MASK_CODE) != string::npos)
                  {
                     m_bSetupSurgeRunupFloodMaskSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_SETUP_SURGE_RUNUP_FLOOD_MASK_CODE);
                  }
                  if (strRH.find(RASTER_WAVE_FLOOD_LINE_CODE) != string::npos)
                  {
                     m_bRasterWaveFloodLineSave = true;
                     strRH = strRemoveSubstr(&strRH, &RASTER_WAVE_FLOOD_LINE_CODE);
                  }

                  // Check to see if all codes have been removed
                  if (! strRH.empty())
                     strErr = "line " + to_string(nLine) + ": unknown code '" + strRH + "' in list of codes for raster GIS output";
               }
            }
            break;

         case 11:
            // Raster GIS output format (note must retain original case). Blank means use same format as input DEM file (if possible)
            m_strRasterGISOutFormat = strTrimLeft(&strRH);
                                    
            // TODO 065 Remove this when GDAL gpkg output is working correctly
            if (strRH.find("gpkg") != string::npos)
               strErr = "GDAL gpkg create() gives floating-point exception. Please choose another output format.";

            break;

         case 12:
            // If needed, scale GIS raster output values
            strRH = strToLower(&strRH);

            m_bScaleRasterOutput = false;
            if (strRH.find("y") != string::npos)
               m_bScaleRasterOutput = true;
            break;

         case 13:
            // If needed, also output GIS raster world file
            strRH = strToLower(&strRH);

            m_bWorldFile = false;
            if (strRH.find("y") != string::npos)
               m_bWorldFile = true;
            break;

         case 14:
            // Elevations for raster slice output, if desired
            if (! strRH.empty())
            {
               m_bSliceSave = true;

               // OK, so find out whether we're dealing with a single seed or more than one: check for a space
               nPos = strRH.find(SPACE);
               if (nPos != string::npos)
               {
                  // There's a space, so we must have more than one number
                  do
                  {
                     // Get LH bit
                     strTmp = strRH.substr(0, nPos);
                     m_VdSliceElev.push_back(strtod(strTmp.c_str(), NULL));

                     // Get the RH bit
                     strRH = strRH.substr(nPos, strRH.size() - nPos);
                     strRH = strTrimLeft(&strRH);

                     // Now look for another space
                     nPos = strRH.find(SPACE);
                  } while (nPos != string::npos);
               }
               // Read either the single number, or the left-over number
               m_VdSliceElev.push_back(strtod(strTmp.c_str(), NULL));
            }
            break;

         case 15:
            // Vector GIS files to output
            if (strRH.empty())
               strErr = "line " + to_string(nLine) + ": must contain '" + VECTOR_ALL_OUTPUT_CODE +"', or '" + VECTOR_USUAL_OUTPUT_CODE + "', or at least one vector GIS output code";
            else
            {
               strRH = strToLower(&strRH);

               if (strRH.find(VECTOR_ALL_OUTPUT_CODE) != string::npos)
               {
                  // Output all vector files
                  m_bCoastSave =
                  m_bWaveAngleAndHeightSave =
                  m_bNormalsSave =
                  m_bInvalidNormalsSave =
                  m_bAvgWaveAngleAndHeightSave =
                  m_bWaveEnergySinceCollapseSave =
                  m_bMeanWaveEnergySave =
                  m_bBreakingWaveHeightSave =
                  m_bCoastCurvatureSave =
                  m_bPolygonNodeSave =
                  m_bPolygonBoundarySave =
                  m_bCliffNotchSave =
                  m_bShadowBoundarySave =
                  m_bShadowDowndriftBoundarySave =
                  m_bDeepWaterWaveAngleAndHeightSave =
                  m_bWaveSetupSave =
                  m_bStormSurgeSave =
                  m_bRunUpSave =
                  m_bVectorWaveFloodLineSave = true;
               }
               else if (strRH.find(VECTOR_USUAL_OUTPUT_CODE) != string::npos)
               {
                  // Output the "usual" collection of vector output files
                  m_bCoastSave =
                  m_bWaveAngleAndHeightSave =
                  m_bNormalsSave =
                  m_bAvgWaveAngleAndHeightSave =
                  m_bWaveEnergySinceCollapseSave =
                  m_bMeanWaveEnergySave =
                  m_bBreakingWaveHeightSave =
                  m_bPolygonBoundarySave =
                  m_bCliffNotchSave =
                  m_bShadowBoundarySave =
                  m_bShadowDowndriftBoundarySave =
                  m_bDeepWaterWaveAngleAndHeightSave = true;
               }
               else
               {
                  // Output only those vector files for which the user specified the code
                  if (strRH.find(VECTOR_COAST_CODE) != string::npos)
                  {
                     m_bCoastSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_COAST_CODE);
                  }

                  if (strRH.find(VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_CODE) != string::npos)
                  {
                     m_bWaveAngleAndHeightSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_CODE);
                  }

                  if (strRH.find(VECTOR_NORMALS_CODE) != string::npos)
                  {
                     m_bNormalsSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_NORMALS_CODE);
                  }

                  if (strRH.find(VECTOR_INVALID_NORMALS_CODE) != string::npos)
                  {
                     m_bInvalidNormalsSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_INVALID_NORMALS_CODE);
                  }

                  if (strRH.find(VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_CODE) != string::npos)
                  {
                     m_bAvgWaveAngleAndHeightSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_CODE);
                  }

                  if (strRH.find(VECTOR_COAST_CURVATURE_CODE) != string::npos)
                  {
                     m_bCoastCurvatureSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_COAST_CURVATURE_CODE);
                  }

                  if (strRH.find(VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_CODE) != string::npos)
                  {
                     m_bWaveEnergySinceCollapseSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_CODE);
                  }

                  if (strRH.find(VECTOR_MEAN_WAVE_ENERGY_CODE) != string::npos)
                  {
                     m_bMeanWaveEnergySave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_MEAN_WAVE_ENERGY_CODE);
                  }

                  if (strRH.find(VECTOR_BREAKING_WAVE_HEIGHT_CODE) != string::npos)
                  {
                     m_bBreakingWaveHeightSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_BREAKING_WAVE_HEIGHT_CODE);
                  }

                  if (strRH.find(VECTOR_POLYGON_NODE_CODE) != string::npos)
                  {
                     m_bPolygonNodeSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_POLYGON_NODE_CODE);
                  }

                  if (strRH.find(VECTOR_POLYGON_BOUNDARY_CODE) != string::npos)
                  {
                     m_bPolygonBoundarySave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_POLYGON_BOUNDARY_CODE);
                  }

                  if (strRH.find(VECTOR_CLIFF_NOTCH_SIZE_CODE) != string::npos)
                  {
                     m_bCliffNotchSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_CLIFF_NOTCH_SIZE_CODE);
                  }

                  if (strRH.find(VECTOR_SHADOW_BOUNDARY_CODE) != string::npos)
                  {
                     m_bShadowBoundarySave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_SHADOW_BOUNDARY_CODE);
                  }

                  if (strRH.find(VECTOR_DOWNDRIFT_BOUNDARY_CODE) != string::npos)
                  {
                     m_bShadowDowndriftBoundarySave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_DOWNDRIFT_BOUNDARY_CODE);
                  }

                  if (strRH.find(VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_CODE) != string::npos)
                  {
                     m_bDeepWaterWaveAngleAndHeightSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_CODE);
                  }

                  if (strRH.find(VECTOR_WAVE_SETUP_CODE) != string::npos)
                  {
                     m_bWaveSetupSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_WAVE_SETUP_CODE);
                  }

                  if (strRH.find(VECTOR_STORM_SURGE_CODE) != string::npos)
                  {
                     m_bStormSurgeSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_STORM_SURGE_CODE);
                  }

                  if (strRH.find(VECTOR_RUN_UP_CODE) != string::npos)
                  {
                     m_bRunUpSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_RUN_UP_CODE);
                  }

                  if (strRH.find(VECTOR_FLOOD_LINE_CODE) != string::npos)
                  {
                     m_bVectorWaveFloodLineSave = true;
                     strRH = strRemoveSubstr(&strRH, &VECTOR_FLOOD_LINE_CODE);
                  }

                  // Check to see if all codes have been removed
                  if (! strRH.empty())
                     strErr = "line " + to_string(nLine) + ": unknown code '" + strRH + "' in list of vector GIS output codes";
               }
            }
            break;

         case 16:
            // Vector GIS output format (note must retain original case)
            m_strVectorGISOutFormat = strRH;

            if (strRH.empty())
               strErr = "line " + to_string(nLine) + ": vector GIS output format";
            break;

         case 17:
            // Time series files to output
            if (! strRH.empty())
            {
               strRH = strToLower(&strRH);

               // First check for "all"
               if (strRH.find(RASTER_ALL_OUTPUT_CODE) != string::npos)
               {
                  m_bSeaAreaTSSave =
                  m_bStillWaterLevelTSSave =
                  m_bActualPlatformErosionTSSave =
                  m_bCliffCollapseErosionTSSave =
                  m_bCliffCollapseDepositionTSSave =
                  m_bCliffCollapseNetTSSave =
                  m_bBeachErosionTSSave =
                  m_bBeachDepositionTSSave =
                  m_bBeachSedimentChangeNetTSSave =
                  m_bSuspSedTSSave =
                  m_bFloodSetupSurgeTSSave =
                  m_bFloodSetupSurgeRunupTSSave = true;
               }
               else
               {
                  if (strRH.find(TIME_SERIES_SEA_AREA_CODE) != string::npos)
                  {
                     m_bSeaAreaTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_SEA_AREA_CODE);
                  }

                  if (strRH.find(TIME_SERIES_STILL_WATER_LEVEL_CODE) != string::npos)
                  {
                     m_bStillWaterLevelTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_STILL_WATER_LEVEL_CODE);
                  }

                  if (strRH.find(TIME_SERIES_PLATFORM_EROSION_CODE) != string::npos)
                  {
                     m_bActualPlatformErosionTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_PLATFORM_EROSION_CODE);
                  }

                  if (strRH.find(TIME_SERIES_CLIFF_COLLAPSE_EROSION_CODE) != string::npos)
                  {
                     m_bCliffCollapseErosionTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_CLIFF_COLLAPSE_EROSION_CODE);
                  }

                  if (strRH.find(TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_CODE) != string::npos)
                  {
                     m_bCliffCollapseDepositionTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_CLIFF_COLLAPSE_DEPOSITION_CODE);
                  }

                  if (strRH.find(TIME_SERIES_CLIFF_COLLAPSE_NET_CODE) != string::npos)
                  {
                     m_bCliffCollapseNetTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_CLIFF_COLLAPSE_NET_CODE);
                  }

                  if (strRH.find(TIME_SERIES_BEACH_EROSION_CODE) != string::npos)
                  {
                     m_bBeachErosionTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_BEACH_EROSION_CODE);
                  }

                  if (strRH.find(TIME_SERIES_BEACH_DEPOSITION_CODE) != string::npos)
                  {
                     m_bBeachDepositionTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_BEACH_DEPOSITION_CODE);
                  }

                  if (strRH.find(TIME_SERIES_BEACH_CHANGE_NET_CODE) != string::npos)
                  {
                     m_bBeachSedimentChangeNetTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_BEACH_CHANGE_NET_CODE);
                  }

                  if (strRH.find(TIME_SERIES_SUSPENDED_SEDIMENT_CODE) != string::npos)
                  {
                     m_bSuspSedTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_SUSPENDED_SEDIMENT_CODE);
                  }

                  if (strRH.find(TIME_SERIES_FLOOD_SETUP_SURGE_CODE) != string::npos)
                  {
                     m_bFloodSetupSurgeTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_FLOOD_SETUP_SURGE_CODE);
                  }

                  if (strRH.find(TIME_SERIES_FLOOD_SETUP_SURGE_RUNUP_CODE) != string::npos)
                  {
                     m_bFloodSetupSurgeRunupTSSave = true;
                     strRH = strRemoveSubstr(&strRH, &TIME_SERIES_FLOOD_SETUP_SURGE_RUNUP_CODE);
                  }

                  // Check to see if all codes have been removed
                  if (! strRH.empty())
                     strErr = "line " + to_string(nLine) + ": unknown code '" + strRH + "' in list of time series output files";
               }
            }
            break;

         case 18:
            // Vector coastline smoothing algorithm: 0 = none, 1 = running mean, 2 = Savitsky-Golay
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for coastline smoothing algorithm '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nCoastSmooth = stoi(strRH);

            if ((m_nCoastSmooth < SMOOTH_NONE) || (m_nCoastSmooth > SMOOTH_SAVITZKY_GOLAY))
               strErr = "line " + to_string(nLine) + ": coastline vector smoothing algorithm";

            break;

         case 19:
            // Size of coastline smoothing window: must be odd
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for coastline smoothing window '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nCoastSmoothWindow = stoi(strRH);

            if ((m_nCoastSmoothWindow <= 0) || !(m_nCoastSmoothWindow % 2))
               strErr = "line " + to_string(nLine) + ": size of coastline vector smoothing window (must be > 0 and odd)";

            break;

         case 20:
            // Order of coastline profile smoothing polynomial for Savitsky-Golay: usually 2 or 4, max is 6
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for Savitsky-Golay polynomial for coastline smoothing '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nSavGolCoastPoly = stoi(strRH);

            if ((m_nSavGolCoastPoly <= 0) || (m_nSavGolCoastPoly > 6))
               strErr = "line " + to_string(nLine) + ": value of Savitsky-Golay polynomial for coastline smoothing (must be <= 6)";

            break;

         case 21:
            // Grid edge(s) to omit when searching for coastline [NSWE]
            strRH = strToLower(&strRH);

            bFound = false;
            if (strRH.find("n") != string::npos)
            {
               m_bOmitSearchNorthEdge =
                   bFound = true;
            }
            if (strRH.find("s") != string::npos)
            {
               m_bOmitSearchSouthEdge =
                   bFound = true;
            }
            if (strRH.find("w") != string::npos)
            {
               m_bOmitSearchWestEdge =
                   bFound = true;
            }
            if (strRH.find("e") != string::npos)
            {
               m_bOmitSearchEastEdge =
                   bFound = true;
            }

            if (! bFound)
               strErr = "line " + to_string(nLine) + ": at least one edge must be omitted from coastline search";

            break;

         case 22:
            // Profile slope running-mean smoothing window size: must be odd
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for size of coastline smoothing window '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nProfileSmoothWindow = stoi(strRH);

            if ((m_nProfileSmoothWindow < 0) || (m_nProfileSmoothWindow > 0 && !(m_nProfileSmoothWindow % 2)))
               strErr = "line " + to_string(nLine) + ": size of profile vector smoothing window (must be >= 0, if > 0 must be odd)";

            break;

         case 23:
            // Max local slope (m/m), first check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for max local slope '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dProfileMaxSlope = strtod(strRH.c_str(), NULL);

            if (m_dProfileMaxSlope <= 0)
               strErr = "line " + to_string(nLine) + ": max local slope must be > 0";
            break;

         case 24:
            // Maximum elevation of beach above SWL, first check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for maximum elevation of beach above SWL '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dMaxBeachElevAboveSWL = strtod(strRH.c_str(), NULL);

            if (m_dMaxBeachElevAboveSWL < 0)
               strErr = "line " + to_string(nLine) + ": maximum elevation of beach above SWL must be >= 0";
            break;

         // ------------------------------------------------- Raster GIS layers ------------------------------------------------
         case 25:
            // Number of sediment layers
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for number of sediment layers '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nLayers = stoi(strRH);

            if (m_nLayers < 1)
            {
               strErr = "line " + to_string(nLine) + ": must be at least one sediment layer";
               break;
            }

            // OK we know the number of layers, so add elements to these vectors
            for (int j = 0; j < m_nLayers; j++)
            {
               m_VstrInitialFineUnconsSedimentFile.push_back("");
               m_VstrInitialSandUnconsSedimentFile.push_back("");
               m_VstrInitialCoarseUnconsSedimentFile.push_back("");
               m_VstrInitialFineConsSedimentFile.push_back("");
               m_VstrInitialSandConsSedimentFile.push_back("");
               m_VstrInitialCoarseConsSedimentFile.push_back("");
               m_VstrGDALIUFDriverCode.push_back("");
               m_VstrGDALIUFDriverDesc.push_back("");
               m_VstrGDALIUFProjection.push_back("");
               m_VstrGDALIUFDataType.push_back("");
               m_VstrGDALIUSDriverCode.push_back("");
               m_VstrGDALIUSDriverDesc.push_back("");
               m_VstrGDALIUSProjection.push_back("");
               m_VstrGDALIUSDataType.push_back("");
               m_VstrGDALIUCDriverCode.push_back("");
               m_VstrGDALIUCDriverDesc.push_back("");
               m_VstrGDALIUCProjection.push_back("");
               m_VstrGDALIUCDataType.push_back("");
               m_VstrGDALICFDriverCode.push_back("");
               m_VstrGDALICFDriverDesc.push_back("");
               m_VstrGDALICFProjection.push_back("");
               m_VstrGDALICFDataType.push_back("");
               m_VstrGDALICSDriverCode.push_back("");
               m_VstrGDALICSDriverDesc.push_back("");
               m_VstrGDALICSProjection.push_back("");
               m_VstrGDALICSDataType.push_back("");
               m_VstrGDALICCDriverCode.push_back("");
               m_VstrGDALICCDriverDesc.push_back("");
               m_VstrGDALICCProjection.push_back("");
               m_VstrGDALICCDataType.push_back("");
            }
            break;

         case 26:
            // Basement DEM file (can be blank)
            if (! strRH.empty())
            {
#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif
               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                  // It has an absolute path, so use it 'as is'
                  m_strInitialBasementDEMFile = strRH;
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInitialBasementDEMFile = m_strCMEDir;
                  m_strInitialBasementDEMFile.append(strRH);
               }
            }
            break;

         case 27:
            // Read 6 x sediment files for each layer
            for (int nLayer = 0; nLayer < m_nLayers; nLayer++)
            {
               for (int j = 1; j <= 6; j++)
               {
                  if (! bFirst)
                  {
                     do
                     {
                        if (! getline(InStream, strRec))
                        {
                           cerr << ERR << "premature end of file in " << m_strDataPathName << endl;
                           return false;
                        }

                        nLine++;

                        // Trim off leading and trailing whitespace
                        strRec = strTrim(&strRec);
                     }
                     // If it is a blank line or a comment then ignore it
                     while (strRec.empty() || (strRec[0] == QUOTE1) || (strRec[0] == QUOTE2));

                     // Not blank or a comment, so find the colon: lines MUST have a colon separating data from leading description portion
                     nPos = strRec.find(COLON);
                     if (nPos == string::npos)
                     {
                        // Error: badly formatted (no colon)
                        cerr << ERR << "on line " << to_string(nLine) << ": badly formatted (no ':') in " << m_strDataPathName << endl << strRec << endl;
                        return false;
                     }

                     // Strip off leading portion (the bit up to and including the colon)
                     strRH = strRec.substr(nPos + 1);
// ERROR                     strRH.resize(nPos);

                     // Remove leading whitespace after the colon
                     strRH = strTrimLeft(&strRH);

                     // Look for a trailing comment, if found then terminate string at that point and trim off any trailing whitespace
                     nPos = strRH.rfind(QUOTE1);
                     if (nPos != string::npos)
                        strRH.resize(nPos);

                     nPos = strRH.rfind(QUOTE2);
                     if (nPos != string::npos)
                        strRH.resize(nPos);

                     // Trim trailing spaces
                     strRH = strTrimRight(&strRH);

#ifdef _WIN32
                     // For Windows, make sure has backslashes, not Unix-style slashes
                     strRH = pstrChangeToBackslash(&strRH);
#endif
                  }

                  bFirst = false;

                  switch (j)
                  {
                  case 1:
                     // Initial fine unconsolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
                        // Set switch
                        m_bHaveFineSediment = true;
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif

                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialFineUnconsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialFineUnconsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialFineUnconsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 2:
                     // Initial sand unconsolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
                        // Set switch
                        m_bHaveSandSediment = true;
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialSandUnconsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialSandUnconsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialSandUnconsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 3:
                     // Initial coarse unconsolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
                        // Set switch
                        m_bHaveCoarseSediment = true;
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialCoarseUnconsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialCoarseUnconsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialCoarseUnconsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 4:
                     // Initial fine consolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
                        // Set switch
                        m_bHaveFineSediment = true;
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialFineConsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialFineConsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialFineConsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 5:
                     // Initial sand consolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
                        // Set switch
                        m_bHaveSandSediment = true;
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialSandConsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialSandConsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialSandConsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;

                  case 6:
                     // Initial coarse consolidated sediment depth GIS file (can be blank)
                     if (! strRH.empty())
                     {
                        // Set switch
                        m_bHaveCoarseSediment = true;
#ifdef _WIN32
                        // For Windows, make sure has backslashes, not Unix-style slashes
                        strRH = pstrChangeToBackslash(&strRH);
#endif
                        // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                        if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                        {
                           // It has an absolute path, so use it 'as is'
                           m_VstrInitialCoarseConsSedimentFile[nLayer] = strRH;
                        }
                        else
                        {
                           // It has a relative path, so prepend the CoastalME dir
                           m_VstrInitialCoarseConsSedimentFile[nLayer] = m_strCMEDir;
                           m_VstrInitialCoarseConsSedimentFile[nLayer].append(strRH);
                        }
                     }
                     break;
                  }

                  // Did an error occur?
                  if (! strErr.empty())
                  {
                     // Error in input to run details file
                     cerr << ERR << "reading " << strErr << " in " << m_strDataPathName << endl
                          << "'" << strRec << "'" << endl;
                     InStream.close();
                     return false;
                  }
               }
            }
            break;

         case 28:
            // Initial suspended sediment depth GIS file (can be blank)
            if (! strRH.empty())
            {
               // Set switch
               m_bHaveFineSediment = true;
#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif
               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strInitialSuspSedimentFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInitialSuspSedimentFile = m_strCMEDir;
                  m_strInitialSuspSedimentFile.append(strRH);
               }
            }
            break;

         case 29:
            // Initial Landform class GIS file (can be blank)
            if (! strRH.empty())
            {
#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif
               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strInitialLandformFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInitialLandformFile = m_strCMEDir;
                  m_strInitialLandformFile.append(strRH);
               }
            }
            break;

         case 30:
            // Initial Intervention class GIS file (can be blank: if so then intervention height file must also be blank)
            if (! strRH.empty())
            {
#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif
               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strInterventionClassFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInterventionClassFile = m_strCMEDir;
                  m_strInterventionClassFile.append(strRH);
               }

               // Set the save switches
               m_bInterventionClassSave = true;
               m_bInterventionHeightSave = true;
            }
            break;

         case 31:
            // Initial Intervention height GIS file (can be blank: if so then intervention class file must also be blank)
            if (strRH.empty())
            {
               if (! m_strInterventionClassFile.empty())
                  strErr = "line " + to_string(nLine) + ": must specify both intervention class and intervention height files";
               break;
            }
            else
            {
               if (m_strInterventionClassFile.empty())
               {
                  strErr = "line " + to_string(nLine) + ": must specify both intervention class and intervention height files";
                  break;
               }

#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif

               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strInterventionHeightFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strInterventionHeightFile = m_strCMEDir;
                  m_strInterventionHeightFile.append(strRH);
               }
            }
            break;

            // ---------------------------------------------------- Hydrology data ------------------------------------------------
         case 32:
            // Wave propagation model [0 = COVE, 1 = CShore]
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for wave propagation model '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nWavePropagationModel = stoi(strRH);

            if ((m_nWavePropagationModel != WAVE_MODEL_COVE) && (m_nWavePropagationModel != WAVE_MODEL_CSHORE))
               strErr = "line " + to_string(nLine) + ": wave propagation model must be 0 or 1";

            break;

         case 33:
            // Density of sea water (kg/m3), first check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for sea water density '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dSeaWaterDensity = strtod(strRH.c_str(), NULL);

            if (m_dSeaWaterDensity <= 0)
               strErr = "line " + to_string(nLine) + ": sea water density must be > 0";
            break;

         case 34:
            // Initial still water level (m), first check that this is a valid double TODO 041 Make this a per-timestep SWL file
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for initial SWL '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dOrigSWL = strtod(strRH.c_str(), NULL);
            break;

         case 35:
            // Final still water level (m) [blank = same as initial SWL]
            if (strRH.empty())
               m_dFinalSWL = m_dOrigSWL;
            else
            {
               // Check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for final SWL '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dFinalSWL = strtod(strRH.c_str(), NULL);
            }
            break;

         case 36:
            // Deep water wave height (m) or a file of point vectors giving deep water wave height (m) and orientation (for units, see below)
            if (strRH.empty())
               strErr = "line " + to_string(nLine) + ": deep water wave height in " + m_strDataPathName + " must be either a number or a filename (filename must not start with a number)";
            else
            {
               if (isdigit(strRH.at(0))) // If this starts with a number then is a single value, otherwise is a filename. Note that filename must not start with number
               {
                  // Just one value of wave height for all deep water cells, first check that this is a valid double
                  if (! bIsStringValidDouble(strRH))
                  {
                     strErr = "line " + to_string(nLine) + ": invalid floating point number for deep water wave height '" + strRH + "' in " + m_strDataPathName;
                     break;
                  }

                  m_bSingleDeepWaterWaveValues = true;
                  m_bHaveWaveStationData = false;

                  m_dAllCellsDeepWaterWaveHeight = strtod(strRH.c_str(), NULL);

                  if (m_dAllCellsDeepWaterWaveHeight <= 0)
                     strErr = "line " + to_string(nLine) + ": deep water wave height must be > 0";
               }
               else
               {
                  // We are reading deep water wave height and deep water wave orientation from two files. This first file is a point shape file with the location of the buoys and integer ID for each one
                  m_bHaveWaveStationData = true;
#ifdef _WIN32
                  // For Windows, make sure has backslashes, not Unix-style slashes
                  strRH = pstrChangeToBackslash(&strRH);
#endif
                  // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                  if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                  {
                     // It has an absolute path, so use it 'as is'
                     m_strDeepWaterWaveStationsShapefile = strRH;
                  }
                  else
                  {
                     // It has a relative path, so prepend the CoastalME dir
                     m_strDeepWaterWaveStationsShapefile = m_strCMEDir;
                     m_strDeepWaterWaveStationsShapefile.append(strRH);
                  }
               }
            }
            break;

         case 37:
            // Deep water wave height time series file
            if (m_bHaveWaveStationData)
            {
               // Only read this if we have a file for wave height points. Each point in m_strDeepWaterWavesTimeSeriesFile is a triad of wave height, orientation and period for each time step
               if (strRH.empty())
               {
                  strErr = "line " + to_string(nLine) + ": filename missing for deep water wave height time series";
                  break;
               }

#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif

               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strDeepWaterWavesTimeSeriesFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strDeepWaterWavesTimeSeriesFile = m_strCMEDir;
                  m_strDeepWaterWavesTimeSeriesFile.append(strRH);
               }
            }
            break;

         case 38:
            // Deep water wave orientation in input CRS: this is the oceanographic convention i.e. direction TOWARDS which the waves move (in degrees clockwise from north)
            if (! m_bHaveWaveStationData)
            {
               // Only read this if we have just a single value of wave height for all deep water cells. Check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for deep water wave orientation '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dAllCellsDeepWaterWaveAngle = strtod(strRH.c_str(), NULL);

               if (m_dAllCellsDeepWaterWaveAngle < 0)
                  strErr = "line " + to_string(nLine) + ": deep water wave orientation must be zero degrees or more";
               else if (m_dAllCellsDeepWaterWaveAngle >= 360)
                  strErr = "line " + to_string(nLine) + ": deep water wave orientation must be less than 360 degrees";
            }
            break;

         case 39:
            // Wave period (sec)
            if (! m_bHaveWaveStationData)
            {
               // Only read this if we also have just a single value of wave height for all deep water cells. Check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for wave period '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dAllCellsDeepWaterWavePeriod = strtod(strRH.c_str(), NULL);

               if (m_dAllCellsDeepWaterWavePeriod <= 0)
                  strErr = "line " + to_string(nLine) + ": wave period must be > 0";
            }
            break;

         case 40:
            // Tide data file (can be blank). This is the change (m) from still water level for each timestep
            if (! strRH.empty())
            {
#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif
               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                  // It has an absolute path, so use it 'as is'
                  m_strTideDataFile = strRH;
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strTideDataFile = m_strCMEDir;
                  m_strTideDataFile.append(strRH);
               }
            }
            break;

         case 41:
            // Breaking wave height-to-depth ratio, check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for breaking wave height to depth ratio '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dBreakingWaveHeightDepthRatio = strtod(strRH.c_str(), NULL);

            if (m_dBreakingWaveHeightDepthRatio <= 0)
               strErr = "line " + to_string(nLine) + ": breaking wave height to depth ratio must be > 0";
            break;

         // ----------------------------------------------------- Sediment data ------------------------------------------------
         case 42:
            // Simulate coast platform erosion?
            strRH = strToLower(&strRH);

            if (strRH.find("y") != string::npos)
               m_bDoShorePlatformErosion = true;
            break;

         case 43:
            // If simulating coast platform erosion, R (coast platform resistance to erosion) values along profile, see Walkden & Hall, 2011
            if (m_bDoShorePlatformErosion)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for R (coast platform resistance to erosion) '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dR = strtod(strRH.c_str(), NULL);

               if (m_dR <= 0)
                  strErr = "line " + to_string(nLine) + ": R (coast platform resistance to erosion) value must be > 0";
            }
            break;

         case 44:
            // Simulate beach sediment transport?
            strRH = strToLower(&strRH);

            m_bDoBeachSedimentTransport = false;
            if (strRH.find("y") != string::npos)
               m_bDoBeachSedimentTransport = true;
            break;

         case 45:
            // If simulating beach sediment transport, beach sediment transport at grid edges [0 = closed, 1 = open, 2 = re-circulate]
            if (m_bDoBeachSedimentTransport)
            {
               if (! bIsStringValidInt(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid integer for beach sediment transport at grid edges '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_nUnconsSedimentHandlingAtGridEdges = stoi(strRH);

               if ((m_nUnconsSedimentHandlingAtGridEdges < GRID_EDGE_CLOSED) || (m_nUnconsSedimentHandlingAtGridEdges > GRID_EDGE_RECIRCULATE))
                  strErr = "line " + to_string(nLine) + ": switch for handling of beach sediment at grid edges must be 0, 1, or 2";
            }
            break;

         case 46:
            // If simulating beach sediment transport, beach erosion/deposition equation [0 = CERC, 1 = Kamphuis]
            if (m_bDoBeachSedimentTransport)
            {
               if (! bIsStringValidInt(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid integer for beach erosion/deposition equation '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_nBeachErosionDepositionEquation = stoi(strRH);

               if ((m_nBeachErosionDepositionEquation != UNCONS_SEDIMENT_EQUATION_CERC) && (m_nBeachErosionDepositionEquation != UNCONS_SEDIMENT_EQUATION_KAMPHUIS))
                  strErr = "line " + to_string(nLine) + ": switch for beach erosion/deposition equation must be 0 or 1";
            }
            break;

         case 47:
            // Median size of fine sediment (mm), always needed [0 = default, only for Kamphuis eqn]. First check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for median particle size of fine sediment '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dD50Fine = strtod(strRH.c_str(), NULL);

            if (m_dD50Fine < 0)
               strErr = "line " + to_string(nLine) + ": median particle size of fine sediment must be > 0";
            else if (bFPIsEqual(m_dD50Fine, 0.0, TOLERANCE))
               // Use default value
               m_dD50Fine = D50_FINE_DEFAULT;

            break;

         case 48:
            // Median size of sand sediment (mm), always needed [0 = default, only for Kamphuis eqn]. First check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for median particle size of sand sediment '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dD50Sand = strtod(strRH.c_str(), NULL);

            if (m_dD50Sand < 0)
               strErr = "line " + to_string(nLine) + ": median particle size of sand sediment must be > 0";
            else if (bFPIsEqual(m_dD50Sand, 0.0, TOLERANCE))
               // Use default value
               m_dD50Sand = D50_SAND_DEFAULT;

            break;

         case 49:
            // Median size of coarse sediment (mm), always needed [0 = default, only for Kamphuis eqn]. First check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for median particle size of coarse sediment '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dD50Coarse = strtod(strRH.c_str(), NULL);

            if (m_dD50Coarse < 0)
               strErr = "line " + to_string(nLine) + ": median particle size of coarse sediment must be > 0";
            else if (bFPIsEqual(m_dD50Coarse, 0.0, TOLERANCE))
               // Use default value
               m_dD50Coarse = D50_COARSE_DEFAULT;

            break;

         case 50:
            // Density of unconsolidated beach sediment (kg/m3)
            if (m_bDoBeachSedimentTransport)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for density of beach sediment '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dBeachSedimentDensity = strtod(strRH.c_str(), NULL);

               if (m_dBeachSedimentDensity <= 0)
                  strErr = "line " + to_string(nLine) + ": density of beach sediment must be > 0";
            }
            break;

         case 51:
            // Beach sediment porosity
            if (m_bDoBeachSedimentTransport)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for porosity of beach sediment '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dBeachSedimentPorosity = strtod(strRH.c_str(), NULL);

               if (m_dBeachSedimentPorosity <= 0)
                  strErr = "line " + to_string(nLine) + ": porosity of beach sediment must be > 0";
            }
            break;

         case 52:
            // Relative erodibility (0 - 1) of fine-sized sediment, always needed. First check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for erodibility of fine-sized sediment '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dFineErodibility = strtod(strRH.c_str(), NULL);

            if ((m_dFineErodibility < 0) || (m_dFineErodibility > 1))
               strErr = "line " + to_string(nLine) + ": relative erodibility of fine-sized sediment must be between 0 and 1";

            break;

         case 53:
            // Relative erodibility (0 - 1) of sand-sized sediment, always needed. First check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for erodibility of sand-sized sediment '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dSandErodibility = strtod(strRH.c_str(), NULL);

            if ((m_dSandErodibility < 0) || (m_dSandErodibility > 1))
               strErr = "line " + to_string(nLine) + ": relative erodibility of sand-sized sediment must be between 0 and 1";

            break;

         case 54:
            // Relative erodibility (0 - 1) of coarse-sized sediment, always needed. First check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for erodibility of coarse-sized sediment '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dCoarseErodibility = strtod(strRH.c_str(), NULL);

            if ((m_dCoarseErodibility < 0) || (m_dCoarseErodibility > 1))
            {
               strErr = "line " + to_string(nLine) + ": relative erodibility of coarse-sized sediment must be between 0 and 1";
               break;
            }

            if ((m_dFineErodibility + m_dSandErodibility + m_dCoarseErodibility) <= 0)
               strErr = "line " + to_string(nLine) + ": must have at least one non-zero erodibility value";

            break;

         case 55:
            // Transport parameter KLS in CERC equation
            if (m_bDoBeachSedimentTransport)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for transport parameter KLS of CERC equation '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dKLS = strtod(strRH.c_str(), NULL);
               // However, many sites do not have transport data available to calibrate K, and for design applications without calibration data the CERC formula provides only order-of-magnitude accuracy (Fowler et al., 1995; Wang et al., 1998). The recommended value of K = 0.39 has been commonly used to represent the potential longshore transport rate. However, Miller (1998) found that the CERC formula sometimes over and sometimes under predicted longshore transport rate for measurements during storms, indicating the value of K also can be higher than 0.39
               // TODO 042 Should this be a user input, or not? The comment above seems inconclusive
               // m_dKLS = tMin(m_dKLS, 0.39);

               if (m_dKLS <= 0)
                  strErr = "line " + to_string(nLine) + ": transport parameter KLS of CERC equation must be > 0";
            }
            break;

         case 56:
            // Transport parameter for Kamphuis equation
            if (m_bDoBeachSedimentTransport)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for transport parameter of Kamphuis equation '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dKamphuis = strtod(strRH.c_str(), NULL);

               if (m_dKamphuis <= 0)
                  strErr = "line " + to_string(nLine) + ": transport parameter of Kamphuis equation must be > 0";
            }
            break;

         case 57:
            // Berm height i.e. height above SWL of start of depositional Dean profile
            if (m_bDoBeachSedimentTransport)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for Dean profile start height above SWL '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dDeanProfileStartAboveSWL = strtod(strRH.c_str(), NULL);

               if (m_dDeanProfileStartAboveSWL < 0)
                  strErr = "line " + to_string(nLine) + ": Berm height (Dean profile start height above SWL) must be >= 0";
            }
            break;

         // ------------------------------------------------ Cliff collapse data -----------------------------------------------
         case 58:
            // Simulate cliff collapse?
            strRH = strToLower(&strRH);

            if (strRH.find("y") != string::npos)
               m_bDoCliffCollapse = true;
            break;

         case 59:
            // Cliff resistance to erosion
            if (m_bDoCliffCollapse)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for cliff resistance to erosion '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dCliffErosionResistance = strtod(strRH.c_str(), NULL);

               if (m_dCliffErosionResistance <= 0)
                  strErr = "line " + to_string(nLine) + ": cliff resistance to erosion must be > 0";
            }
            break;

         case 60:
            // Notch overhang at collapse (m)
            if (m_bDoCliffCollapse)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for cliff notch overhang at collapse '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dNotchDepthAtCollapse = strtod(strRH.c_str(), NULL);

               if (m_dNotchDepthAtCollapse <= 0)
                  strErr = "line " + to_string(nLine) + ": cliff notch overhang at collapse must be > 0";
            }
            break;

         case 61:
            // Notch base below still water level (m)
            if (m_bDoCliffCollapse)
            {
               m_dNotchBaseBelowSWL = strtod(strRH.c_str(), NULL);
               if (m_dNotchBaseBelowSWL < 0)
                  strErr = "line " + to_string(nLine) + ": cliff notch base below still water level must be > 0";
            }
            break;

         case 62:
            // Scale parameter A for cliff deposition (m^(1/3)) [0 = auto]
            if (m_bDoCliffCollapse)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for scale parameter A for cliff deposition '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dCliffDepositionA = strtod(strRH.c_str(), NULL);

               if (m_dCliffDepositionA < 0)
                  strErr = "line " + to_string(nLine) + ": scale parameter A for cliff deposition must be 0 [= auto] or greater";
            }
            break;

         case 63:
            // Approximate planview width of cliff collapse talus (in m)
            if (m_bDoCliffCollapse)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for width of cliff collapse talus '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dCliffDepositionPlanviewWidth = strtod(strRH.c_str(), NULL);

               if (m_dCliffDepositionPlanviewWidth <= 0)
                  strErr = "line " + to_string(nLine) + ": planview width of cliff deposition must be > 0";

            }
            break;

         case 64:
            // Planview length of cliff deposition talus (m)
            if (m_bDoCliffCollapse)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for planview length of cliff deposition '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dCliffTalusMinDepositionLength = strtod(strRH.c_str(), NULL);

               if (m_dCliffTalusMinDepositionLength <= 0)
                  strErr = "line " + to_string(nLine) + ": planview length of cliff deposition must be > 0";
            }
            break;

         case 65:
            // Height of landward end of talus, as a fraction of cliff elevation
            if (m_bDoCliffCollapse)
            {
               // First check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for height of cliff collapse (as a fraction of cliff elevation) '" + strRH + "' in " + m_strDataPathName;
                  break;
               }

               m_dMinCliffTalusHeightFrac = strtod(strRH.c_str(), NULL);

               if (m_dMinCliffTalusHeightFrac < 0)
                  strErr = "line " + to_string(nLine) + ": height of cliff collapse (as a fraction of cliff elevation) must be >= 0";
            }
            break;

         // -------------------------------------------------- Input events data -----------------------------------------------
         case 66:
            // Simulate riverine flooding?
            strRH = strToLower(&strRH);

            if (strRH.find("y") != string::npos)
            {
               m_bDoRiverineFlooding = true;
               m_bSetupSurgeFloodMaskSave = true;
               m_bSetupSurgeRunupFloodMaskSave = true;
               m_bRasterWaveFloodLineSave = true;
            }
            break;

         case 67:
            // Output riverine flooding vector files
            if (m_bDoRiverineFlooding)
            {
               if (! strRH.empty())
               {
                  // Convert to lower case
                  strRH = strToLower(&strRH);

                  // First look for "all"
                  if (strRH.find(VECTOR_ALL_RIVER_FLOOD_OUTPUT_CODE) != string::npos)
                  {
                     m_bFloodSWLSetupLine =
                     m_bFloodSWLSetupSurgeLine =
                     m_bFloodSWLSetupSurgeRunupLine = true;
                  }
                  else
                  {
                     // We are not outputting all vector flood GIS files, so set switches (and remove strings) for those optional files for which the user specified the code
                     if (strRH.find(VECTOR_FLOOD_SWL_SETUP_LINE_CODE) != string::npos)
                     {
                        m_bFloodSWLSetupLine = true;
                        strRH = strRemoveSubstr(&strRH, &VECTOR_FLOOD_SWL_SETUP_LINE_CODE);
                     }

                     if (strRH.find(VECTOR_FLOOD_SWL_SETUP_SURGE_LINE_CODE) != string::npos)
                     {
                        m_bFloodSWLSetupSurgeLine = true;
                        strRH = strRemoveSubstr(&strRH, &VECTOR_FLOOD_SWL_SETUP_SURGE_LINE_CODE);
                     }

                     if (strRH.find(VECTOR_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE_CODE) != string::npos)
                     {
                        m_bFloodSWLSetupSurgeRunupLine = true;
                        strRH = strRemoveSubstr(&strRH, &VECTOR_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE_CODE);
                     }

                     // Check to see if all codes have been removed
                     if (! strRH.empty())
                        strErr = "line " + to_string(nLine) + ": unknown code '" + strRH + "' in list of riverine flooding output codes";
                  }
               }
               else
                  strErr = "line " + to_string(nLine) + ": if simulating riverine flooding, must contain '" + VECTOR_ALL_RIVER_FLOOD_OUTPUT_CODE + "' or at least one vector GIS output code for riverine flooding";
            }
            break;

         case 68:
            if (m_bDoRiverineFlooding)
            {
               // Run-up equation?
               if (bIsStringValidInt(strRH))
               {
                  m_nRunUpEquation = stoi(strRH);
               }
               else
                  strErr = "line " + to_string(nLine) + ": invalid code for run-up equation used in simulating floods";
            }
            break;

         case 69:
            if (m_bDoRiverineFlooding && m_bVectorWaveFloodLineSave)
            {
               // Characteristic locations for flood?
               strRH = strToLower(&strRH);

               m_bFloodLocation = false;
               if (strRH.find("y") != string::npos)
               {
                  m_bFloodLocation = true;
               }
            }
            break;

         case 70:
            if (m_bDoRiverineFlooding)
            {
               // Path of location points file
               if (! strRH.empty())
               {
#ifdef _WIN32
                  // For Windows, make sure has backslashes, not Unix-style slashes
                  strRH = pstrChangeToBackslash(&strRH);
#endif
                  // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                  if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                     // It has an absolute path, so use it 'as is'
                     m_strFloodLocationShapefile = strRH;
                  else
                  {
                     // It has a relative path, so prepend the CoastalME dir
                     m_strFloodLocationShapefile = m_strCMEDir;
                     m_strFloodLocationShapefile.append(strRH);
                  }

                  // Set the switch
                  m_bFloodLocation = true;
               }
               else
                  strErr = "line " + to_string(nLine) + ": path of location points file must not be empty if simulating floods";
            }
            break;

         case 71:
            // Simulate sediment input?
            strRH = strToLower(&strRH);

            if (strRH.find("y") != string::npos)
            {
               m_bSedimentInput = true;
               m_bSedimentInputEventSave = true;
            }
            break;

         case 72:
            // Sediment input location (point or line shapefile)
            if (m_bSedimentInput)
            {
               if (! strRH.empty())
               {
#ifdef _WIN32
                  // For Windows, make sure has backslashes, not Unix-style slashes
                  strRH = pstrChangeToBackslash(&strRH);
#endif
                  // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
                  if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
                     // It has an absolute path, so use it 'as is'
                     m_strSedimentInputEventShapefile = strRH;
                  else
                  {
                     // It has a relative path, so prepend the CoastalME dir
                     m_strSedimentInputEventShapefile = m_strCMEDir;
                     m_strSedimentInputEventShapefile.append(strRH);
                  }
               }
            }
            break;

         case 73:
            // Sediment input type: required if have shapefile [P = Point, C = coast block, L = line]
            if (m_bSedimentInput)
            {
               strRH = strToLower(&strRH);

               if (strRH.find("p") != string::npos)
                  m_bSedimentInputAtPoint = true;
               else if (strRH.find("c") != string::npos)
                  m_bSedimentInputAtCoast = true;
               else if (strRH.find("l") != string::npos)
                  m_bSedimentInputAlongLine = true;
               else
                  strErr = "line " + to_string(nLine) + ": Sediment input type must be P, C, or L";
            }
            break;

         case 74:
            // Sediment input details file (required if have shapefile)
            if (m_bSedimentInput)
            {
               if (strRH.empty())
               {
                  strErr = "line " + to_string(nLine) + ": filename missing for sediment input time series";
                  break;
               }

#ifdef _WIN32
               // For Windows, make sure has backslashes, not Unix-style slashes
               strRH = pstrChangeToBackslash(&strRH);
#endif

               // Now check for leading slash, or leading Unix home dir symbol, or occurrence of a drive letter
               if ((strRH[0] == PATH_SEPARATOR) || (strRH[0] == TILDE) || (strRH[1] == COLON))
               {
                  // It has an absolute path, so use it 'as is'
                  m_strSedimentInputEventTimeSeriesFile = strRH;
               }
               else
               {
                  // It has a relative path, so prepend the CoastalME dir
                  m_strSedimentInputEventTimeSeriesFile = m_strCMEDir;
                  m_strSedimentInputEventTimeSeriesFile.append(strRH);
               }
            }
            break;

         // ------------------------------------------------------ Other data --------------------------------------------------
         case 75:
            // Gravitational acceleration (m2/s). First check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for gravitational acceleration '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dG = strtod(strRH.c_str(), NULL);

            if (m_dG <= 0)
               strErr = "line " + to_string(nLine) + ": gravitational acceleration must be > 0";
            break;

         case 76:
            // Spacing of coastline normals (m)
            m_dCoastNormalAvgSpacing = strtod(strRH.c_str(), NULL);

            if (bFPIsEqual(m_dCoastNormalAvgSpacing, 0.0, TOLERANCE))
               m_nCoastNormalAvgSpacing = MIN_PROFILE_SPACING; // In cells, we will set m_dCoastNormalAvgSpacing later when we know m_dCellSide
            else if (m_dCoastNormalAvgSpacing < 0)
               strErr = "line " + to_string(nLine) + ": spacing of coastline normals must be > 0";
            break;

         case 77:
            // Random factor for spacing of normals  [0 to 1, 0 = deterministic], check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for random factor for spacing of coastline normals '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dCoastNormalRandSpacingFactor = strtod(strRH.c_str(), NULL);

            if (m_dCoastNormalRandSpacingFactor < 0)
               strErr = "line " + to_string(nLine) + ": random factor for spacing of coastline normals must be >= 0";
            else if (m_dCoastNormalRandSpacingFactor > 1)
               strErr = "line " + to_string(nLine) + ": random factor for spacing of coastline normals must be < 1";
            break;

         case 78:
            // Length of coastline normals (m), check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for length of coastline normals '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dCoastNormalLength = strtod(strRH.c_str(), NULL);

            if (m_dCoastNormalLength <= 0)
               strErr = "line " + to_string(nLine) + ": length of coastline normals must be > 0";
            break;

         case 79:
            // Maximum number of 'cape' normals
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for maximum number of 'cape' normals '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nNaturalCapeNormals = stoi(strRH);

            if (m_nNaturalCapeNormals < 0)
               strErr = "line " + to_string(nLine) + ": number of 'cape' normals must be >= 0";

            break;

         case 80:
            // Start depth for wave calcs (ratio to deep water wave height), check that this is a valid double
            if (! bIsStringValidDouble(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number for start depth for wave calcs '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_dWaveDepthRatioForWaveCalcs = strtod(strRH.c_str(), NULL);

            if (m_dWaveDepthRatioForWaveCalcs <= 0)
               strErr = "line " + to_string(nLine) + ": start depth for wave calcs must be > 0";
            break;


         // ----------------------------------------------------- Testing only -------------------------------------------------
         case 81:
            // Output profile data?
            strRH = strToLower(&strRH);

            m_bOutputProfileData = false;
            if (strRH.find("y") != string::npos)
               m_bOutputProfileData = true;

            // TODO 043 What about randomness of profile spacing, since profile location is determined by curvature?

            break;

         case 82:
            // Numbers of profiles to be saved
            if (m_bOutputProfileData)
            {
               VstrTmp = VstrSplit(&strRH, SPACE);
               for (unsigned int j = 0; j < VstrTmp.size(); j++)
               {
                  VstrTmp[j] = strTrim(&VstrTmp[j]);

                  if (! bIsStringValidInt(VstrTmp[j]))
                  {
                     strErr = "line " + to_string(nLine) + ": invalid integer for profile to be saved '" + VstrTmp[j] + "' in " + m_strDataPathName;
                     break;
                  }

                  int nTmp = stoi(VstrTmp[j]);

                  if (nTmp < 0)
                  {
                     strErr = "line " + to_string(nLine) + ": Profile number for saving must be >= 0";
                     break;
                  }

                  m_VnProfileToSave.push_back(nTmp);
               }
            }
            break;

         case 83:
            // Timesteps to save profiles
            if (m_bOutputProfileData)
            {
               VstrTmp = VstrSplit(&strRH, SPACE);
               for (unsigned int j = 0; j < VstrTmp.size(); j++)
               {
                  VstrTmp[j] = strTrim(&VstrTmp[j]);
                  unsigned long ulTmp = atol(VstrTmp[j].c_str());
                  if (ulTmp < 1)
                  {
                     strErr = "line " + to_string(nLine) + ": Timestep for profile saves must >= 1";
                     break;
                  }

                  m_VulProfileTimestep.push_back(ulTmp);
               }
            }
            break;

         case 84:
            // Output parallel profile data?
            strRH = strToLower(&strRH);

            m_bOutputParallelProfileData = false;
            if (strRH.find("y") != string::npos)
               m_bOutputParallelProfileData = true;
            break;

         case 85:
            // Output erosion potential look-up data?
            strRH = strToLower(&strRH);

            m_bOutputErosionPotentialData = false;
            if (strRH.find("y") != string::npos)
               m_bOutputErosionPotentialData = true;
            break;

         case 86:
            // Erode shore platform in alternate direction each timestep?
            strRH = strToLower(&strRH);

            m_bErodeShorePlatformAlternateDirection = false;
            if (strRH.find("y") != string::npos)
               m_bErodeShorePlatformAlternateDirection = true;
            break;

         case 87:
            // Size of moving window for coastline curvature calculation (must be odd)
            if (! bIsStringValidInt(strRH))
            {
               strErr = "line " + to_string(nLine) + ": invalid integer for size of moving window for coastline curvature calculation '" + strRH + "' in " + m_strDataPathName;
               break;
            }

            m_nCoastCurvatureMovingWindowSize = stoi(strRH);

            if ((m_nCoastCurvatureMovingWindowSize <= 0) || ! (m_nCoastCurvatureMovingWindowSize % 2))
               strErr = "line " + to_string(nLine) + ": size of moving window for coastline curvature calculation (must be > 0 and odd)";

            break;
         }

         // Did an error occur?
         if (! strErr.empty())
         {
            // Error in input to run details file
            cerr << endl
                 << ERR << strErr << ".\nPlease edit " << m_strDataPathName << " and change this line:" << endl
                 << "'" << strRec << "'" << endl
                 << endl;
            InStream.close();
            return false;
         }
      }
   }

   // Close file
   InStream.close();

   // Finally, need to check that we have at least one raster file, so that we know the grid size and units (and preferably also the projection)
   bool bNoRasterFiles = true;
   if ((! m_strInitialBasementDEMFile.empty()) || (! m_strInitialSuspSedimentFile.empty()) || (! m_strInitialLandformFile.empty()) || (! m_strInterventionHeightFile.empty()))
      bNoRasterFiles = false;

   for (int j = 0; j < m_nLayers; j++)
   {
      if ((! m_VstrInitialFineUnconsSedimentFile[j].empty()) || (! m_VstrInitialSandUnconsSedimentFile[j].empty()) || (! m_VstrInitialCoarseUnconsSedimentFile[j].empty()) || (! m_VstrInitialFineConsSedimentFile[j].empty()) || (! m_VstrInitialSandConsSedimentFile[j].empty()) || (! m_VstrInitialCoarseConsSedimentFile[j].empty()))
         bNoRasterFiles = false;
   }

   if (bNoRasterFiles)
   {
      // No raster files
      cerr << ERR << "at least one raster GIS file is needed" << endl;
      return false;
   }

   return true;
}

//===============================================================================================================================
//! Reads the tide time series data
//===============================================================================================================================
int CSimulation::nReadTideDataFile()
{
   // Create an ifstream object
   ifstream InStream;

   // Try to open the file for input
   InStream.open(m_strTideDataFile.c_str(), ios::in);

   // Did it open OK?
   if (!InStream.is_open())
   {
      // Error: cannot open tide data file for input
      cerr << ERR << "cannot open " << m_strTideDataFile << " for input" << endl;
      return RTN_ERR_TIDEDATAFILE;
   }

   // Opened OK
   int nLine = 0;
   string strRec;

   // Now read the data from the file
   while (getline(InStream, strRec))
   {
      nLine++;

      // Trim off leading and trailing whitespace
      strRec = strTrim(&strRec);

      // If it is a blank line or a comment then ignore it
      if ((strRec.empty()) || (strRec[0] == QUOTE1) || (strRec[0] == QUOTE2))
         continue;

      // Check that this is a valid double
      if (! bIsStringValidDouble(strRec))
      {
         cerr << ERR << "invalid floating point number for tide data '" << strRec << "' in " << m_strTideDataFile << endl;
         return RTN_ERR_TIDEDATAFILE;
      }

      // Convert to double then append the value to the vector
      m_VdTideData.push_back(strtod(strRec.c_str(), NULL));
   }

   // Close file
   InStream.close();

   return RTN_OK;
}

//===============================================================================================================================
//! Reads the shape of the erosion potential distribution (see shape function in Walkden & Hall, 2005)
//===============================================================================================================================
int CSimulation::nReadShapeFunctionFile()
{
   // Sort out the path and filename
   m_strSCAPEShapeFunctionFile = m_strCMEDir;
   m_strSCAPEShapeFunctionFile.append(SCAPE_DIR);
   m_strSCAPEShapeFunctionFile.append(SCAPE_SHAPE_FUNCTION_FILE);

   // Create an ifstream object
   ifstream InStream;

   // Try to open the file for input
   InStream.open(m_strSCAPEShapeFunctionFile.c_str(), ios::in);

   // Did it open OK?
   if (!InStream.is_open())
   {
      // Error: cannot open shape function file for input
      cerr << ERR << "cannot open " << m_strSCAPEShapeFunctionFile << " for input" << endl;
      return RTN_ERR_SCAPE_SHAPE_FUNCTION_FILE;
   }

   // Opened OK
   int nLine = 0;
   int nExpected = 0, nRead = 0;
   string strRec;

   // Read in the number of data lines expected
   InStream >> nExpected;

   // Set up the vectors to hold the input data
   vector<double>
       VdDepthOverDB,
       VdErosionPotential,
       VdErosionPotentialFirstDeriv;

   // Now read the rest of the data from the file to get the Erosion Potential Shape function
   while (getline(InStream, strRec))
   {
      nLine++;
   
      // Trim off leading and trailing whitespace
      strRec = strTrim(&strRec);

      // If it is a blank line or a comment then ignore it
      if ((strRec.empty()) || (strRec[0] == QUOTE1) || (strRec[0] == QUOTE2))
         continue;

      // It isn't so increment counter
      nRead++;

      // Split the string, and remove whitespace
      vector<string> strTmp = VstrSplit(&strRec, SPACE);
      for (unsigned int i = 0; i < strTmp.size(); i++)
      {
         // Remove leading and trailing whitespace
         strTmp[i] = strTrim(&strTmp[i]);

         // Check that this is a valid double
         if (! bIsStringValidDouble(strTmp[i]))
         {
            cerr << ERR << "on line " + to_string(nLine) + "  invalid floating point number for Erosion Potential Shape data '" << strTmp[i] << "' in " << m_strSCAPEShapeFunctionFile << endl;
            return RTN_ERR_SCAPE_SHAPE_FUNCTION_FILE;
         }
      }

      // Convert to doubles then append the values to the vectors
      VdDepthOverDB.push_back(strtod(strTmp[0].c_str(), NULL));
      VdErosionPotential.push_back(strtod(strTmp[1].c_str(), NULL));
      VdErosionPotentialFirstDeriv.push_back(strtod(strTmp[2].c_str(), NULL));
   }
   // Now create the look up table values
   m_VdErosionPotential = VdErosionPotential;
   m_VdDepthOverDB = VdDepthOverDB;
   m_dDepthOverDBMax = VdDepthOverDB[14];

   // Close file
   InStream.close();

   // Did we read in what we expected?
   if (nExpected != nRead)
   {
      cout << ERR << "read in " << nRead << " lines from " << m_strSCAPEShapeFunctionFile << " but " << nExpected << " lines expected" << endl;
      return RTN_ERR_SCAPE_SHAPE_FUNCTION_FILE;
   }

   // Is the shape funcion well defined? i.e. it must be -ve or 0.0 for all values
   for (unsigned int j = 0; j < m_VdErosionPotential.size(); j++)
   {
      if (m_VdErosionPotential[j] > 0)
      {
         cout << ERR << " in " << m_strSCAPEShapeFunctionFile << ", erosion potential function cannot be positive" << endl;
         return RTN_ERR_SCAPE_SHAPE_FUNCTION_FILE;
      }
   }

   // OK, now use this data to create a look-up table to be used for the rest of the simulation
   if (! bCreateErosionPotentialLookUp(&VdDepthOverDB, &VdErosionPotential, &VdErosionPotentialFirstDeriv))
   {
      cout << ERR << "line " + to_string(nLine) + "  in " << m_strSCAPEShapeFunctionFile << ": erosion potential function is unbounded for high values of depth over DB" << endl;
      return RTN_ERR_SCAPE_SHAPE_FUNCTION_FILE;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Reads the deep water wave station time series data
//===============================================================================================================================
int CSimulation::nReadWaveStationTimeSeriesFile(int const nWaveStations)
{
   // Create an ifstream object
   ifstream InStream;

   // Try to open the file for input
   InStream.open(m_strDeepWaterWavesTimeSeriesFile.c_str(), ios::in);

   // Did it open OK?
   if (!InStream.is_open())
   {
      // Error: cannot open time series file for input
      cerr << ERR << "cannot open " << m_strDeepWaterWavesTimeSeriesFile << " for input" << endl;
      return RTN_ERR_READING_SEDIMENT_INPUT_EVENT;
   }

   // Opened OK
   int
       nLine = 0,
       nExpectedStations = 0,
       nRead = 0,
       nTimeStepsRead = 0;
   string strRec, strErr;

   // Read each line, ignoring comment lines
   while (getline(InStream, strRec))
   {
      nLine++;

      // Trim off leading and trailing whitespace
      strRec = strTrim(&strRec);

      // If it is a blank line or a comment then ignore it
      if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2))
      {
         // It isn't so increment counter
         nRead++;

         // The header lines (the first four lines of the file) contains leading description separated by a colon from the data
         if (nRead < 5)
         {
            // Find the colon: note that lines MUST have a colon separating data from leading description portion
            size_t nPos = strRec.find(COLON);
            if (nPos == string::npos)
            {
               // Error: badly formatted (no colon)
               cerr << ERR << "on line " << to_string(nLine) << "badly formatted (no ':') in " << m_strDeepWaterWavesTimeSeriesFile << endl
                    << "'" << strRec << "'" << endl;
               return RTN_ERR_READING_DEEP_WATER_WAVE_DATA;
            }

            if (nPos == strRec.size() - 1)
            {
               // Error: badly formatted (colon with nothing following)
               cerr << ERR << "on line " << to_string(nLine) << "badly formatted (nothing following ':') in " << m_strDeepWaterWavesTimeSeriesFile << endl
                    << "'" << strRec << "'" << endl;
               return RTN_ERR_READING_DEEP_WATER_WAVE_DATA;
            }

            // Strip off leading portion (the bit up to and including the colon)
            string strRH = strRec.substr(nPos + 1);

            // Remove leading whitespace
            strRH = strTrimLeft(&strRH);

            // Look for a trailing comment, if found then terminate string at that point and trim off any trailing whitespace
            nPos = strRH.rfind(QUOTE1);
            if (nPos != string::npos)
               strRH.resize(nPos);

            nPos = strRH.rfind(QUOTE2);
            if (nPos != string::npos)
               strRH.resize(nPos);

            // Remove trailing whitespace
            strRH = strTrimRight(&strRH);

            int
                nSec = 0,
                nMin = 0,
                nHour = 0,
                nDay = 0,
                nMonth = 0,
                nYear = 0;
            double
                dMult,
                dThisIter;
            vector<string> VstrTmp;

            switch (nRead)
            {
            case 1:
               // Get the start date/time for this data, format is [hh-mm-ss dd/mm/yyyy]
               VstrTmp = VstrSplit(&strRH, SPACE);

               // Both date and time here?
               if (VstrTmp.size() < 2)
               {
                  strErr = "line " + to_string(nLine) + ": must have both date and time for start of data in";
                  break;
               }

               // OK, first sort out the time
               if (! bParseTime(&VstrTmp[0], nHour, nMin, nSec))
               {
                  strErr = "line " + to_string(nLine) + ": could not understand start time for data";
                  break;
               }

               // Next sort out the date
               if (! bParseDate(&VstrTmp[1], nDay, nMonth, nYear))
               {
                  strErr = "line " + to_string(nLine) + ": could not understand start date for data";
                  break;
               }

               // Compare time and date with simulation time and date
               if ((nSec != m_nSimStartSec) ||
                   (nMin != m_nSimStartMin) ||
                   (nHour != m_nSimStartHour) ||
                   (nDay != m_nSimStartDay) ||
                   (nMonth != m_nSimStartMonth) ||
                   (nYear != m_nSimStartYear))
               {
                  strErr = "line " + to_string(nLine) + ": start time and date for wave time series data differs from simulation start time and date,";
                  break;
               }

               break;

            case 2:
               // Get the timestep of this data (in hours or days)
               strRH = strToLower(&strRH);

               dMult = dGetTimeMultiplier(&strRH);
               if (static_cast<int>(dMult) == TIME_UNKNOWN)
               {
                  strErr = "line " + to_string(nLine) + ": unknown units for timestep";
                  break;
               }

               // We have the multiplier, now calculate the timestep in hours: look for the whitespace between the number and unit
               nPos = strRH.rfind(SPACE);
               if (nPos == string::npos)
               {
                  strErr = "line " + to_string(nLine) + ": format of timestep line";
                  break;
               }

               // Cut off rh bit of string
               strRH.resize(nPos);

               // Remove trailing spaces
               strRH = strTrimRight(&strRH);

               // Check that this is a valid double
               if (! bIsStringValidDouble(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for timestep";
                  break;
               }

               dThisIter = strtod(strRH.c_str(), NULL) * dMult; // in hours

               if (dThisIter <= 0)
                  strErr = "line " + to_string(nLine) + ": timestep must be > 0";

               if (! bFPIsEqual(dThisIter, m_dTimeStep, TOLERANCE))
                  strErr = "line " + to_string(nLine) + ": timestep must be the same as the simulation timestep";

               break;

            case 3:
               // Read the number of stations
               if (! bIsStringValidInt(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid integer for number of wave stations '" + strRH + "' in " + m_strDeepWaterWavesTimeSeriesFile;
                  break;
               }

               nExpectedStations = stoi(strRH);

               // Check that the number of expected stations is equal to the number of stations on the point shape file
               if (nExpectedStations != nWaveStations)
               {
                  // Error: number of points on shape file does not match the number of stations on the wave time series file
                  strErr = "line " + to_string(nLine) + ": number of wave stations in " + m_strDeepWaterWaveStationsShapefile + " is " + to_string(nWaveStations) + " but we have " + to_string(nExpectedStations) + " stations";

                  break;
               }

               break;

            case 4:
               // Read the expected number of time steps in the file
               if (! bIsStringValidInt(strRH))
               {
                  strErr = "line " + to_string(nLine) + ": invalid integer for expected number of time steps '" + strRH + "' in " + m_strDeepWaterWaveStationsShapefile;
                  break;
               }

               m_nDeepWaterWaveDataNumTimeSteps = stoi(strRH);

               if (m_nDeepWaterWaveDataNumTimeSteps < 1)
               {
                  // Error: must have value(s) for at least one timestep
                  strErr = "line " + to_string(nLine) + ": must have values for at least one timestep";
                  break;
               }

               break;
            }
         }
         else
         {
            // This is not a header line
            nTimeStepsRead++;

            // Read in each wave attribute for each time step and station: split the string, and remove whitespace
            vector<string> VstrTmp = VstrSplit(&strRec, COMMA);
            for (unsigned int i = 0; i < VstrTmp.size(); i++) // VstrTmp.size() should be 3 x nExpectedStations
            {
               // Remove leading and trailing whitespace
               VstrTmp[i] = strTrim(&VstrTmp[i]);

               // Check that this is a valid double
               if (! bIsStringValidDouble(VstrTmp[i]))
               {
                  strErr = "line " + to_string(nLine) + ": invalid floating point number for deep water wave value '" + VstrTmp[i] + "' in " + m_strDeepWaterWavesTimeSeriesFile;
                  break;
               }
            }

            // Convert to doubles then append the values to the vectors
            int n = 0;
            for (int i = 0; i < nExpectedStations; i++)
            {
               m_VdTSDeepWaterWaveStationHeight.push_back(strtod(VstrTmp[n + 0].c_str(), NULL));
               m_VdTSDeepWaterWaveStationAngle.push_back(strtod(VstrTmp[n + 1].c_str(), NULL));
               m_VdTSDeepWaterWaveStationPeriod.push_back(strtod(VstrTmp[n + 2].c_str(), NULL));

               // Check some simple wave input stats
               if (m_VdTSDeepWaterWaveStationHeight.back() > m_dMaxUserInputWaveHeight)
                  m_dMaxUserInputWaveHeight = m_VdTSDeepWaterWaveStationHeight.back();

               if (m_VdTSDeepWaterWaveStationPeriod.back() > m_dMaxUserInputWavePeriod)
                  m_dMaxUserInputWavePeriod = m_VdTSDeepWaterWaveStationPeriod.back();

               n += 3;
            }
         }
      }

      // Did an error occur?
      if (! strErr.empty())
      {
         // Error in input to initialisation file
         cerr << ERR << strErr << " in deep water wave time series file " << m_strDeepWaterWavesTimeSeriesFile << endl
              << "'" << strRec << "'" << endl;
         InStream.close();

         return RTN_ERR_READING_DEEP_WATER_WAVE_DATA;
      }
   }

   if (nTimeStepsRead != m_nDeepWaterWaveDataNumTimeSteps)
   {
      // Error: number of timesteps read does not match the number given in the file's header
      cerr << ERR << "in " << m_strDeepWaterWavesTimeSeriesFile << ", data for " << nTimeStepsRead << " timesteps was read, but " << m_nDeepWaterWaveDataNumTimeSteps << " timesteps were specified in the file's header" << endl;

      return RTN_ERR_READING_DEEP_WATER_WAVE_DATA;
   }

   // Close file
   InStream.close();

   // Did we read in what we expected?
   unsigned int nTotExpected = nExpectedStations * m_nDeepWaterWaveDataNumTimeSteps;
   if (m_VdTSDeepWaterWaveStationHeight.size() != nTotExpected)
   {
      cout << ERR << "read in " << m_VdTSDeepWaterWaveStationHeight.size() << " lines from " << m_strDeepWaterWavesTimeSeriesFile << " but " << nTotExpected << " values expected" << endl;

      return RTN_ERR_READING_DEEP_WATER_WAVE_DATA;
   }

   if (m_VdTSDeepWaterWaveStationAngle.size() != nTotExpected)
   {
      cout << ERR << "read in " << m_VdTSDeepWaterWaveStationAngle.size() << " lines from " << m_strDeepWaterWavesTimeSeriesFile << " but " << nTotExpected << " values expected" << endl;

      return RTN_ERR_READING_DEEP_WATER_WAVE_DATA;
   }

   if (m_VdTSDeepWaterWaveStationPeriod.size() != nTotExpected)
   {
      cout << ERR << "read in " << m_VdTSDeepWaterWaveStationPeriod.size() << " lines from " << m_strDeepWaterWavesTimeSeriesFile << " but " << nTotExpected << " values expected" << endl;

      return RTN_ERR_READING_DEEP_WATER_WAVE_DATA;
   }

   // All is OK, so we can now initialize the vectors that will store this timestep's deep water wave values
   for (int j = 0; j < nExpectedStations; j++)
   {
      m_VdThisIterDeepWaterWaveStationHeight.push_back(DBL_NODATA);
      m_VdThisIterDeepWaterWaveStationAngle.push_back(DBL_NODATA);
      m_VdThisIterDeepWaterWaveStationPeriod.push_back(DBL_NODATA);
   }

   // Finally, check whether the wave data will 'wrap' i.e. whether the number of timesteps is less than the total number of timesteps in the simulation
   int nSimulationTimeSteps = static_cast<int>(floor(m_dSimDuration / m_dTimeStep));
   if (m_nDeepWaterWaveDataNumTimeSteps < nSimulationTimeSteps)
   {
      m_dWaveDataWrapHours = m_nDeepWaterWaveDataNumTimeSteps * m_dTimeStep;
      string strTmp = "Deep water wave data will wrap every " + (m_nDeepWaterWaveDataNumTimeSteps > 1 ? to_string(m_nDeepWaterWaveDataNumTimeSteps) + " " : "") + "time step" + (m_nDeepWaterWaveDataNumTimeSteps > 1 ? "s" : "") + " (every " + to_string(m_dWaveDataWrapHours) + " hours)\n";

      cout << NOTE << strTmp;
   }

   return RTN_OK;
}

//===============================================================================================================================
//! Reads the sediment input events time series data
//===============================================================================================================================
int CSimulation::nReadSedimentInputEventTimeSeriesFile(void)
{
   // Create an ifstream object
   ifstream InStream;

   // Try to open the file for input
   InStream.open(m_strSedimentInputEventTimeSeriesFile.c_str(), ios::in);

   // Did it open OK?
   if (!InStream.is_open())
   {
      // Error: cannot open time series file for input
      cerr << ERR << "cannot open " << m_strSedimentInputEventTimeSeriesFile << " for input" << endl;
      return RTN_ERR_READING_SEDIMENT_INPUT_EVENT;
   }

   // Opened OK
   int
      nLine = 0,
      nRead = 0;
   string strRec, strErr;

   // Read each line, ignoring comment lines
   while (getline(InStream, strRec))
   {
      nLine++;

      // Trim off leading and trailing whitespace
      strRec = strTrim(&strRec);

      // If it is a blank line or a comment then ignore it
      if ((! strRec.empty()) && (strRec[0] != QUOTE1) && (strRec[0] != QUOTE2))
      {
         // It isn't so increment counter
         nRead++;

         // Split at commas
         vector<string> VstrTmp = VstrSplit(&strRec, COMMA);

         // Check that have all we need
         unsigned int nTarget = 7;
         if (m_bSedimentInputAlongLine || m_bSedimentInputAtPoint)
            nTarget = 5;
         if (VstrTmp.size() < nTarget)
         {
            strErr = "line " + to_string(nLine) + ": too few data items on data line '" + to_string(nRead) + "' in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         // First item is the Location ID of the sediment input event (same as the ID in the shapefile)
         if (! bIsStringValidInt(VstrTmp[0]))
         {
            strErr = "line " + to_string(nLine) + ": invalid integer for Location ID of sediment input event '" + VstrTmp[0] + "' in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         int nID = stoi(strTrim(&VstrTmp[0]));

         // OK, check the ID against IDs read in from the shapefile
         auto result = find(m_VnSedimentInputLocationID.begin(), m_VnSedimentInputLocationID.end(), nID);
         if (result == m_VnSedimentInputLocationID.end())
         {
            strErr = "line " + to_string(nLine) + ": invalid Location ID '" + to_string(nID) + "' for sediment input event location event in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         // Next get the timestep at which the sediment input event occurs. This may be specified either as a relative time (i.e. a number of hours or days after the simulation start) or as an absolute time (i.e. a time/date in the format hh-mm-ss dd/mm/yyyy)
         unsigned long ulEventTimeStep = ulConvertToTimestep(&VstrTmp[1]);
         if (ulEventTimeStep == SEDIMENT_INPUT_EVENT_ERROR)
         {
            strErr = "line " + to_string(nLine) + ": invalid time and/or date '" + VstrTmp[1] + "' for sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         // Then the volume (m3) of fine sediment, first check that this is a valid double
         if (! bIsStringValidDouble(VstrTmp[2]))
         {
            strErr = "line " + to_string(nLine) + ": invalid floating point number '" + VstrTmp[2] + "' for fine sediment volume for sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         double dFineSedVol = stod(strTrim(&VstrTmp[2]));
         if (dFineSedVol < 0)
         {
            strErr = "line " + to_string(nLine) + ": negative number '" + to_string(dFineSedVol) + "' for fine sediment volume for sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         if (dFineSedVol > 0)
            m_bHaveFineSediment = true;

         // Then the volume (m3) of sand sediment, first check that this is a valid double
         if (! bIsStringValidDouble(VstrTmp[3]))
         {
            strErr = "line " + to_string(nLine) + ": invalid floating point number '" + VstrTmp[3] + "' for sand-sized sediment volume for sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         double dSandSedVol = stod(strTrim(&VstrTmp[3]));
         if (dSandSedVol < 0)
         {
            strErr = "line " + to_string(nLine) + ": negative number '" + to_string(dSandSedVol) + "' for sand-sized sediment volume for sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         if (dSandSedVol > 0)
            m_bHaveSandSediment = true;

         // Then the volume (m3) of coarse sediment, first check that this is a valid double
         if (! bIsStringValidDouble(VstrTmp[4]))
         {
            strErr = "line " + to_string(nLine) + ": invalid floating point number '" + VstrTmp[4] + "' for coarse sediment volume for sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         double dCoarseSedVol = stod(strTrim(&VstrTmp[4]));
         if (dCoarseSedVol < 0)
         {
            strErr = "line " + to_string(nLine) + ": negative number '" + to_string(dCoarseSedVol) + "' for coarse sediment volume of sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
            break;
         }

         if (dCoarseSedVol > 0)
            m_bHaveCoarseSediment = true;

         // Only read the last two items if we have on-coast sediment block sediment input
         double
             dLen = 0,
             dWidth = 0;
         //  dThick = 0;
         if (m_bSedimentInputAtCoast)
         {

            // The coast-normal length (m) of the sediment block, first check that this is a valid double
            if (! bIsStringValidDouble(VstrTmp[5]))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number '" + VstrTmp[5] + "' for coast-normal length of sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
               break;
            }

            dLen = stod(strTrim(&VstrTmp[5]));
            if (dLen <= 0)
            {
               strErr = "line " + to_string(nLine) + ": coast-normal length L of the sediment block '" + to_string(dLen) + "' must be > 0 in " + m_strSedimentInputEventTimeSeriesFile;
               break;
            }

            // The along-coast width (m) of the sediment block, first check that this is a valid double
            if (! bIsStringValidDouble(VstrTmp[6]))
            {
               strErr = "line " + to_string(nLine) + ": invalid floating point number '" + VstrTmp[6] + "' for along-coast width of sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
               break;
            }

            dWidth = stod(strTrim(&VstrTmp[6]));
            if ((! m_bSedimentInputAtCoast) && (dWidth <= 0))
            {
               strErr = "line " + to_string(nLine) + ": along-coast width W (m) of the sediment block '" + to_string(dWidth) + "' must be > 0 in " + m_strSedimentInputEventTimeSeriesFile;
               break;
            }

            // The along-coast thickness (m) of the sediment block, first check that this is a valid double
            // if (! bIsStringValidDouble(VstrTmp[7]))
            // {
            //    strErr = "line " + to_string(nLine) + ": invalid floating point number '" + VstrTmp[7] + "' for along-coast thickness of sediment input event in " + m_strSedimentInputEventTimeSeriesFile;
            //    break;
            // }

            // dThick = stod(strTrim(&VstrTmp[7]));
            // if ((! m_bSedimentInputAtCoast) && (dWidth <= 0))
            // {
            //    strErr = "line " + to_string(nLine) + ": along-coast thickness E (m) of the sediment block '" + to_string(dThick) + "' must be > 0 in " + m_strSedimentInputEventTimeSeriesFile;
            //    break;
            // }
         }

         // Create the CSedInputEvent object
         CSedInputEvent* pEvent = new CSedInputEvent(nID, ulEventTimeStep, dFineSedVol, dSandSedVol, dCoarseSedVol, dLen, dWidth); //, dThick);

         // And store it in the m_pVSedInputEvent 2D vector
         m_pVSedInputEvent.push_back(pEvent);
      }
   }

   // Did an error occur?
   if (! strErr.empty())
   {
      // Error in input to initialisation file
      cerr << ERR << strErr << " in sediment input event file " << m_strSedimentInputEventTimeSeriesFile << endl
           << "'" << strRec << "'" << endl;
      InStream.close();

      return RTN_ERR_READING_SEDIMENT_INPUT_EVENT;
   }

   // Close file
   InStream.close();

   return RTN_OK;
}
