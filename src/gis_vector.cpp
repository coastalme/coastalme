/*!
 *
 * \file gis_vector.cpp
 * \brief These functions use GDAL to read and write vector GIS files in several formats. This version will build with GDAL version 2
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
#include <cfloat>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::ios;

#include <sstream>
using std::stringstream;

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "cme.h"
#include "simulation.h"
#include "coast.h"
#include "cliff.h"

//===============================================================================================================================
//! Reads vector GIS datafiles using OGR
//===============================================================================================================================
int CSimulation::nReadVectorGISFile(int const nDataItem)
{
   int
       nMaxLayer = 0,
       nNeedGeometry = 0;

   string
       strGISFile,
       strGeometry;

   // Set up file name and constraints
   switch (nDataItem)
   {
      case (DEEP_WATER_WAVE_STATIONS_VEC):
         strGISFile = m_strDeepWaterWaveStationsShapefile;
         nMaxLayer = DEEP_WATER_WAVE_STATIONS_MAX_LAYER;
         nNeedGeometry = DEEP_WATER_WAVE_STATIONS_POINT_GEOMETRY;
         break;

      case (SEDIMENT_INPUT_EVENT_LOCATION_VEC):
         strGISFile = m_strSedimentInputEventShapefile;
         nMaxLayer = SEDIMENT_INPUT_EVENT_LOCATION_MAX_LAYER;
         if (m_bSedimentInputAtPoint || m_bSedimentInputAtCoast)
            nNeedGeometry = SEDIMENT_INPUT_EVENT_LOCATION_POINT_GEOMETRY;
         else if (m_bSedimentInputAlongLine)
            nNeedGeometry = SEDIMENT_INPUT_EVENT_LOCATION_LINE_GEOMETRY;
         break;

      case (FLOOD_LOCATION_VEC):
         strGISFile = m_strFloodLocationShapefile;
         nMaxLayer = FLOOD_LOCATION_MAX_LAYER;
         nNeedGeometry = FLOOD_LOCATION_POINT_GEOMETRY;
         break;
   }

   // Open the GDAL/OGR datasource
   GDALDataset* pOGRDataSource = static_cast<GDALDataset*>(GDALOpenEx(strGISFile.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL));
   if (pOGRDataSource == NULL)
   {
      // Can't open file (note will already have sent GDAL error message to stdout)
      cerr << ERR << "cannot open " << strGISFile << " for input: " << CPLGetLastErrorMsg() << endl;
      return RTN_ERR_VECTOR_FILE_READ;
   }

   // Find out number of layers, and compare with the required number
   int nLayer = pOGRDataSource->GetLayerCount();
   if (nLayer > nMaxLayer)
      LogStream << WARN << "need " << nMaxLayer << (nMaxLayer > 1 ? "layers" : "layer") << " in " << strGISFile << ", " << nLayer << " found. Only the first " << nMaxLayer << (nMaxLayer > 1 ? "layers" : "layer") << " will be read." << endl;

   for (int n = 0; n < nMaxLayer; n++)
   {
      // Open this layer
      OGRLayer* pOGRLayer;
      pOGRLayer = pOGRDataSource->GetLayer(n);

      // Get features from the layer
      OGRFeature* pOGRFeature;

      // Make sure we are at the beginning of the layer
      pOGRLayer->ResetReading();

      // Now iterate for all features in the layer
      while ((pOGRFeature = pOGRLayer->GetNextFeature()) != NULL)
      {
         // First get the geometry for this feature
         OGRGeometry* pOGRGeometry;
         pOGRGeometry = pOGRFeature->GetGeometryRef();
         if (pOGRGeometry == NULL)
         {
            cerr << ERR << " null geometry in " << strGISFile << "." << endl;
            return RTN_ERR_VECTOR_FILE_READ;
         }

         // Now get the geometry type
         int
             nGeometry = wkbFlatten(pOGRGeometry->getGeometryType()),
             nThisGeometry = 0;

         switch (nGeometry)
         {
            case wkbPoint:
               nThisGeometry = VEC_GEOMETRY_POINT;
               strGeometry = "point";
               break;

            case wkbLineString:
               nThisGeometry = VEC_GEOMETRY_LINE;
               strGeometry = "line";
               break;

            case wkbPolygon:
               nThisGeometry = VEC_GEOMETRY_POLYGON;
               strGeometry = "polygon";
               break;

            default:
               nThisGeometry = VEC_GEOMETRY_OTHER;
               strGeometry = "other";
               break;
         }

         // Have we got the expected geometry type?
         if (nThisGeometry != nNeedGeometry)
         {
            // Error, we do not have the desired geometry
            string strNeedGeometry;
            switch (nNeedGeometry)
            {
               case VEC_GEOMETRY_POINT:
                  strNeedGeometry = "point";
                  break;

               case VEC_GEOMETRY_LINE:
                  strNeedGeometry = "line";
                  break;

               case VEC_GEOMETRY_POLYGON:
                  strNeedGeometry = "polygon";
                  break;

               case VEC_GEOMETRY_OTHER:
                  strNeedGeometry = "other";
                  break;
            }

            cerr << strGeometry << " data found in " << strGISFile << ", but " << strNeedGeometry << " data is needed" << endl;
            return RTN_ERR_VECTOR_FILE_READ;
         }

         // The geometry type is OK, so process the geometry data
         int nPoints = 0;
         OGRPoint* pOGRPoint;
         OGRLineString* pOGRLineString;
         switch (nDataItem)
         {
            case (DEEP_WATER_WAVE_STATIONS_VEC):
               // Point data
               pOGRPoint = static_cast<OGRPoint*>(pOGRGeometry);

               // Convert the wave station co-ordinates to grid CRS and store them: we will use these in the spatial interpolation of deep water waves
               m_VdDeepWaterWaveStationX.push_back(dExtCRSXToGridX(pOGRPoint->getX()));
               m_VdDeepWaterWaveStationY.push_back(dExtCRSYToGridY(pOGRPoint->getY()));
               break;

            case (SEDIMENT_INPUT_EVENT_LOCATION_VEC):
               if (m_bSedimentInputAtPoint || m_bSedimentInputAtCoast)
               {
                  // Point data
                  pOGRPoint = static_cast<OGRPoint*>(pOGRGeometry);

                  int
                     nPointGridX = nRound(dExtCRSXToGridX(pOGRPoint->getX())),
                     nPointGridY = nRound(dExtCRSYToGridY(pOGRPoint->getY()));

                  // Check point data is inside the mesh
                  if ((nPointGridX < 0) || (nPointGridX > m_nXGridMax))
                     return RTN_ERR_SEDIMENT_INPUT_EVENT_LOCATION;

                  if ((nPointGridY < 0) || (nPointGridY > m_nYGridMax))
                     return RTN_ERR_SEDIMENT_INPUT_EVENT_LOCATION;

                  // Convert the point co-ordinates to grid CRS and store them
                  m_VdSedimentInputLocationX.push_back(nPointGridX);
                  m_VdSedimentInputLocationY.push_back(nPointGridY);
               }
               else if (m_bSedimentInputAlongLine)
               {
                  // Line data
                  pOGRLineString = static_cast<OGRLineString*>(pOGRGeometry);

                  nPoints = pOGRLineString->getNumPoints();
                  for (int i = 0; i < nPoints; i++)
                  {
                     // Convert the co-ordinates for each point in the line to grid CRS and store them
                     m_VdSedimentInputLocationX.push_back(dExtCRSXToGridX(pOGRLineString->getX(i)));
                     m_VdSedimentInputLocationY.push_back(dExtCRSYToGridY(pOGRLineString->getY(i)));
                  }
               }
               break;

            case (FLOOD_LOCATION_VEC):

               // Point data
               pOGRPoint = static_cast<OGRPoint*>(pOGRGeometry);

               // Convert the wave station co-ordinates to grid CRS and store them: we will use these in the spatial interpolation of deep water waves
               m_VdFloodLocationX.push_back(pOGRPoint->getX());
               m_VdFloodLocationY.push_back(pOGRPoint->getY());

               break;
         }

         // Now get the attributes of this feature
         OGRFeatureDefn* pOGRFeatureDefn = pOGRLayer->GetLayerDefn();

         int
             nFieldIndex = -1,
             nID = -1;

         switch (nDataItem)
         {
            case (DEEP_WATER_WAVE_STATIONS_VEC):
               // First get the station ID
               nFieldIndex = pOGRFeatureDefn->GetFieldIndex(DEEP_WATER_WAVE_STATION_ID.c_str());
               if (nFieldIndex == -1)
               {
                  // Can't find this field in the vector file
                  cerr << ERR << "cannot find " << DEEP_WATER_WAVE_STATION_ID << " field in " << strGISFile << ": " << CPLGetLastErrorMsg() << endl;
                  return RTN_ERR_VECTOR_FILE_READ;
               }

               // Get the Station ID for this point
               nID = pOGRFeature->GetFieldAsInteger(nFieldIndex);
               m_VnDeepWaterWaveStationID.push_back(nID);

               break;

            case (SEDIMENT_INPUT_EVENT_LOCATION_VEC):
               // First get the station ID
               nFieldIndex = pOGRFeatureDefn->GetFieldIndex(SEDIMENT_INPUT_EVENT_LOCATION_ID.c_str());
               if (nFieldIndex == -1)
               {
                  // Can't find this field in the vector file
                  cerr << ERR << "cannot find " << SEDIMENT_INPUT_EVENT_LOCATION_ID << " field in " << strGISFile << ": " << CPLGetLastErrorMsg() << endl;
                  return RTN_ERR_VECTOR_FILE_READ;
               }

               // Get the Station ID for this point (note that we read it as an integer, not caring what type of field is actually in the shapefile)
               nID = pOGRFeature->GetFieldAsInteger(nFieldIndex);

               if (m_bSedimentInputAtPoint || m_bSedimentInputAtCoast)
                  // Save the ID for the point
                  m_VnSedimentInputLocationID.push_back(nID);

               else if (m_bSedimentInputAlongLine)
               {
                  // Save the ID for every point in the line
                  for (int i = 0; i < nPoints; i++)
                     m_VnSedimentInputLocationID.push_back(nID);
               }

               break;

            case (FLOOD_LOCATION_VEC):
               // First get the station ID
               nFieldIndex = pOGRFeatureDefn->GetFieldIndex(FLOOD_LOCATION_ID.c_str());
               if (nFieldIndex == -1)
               {
                  // Can't find this field in the vector file
                  cerr << ERR << "cannot find " << FLOOD_LOCATION_ID << " field in " << strGISFile << ": " << CPLGetLastErrorMsg() << endl;
                  return RTN_ERR_VECTOR_FILE_READ;
               }

               // Get the Station ID for this point (note that we read it as an integer, not caring what type of field is actually in the shapefile)
               nID = pOGRFeature->GetFieldAsInteger(nFieldIndex);

               // Save the ID for the point
               m_VnFloodLocationID.push_back(nID);
               break;
         }
      }

      // Get rid of the Feature object
      OGRFeature::DestroyFeature(pOGRFeature);
   }

   // Save some info, to be shown in the text output
   switch (nDataItem)
   {
      case (DEEP_WATER_WAVE_STATIONS_VEC):
         m_strOGRDWWVDriverCode = pOGRDataSource->GetDriverName();
         m_strOGRDWWVDataType = "integer";
         m_strOGRDWWVGeometry = strGeometry;
         break;

      case (SEDIMENT_INPUT_EVENT_LOCATION_VEC):
         m_strOGRSedInputDriverCode = pOGRDataSource->GetDriverName();
         m_strOGRSedInputGeometry = "integer";
         m_strOGRSedInputDataType = strGeometry;
         break;

      case (FLOOD_LOCATION_VEC):
         m_strOGRFloodDriverCode = pOGRDataSource->GetDriverName();
         m_strOGRFloodGeometry = "integer";
         m_strOGRFloodDataType = strGeometry;
         break;
   }

   // Clean up: get rid of the data source object
   GDALClose(pOGRDataSource);

   return RTN_OK;
}

//===============================================================================================================================
//! Writes vector GIS files using OGR
//===============================================================================================================================
bool CSimulation::bWriteVectorGISFile(int const nDataItem, string const *strPlotTitle)
{
   // Begin constructing the file name for this save
   string strFilePathName(m_strOutPath);
   stringstream strstrFileName;

   switch (nDataItem)
   {
      case (VECTOR_PLOT_COAST):
         strFilePathName.append(VECTOR_COAST_NAME);
         strstrFileName << VECTOR_COAST_NAME;
         break;

      case (VECTOR_PLOT_NORMALS):
         strFilePathName.append(VECTOR_NORMALS_NAME);
         strstrFileName << VECTOR_NORMALS_NAME;
         break;

      case (VECTOR_PLOT_INVALID_NORMALS):
         strFilePathName.append(VECTOR_INVALID_NORMALS_NAME);
         break;

      case (VECTOR_PLOT_COAST_CURVATURE):
         strFilePathName.append(VECTOR_COAST_CURVATURE_NAME);
         break;

      case (VECTOR_PLOT_WAVE_ANGLE_AND_HEIGHT):
         strFilePathName.append(VECTOR_WAVE_ANGLE_AND_HEIGHT_NAME);
         break;

      case (VECTOR_PLOT_AVG_WAVE_ANGLE_AND_HEIGHT):
         strFilePathName.append(VECTOR_AVG_WAVE_ANGLE_AND_HEIGHT_NAME);
         break;

      case (VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE):
         strFilePathName.append(VECTOR_WAVE_ENERGY_SINCE_COLLAPSE_NAME);
         break;

      case (VECTOR_PLOT_MEAN_WAVE_ENERGY):
         strFilePathName.append(VECTOR_MEAN_WAVE_ENERGY_NAME);
         break;

      case (VECTOR_PLOT_BREAKING_WAVE_HEIGHT):
         strFilePathName.append(VECTOR_BREAKING_WAVE_HEIGHT_NAME);
         break;

      case (VECTOR_PLOT_POLYGON_NODES):
         strFilePathName.append(VECTOR_POLYGON_NODE_NAME);
         break;

      case (VECTOR_PLOT_POLYGON_BOUNDARY):
         strFilePathName.append(VECTOR_POLYGON_BOUNDARY_NAME);
         break;

      case (VECTOR_PLOT_CLIFF_NOTCH_SIZE):
         strFilePathName.append(VECTOR_CLIFF_NOTCH_SIZE_NAME);
         break;

      case (VECTOR_PLOT_SHADOW_BOUNDARY):
         strFilePathName.append(VECTOR_SHADOW_BOUNDARY_NAME);
         break;

      case (VECTOR_PLOT_DOWNDRIFT_BOUNDARY):
         strFilePathName.append(VECTOR_DOWNDRIFT_BOUNDARY_NAME);
         break;

      case (VECTOR_PLOT_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT):
         strFilePathName.append(VECTOR_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT_NAME);
         break;

      case (VECTOR_PLOT_WAVE_SETUP):
         strFilePathName.append(VECTOR_WAVE_SETUP_NAME);
         break;

      case (VECTOR_PLOT_STORM_SURGE):
         strFilePathName.append(VECTOR_STORM_SURGE_NAME);
         break;

      case (VECTOR_PLOT_RUN_UP):
         strFilePathName.append(VECTOR_RUN_UP_NAME);
         break;

      case (VECTOR_PLOT_FLOOD_LINE):
         strFilePathName.append(VECTOR_FLOOD_LINE_NAME);
         strstrFileName << VECTOR_FLOOD_LINE_NAME;
         break;

         // case (VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_LINE):
         //    strFilePathName.append(VECTOR_FLOOD_SWL_SETUP_SURGE_LINE_NAME);
         //    strstrFileName << VECTOR_FLOOD_SWL_SETUP_SURGE_LINE_NAME;
         //    break;

         // case (VECTOR_PLOT_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE):
         //    strFilePathName.append(VECTOR_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE_NAME);
         //    strstrFileName << VECTOR_FLOOD_SWL_SETUP_SURGE_RUNUP_LINE_NAME;
         //    break;
   }

   // Append the 'save number' to the filename, and prepend zeros to the save number
   strFilePathName.append("_");
   stringstream ststrTmp;
   if (m_bGISSaveDigitsSequential)
   {
      // Save number is sequential
      ststrTmp << FillToWidth('0', m_nGISMaxSaveDigits) << m_nGISSave;
   }
   else
   {
      // Save number is iteration
      ststrTmp << FillToWidth('0', m_nGISMaxSaveDigits) << m_ulIter;
   }
   strFilePathName.append(ststrTmp.str());
   strstrFileName << ststrTmp.str();

   // Make a copy of the filename without any extension
   // string strFilePathNameNoExt = strFilePathName;

   // If desired, append an extension
   if (! m_strOGRVectorOutputExtension.empty())
      strFilePathName.append(m_strOGRVectorOutputExtension);

   // Set up the vector driver
   GDALDriver* pGDALDriver = GetGDALDriverManager()->GetDriverByName(m_strVectorGISOutFormat.c_str());
   if (pGDALDriver == NULL)
   {
      cerr << ERR << "vector GIS output driver " << m_strVectorGISOutFormat << CPLGetLastErrorMsg() << endl;
      return false;
   }

   // Now create the dataset
   GDALDataset* pGDALDataSet = pGDALDriver->Create(strFilePathName.c_str(), 0, 0, 0, GDT_Unknown, m_papszGDALVectorOptions);
   if (pGDALDataSet == NULL)
   {
      cerr << ERR << "cannot create " << m_strVectorGISOutFormat << " named " << strFilePathName << "\n"
           << CPLGetLastErrorMsg() << endl;
      return false;
   }

   // Create a spatial reference object
   OGRSpatialReference OGRSpatialRef;

   // And tell it about the co-ordinate system used by the basement raster layer
   if (m_strGDALBasementDEMProjection.empty())
   {
      OGRSpatialRef.importFromWkt(m_strGDALBasementDEMProjection.c_str());
   }
   else
   {
      OGRSpatialRef.importFromEPSG(25830); // TODO 035 Change for any EPSG
   }

   OGRwkbGeometryType eGType = wkbUnknown;
   string strType = "unknown";

   // Now create the output layer
   // // OGRLayer* pOGRLayer = pGDALDataSet->CreateLayer(strFilePathNameNoExt.c_str(), &OGRSpatialRef, eGType, m_papszGDALVectorOptions);
   // if (EQUAL(m_strVectorGISOutFormat.c_str(), "geojson"))
   // {
   //    CPLSetConfigOption("GDAL_VALIDATE_CREATION_OPTIONS", "NO");
   //    m_papszGDALVectorOptions = CSLSetNameValue(m_papszGDALVectorOptions, "COORDINATE_PRECISION", "2");
   // }

   OGRLayer* pOGRLayer = pGDALDataSet->CreateLayer(strstrFileName.str().c_str(), &OGRSpatialRef, eGType, m_papszGDALVectorOptions);

   if (pOGRLayer == NULL)
   {
      cerr << ERR << "cannot create '" << strType << "' layer in " << strFilePathName << "\n"
           << CPLGetLastErrorMsg() << endl;
      return false;
   }

   // Need to switch off OGR error messages here, since real (floating point) fields in ESRI shapefiles are treated as width 24 with 15 decimal places of precision (unless an explicit width is given). If fields exceed this, then a "not successfully written. Possibly due to too larger number with respect to field width" error message is shown. Get this to fail silently
   CPLPushErrorHandler(CPLQuietErrorHandler); 
   
   switch (nDataItem)
   {
      case (VECTOR_PLOT_COAST):
      {
         eGType = wkbLineString;
         strType = "line";

         // The layer has been created, so create an integer-numbered value (the number of the coast object) for the multi-line
         string strFieldValue1 = "Coast";
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now do features
         OGRLineString OGRls;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            // Create a feature object, one per coast
            OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

            // Set the feature's attribute (the coast number)
            pOGRFeature->SetField(strFieldValue1.c_str(), i);

            // Now attach a geometry to the feature object
            for (int j = 0; j < m_VCoast[i].pLGetCoastlineExtCRS()->nGetSize(); j++)
               //  In external CRS
               OGRls.addPoint(m_VCoast[i].pPtGetCoastlinePointExtCRS(j)->dGetX(), m_VCoast[i].pPtGetCoastlinePointExtCRS(j)->dGetY());

            pOGRFeature->SetGeometry(&OGRls);

            // Create the feature in the output layer
            if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
            {
               cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << " for coast " << i << " in " << strFilePathName << "\n"
                  << CPLGetLastErrorMsg() << endl;
               return false;
            }

            // Tidy up: empty the line string and get rid of the feature object
            OGRls.empty();
            OGRFeature::DestroyFeature(pOGRFeature);
         }

         break;
      }
   
      case (VECTOR_PLOT_FLOOD_LINE):
      {
         eGType = wkbLineString;
         strType = "line";

         // The layer has been created, so create an integer-numbered value (the number of the coast object) for the multi-line
         string strFieldValue1 = "NMR";
         string strFieldValue2 = "tiempo";
         string strFieldValue3 = "surge_mm";
         // string strFieldValue4 = "eta-surge(mm)";
         string strFieldValue4 = "runup_mm";

         // Create a feature with general properties
         // OGRLineString OGRls;
         // OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

         // Testing coordinate system
         // if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         // {
         //    cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
         //         << CPLGetLastErrorMsg() << endl;
         //    return false;
         // }
         // if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         // {
         //    cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n"
         //         << CPLGetLastErrorMsg() << endl;
         //    return false;
         // }

         // pOGRFeature->SetGeometry(&OGRls);
         // Create the feature in the output layer
         // if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
         // {
         //    cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << " in " << strFilePathName << "\n"
         //         << CPLGetLastErrorMsg() << endl;
         //    return false;
         // }
         // OGRFeature::DestroyFeature(pOGRFeature);
         // Create a feature object, one per coast

         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTReal);
         OGRFieldDefn OGRField2(strFieldValue2.c_str(), OFTReal);
         OGRFieldDefn OGRField3(strFieldValue3.c_str(), OFTInteger64);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField3) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 3 '" << strFieldValue3 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now do features
         OGRLineString OGR2ls;
         // for (int i = 0; i < static_cast<int>(m_VFloodWaveSetupSurge.size()); i++)
         // {
         //    OGRFeature* pOGR2Feature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());
         //    pOGR2Feature->SetField(strFieldValue1.c_str(), m_dThisIterSWL);
         //    pOGR2Feature->SetField(strFieldValue2.c_str(), m_ulIter);
         //    int setup_level = int(m_dThisIterDiffWaveSetupWaterLevel * 1000);
         //    pOGR2Feature->SetField(strFieldValue3.c_str(), setup_level);
         //    // Set the feature's attribute (the coast number)
         //    // Now attach a geometry to the feature object
         //    for (int j = 0; j < m_VFloodWaveSetup[i].pLGetCoastlineExtCRS()->nGetSize(); j++)
         //    {
         //       //  In external CRS
         //       // Add SWL + wave setup line
         //       OGR2ls.addPoint(m_VFloodWaveSetup[i].pPtGetCoastlinePointExtCRS(j)->dGetX(), m_VFloodWaveSetup[i].pPtGetCoastlinePointExtCRS(j)->dGetY());
         //    }

         //    pOGR2Feature->SetGeometry(&OGR2ls);
         //    OGR2ls.empty();

         //    // Create the feature in the output layer
         //    if (pOGRLayer->CreateFeature(pOGR2Feature) != OGRERR_NONE)
         //    {
         //       cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << " for coast " << i << " in " << strFilePathName << "\n"
         //            << CPLGetLastErrorMsg() << endl;
         //       return false;
         //    }

         //    // Tidy up: empty the line string and get rid of the feature object
         //    // OGRls.empty();
         //    OGRFeature::DestroyFeature(pOGR2Feature);
         // }
         if (m_bFloodSWLSetupSurgeLine)
         {
            // Create a feature object, one per coast
            OGRFeature* pOGR3Feature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());
            OGRFieldDefn OGRField6(strFieldValue1.c_str(), OFTReal);
            OGRFieldDefn OGRField7(strFieldValue2.c_str(), OFTReal);
            OGRFieldDefn OGRField4(strFieldValue3.c_str(), OFTInteger64);
            if (pOGRLayer->CreateField(&OGRField6) != OGRERR_NONE)
            {
               cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
                  << CPLGetLastErrorMsg() << endl;
               return false;
            }
            if (pOGRLayer->CreateField(&OGRField7) != OGRERR_NONE)
            {
               cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n"
                  << CPLGetLastErrorMsg() << endl;
               return false;
            }
            if (pOGRLayer->CreateField(&OGRField4) != OGRERR_NONE)
            {
               cerr << ERR << "cannot create " << strType << " attribute field 4 '" << strFieldValue4 << "' in " << strFilePathName << "\n"
                  << CPLGetLastErrorMsg() << endl;
               return false;
            }

            // OK, now do features
            OGRLineString OGR3ls;
            for (int i = 0; i < static_cast<int>(m_VFloodWaveSetupSurge.size()); i++)
            {

               pOGR3Feature->SetField(strFieldValue1.c_str(), m_dThisIterSWL);
               pOGR3Feature->SetField(strFieldValue2.c_str(), static_cast<double>(m_bGISSaveDigitsSequential ? m_nGISSave : m_ulIter));
               int surge_level = int(m_dThisIterDiffWaveSetupSurgeWaterLevel * 1000);
               pOGR3Feature->SetField(strFieldValue3.c_str(), surge_level);
               // Set the feature's attribute (the coast number)
               // Now attach a geometry to the feature object
               for (int j = 0; j < m_VFloodWaveSetupSurge[i].pLGetCoastlineExtCRS()->nGetSize(); j++)
               {
                  //  In external CRS
                  // Add SWL + wave setup + storm surge line
                  OGR3ls.addPoint(m_VFloodWaveSetupSurge[i].pPtGetCoastlinePointExtCRS(j)->dGetX(), m_VFloodWaveSetupSurge[i].pPtGetCoastlinePointExtCRS(j)->dGetY());
               }

               pOGR3Feature->SetGeometry(&OGR3ls);
               OGR3ls.empty();

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGR3Feature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << " for coast " << i << " in " << strFilePathName << "\n"
                     << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Tidy up: empty the line string and get rid of the feature object
               // OGRls.empty();
               OGRFeature::DestroyFeature(pOGR3Feature);
            }
         }
         if (m_bFloodSWLSetupSurgeRunupLine)
         {
            // Create a feature object, one per coast
            OGRFeature* pOGR4Feature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());
            OGRFieldDefn OGRField8(strFieldValue1.c_str(), OFTReal);
            OGRFieldDefn OGRField9(strFieldValue2.c_str(), OFTReal);
            OGRFieldDefn OGRField5(strFieldValue4.c_str(), OFTInteger64);
            if (pOGRLayer->CreateField(&OGRField8) != OGRERR_NONE)
            {
               cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
                  << CPLGetLastErrorMsg() << endl;
               return false;
            }
            if (pOGRLayer->CreateField(&OGRField9) != OGRERR_NONE)
            {
               cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n"
                  << CPLGetLastErrorMsg() << endl;
               return false;
            }
            if (pOGRLayer->CreateField(&OGRField5) != OGRERR_NONE)
            {
               cerr << ERR << "cannot create " << strType << " attribute field 5 '" << strFieldValue4 << "' in " << strFilePathName << "\n"
                  << CPLGetLastErrorMsg() << endl;
               return false;
            }
            // OK, now do features
            OGRLineString OGR4ls;
            for (int i = 0; i < static_cast<int>(m_VFloodWaveSetupSurgeRunup.size()); i++)
            {
               pOGR4Feature->SetField(strFieldValue1.c_str(), m_dThisIterSWL);
               pOGR4Feature->SetField(strFieldValue2.c_str(), static_cast<double>(m_bGISSaveDigitsSequential ? m_nGISSave : m_ulIter));
               int runup_level = int(m_dThisIterDiffWaveSetupSurgeRunupWaterLevel * 1000);
               pOGR4Feature->SetField(strFieldValue4.c_str(), runup_level);

               // Set the feature's attribute (the coast number)
               // Now attach a geometry to the feature object
               for (int j = 0; j < m_VFloodWaveSetupSurgeRunup[i].pLGetCoastlineExtCRS()->nGetSize(); j++)
               {
                  //  In external CRS
                  // Add SWL + wave setup + surge + runup line
                  OGR4ls.addPoint(m_VFloodWaveSetupSurgeRunup[i].pPtGetCoastlinePointExtCRS(j)->dGetX(), m_VFloodWaveSetupSurgeRunup[i].pPtGetCoastlinePointExtCRS(j)->dGetY());
               }

               pOGR4Feature->SetGeometry(&OGR4ls);
               OGR4ls.empty();

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGR4Feature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << " for coast " << i << " in " << strFilePathName << "\n"
                     << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Tidy up: empty the line string and get rid of the feature object
               // OGRls.empty();
               OGRFeature::DestroyFeature(pOGR4Feature);
            }
         }
         // OGRls.empty();

         break;
      }

      case (VECTOR_PLOT_NORMALS):
      case (VECTOR_PLOT_INVALID_NORMALS):
      {
         eGType = wkbLineString;
         strType = "line";

         // The layer has been created, so create an integer-numbered value (the number of the normal) associated with the line
         string strFieldValue1 = "Normal";
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // Also create other integer-numbered values for the category codes of the coastline-normalprofile
         string
            strFieldValue2 = "StartCoast",
            strFieldValue3 = "EndCoast",
            strFieldValue4 = "HitLand",
            strFieldValue5 = "HitCoast",
            strFieldValue6 = "HitNormal";
         OGRFieldDefn
            OGRField2(strFieldValue2.c_str(), OFTInteger),
            OGRField3(strFieldValue3.c_str(), OFTInteger),
            OGRField4(strFieldValue4.c_str(), OFTInteger),
            OGRField5(strFieldValue5.c_str(), OFTInteger),
            OGRField6(strFieldValue6.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField3) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 3 '" << strFieldValue3 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField4) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 4 '" << strFieldValue4 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField5) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 5 '" << strFieldValue5 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }
         if (pOGRLayer->CreateField(&OGRField6) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 6 '" << strFieldValue6 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now create features
         OGRLineString OGRls;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            for (int j = 0; j < m_VCoast[i].nGetNumProfiles(); j++)
            {
               CGeomProfile* pProfile = m_VCoast[i].pGetProfile(j);

               if (((nDataItem == VECTOR_PLOT_NORMALS) && (pProfile->bOKIncStartAndEndOfCoast())) || ((nDataItem == VECTOR_PLOT_INVALID_NORMALS) && (!pProfile->bOKIncStartAndEndOfCoast())))
               {
                  // Create a feature object, one per profile
                  OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

                  // Set the feature's attributes
                  pOGRFeature->SetField(strFieldValue1.c_str(), j);
                  pOGRFeature->SetField(strFieldValue2.c_str(), 0);
                  pOGRFeature->SetField(strFieldValue3.c_str(), 0);
                  pOGRFeature->SetField(strFieldValue4.c_str(), 0);
                  pOGRFeature->SetField(strFieldValue5.c_str(), 0);
                  pOGRFeature->SetField(strFieldValue6.c_str(), 0);
                  if (pProfile->bStartOfCoast())
                     pOGRFeature->SetField(strFieldValue2.c_str(), 1);
                  if (pProfile->bEndOfCoast())
                     pOGRFeature->SetField(strFieldValue3.c_str(), 1);
                  if (pProfile->bHitLand())
                     pOGRFeature->SetField(strFieldValue4.c_str(), 1);
                  if (pProfile->bHitCoast())
                     pOGRFeature->SetField(strFieldValue5.c_str(), 1);
                  if (pProfile->bHitAnotherProfile())
                     pOGRFeature->SetField(strFieldValue6.c_str(), 1);

                  // Now attach a geometry to the feature object
                  for (int k = 0; k < pProfile->nGetProfileSize(); k++)
                     OGRls.addPoint(pProfile->pPtGetPointInProfile(k)->dGetX(), pProfile->pPtGetPointInProfile(k)->dGetY());

                  pOGRFeature->SetGeometry(&OGRls);
                  OGRls.empty();

                  // Create the feature in the output layer
                  if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
                  {
                     cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << " for coast " << i << " and profile " << j << " in " << strFilePathName << "\n"
                        << CPLGetLastErrorMsg() << endl;
                     return false;
                  }

                  // Tidy up: get rid of the feature object
                  OGRFeature::DestroyFeature(pOGRFeature);
               }
            }
         }

         break;
      }

      case (VECTOR_PLOT_COAST_CURVATURE):
      case (VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE):
      case (VECTOR_PLOT_MEAN_WAVE_ENERGY):
      case (VECTOR_PLOT_BREAKING_WAVE_HEIGHT):
      case (VECTOR_PLOT_POLYGON_NODES):
      case (VECTOR_PLOT_WAVE_SETUP):
      case (VECTOR_PLOT_STORM_SURGE):
      case (VECTOR_PLOT_RUN_UP):
      case (VECTOR_PLOT_CLIFF_NOTCH_SIZE):
      {
         eGType = wkbPoint;
         strType = "point";

         // The layer has been created, so create a real-numbered value associated with each point
         string strFieldValue1;
         if (nDataItem == VECTOR_PLOT_COAST_CURVATURE)
            strFieldValue1 = "Curve";
         else if (nDataItem == VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE)
            strFieldValue1 = "SC_Energy";
         else if (nDataItem == VECTOR_PLOT_MEAN_WAVE_ENERGY)
            strFieldValue1 = "MeanEnergy";
         else if (nDataItem == VECTOR_PLOT_BREAKING_WAVE_HEIGHT)
            strFieldValue1 = "Height";
         else if (nDataItem == VECTOR_PLOT_POLYGON_NODES)
            strFieldValue1 = "Node";
         else if (nDataItem == VECTOR_PLOT_CLIFF_NOTCH_SIZE)
            strFieldValue1 = "Notch";
         else if (nDataItem == VECTOR_PLOT_WAVE_SETUP)
            strFieldValue1 = "Wavesetup";
         else if (nDataItem == VECTOR_PLOT_STORM_SURGE)
            strFieldValue1 = "Stormsurge";
         else if (nDataItem == VECTOR_PLOT_RUN_UP)
            strFieldValue1 = "Runup";

         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now create features
         OGRLineString OGRls;
         OGRMultiLineString OGRmls;
         OGRPoint OGRPt;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            for (int j = 0; j < m_VCoast[i].pLGetCoastlineExtCRS()->nGetSize(); j++)
            {
               // Create a feature object, one per coastline point
               OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

               // Set the feature's geometry (in external CRS)
               OGRPt.setX(m_VCoast[i].pPtGetCoastlinePointExtCRS(j)->dGetX());
               OGRPt.setY(m_VCoast[i].pPtGetCoastlinePointExtCRS(j)->dGetY());
               pOGRFeature->SetGeometry(&OGRPt);

               if (nDataItem == VECTOR_PLOT_COAST_CURVATURE)
               {
                  double dCurvature = m_VCoast[i].dGetDetailedCurvature(j);
                  if (bFPIsEqual(dCurvature, DBL_NODATA, TOLERANCE))
                     continue;

                  // Set the feature's attribute
                  pOGRFeature->SetField(strFieldValue1.c_str(), dCurvature);
               }
               else if (nDataItem == VECTOR_PLOT_WAVE_ENERGY_SINCE_COLLAPSE)
               {
                  // Set the feature's attribute
                  if (m_VCoast[i].pGetCoastLandform(j) == NULL)
                     pOGRFeature->SetField(strFieldValue1.c_str(), DBL_NODATA);
                  else
                     pOGRFeature->SetField(strFieldValue1.c_str(), m_VCoast[i].pGetCoastLandform(j)->dGetTotAccumWaveEnergy());
               }
               else if (nDataItem == VECTOR_PLOT_MEAN_WAVE_ENERGY)
               {
                  // Set the feature's attribute
                  if (m_VCoast[i].pGetCoastLandform(j) == NULL)
                     pOGRFeature->SetField(strFieldValue1.c_str(), DBL_NODATA);
                  else
                  {
                     double dEnergy = m_VCoast[i].pGetCoastLandform(j)->dGetTotAccumWaveEnergy();
                     dEnergy *= 24;
                     dEnergy /= m_dSimElapsed; // Is in energy units per day

                     pOGRFeature->SetField(strFieldValue1.c_str(), dEnergy);
                  }
               }
               else if (nDataItem == VECTOR_PLOT_BREAKING_WAVE_HEIGHT)
               {
                  // Set the feature's attribute
                  double dHeight = m_VCoast[i].dGetBreakingWaveHeight(j);
                  pOGRFeature->SetField(strFieldValue1.c_str(), dHeight);
               }
               else if (nDataItem == VECTOR_PLOT_WAVE_SETUP)
               {
                  // Set the feature's attribute
                  double dWaveSetupSurge = m_VCoast[i].dGetWaveSetupSurge(j);
                  pOGRFeature->SetField(strFieldValue1.c_str(), dWaveSetupSurge);
               }
               // else if (nDataItem == VECTOR_PLOT_STORM_SURGE)
               // {
               //    // Set the feature's attribute
               //    double dStormSurge = m_VCoast[i].dGetStormSurge(j);
               //    pOGRFeature->SetField(strFieldValue1.c_str(), dStormSurge);
               // }
               else if (nDataItem == VECTOR_PLOT_RUN_UP)
               {
                  // Set the feature's attribute
                  double dRunUp = m_VCoast[i].dGetRunUp(j);
                  pOGRFeature->SetField(strFieldValue1.c_str(), dRunUp);
               }
               else if (nDataItem == VECTOR_PLOT_POLYGON_NODES)
               {
                  int nNode = m_VCoast[i].nGetPolygonNode(j);
                  if (nNode == INT_NODATA)
                     continue;

                  // Set the feature's attribute
                  pOGRFeature->SetField(strFieldValue1.c_str(), nNode);
               }
               else if (nDataItem == VECTOR_PLOT_CLIFF_NOTCH_SIZE)
               {
                  CACoastLandform* pCoastLandform = m_VCoast[i].pGetCoastLandform(j);
                  if (pCoastLandform == NULL)
                     pOGRFeature->SetField(strFieldValue1.c_str(), DBL_NODATA);
                  else
                  {
                     int nCategory = pCoastLandform->nGetLandFormCategory();
                     double dNotchDepth = 0.0;

                     if (nCategory == LF_CAT_CLIFF)
                     {
                        CRWCliff const* pCliff = reinterpret_cast<CRWCliff*>(pCoastLandform);

                        // Get attribute values from the cliff object
                        dNotchDepth = pCliff->dGetNotchDepth();
                     }

                     // Set the feature's attribute
                     pOGRFeature->SetField(strFieldValue1.c_str(), dNotchDepth);
                  }
               }

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create " << strType << " feature " << strPlotTitle << " for coast " << i << " point " << j << " in " << strFilePathName << "\n"
                     << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Get rid of the feature object
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }

         break;
      }

      case (VECTOR_PLOT_WAVE_ANGLE_AND_HEIGHT):
      {
         eGType = wkbPoint;
         strType = "point";

         // The layer has been created, so create real-numbered values associated with each point
         string
            strFieldValue1 = "Angle",
            strFieldValue2 = "Height";

         // Create the first field
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // Create the second field
         OGRFieldDefn OGRField2(strFieldValue2.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now create features
         OGRLineString OGRls;
         OGRMultiLineString OGRmls;
         OGRPoint OGRPt;

         for (int nX = 0; nX < m_nXGridMax; nX++)
         {
            for (int nY = 0; nY < m_nYGridMax; nY++)
            {
               // Only output a value if the cell is a sea cell which is not in the active zone (wave height and angle values are meaningless if in the active zone)
               if ((m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()) && (! m_pRasterGrid->m_Cell[nX][nY].bIsInActiveZone()))
               {
                  // Create a feature object, one per sea cell
                  OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

                  // Set the feature's geometry (in external CRS)
                  OGRPt.setX(dGridCentroidXToExtCRSX(nX));
                  OGRPt.setY(dGridCentroidYToExtCRSY(nY));
                  pOGRFeature->SetGeometry(&OGRPt);

                  double
                     dOrientation = m_pRasterGrid->m_Cell[nX][nY].dGetWaveAngle(),
                     dHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight();

                  if (bFPIsEqual(dHeight, DBL_NODATA, TOLERANCE) || bFPIsEqual(dOrientation, DBL_NODATA, TOLERANCE))
                     continue;

                  // Set the feature's attributes
                  pOGRFeature->SetField(strFieldValue1.c_str(), dOrientation);
                  pOGRFeature->SetField(strFieldValue2.c_str(), dHeight);

                  // Create the feature in the output layer
                  if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
                  {
                     cerr << ERR << "cannot create " << strType << " feature " << strPlotTitle << " for cell [" << nX << "][" << nY << "] in " << strFilePathName << "\n"
                        << CPLGetLastErrorMsg() << endl;
                     return false;
                  }

                  // Get rid of the feature object
                  OGRFeature::DestroyFeature(pOGRFeature);
               }
            }
         }

         break;
      }

      case (VECTOR_PLOT_AVG_WAVE_ANGLE_AND_HEIGHT):
      {
         eGType = wkbPoint;
         strType = "point";

         // The layer has been created, so create real-numbered values associated with each point
         string
            strFieldValue1 = "Angle",
            strFieldValue2 = "Height";

         // Create the first field
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // Create the second field
         OGRFieldDefn OGRField2(strFieldValue2.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now create features
         OGRLineString OGRls;
         OGRMultiLineString OGRmls;
         OGRPoint OGRPt;

         for (int nX = 0; nX < m_nXGridMax; nX++)
         {
            for (int nY = 0; nY < m_nYGridMax; nY++)
            {
               // Create a feature object, one per sea cell
               OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

               // Set the feature's geometry (in external CRS)
               OGRPt.setX(dGridCentroidXToExtCRSX(nX));
               OGRPt.setY(dGridCentroidYToExtCRSY(nY));
               pOGRFeature->SetGeometry(&OGRPt);

               double
                  dOrientation = m_pRasterGrid->m_Cell[nX][nY].dGetTotWaveAngle() / static_cast<double>(m_ulIter),
                  dHeight = m_pRasterGrid->m_Cell[nX][nY].dGetWaveHeight() / static_cast<double>(m_ulIter);

               if (bFPIsEqual(dHeight, DBL_NODATA, TOLERANCE) || bFPIsEqual(dOrientation, DBL_NODATA, TOLERANCE))
                  continue;

               // Set the feature's attributes
               pOGRFeature->SetField(strFieldValue1.c_str(), dOrientation);
               pOGRFeature->SetField(strFieldValue2.c_str(), dHeight);

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create " << strType << " feature " << strPlotTitle << " for cell [" << nX << "][" << nY << "] in " << strFilePathName << "\n"
                     << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Get rid of the feature object
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }

         break;
      }

      case (VECTOR_PLOT_POLYGON_BOUNDARY):
      {
         eGType = wkbPolygon;
         strType = "polygon";

         // The layer has been created, so create two integer-numbered values (the number of the polygon object, and the number of the coast point which is the polygon's node) for the polygon
         string
            strFieldValue1 = "Polygon",
            strFieldValue2 = "CoastNode",
            strFieldValue3 = "TotSedChng",
            strFieldValue4 = "FinSedChng",
            strFieldValue5 = "SndSedChng",
            strFieldValue6 = "CrsSedChng";

         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField2(strFieldValue2.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField3(strFieldValue3.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField3) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 3 '" << strFieldValue3 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField4(strFieldValue4.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField4) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 4 '" << strFieldValue4 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField5(strFieldValue5.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField5) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 5 '" << strFieldValue5 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         OGRFieldDefn OGRField6(strFieldValue6.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField6) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 6 '" << strFieldValue6 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now do features
         OGRLineString OGRls;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            for (int j = 0; j < m_VCoast[i].nGetNumPolygons(); j++)
            {
               // Create a feature object, one per polygon
               OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

               CGeomCoastPolygon* pPolygon = m_VCoast[i].pGetPolygon(j);

               // Set the feature's attributes
               pOGRFeature->SetField(strFieldValue1.c_str(), j);
               pOGRFeature->SetField(strFieldValue2.c_str(), pPolygon->nGetNodeCoastPoint());
               pOGRFeature->SetField(strFieldValue3.c_str(), pPolygon->dGetBeachDepositionAndSuspensionAllUncons());
               pOGRFeature->SetField(strFieldValue4.c_str(), pPolygon->dGetSuspensionUnconsFine());
               pOGRFeature->SetField(strFieldValue5.c_str(), pPolygon->dGetBeachDepositionUnconsSand());
               pOGRFeature->SetField(strFieldValue6.c_str(), pPolygon->dGetBeachDepositionUnconsCoarse());

               // Now attach a geometry to the feature object
               for (int n = 0; n < pPolygon->nGetBoundarySize(); n++)
                  //  In external CRS
                  OGRls.addPoint(pPolygon->pPtGetBoundaryPoint(n)->dGetX(), pPolygon->pPtGetBoundaryPoint(n)->dGetY());

               pOGRFeature->SetGeometry(&OGRls);

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create " << strType << " feature " << strPlotTitle << " for coast " << i << " polygon " << j << " in " << strFilePathName << "\n"
                     << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Tidy up: empty the line string and get rid of the feature object
               OGRls.empty();
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }

         break;
      }

      case (VECTOR_PLOT_SHADOW_BOUNDARY):
      {
         eGType = wkbLineString;
         strType = "line";

         // Create an integer-numbered value (the number of the shadow boundary line object) for the multi-line
         string strFieldValue1 = "ShadowLine";
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now do features
         OGRLineString OGRls;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            for (int j = 0; j < m_VCoast[i].nGetNumShadowBoundaries(); j++)
            {
               // Create a feature object, one per coast
               OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

               // Set the feature's attribute (the shadow boundary line number)
               pOGRFeature->SetField(strFieldValue1.c_str(), j);

               // Now attach a geometry to the feature object
               CGeomLine LShadow = *m_VCoast[i].pGetShadowBoundary(j);
               for (int nn = 0; nn < LShadow.nGetSize(); nn++)
                  OGRls.addPoint(LShadow.dGetXAt(nn), LShadow.dGetYAt(nn));

               pOGRFeature->SetGeometry(&OGRls);

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << j << " for coast " << i << " in " << strFilePathName << "\n"
                     << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Tidy up: empty the line string and get rid of the feature object
               OGRls.empty();
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }

         break;
      }

      case (VECTOR_PLOT_DOWNDRIFT_BOUNDARY):
      {
         eGType = wkbLineString;
         strType = "line";

         // Create an integer-numbered value (the number of the downdrift boundary line object) for the multi-line
         string strFieldValue1 = "DdriftLine";
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTInteger);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now do features
         OGRLineString OGRls;

         for (int i = 0; i < static_cast<int>(m_VCoast.size()); i++)
         {
            for (int j = 0; j < m_VCoast[i].nGetNumShadowDowndriftBoundaries(); j++)
            {
               // Create a feature object, one per coast
               OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

               // Set the feature's attribute (the downdrift boundary line number)
               pOGRFeature->SetField(strFieldValue1.c_str(), j);

               // Now attach a geometry to the feature object
               CGeomLine LDowndrift = *m_VCoast[i].pGetShadowDowndriftBoundary(j);
               for (int nn = 0; nn < LDowndrift.nGetSize(); nn++)
                  OGRls.addPoint(LDowndrift.dGetXAt(nn), LDowndrift.dGetYAt(nn));

               pOGRFeature->SetGeometry(&OGRls);

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create  " << strType << " feature " << strPlotTitle << j << " for coast " << i << " in " << strFilePathName << "\n"
                     << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Tidy up: empty the line string and get rid of the feature object
               OGRls.empty();
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }

         break;
      }

      case (VECTOR_PLOT_DEEP_WATER_WAVE_ANGLE_AND_HEIGHT):
      {
         eGType = wkbPoint;
         strType = "point";

         // The layer has been created, so create real-numbered values associated with each point
         string
            strFieldValue1 = "Angle",
            strFieldValue2 = "Height";

         // Create the first field
         OGRFieldDefn OGRField1(strFieldValue1.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField1) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 1 '" << strFieldValue1 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // Create the second field
         OGRFieldDefn OGRField2(strFieldValue2.c_str(), OFTReal);
         if (pOGRLayer->CreateField(&OGRField2) != OGRERR_NONE)
         {
            cerr << ERR << "cannot create " << strType << " attribute field 2 '" << strFieldValue2 << "' in " << strFilePathName << "\n"
               << CPLGetLastErrorMsg() << endl;
            return false;
         }

         // OK, now create features
         OGRLineString OGRls;
         OGRMultiLineString OGRmls;
         OGRPoint OGRPt;

         for (int nX = 0; nX < m_nXGridMax; nX++)
         {
            for (int nY = 0; nY < m_nYGridMax; nY++)
            {
               // Create a feature object, one per cell (does this whether the sea is a sea cell or a land cell)
               OGRFeature* pOGRFeature = OGRFeature::CreateFeature(pOGRLayer->GetLayerDefn());

               // Set the feature's geometry (in external CRS)
               OGRPt.setX(dGridCentroidXToExtCRSX(nX));
               OGRPt.setY(dGridCentroidYToExtCRSY(nY));
               pOGRFeature->SetGeometry(&OGRPt);

               double
                  dOrientation = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveAngle(),
                  dHeight = m_pRasterGrid->m_Cell[nX][nY].dGetCellDeepWaterWaveHeight();

               if (bFPIsEqual(dHeight, DBL_NODATA, TOLERANCE) || bFPIsEqual(dOrientation, DBL_NODATA, TOLERANCE) || (! m_pRasterGrid->m_Cell[nX][nY].bIsInContiguousSea()))
                  continue;

               // Set the feature's attributes
               pOGRFeature->SetField(strFieldValue1.c_str(), dOrientation);
               pOGRFeature->SetField(strFieldValue2.c_str(), dHeight);

               // Create the feature in the output layer
               if (pOGRLayer->CreateFeature(pOGRFeature) != OGRERR_NONE)
               {
                  cerr << ERR << "cannot create " << strType << " feature " << strPlotTitle << " for cell [" << nX << "][" << nY << "] in " << strFilePathName << "\n"
                     << CPLGetLastErrorMsg() << endl;
                  return false;
               }

               // Get rid of the feature object
               OGRFeature::DestroyFeature(pOGRFeature);
            }
         }
         break;
      }
   }

   CPLPopErrorHandler();
   
   // Get rid of the dataset object
   GDALClose(pGDALDataSet);

   return true;
}
