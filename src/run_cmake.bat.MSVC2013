@echo off

rem Change this to change build type
set buildtype=Debug
rem set buildtype=Release
rem set buildtype=RelWithDebInfo
rem set buildtype=MinSizeRel

echo CoastalME: starting CMake for 64-bit Windows (using MS Visual Studio 2013, %buildtype% build)

del CMakeCache.txt
rem cmake -G "Visual Studio 12 2013 Win64" -DCMAKE_BUILD_TYPE=%buildtype% -DGDAL_LIBRARY=$ENV{GDAL_DIR} -DGDAL_INCLUDE_DIR=$ENV{GDAL_DIR} .
cmake -G "Visual Studio 12 2013 Win64" -DCMAKE_BUILD_TYPE=%buildtype% .
