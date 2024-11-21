#!/bin/sh

# Change this to change build type
buildtype=Debug
#buildtype=Release
#buildtype=RelWithDebInfo
#buildtype=MinSizeRel

echo "CoastalME: starting CMake for Windows (using Cygwin, $buildtype build)"

rm -f CMakeCache.txt
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype .
