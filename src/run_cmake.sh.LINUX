#!/bin/sh

# Change this to change build type
buildtype=Debug
#buildtype=Release
#buildtype=RelWithDebInfo
#buildtype=MinSizeRel
#buildtype=gcov
#buildtype=Callgrind

echo "CoastalME: starting CMake for Linux (using gcc, $buildtype build)"

rm -f CMakeCache.txt
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype .
