#!/bin/sh

# Change this to change build type
buildtype=Debug
#buildtype=Release
#buildtype=RelWithDebInfo        # Not yet implemented in CMakeLists.txt
#buildtype=MinSizeRel            # Not yet implemented in CMakeLists.txt
#buildtype=gcov
#buildtype=Callgrind

# Change this to select the CShore library type
#cshorelibrary=STATIC
cshorelibrary=SHARED

# Change this to select CShore input/output method
#cshoreinout=FILE
cshoreinout=ARG
#cshoreinout=BOTH

# Always build CShore
echo "Building all versions of the CShore library"
echo ""
cd cshore
./make_cshore_lib.sh
cd ..

# Now run CMake for CoastalME
echo ""
echo "================================================================="
echo ""
echo "CoastalME: starting CMake for Linux (using gcc, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
echo ""

rm -f CMakeCache.txt
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout .
#cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout -DCMAKE_VERBOSE_MAKEFILE=ON .
#cmake -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout . -G"CodeBlocks - Unix Makefiles"

echo ""
echo "================================================================="
echo ""

echo "Finished CMake for Linux (using gcc, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
echo ""
echo "================================================================="

