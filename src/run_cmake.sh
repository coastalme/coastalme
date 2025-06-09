#!/bin/sh

# Change this to change build type
#buildtype=Debug
buildtype=Release
#buildtype=Prerelease
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
echo "CoastalME: starting CMake for macOS (using gcc-15 with OpenMP, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
echo ""

rm -f CMakeCache.txt
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout 
# CC=/opt/homebrew/bin/gcc-15 CXX=/opt/homebrew/bin/g++-15 cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 .
#cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout -DCMAKE_VERBOSE_MAKEFILE=ON .
#cmake -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout . -G"CodeBlocks - Unix Makefiles"

echo ""
echo "================================================================="
echo ""

echo "Finished CMake for macOS (using gcc-15 with OpenMP, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
echo ""
echo "================================================================="

