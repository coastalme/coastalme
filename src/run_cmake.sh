#!/bin/sh

# Change this to change build type
buildtype=Debug
#buildtype=Release
#buildtype=Prerelease             # DO NOT run using debugger e.g. gdb
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
echo ""
rm -f ./lib/*
cd cshore

if [[ $OSTYPE == 'darwin'* ]]; then
   cshorelibrary=STATIC
   ./make_cshore_lib.sh $cshorelibrary $buildtype $cshoreinout
else
   ./make_cshore_lib.sh
fi
cd ..
# Note: The cshore Makefile now correctly names libraries for MacOS automatically
echo ""

# Now run CMake for CoastalME
echo ""
echo "================================================================="
echo ""

# On Mac, switch to gcc15 to ensure openMP compatability
if [ "$OSTYPE" = "darwin"* ]; then
   export CC=gcc-15
   export CXX=g++-15
   CMAKE_COMPILER_ARGS="-DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15"
else
   echo "Using default compiler, OpenMP may not be available"
   echo "For OpenMP support, install GCC or install libomp for Clang"
   CMAKE_COMPILER_ARGS=""
fi
echo ""

rm -f CMakeCache.txt
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout $CMAKE_COMPILER_ARGS .
#cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout -DCMAKE_VERBOSE_MAKEFILE=ON $CMAKE_COMPILER_ARGS .
#cmake -DCMAKE_BUILD_TYPE=$buildtype -DCSHORE_LIBRARY=$cshorelibrary -DCSHORE_INOUT=$cshoreinout $CMAKE_COMPILER_ARGS . -G"CodeBlocks - Unix Makefiles"

echo ""
echo "================================================================="
echo ""

echo "Finished CMake ($buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
echo ""
echo "================================================================="

