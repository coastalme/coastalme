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

# Check for available compilers and choose the best option for OpenMP
if [ -f "/opt/homebrew/bin/gcc-15" ] && [ -f "/opt/homebrew/bin/g++-15" ]; then
    echo "CoastalME: starting CMake for macOS (using gcc-15 with OpenMP, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
    export CC=/opt/homebrew/bin/gcc-15
    export CXX=/opt/homebrew/bin/g++-15
    CMAKE_COMPILER_ARGS="-DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15"
elif command -v gcc-15 >/dev/null 2>&1 && command -v g++-15 >/dev/null 2>&1; then
    echo "CoastalME: starting CMake for macOS (using system gcc-15 with OpenMP, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
    export CC=gcc-15
    export CXX=g++-15
    CMAKE_COMPILER_ARGS="-DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15"
elif brew list libomp >/dev/null 2>&1; then
    echo "CoastalME: starting CMake for macOS (using Clang with libomp, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
    CMAKE_COMPILER_ARGS=""
else
    echo "CoastalME: starting CMake for macOS (using default compiler - OpenMP may not be available, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
    echo "For OpenMP support, install GCC: brew install gcc"
    echo "Or install libomp for Clang: brew install libomp"
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

echo "Finished CMake for macOS (using gcc-15 with OpenMP, $buildtype build, $cshorelibrary CShore library, CShore input/output method=$cshoreinout)"
echo ""
echo "================================================================="

