# ##############################################################################
#
# CMakeLists.txt for CoastalME (Coastal Modelling Environment)
#
# ##############################################################################
# Set CMake required version
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# This is needed if we want to use relative paths for e.g. the install dir. Such
# as ${CMAKE_SOURCE_DIR}/.. This was deprecated in CMake 3.31
if(${CMAKE_VERSION} VERSION_GREATER "3.31")
  cmake_policy(SET CMP0177 OLD)
endif()

# Create a .json file for use with clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ##############################################################################
# Detect the build environment
if(UNIX)
  if(APPLE)
    set(MACOS TRUE)
    message("Creating macOS-Unix Makefile for CoastalME")
  elseif(CYGWIN)
    message("Creating Cygwin Makefile for CoastalME")
  else()
    set(LINUX TRUE)
    message("Creating Linux-Unix Makefile for CoastalME")
  endif()

elseif(WIN32)
  message(
    FATAL_ERROR
      "CoastalME cannot be built on Windows, using the Visual Studio compiler... yet"
  )
elseif(MINGW)
  message(
    FATAL_ERROR
      "CoastalME cannot be built on Windows, using the MinGW compiler... yet")
endif()

message("")

# ##############################################################################
# First set up project-specific defaults
project(CME)

set(CMAKE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}") # Assume that we are in the
                                                    # src folder
set(CMAKE_INSTALL_DIR "${CMAKE_SOURCE_DIR}/..") # Install in folder above

file(GLOB CME_SOURCE_FILES *.cpp) # Find all cpp files in src directory
set(CME_EXECUTABLE cme) # Set executale name

# Set the path for CoastalME's own cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# ##############################################################################
# Detect the compiler and set compiler-specific flags etc.
message("Checking compiler")

# We are using C++
enable_language(CXX)

# We need c++17 (for filesystem stuff, not that important howerver), so first
# set the flag (note that later on we check the compiler for c++17 capability)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# find system CPU architecture
execute_process(
  COMMAND uname -m
  COMMAND tr -d '\n'
  OUTPUT_VARIABLE ARCHITECTURE)
# Set CPU flag
if(${ARCHITECTURE} STREQUAL "arm64")
  message("Buiding for arm64")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCPU -Dxarm64")
elseif(${ARCHITECTURE} STREQUAL "x86_64")
  message("Buiding for x86 64bit")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCPU -Dx86")
else()
  message(FATAL_ERROR "CoastalME cannot be built on 32bit systems")
endif()

# Use 64-bit libraries
set(FIND_LIBRARY_USE_LIB64_PATHS true)

# Convert build type string to lower case
string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

# Compiler specific settings
if(CMAKE_COMPILER_IS_GNUCC)
  message(
    "C++ compiler is ${CMAKE_CXX_COMPILER_ID} g++ ${CMAKE_CXX_COMPILER_VERSION}"
  )
  # GCC compiler defaults
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

  # Build type specific settings
  if(CMAKE_BUILD_TYPE_LC STREQUAL "debug")
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} \
        -W -Wall -Wextra \
        -Wpointer-arith -Wcast-qual -Wcast-align \
        -Wconversion \
        -Wwrite-strings -Wredundant-decls \
        -Wshadow -Wuninitialized -Wnull-dereference \
        -Wformat -Wuseless-cast -Wempty-body -Wfloat-equal \
        -Wformat-overflow -Wformat-signedness \
        -Wno-strict-overflow \
        ")

    # Also use these with debug builds to check for vulnerabilities at compile
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "13.0")
      set(CMAKE_CXX_FLAGS
          "${CMAKE_CXX_FLAGS} \
          -Wformat=2 \
          -Wimplicit-fallthrough \
          -Werror=format-security \
          -U_FORTIFY_SOURCE \
          -D_FORTIFY_SOURCE=3 \
          -D_GLIBCXX_ASSERTIONS \
          -fstrict-flex-arrays=3 \
          -fstack-clash-protection \
          -fstack-protector-strong \
          -Wduplicated-cond \
          -Wduplicated-branches\
          -Wlogical-op \
          -Wrestrict\
          -Wnull-dereference \
          -Wold-style-cast \
          -Woverloaded-virtual\
          -Wsuggest-override \
          -Wmisleading-indentation \
          ")
      # -Wtrampolines -Wbidi-chars=any -fPIE -pie (for x86_64)
      # -fcf-protection=full
      if(LINUX)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} \
                    -Wl,-z,nodlopen -Wl,-z,noexecstack -Wl,-z,relro \
                    -Wl,-z,now -Wl,--as-needed -Wl,--no-copy-dt-needed-entries \
                    ")
      endif()
    endif()
    # ===============================================================================================================

    # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
  endif()

  if(CMAKE_BUILD_TYPE_LC STREQUAL "prerelease")
    # Pre-release build using gcc, do this before a release to check for memory
    # access problems. Note that -fsanitize=address (and the other sanitize
    # options?) DO NOT work if run under gdb #
    # ===============================================================================================================
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
      set(CMAKE_CXX_FLAGS
          "${CMAKE_CXX_FLAGS}\
          -fsanitize=undefined \
          -fsanitize=shift \
          -fsanitize=integer-divide-by-zero \
          -fsanitize=unreachable \
          -fsanitize=null \
          -fsanitize=return \
          -fsanitize=vla-bound \
          -fsanitize=float-divide-by-zero \
          -fsanitize=float-cast-overflow \
          -fsanitize=bounds \
          -fsanitize=alignment \
          -fsanitize=object-size \
          -fsanitize=vptr \
          -fsanitize=address \
          ")
      # Note  -fsanitize-ignorelist=sanitize_ignorelist.txt not yet working for
      # gcc, need this to hide GDAL errors
    endif()
    # #
    # ===============================================================================================================
  endif()

  if(CMAKE_BUILD_TYPE_LC STREQUAL "release")
    # Production build using gcc #
    # ===============================================================================================================
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} \
        -O3 \
        -Wall \
        -fno-delete-null-pointer-checks \
        -fno-strict-overflow \
        -fno-strict-aliasing \
        -ftrivial-auto-var-init=zero \
        -flto=auto \
        -ffast-math \
        -funroll-loops \
        -DNDEBUG \
        ")
    #
    # ===============================================================================================================
  endif()

  # Finally set some gcc-specific linkerrflags, for all build types
  if(LINUX)
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed -Wl,--no-undefined")
  endif()

elseif(MSVC)
  message("Compiler is Visual Studio")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
endif()

message("")

# ##############################################################################
# If not specified, assume a Debug build
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE_LC debug)
endif()

# Sort out compiler flags for various build types
if(CMAKE_BUILD_TYPE_LC STREQUAL "debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${CME_BINARY_DIR}/Debug
      CACHE PATH "Directory for Debug builds" FORCE)
  # For Debug builds on any platform, add the model's "internal" debug flag
  # (used for extra logging) etc.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_DEBUG")
endif()

if(CMAKE_BUILD_TYPE_LC STREQUAL "release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${CME_BINARY_DIR}/Release
      CACHE PATH "Directory for Release builds" FORCE)
endif()

if(CMAKE_BUILD_TYPE_LC STREQUAL "gcov")
  # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 --coverage")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${CME_BINARY_DIR}/gcov
      CACHE PATH "Directory for gcov builds" FORCE)
endif()

if(CMAKE_BUILD_TYPE_LC STREQUAL "valgrind")
  # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${CME_BINARY_DIR}/Callgrind
      CACHE PATH "Directory for Callgrind builds" FORCE)
endif()

# ##############################################################################
# Next, find libraries
message("Finding libraries")

# ##############################################################################
# Find OpenMP for parallelization First try to find OpenMP using the standard
# CMake module
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found via find_package: ${OpenMP_CXX_VERSION}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(LIBS ${LIBS} OpenMP::OpenMP_CXX)
  message(STATUS "OpenMP enabled for parallel grid operations")
else()
  # Standard detection failed, try manual detection based on compiler
  message(
    STATUS "Standard OpenMP detection failed, attempting manual detection...")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Test if GCC supports OpenMP
    message(STATUS "Testing GCC OpenMP support...")
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_FLAGS "-fopenmp")
    set(CMAKE_REQUIRED_LIBRARIES "-fopenmp")
    check_cxx_source_compiles(
      "
         #include <omp.h>
         int main() {
            #pragma omp parallel
            {
               int tid = omp_get_thread_num();
            }
            return 0;
         }"
      OPENMP_TEST_COMPILE)

    if(OPENMP_TEST_COMPILE)
      message(STATUS "GCC OpenMP support confirmed - adding -fopenmp flags")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
      set(OPENMP_ENABLED TRUE)
    else()
      message(STATUS "GCC does not support OpenMP")
      set(OPENMP_ENABLED FALSE)
    endif()

  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Test Clang OpenMP support (may need libomp)
    message(STATUS "Testing Clang OpenMP support...")
    include(CheckCXXSourceCompiles)

    # Try with -fopenmp first
    set(CMAKE_REQUIRED_FLAGS "-fopenmp")
    set(CMAKE_REQUIRED_LIBRARIES "-fopenmp")
    check_cxx_source_compiles(
      "
         #include <omp.h>
         int main() {
            #pragma omp parallel
            {
               int tid = omp_get_thread_num();
            }
            return 0;
         }"
      CLANG_OPENMP_FOPENMP)

    if(CLANG_OPENMP_FOPENMP)
      message(STATUS "Clang OpenMP support confirmed with -fopenmp")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
      set(OPENMP_ENABLED TRUE)
    else()
      # Try with -fopenmp=libomp (for Homebrew clang)
      set(CMAKE_REQUIRED_FLAGS "-fopenmp=libomp")
      set(CMAKE_REQUIRED_LIBRARIES "-fopenmp=libomp")
      check_cxx_source_compiles(
        "
            #include <omp.h>
            int main() {
               #pragma omp parallel
               {
                  int tid = omp_get_thread_num();
               }
               return 0;
            }"
        CLANG_OPENMP_LIBOMP)

      if(CLANG_OPENMP_LIBOMP)
        message(STATUS "Clang OpenMP support confirmed with -fopenmp=libomp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp=libomp")
        set(OPENMP_ENABLED TRUE)
      else()
        message(STATUS "Clang does not support OpenMP")
        set(OPENMP_ENABLED FALSE)
      endif()
    endif()

  else()
    message(STATUS "Unknown compiler - attempting generic OpenMP detection")
    set(OPENMP_ENABLED FALSE)
  endif()

  # Reset CMAKE_REQUIRED variables
  unset(CMAKE_REQUIRED_FLAGS)
  unset(CMAKE_REQUIRED_LIBRARIES)

  if(NOT OPENMP_ENABLED)
    message(WARNING "OpenMP not found - grid operations will run sequentially")
    message(STATUS "To enable OpenMP on macOS:")
    message(STATUS "  1. Install GCC: brew install gcc")
    message(
      STATUS "  2. Use GCC for compilation: CC=gcc-15 CXX=g++-15 cmake ...")
    message(STATUS "  3. Or install libomp for Clang: brew install libomp")
  else()
    message(STATUS "OpenMP enabled for parallel grid operations")
  endif()
endif()

# ##############################################################################
# CoastalME cannot be built without GDAL, which is linked dynamically
find_package(GDAL REQUIRED)
if(GDAL_FOUND)
  # OK, we found GDAL
  message(STATUS "GDAL_INCLUDE_DIRS=${GDAL_INCLUDE_DIRS}")
  message(STATUS "GDAL_LIBRARIES=${GDAL_LIBRARIES}")
  message(STATUS "GDAL_VERSION=${GDAL_VERSION}")

  # Now check the GDAL version: CoastalME can only be built with GDAL 2.1 or
  # above
  if(GDAL_VERSION VERSION_LESS "2.1")
    message(
      FATAL_ERROR
        "CoastalME requires GDAL 2.1 or later. You have version ${GDAL_VERSION}. Please upgrade"
    )
  endif()

  set(LIBS ${LIBS} ${GDAL_LIBRARIES})
  set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GDAL_INCLUDE_DIRS})
  message(STATUS "LIBS=${LIBS}")
  message(STATUS "CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")
endif()

# ##############################################################################
# CoastalME also requires the CShore library. This can be linked statically or
# dynamically (i.e. shared). If not specified, assume a shared library
if(NOT CSHORE_LIBRARY)
  set(CSHORE_LIBRARY_LC "shared")
else()
  # The library link type was specified, so create a lower-case version
  string(TOLOWER "${CSHORE_LIBRARY}" CSHORE_LIBRARY_LC)
endif()

# If not specified, assume that the CShore library communicates with CoastalME
# via argument passing
if(NOT CSHORE_INOUT)
  set(CSHORE_INOUT_LC "arg")
else()
  # The communication method was specified, so make a lower-case version
  string(TOLOWER "${CSHORE_INOUT}" CSHORE_INOUT_LC)
endif()

if(CSHORE_INOUT_LC STREQUAL "file")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCSHORE_FILE_INOUT")
elseif(CSHORE_INOUT_LC STREQUAL "arg")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCSHORE_ARG_INOUT")
elseif(CSHORE_INOUT_LC STREQUAL "both")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCSHORE_BOTH")
else()
  message(
    FATAL_ERROR
      "Invalid value specified for communication with the CShore library: ${CSHORE_INOUT}"
  )
endif()

if(UNIX)
  # Put the correct version of the CShore library into
  # ${CMAKE_SOURCE_DIR}/lib/libcshore.a"
  if(UNIX
     AND NOT APPLE
     AND NOT CYGWIN)
    if(CSHORE_INOUT_LC STREQUAL "file")
      if((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL
                                                    "valgrind"))
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.LINUX.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      else()
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.LINUX.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      endif()

    elseif(CSHORE_INOUT_LC STREQUAL "arg")
      if((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL
                                                    "valgrind"))
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.LINUX.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      else()
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.LINUX.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      endif()

    elseif(CSHORE_INOUT_LC STREQUAL "both")
      if((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL
                                                    "valgrind"))
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.LINUX.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      else()
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.LINUX.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      endif()
    endif()
  endif()

  if(CYGWIN)
    if(CSHORE_INOUT_LC STREQUAL "file")
      if(CMAKE_BUILD_TYPE_LC STREQUAL "debug")
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.CYGWIN.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.CYGWIN.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      else()
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.CYGWIN.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.CYGWIN.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      endif()

    elseif(CSHORE_INOUT_LC STREQUAL "arg")
      if(CMAKE_BUILD_TYPE_LC STREQUAL "debug")
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.CYGWIN.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.CYGWIN.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      else()
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.CYGWIN.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.CYGWIN.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      endif()

    elseif(CSHORE_INOUT_LC STREQUAL "both")
      if(CMAKE_BUILD_TYPE_LC STREQUAL "debug")
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.CYGWIN.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.CYGWIN.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      else()
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.CYGWIN.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.CYGWIN.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      endif()
    endif()
  endif()

  if(APPLE)
    if(CSHORE_INOUT_LC STREQUAL "file")
      if((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL
                                                    "valgrind"))
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.MACOS.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.MACOS.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      else()
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.MACOS.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.MACOS.FILEINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      endif()

    elseif(CSHORE_INOUT_LC STREQUAL "arg")
      if((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL
                                                    "valgrind"))
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.MACOS.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.MACOS.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      else()
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.MACOS.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.MACOS.ARGINOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      endif()

    elseif(CSHORE_INOUT_LC STREQUAL "both")
      if((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL
                                                    "valgrind"))
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.MACOS.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.MACOS.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      else()
        if(CSHORE_LIBRARY_LC STREQUAL "shared")
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.MACOS.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
        else()
          execute_process(
            COMMAND
              ln -s -f
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.MACOS.ARGINBOTHOUT"
              "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
        endif()
      endif()
    endif()
  endif()
endif()

if(WIN32)
  # TODO
endif()

# OK, now find the CShore library (and GFortran, if the CShore linkage is
# static)
find_package(CShoreLib REQUIRED)

# message (STATUS "CSHORELIB_DIR=${CSHORELIB_DIR}")
set(LIBS ${LIBS} ${CSHORELIB_DIR})
message(STATUS "CSHORELIB_INCLUDE_DIR=${CSHORELIB_INCLUDE_DIR}")
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CSHORELIB_INCLUDE_DIR})

if(CSHORE_LIBRARY_LC STREQUAL "static")
  set(LIBS ${LIBS} ${LIBGFORTRAN_LIBRARIES} ${LIBQUADMATH_LIBRARIES})
  set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GFORTRAN_INCLUDE_DIR})
  # include_directories (${GFORTRAN_INCLUDE_DIR}) target_link_libraries
  # (${CME_EXECUTABLE} ${LIBS} ${LIBGFORTRAN_LIBRARIES}) target_link_libraries
  # (${CME_EXECUTABLE} ${LIBS} ${LIBQUADMATH_LIBRARIES})

else()
  # Assume that libcshore.a is either in one of the system library dirs, or is
  # in a location specified by export LD_LIBRARY_PATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # We have the absolute path to the CShore library, now get the relative path
  # from the install dir
  file(RELATIVE_PATH CSHORELIB_REL_FILEPATH ${CMAKE_INSTALL_DIR}
       ${CSHORELIB_DIR})

  # Remove the file name
  get_filename_component(CSHORELIB_REL_PATH ${CSHORELIB_REL_FILEPATH} DIRECTORY)

  # And get the absolute path (without the file name)
  get_filename_component(CSHORELIB_ABS_PATH ${CSHORELIB_DIR} DIRECTORY)

  # Get GDAL library directory
  get_filename_component(GDAL_LIB_DIR ${GDAL_LIBRARIES} DIRECTORY)

  # Set platform-specific RPATH configuration
  if(APPLE)
    # macOS uses @loader_path instead of $ORIGIN
    set(CMAKE_INSTALL_RPATH
        "@loader_path/${CSHORELIB_REL_PATH}:${CSHORELIB_ABS_PATH}:${GDAL_LIB_DIR}:/opt/homebrew/lib:/usr/local/lib"
    )
  else()
    # Linux and other Unix systems use $ORIGIN
    set(CMAKE_INSTALL_RPATH
        "$ORIGIN/${CSHORELIB_REL_PATH}:$$ORIGIN/${CSHORELIB_REL_PATH}:${CSHORELIB_ABS_PATH}:${GDAL_LIB_DIR}:/usr/local/lib"
    )
  endif()
endif()

# ##############################################################################
# The important bits
include_directories(SYSTEM ${CMAKE_INCLUDE_PATH})
include_directories(${CME_SOURCE_DIR})
add_executable(${CME_EXECUTABLE} ${CME_SOURCE_FILES})
target_link_libraries(${CME_EXECUTABLE} ${LIBS})

install(TARGETS ${CME_EXECUTABLE} RUNTIME DESTINATION ${CMAKE_INSTALL_DIR})

# Check to see if the compiler supports c++17
set_property(TARGET ${CME_EXECUTABLE} PROPERTY CXX_STANDARD 17)
set_property(TARGET ${CME_EXECUTABLE} PROPERTY CXX_STANDARD_REQUIRED ON)

# ##############################################################################

# ##############################################################################
# Tell the user what we have found
message("")
message("Writing Makefile")
message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_INSTALL_DIR = ${CMAKE_INSTALL_DIR}")
# message (STATUS
# "CMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")
message(STATUS "LIBS=${LIBS}")
# message (STATUS "CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}") message (STATUS
# "CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
message(STATUS "CMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH}")

# ##############################################################################
# Some extra messages
if(CMAKE_BUILD_TYPE_LC STREQUAL "valgrind")
  message("")
  message("When the build has finished, use valgrind/callgrind as follows.")
  message("")
  message("To check for memory leaks:")
  message(
    "valgrind --leak-check=yes --suppressions=system-libs.supp --track-origins=yes ./cme &> valgrind.txt"
  )
  message("Then look at valgrind.txt")
  message("")
  message("Or to check coverage:")
  message("valgrind --tool=callgrind ./cme")
  message("Then run:")
  message(
    "callgrind_annotate --auto=yes callgrind.out.XXXXX > ./profile/callgrind/callgrind.txt"
  )
  message(
    "where XXXXX is the number of the callgrind.out.XXXXX that was produced by valgrind. Then look at ./profile/callgrind.txt"
  )
  message("")
endif()

if(CMAKE_BUILD_TYPE_LC STREQUAL "gcov")
  message("")
  message("When the build has finished, use gcov/lcov as follows:")
  message("")
  message("./cme")
  message(
    "lcov --capture --directory ./src/CMakeFiles/cme.dir/ --output-file ./profile/lcov_output/coverage.info"
  )
  message("cd ./profile/lcov_output")
  message("genhtml coverage.info")
  message("")
  message("Then look at index.html in your browser")
  message("")
endif()

# ##############################################################################
