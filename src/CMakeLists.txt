#########################################################################################
#
# CMakeLists.txt for CoastalME (Coastal Modelling Environment)
#
#########################################################################################
# Set CMake required version
cmake_minimum_required (VERSION 3.5.0 FATAL_ERROR)

# Create a .json file for use with clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#########################################################################################
# Detect the build environment
if (UNIX)
   if (UNIX AND NOT APPLE AND NOT CYGWIN)
      set (LINUX TRUE)
      message ("Creating Linux-Unix Makefile for CoastalME")
   endif ()

   if (CYGWIN)
      message ("Creating Cygwin Makefile for CoastalME")
   endif ()

   if (APPLE)
      set (LINUX TRUE)
      message ("Creating OSX-Unix Makefile for CoastalME")
   endif ()
endif ()

if (WIN32)
   message (FATAL_ERROR "CoastalME cannot be built on Windows, using the Visual Studio compiler... yet")
endif ()

if (MINGW)
   message (FATAL_ERROR "CoastalME cannot be built on Windows, using the MinGW compiler... yet")
endif ()

message ("")

#########################################################################################
# First set up project-specific defaults
project (CME)

set (CMAKE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")           # Assume that we are in the src folder
set (CMAKE_INSTALL_DIR "${CMAKE_SOURCE_DIR}/..")

# set (CMAKE_SOURCE_DIR "$ENV{HOME}/coast/CoastalME/src")      # May be changed by the user
# set (CMAKE_INSTALL_DIR "$ENV{HOME}/coast/CoastalME")         # May be changed by the user

file (GLOB CME_SOURCE_FILES *.cpp)
set (CME_EXECUTABLE cme)

# Set the path for CoastalME's own cmake modules
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#########################################################################################
# Detect the compiler and set compiler-specific flags etc.
message ("Checking compiler")

# We are using c++
enable_language (CXX)

# We need c++11, first set the c++11 flag (note that later on we check the compiler for c++11 capability)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Set the CoastalME CPU flag TODO Need to change this for different CPU architectures
if (APPLE)
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCPU -Dxarm64")
else()
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCPU -Dx86")
endif()

# Use 64-bit libraries
set (FIND_LIBRARY_USE_LIB64_PATHS true)

# Convert build type string to lower case
string (TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LC)

if (CMAKE_COMPILER_IS_GNUCC)
   message ("C++ compiler is ${CMAKE_CXX_COMPILER_ID} g++ ${CMAKE_CXX_COMPILER_VERSION}")

   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

   # For Debug builds under gcc
   if (CMAKE_BUILD_TYPE_LC STREQUAL "debug")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wpointer-arith -Wconversion -Wcast-qual -Wcast-align -Wwrite-strings -Wredundant-decls -Wno-strict-overflow -Wshadow -Wuninitialized -Wnull-dereference -Wformat -Wformat-overflow -Wformat-signedness -Wuseless-cast -Wempty-body -Wfloat-equal")

      # # Use these before a release ====================================================================================
      # if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
      #    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -fsanitize=shift -fsanitize=integer-divide-by-zero -fsanitize=unreachable -fsanitize=null -fsanitize=return -fsanitize=vla-bound -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fsanitize=bounds -fsanitize=alignment -fsanitize=object-size -fsanitize=vptr")
      #    #  -fsanitize=address
      # endif ()
      # # END Use these before a release ================================================================================

#      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
   endif ()

   # Finally set some gcc-specific linker flags
   set (CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed -Wl,--no-undefined")
endif ()

if (MSVC)
   message ("Compiler is Visual Studio")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
endif ()

message ("")

#########################################################################################
# If not specified, assume a Debug build
if (NOT CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE_LC debug)
endif ()

# Sort out compiler flags for various build types
if (CMAKE_BUILD_TYPE_LC STREQUAL "debug")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
   set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g")
   set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Debug CACHE PATH "Directory for Debug builds" FORCE)
   # For Debug builds on any platform, add the model's "internal" debug flag (used for extra logging) etc.
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_DEBUG")
endif ()

if (CMAKE_BUILD_TYPE_LC STREQUAL "release")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
   set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Release CACHE PATH "Directory for Release builds" FORCE)
endif ()

if (CMAKE_BUILD_TYPE_LC STREQUAL "gcov")
#   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 --coverage")
   set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/gcov CACHE PATH "Directory for gcov builds" FORCE)
endif ()

if (CMAKE_BUILD_TYPE_LC STREQUAL "valgrind")
#   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
   set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
   set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
   set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CME_BINARY_DIR}/Callgrind CACHE PATH "Directory for Callgrind builds" FORCE)
endif ()

#########################################################################################
# Next, find libraries
message ("Finding libraries")

#########################################################################################
# CoastalME cannot be built without GDAL, which is linked dynamically
find_package(GDAL REQUIRED)
if (GDAL_FOUND)
   # OK, we found GDAL
   message (STATUS "GDAL_INCLUDE_DIRS=${GDAL_INCLUDE_DIRS}")
   message (STATUS "GDAL_LIBRARIES=${GDAL_LIBRARIES}")
   message (STATUS "GDAL_VERSION=${GDAL_VERSION}")

   # Now check the GDAL version: CoastalME can only be built with GDAL 2.1 or above
   if (GDAL_VERSION VERSION_LESS "2.1")
      message (FATAL_ERROR "CoastalME requires GDAL 2.1 or later. You have version ${GDAL_VERSION}. Please upgrade")
   endif ()

   set (LIBS ${LIBS} ${GDAL_LIBRARIES})
   set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GDAL_INCLUDE_DIRS})
   message (STATUS "LIBS=${LIBS}")
   message (STATUS "CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")
endif ()

#########################################################################################
# CoastalME also requires the CShore library. This can be linked statically or dynamically (i.e. shared). If not specified, assume a shared library
if (NOT CSHORE_LIBRARY)
   set (CSHORE_LIBRARY_LC "shared")
else ()
   # The library link type was specified, so create a lower-case version
   string (TOLOWER "${CSHORE_LIBRARY}" CSHORE_LIBRARY_LC )
endif ()

# If not specified, assume that the CShore library communicates with CoastalME via argument passing
if (NOT CSHORE_INOUT)
   set (CSHORE_INOUT_LC "arg")
else ()
   # The communication method was specified, so make a lower-case version
   string (TOLOWER "${CSHORE_INOUT}" CSHORE_INOUT_LC )
endif ()

if (CSHORE_INOUT_LC STREQUAL "file")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCSHORE_FILE_INOUT")
elseif (CSHORE_INOUT_LC STREQUAL "arg")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCSHORE_ARG_INOUT")
elseif (CSHORE_INOUT_LC STREQUAL "both")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCSHORE_BOTH")
else ()
   message (FATAL_ERROR "Invalid value specified for communication with the CShore library: ${CSHORE_INOUT}")
endif ()

if (UNIX)
   # Put the correct version of the CShore library into ${CMAKE_SOURCE_DIR}/lib/libcshore.a"
   if (UNIX AND NOT APPLE AND NOT CYGWIN)
      if (CSHORE_INOUT_LC STREQUAL "file")
         if ((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL "valgrind"))
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.LINUX.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         else ()
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.LINUX.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         endif ()

      elseif (CSHORE_INOUT_LC STREQUAL "arg")
         if ((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL "valgrind"))
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.LINUX.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         else ()
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.LINUX.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         endif ()

      elseif (CSHORE_INOUT_LC STREQUAL "both")
         if ((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL "valgrind"))
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.LINUX.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         else ()
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.LINUX.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         endif ()
      endif ()
   endif ()

   if (CYGWIN)
      if (CSHORE_INOUT_LC STREQUAL "file")
         if (CMAKE_BUILD_TYPE_LC STREQUAL "debug")
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.CYGWIN.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.CYGWIN.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         else ()
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.CYGWIN.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.CYGWIN.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         endif ()

      elseif (CSHORE_INOUT_LC STREQUAL "arg")
         if (CMAKE_BUILD_TYPE_LC STREQUAL "debug")
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.CYGWIN.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.CYGWIN.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         else ()
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.CYGWIN.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.CYGWIN.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         endif ()

      elseif (CSHORE_INOUT_LC STREQUAL "both")
         if (CMAKE_BUILD_TYPE_LC STREQUAL "debug")
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.CYGWIN.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.CYGWIN.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         else ()
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.CYGWIN.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.CYGWIN.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         endif ()
      endif ()
   endif ()

   if (APPLE)
      if (CSHORE_INOUT_LC STREQUAL "file")
         if ((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL "valgrind"))
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.LINUX.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         else ()
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.LINUX.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX.FILEINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         endif ()

      elseif (CSHORE_INOUT_LC STREQUAL "arg")
         if ((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL "valgrind"))
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.LINUX.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         else ()
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.LINUX.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX.ARGINOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         endif ()

      elseif (CSHORE_INOUT_LC STREQUAL "both")
         if ((CMAKE_BUILD_TYPE_LC STREQUAL "debug") OR (CMAKE_BUILD_TYPE_LC STREQUAL "valgrind"))
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.DEBUG.LINUX.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.DEBUG.LINUX.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         else ()
            if (CSHORE_LIBRARY_LC STREQUAL "shared")
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.so.RELEASE.LINUX.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.so")
            else ()
               execute_process (COMMAND ln -s -f "${CMAKE_SOURCE_DIR}/lib/libcshore.a.RELEASE.LINUX.ARGINBOTHOUT" "${CMAKE_SOURCE_DIR}/lib/libcshore.a")
            endif ()
         endif ()
      endif ()
   endif ()
endif ()

if (WIN32)
   # TODO
endif ()

# OK, now find the CShore library (and GFortran, if the CShore linkage is static)
find_package (CShoreLib REQUIRED)

#message (STATUS "CSHORELIB_DIR=${CSHORELIB_DIR}")
set (LIBS ${LIBS} ${CSHORELIB_DIR})
message (STATUS "CSHORELIB_INCLUDE_DIR=${CSHORELIB_INCLUDE_DIR}")
set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CSHORELIB_INCLUDE_DIR})

if (CSHORE_LIBRARY_LC STREQUAL "static")
   set (LIBS ${LIBS} ${LIBGFORTRAN_LIBRARIES} ${LIBQUADMATH_LIBRARIES})
   set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${GFORTRAN_INCLUDE_DIR})
   #    include_directories (${GFORTRAN_INCLUDE_DIR})
   #    target_link_libraries (${CME_EXECUTABLE} ${LIBS} ${LIBGFORTRAN_LIBRARIES})
   #    target_link_libraries (${CME_EXECUTABLE} ${LIBS} ${LIBQUADMATH_LIBRARIES})

else ()
   # Assume that libcshore.a is either in one of the system library dirs, or is in a location specified by export LD_LIBRARY_PATH
   set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

   # We have the absolute path to the CShore library, now get the relative path from the install dir
   file (RELATIVE_PATH CSHORELIB_REL_FILEPATH ${CMAKE_INSTALL_DIR} ${CSHORELIB_DIR})

   # Remove the file name
   get_filename_component(CSHORELIB_REL_PATH ${CSHORELIB_REL_FILEPATH} DIRECTORY)

   # And get the absolute path (without the file name)
   get_filename_component(CSHORELIB_ABS_PATH ${CSHORELIB_DIR} DIRECTORY)

   # And set the installed executable's RPATH/RUNPATH to point to the CShore library: first as a relative path (two versions of this, $ORIGIN and $$ORIGIN, dunno why), next as as absolute path, and finally a fallback to a system library folder
   #set (CMAKE_INSTALL_RPATH "${CSHORELIB_DIR}:$ORIGIN/src/lib/:$$ORIGIN/src/lib/")
   set (CMAKE_INSTALL_RPATH "$ORIGIN/${CSHORELIB_REL_PATH}:$$ORIGIN/${CSHORELIB_REL_PATH}:${CSHORELIB_ABS_PATH}:/usr/local/lib")
endif ()

#########################################################################################
# The important bits
include_directories (SYSTEM ${CMAKE_INCLUDE_PATH})
include_directories (${CME_SOURCE_DIR})
add_executable (${CME_EXECUTABLE} ${CME_SOURCE_FILES})
target_link_libraries (${CME_EXECUTABLE} ${LIBS})

install (TARGETS ${CME_EXECUTABLE} RUNTIME DESTINATION ${CMAKE_INSTALL_DIR})

# Check to see if the compiler supports c++11
set_property (TARGET ${CME_EXECUTABLE} PROPERTY CXX_STANDARD 11)
set_property (TARGET ${CME_EXECUTABLE} PROPERTY CXX_STANDARD_REQUIRED ON)

#########################################################################################

#########################################################################################
# Tell the user what we have found
message ("")
message ("Writing Makefile")
message (STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message (STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}")
message (STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}")
message (STATUS "CMAKE_INSTALL_DIR = ${CMAKE_INSTALL_DIR}")
# message (STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message (STATUS "CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH}")
message (STATUS "LIBS=${LIBS}")
# message (STATUS "CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
# message (STATUS "CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
message (STATUS "CMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH}")

#########################################################################################
# Some extra messages
if (CMAKE_BUILD_TYPE_LC STREQUAL "valgrind")
   message ("")
   message ("When the build has finished, use valgrind/callgrind as follows.")
   message ("")
   message ("To check for memory leaks:")
   message ("valgrind --leak-check=yes --suppressions=system-libs.supp --track-origins=yes ./cme &> valgrind.txt")
   message ("Then look at valgrind.txt")
   message ("")
   message ("Or to check coverage:")
   message ("valgrind --tool=callgrind ./cme")
   message ("Then run:")
   message ("callgrind_annotate --auto=yes callgrind.out.XXXXX > ./profile/callgrind/callgrind.txt")
   message ("where XXXXX is the number of the callgrind.out.XXXXX that was produced by valgrind. Then look at ./profile/callgrind.txt")
   message ("")
endif ()

if (CMAKE_BUILD_TYPE_LC STREQUAL "gcov")
   message ("")
   message ("When the build has finished, use gcov/lcov as follows:")
   message ("")
   message ("./cme")
   message ("lcov --capture --directory ./src/CMakeFiles/cme.dir/ --output-file ./profile/lcov_output/coverage.info")
   message ("cd ./profile/lcov_output")
   message ("genhtml coverage.info")
   message ("")
   message ("Then look at index.html in your browser")
   message ("")
endif ()


#########################################################################################
