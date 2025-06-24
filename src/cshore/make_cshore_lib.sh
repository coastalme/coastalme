#!/bin/sh

# test if all three command line args provided
if [ $# -eq 3 ];then
          echo "building only specified cshore version"
          lib_type="${1}LIB"
          build_type=${2}
          com_type="${3}INOUT"

          # make clean
          # touch ../lib/libcshore.so
          make -B TO_BUILD=${lib_type} BUILD_VERSION=${build_type} INPUT_AND_OUTPUT=${com_type}
          # echo "make -B TO_BUILD=${lib_type} BUILD_VERSION=${build_type} INPUT_AND_OUTPUT=${com_type}"
else 
          echo "building all cshore versions"
          # Do ARGINOUT versions
          make TO_BUILD=STATICLIB BUILD_VERSION=DEBUG INPUT_AND_OUTPUT=ARGINOUT
          make TO_BUILD=STATICLIB BUILD_VERSION=RELEASE INPUT_AND_OUTPUT=ARGINOUT
          make TO_BUILD=SHAREDLIB BUILD_VERSION=DEBUG INPUT_AND_OUTPUT=ARGINOUT
          make TO_BUILD=SHAREDLIB BUILD_VERSION=RELEASE INPUT_AND_OUTPUT=ARGINOUT
          #
          # # Do ARGINBOTHOUT versions
          make TO_BUILD=STATICLIB BUILD_VERSION=DEBUG INPUT_AND_OUTPUT=ARGINBOTHOUT
          make TO_BUILD=STATICLIB BUILD_VERSION=RELEASE INPUT_AND_OUTPUT=ARGINBOTHOUT
          make TO_BUILD=SHAREDLIB BUILD_VERSION=DEBUG INPUT_AND_OUTPUT=ARGINBOTHOUT
          make TO_BUILD=SHAREDLIB BUILD_VERSION=RELEASE INPUT_AND_OUTPUT=ARGINBOTHOUT
          #
          # # Do FILEINOUT versions
          make TO_BUILD=STATICLIB BUILD_VERSION=DEBUG INPUT_AND_OUTPUT=FILEINOUT
          make TO_BUILD=STATICLIB BUILD_VERSION=RELEASE INPUT_AND_OUTPUT=FILEINOUT
          make TO_BUILD=SHAREDLIB BUILD_VERSION=DEBUG INPUT_AND_OUTPUT=FILEINOUT
          make TO_BUILD=SHAREDLIB BUILD_VERSION=RELEASE INPUT_AND_OUTPUT=FILEINOUT
fi
