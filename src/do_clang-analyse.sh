#!/bin/bash

# clang-tidy *.cpp  -checks=-*,clang-analyzer-*,-clang-analyzer-cplusplus*,-header-filter=.*  > 000_clang-analyze_advice.txt

clang-tidy *.cpp  -checks=-*,clang-analyzer-*,bugprone-*  > 000_clang-analyze_advice.txt
clang-tidy *.h  -checks=-*,clang-analyzer-*,bugprone-*  >> 000_clang-analyze_advice.txt

# cppcoreguidelines-*,misc-*
# modernize-*,openmp-*,performance-*,readability-*,portability-*
