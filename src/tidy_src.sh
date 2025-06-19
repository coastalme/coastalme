#!/bin/bash

for i in *.h; do echo "Formatting $i" && clang-format -i "$i"; done
for i in *.cpp; do echo "Formatting $i" && clang-format -i "$i"; done

