#!/bin/sh

CTEST_BINARY_DIRECTORY=$1

grep -v '//' $CTEST_BINARY_DIRECTORY/CMakeCache.txt | grep -v '^#' | grep -v '^$' > $CTEST_BINARY_DIRECTORY/CMakeCache.clean.txt
