#!/usr/bin/env bash

echo
echo "Starting nightly Trilinos development testing on p90n03: $(date)"
echo

selfdir=$(cd "$(dirname "$0")";pwd)

export TDD_CTEST_TEST_TYPE=Experimental
export CTEST_TEST_TIMEOUT=28800

ctest -V -S $selfdir/../TrilinosDriverDashboard.cmake

echo
echo "Ending nightly Trilinos development testing on p90n03: $(date)"
echo
