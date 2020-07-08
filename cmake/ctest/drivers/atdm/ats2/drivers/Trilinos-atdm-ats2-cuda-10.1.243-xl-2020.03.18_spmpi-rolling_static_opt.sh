#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Experimental
fi

if [[ "${Trilinos_INNER_ENABLE_TESTS}" == "" ]]; then
  # For XL, do not build tests and examples
  export Trilinos_INNER_ENABLE_TESTS=OFF
fi

if [[ "${Trilinos_SKIP_CTEST_ADD_TEST}" == "" ]] ; then
  # For XL, do not run tests
  export Trilinos_SKIP_CTEST_ADD_TEST=TRUE
fi

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats2/local-driver.sh
