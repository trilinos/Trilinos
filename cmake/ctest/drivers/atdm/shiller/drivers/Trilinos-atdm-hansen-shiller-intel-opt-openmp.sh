#!/bin/bash
export CTEST_TEST_TYPE=Experimental
export Trilinos_TRACK=Experimental
export Trilinos_REPOSITORY_LOCATION=https://github.com/trilinos/Trilinos.git
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/jenkins-driver.sh
