#!/bin/bash
export Trilinos_TRACK=ATDM
export Trilinos_REPOSITORY_LOCATION=https://github.com/trilinos/Trilinos.git
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/jenkins-driver.sh
