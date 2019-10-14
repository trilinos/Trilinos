#!/bin/bash

cd ~/development/TEST/rol-trilinos
FAILED=0

./Trilinos/cmake/ctest/drivers/rol-sems/rol-gcc-default-debug.sh
FILE=`ls /home/dridzal/development/TEST/rol-trilinos/BUILD/Testing/Temporary/LastTest_* 2>/dev/null`
if [ -f FILE ]; then
  if grep -q "End Result: TEST FAILED" ${FILE};then
    FAILED=1
  fi
fi

./Trilinos/cmake/ctest/drivers/rol-sems/rol-gcc-default-release.sh
FILE=`ls /home/dridzal/development/TEST/rol-trilinos/BUILD/Testing/Temporary/LastTest_* 2>/dev/null`
if [ -f FILE ]; then
  if grep -q "End Result: TEST FAILED" ${FILE};then
    FAILED=1
  fi
fi

./Trilinos/cmake/ctest/drivers/rol-sems/rol-downstream-gcc-default-release.sh
FILE=`ls /home/dridzal/development/TEST/rol-trilinos/BUILD/Testing/Temporary/LastTest_* 2>/dev/null`
if [ -f FILE ]; then
  if grep -q "End Result: TEST FAILED" ${FILE};then
    FAILED=1
  fi
fi

./Trilinos/cmake/ctest/drivers/rol-sems/rol-intel-latest-release.sh
FILE=`ls /home/dridzal/development/TEST/rol-trilinos/BUILD/Testing/Temporary/LastTest_* 2>/dev/null`
if [ -f FILE ]; then
  if grep -q "End Result: TEST FAILED" ${FILE};then
    FAILED=1
  fi
fi

./Trilinos/cmake/ctest/drivers/rol-sems/rol-gcc-latest-release.sh
FILE=`ls /home/dridzal/development/TEST/rol-trilinos/BUILD/Testing/Temporary/LastTest_* 2>/dev/null`
if [ -f FILE ]; then
  if grep -q "End Result: TEST FAILED" ${FILE};then
    FAILED=1
  fi
fi

./Trilinos/cmake/ctest/drivers/rol-sems/sendMail.sh ${FAILED}
