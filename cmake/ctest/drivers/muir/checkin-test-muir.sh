
#!/bin/bash
if [ "$TRILINOS_BASE_DIR" == "" ] ; then
  TRILINOS_BASE_DIR=../..
fi

TRILINOS_BASE_DIR_ABS=$(readlink -f $TRILINOS_BASE_DIR)

DRIVERS_BASE_DIR="$TRILINOS_BASE_DIR_ABS/Trilinos/cmake/ctest/drivers/fissile4"

echo "
-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON
" > COMMON.config

#
# Built-in Primary Tested (PT) --default-builds (DO NOT MODIFY)
#

echo "
" > MPI_DEBUG.config

echo "
-DCMAKE_C_COMPILER=gcc
" > SERIAL_RELEASE.config

#
# Standard Secondary Tested (ST) --st-extra-builds (DO NOT MODIFY)
#

echo "
-DTrilinos_ENABLE_DEBUG=ON
-DTPL_ENABLE_MPI=ON
" > MPI_DEBUG_ST.config

echo "
-DCMAKE_C_COMPILER=gcc
-DCMAKE_BUILD_TYPE=RELEASE
-DTrilinos_ENABLE_DEBUG=OFF
-DTPL_ENABLE_MPI=OFF
" > SERIAL_RELEASE_ST.config

#
# Create local defaults file if one does not exist
#

_LOCAL_CHECKIN_TEST_DEFAULTS=local-checkin-test-defaults.py
if [ -f $_LOCAL_CHECKIN_TEST_DEFAULTS ] ; then
  echo "File $_LOCAL_CHECKIN_TEST_DEFAULTS already exists, leaving it!"
else
  echo "Creating default file $_LOCAL_CHECKIN_TEST_DEFAULTS!"
  echo "
defaults = [
  \"-j12\",
  \"--ctest-timeout=180\",
  \"--skip-case-no-email\",
  ]
  " > $_LOCAL_CHECKIN_TEST_DEFAULTS
fi

#
# Invocation
#

$TRILINOS_BASE_DIR/Trilinos/checkin-test.py \
"$@"
