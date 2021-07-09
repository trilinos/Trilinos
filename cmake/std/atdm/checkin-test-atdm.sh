#!/bin/bash

echo
echo "***"
echo "*** $0 " "$@"
echo "***"
echo

#
# Get the location of the base Trilinos directory
#

if [ "$ATDM_TRILINOS_DIR" == "" ] ; then
  # Grab from the symlink (only works on Linux)
  _ABS_FILE_PATH=`readlink -f $0` || \
   echo "Could not follow symlink to set TRILINOS_DIR!"
  if [ "$_ABS_FILE_PATH" != "" ] ; then
    export STD_ATDM_DIR=`dirname $_ABS_FILE_PATH`
    export ATDM_TRILINOS_DIR=$STD_ATDM_DIR/../../..
  fi
fi
echo "ATDM_TRILINOS_DIR = '$ATDM_TRILINOS_DIR'"

if [ "$ATDM_TRILINOS_DIR" == "" ] ; then
  echo "ERROR: Cannot determine TRILINOS_DIR (you must be on a non-Linux system or you must have copied the script instead of symlinking it as per instructions)."
  exit 1
fi

if [ "$ATDM_TRIBITS_DIR" == "" ] ; then
  export ATDM_TRIBITS_DIR=$ATDM_TRILINOS_DIR/cmake/tribits
fi
echo "ATDM_TRIBITS_DIR = '${ATDM_TRIBITS_DIR}'"

#
# Load a default env for the system
#

if [ "$ATDM_CHT_DEFAULT_ENV" == "" ] ; then
  ATDM_CHT_DEFAULT_ENV=default
fi
#echo "ATDM_CHT_DEFAULT_ENV = ${ATDM_CHT_DEFAULT_ENV}"

echo
echo "Load some env to get python, cmake, etc ..."
echo
source $STD_ATDM_DIR/load-env.sh ${ATDM_CHT_DEFAULT_ENV}
# NOTE: Above, it does not matter which env you load.  Any of them will
# provide the right python, cmake, etc.

ATDM_CHT_HELP_STR="
usage: ./checkin-test-atdm.sh \\
         <build-name-keys0> <build-name-keys1> ... <build-name-keysn> \\
         [other checkin-test options]

This script drives local configure, build and testing for each of the
<build-name-keyi> builds.

If just 'all' is passsed in for the <build-name-keys> list, then the list of all
of the supported jobs for that current system will be loaded from
<system_name>/all_supported_builds.sh.

To use this script, just symlink this into any desired directory like:

   cd <some-base-dir>/
   ln -s $TRILNOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .

Then, for example, to locally test a few builds for just Kokkos use:

  ./checkin-test-atdm.sh \\
     gnu-opt-openmp intel-debug-serial \\
     --no-enable-fwd-packages --enable-packages=Kokkos \\
     --local-do-all

This will only send email for the final check of all of the builds specified
(email can be turned off by passing in --send-email-to=).

Note that this will auatomatically use the full number processors for building
and running tests as specified in the <system_name>/environment.sh file.
Therefore, be caseful to only run this on a dedicated compute node for the
system and not the login node or it will completely fill up the login node.

To reproduce any build just do:

  cd <build-name-keys>/
  source load-env.sh
  ./do-configure
  make -j16
  ctest -j16

or any combination of configure, build, and test commands you want.

NOTE: To see checkin-test.py --help, run:

  $STD_ATDM_DIR/../../../checkin-test.py --help
"

#
# A) Parse the arguments
#

# A.1) Pull off the initial <build-name-keysi> arguments

ATDM_BUILD_NAME_KEYS_LIST=
while [[ ! "$1" == "--"* ]] && [[ ! "$1" == "" ]] ; do
  if [[ "$ATDM_BUILD_NAME_KEYS_LIST" == "" ]] ; then
    ATDM_BUILD_NAME_KEYS_LIST="$1"
  else
    ATDM_BUILD_NAME_KEYS_LIST="$ATDM_BUILD_NAME_KEYS_LIST $1"
  fi
  shift
done

# A.2) Look for --pull and --push arguments and pull them out

ATDM_CHT_FOUND_HELP=0
ATDM_CHT_FOUND_PULL=0
ATDM_CHT_FOUND_PUSH=0
ATDM_CHT_SEND_EMAIL_TO_ARG=
ATDM_CHT_ENABLE_PACKAGES_ARG=

for ATDM_CHT_CURENT_ARG in "$@" ; do
  if [[ "$ATDM_CHT_CURENT_ARG" == "--help" ]] ; then
    #echo "Found --help"
    ATDM_CHT_FOUND_HELP=1
  elif [[ "$ATDM_CHT_CURENT_ARG" == "-h" ]] ; then
    #echo "Found -h"
    ATDM_CHT_FOUND_HELP=1
  elif [[ "$ATDM_CHT_CURENT_ARG" == "--pull" ]] ; then
    #echo "Found --pull"
    ATDM_CHT_FOUND_PULL=1
  elif [[ "$ATDM_CHT_CURENT_ARG" == "--push" ]] ; then
    #echo "Found --push"
    ATDM_CHT_FOUND_PUSH=1
  elif [[ "$ATDM_CHT_CURENT_ARG" == "--send-email-to"* ]] ; then
    ATDM_CHT_SEND_EMAIL_TO_ARG="$ATDM_CHT_CURENT_ARG"
  elif [[ "$ATDM_CHT_CURENT_ARG" == "--enable-packages"* ]] ; then
    ATDM_CHT_ENABLE_PACKAGES_ARG="$ATDM_CHT_CURENT_ARG"
  fi
done

if [[ "$ATDM_CHT_FOUND_HELP" == "1" ]] ; then
  echo "$ATDM_CHT_HELP_STR"
  exit 0
fi

# A.3) Inspect the <build-name-keys> args and deal with 'all'

if [[ "$ATDM_BUILD_NAME_KEYS_LIST" == "" ]] ; then
  echo
  echo "Error, at least one <build-name-keys> (e.g. gnu-opt-openmp) argument is required!"
  exit 1
fi

if [[ "$ATDM_BUILD_NAME_KEYS_LIST" == "all" ]] ; then
  export ATDM_CONFIG_ALL_SUPPORTED_BUILDS=
  source $STD_ATDM_DIR/$ATDM_CONFIG_SYSTEM_NAME/all_supported_builds.sh
  ATDM_BUILD_NAME_KEYS_LIST="${ATDM_CONFIG_ALL_SUPPORTED_BUILDS[@]}"
fi

ATDM_BUILD_NAME_KEYS_COMMA_LIST=
ATDM_NUM_BULDS=0
for ATDM_BUILD_NAME_KEYS in $ATDM_BUILD_NAME_KEYS_LIST ; do
  if [[ "$ATDM_BUILD_NAME_KEYS_COMMA_LIST" == "" ]] ; then
    ATDM_BUILD_NAME_KEYS_COMMA_LIST="$ATDM_BUILD_NAME_KEYS"
  else
    ATDM_BUILD_NAME_KEYS_COMMA_LIST="$ATDM_BUILD_NAME_KEYS_COMMA_LIST,$ATDM_BUILD_NAME_KEYS"
  fi
  ATDM_NUM_BULDS=$((ATDM_NUM_BULDS+1))
done

#
# B) Remove log files
#

if [ -f checkin-test.final.out ] ; then
  rm checkin-test.final.out
fi

for ATDM_BUILD_NAME_KEYS in $ATDM_BUILD_NAME_KEYS_LIST ; do
  if [ -f checkin-test.$ATDM_BUILD_NAME_KEYS.out ] ; then
    rm checkin-test.$ATDM_BUILD_NAME_KEYS.out
  fi
done

#
# C) Do an initial pull
#

# ToDo: Implement

#
# D) Write default local-checkin-test-defaults.py file
#

_LOCAL_CHECKIN_TEST_DEFAULTS=local-checkin-test-defaults.py
if [ -f $_LOCAL_CHECKIN_TEST_DEFAULTS ] ; then
  echo "File $_LOCAL_CHECKIN_TEST_DEFAULTS already exists, leaving it!"
else
  echo "Creating default file $_LOCAL_CHECKIN_TEST_DEFAULTS!"
  echo "
defaults = [
  \"--enable-all-package=off\",
  \"--no-enable-fwd-packages\",
  ]
  " > $_LOCAL_CHECKIN_TEST_DEFAULTS
fi

#
# E) Loop over individual builds and run them
#

echo
echo "Running configure, build, and/or testing for $ATDM_NUM_BULDS builds:"
for build_name in $ATDM_BUILD_NAME_KEYS_LIST ; do
  echo "    ${build_name}"
done
echo

ATDM_CHT_BUILD_CASE_IDX=0
for ATDM_BUILD_NAME_KEYS in $ATDM_BUILD_NAME_KEYS_LIST ; do
  echo
  echo "***"
  echo "*** $ATDM_CHT_BUILD_CASE_IDX) Process build case $ATDM_BUILD_NAME_KEYS"
  echo "***"
  echo
  $STD_ATDM_DIR/utils/checkin-test-atdm-single.sh $ATDM_BUILD_NAME_KEYS \
    --default-builds= --allow-no-pull --send-email-to="" \
    --test-categories=NIGHTLY --ctest-timeout=600 \
    "$@"
  if [[ "$?" == "0" ]] ; then
    echo "$ATDM_BUILD_NAME_KEYS: PASSED!"
  else
    echo "$ATDM_BUILD_NAME_KEYS: FAILED!"
  fi
  ATDM_CHT_BUILD_CASE_IDX=$((ATDM_CHT_BUILD_CASE_IDX+1))
done

#
# F) Collect the results from all the builds
#

echo
echo "Collect and report final results:"
echo
echo "  ==> See output file checkin-test.final.out" 
echo

${ATDM_TRIBITS_DIR}/ci_support/checkin-test.py \
  --src-dir=$ATDM_TRILINOS_DIR \
  --default-builds= --st-extra-builds=$ATDM_BUILD_NAME_KEYS_COMMA_LIST \
  --allow-no-pull "$ATDM_CHT_ENABLE_PACKAGES_ARG" \
  $ATDM_CHT_SEND_EMAIL_TO_ARG \
  --log-file=checkin-test.final.out \
  &> /dev/null

ATDM_CHT_RETURN_CODE=$?

# NOTE The return value will be 0 if everything passed!

#
# G) Print final status
#

echo
grep -A 1000 "Commit status email being sent" checkin-test.final.out \
  | grep -B 1000 "Commits for repo" \
  | grep -v "Commit status email being sent" \
  | grep -v "Commits for repo"
echo
grep "REQUESTED ACTIONS" checkin-test.final.out

# ToDo: Add logic to --push if requested!

exit $ATDM_CHT_RETURN_CODE
