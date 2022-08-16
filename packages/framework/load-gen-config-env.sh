#
# Source this script to load a GenConfig env and generate the CMake fragment file
#
# Usage:
#
#   cd <any-build-dir>/
#
#   source <trilinosDir>/packages/framework/load-gen-config-env.sh \
#     <full-build-name> [--ci-mode]
#
# where:
#
#   <full-build-name>: Full build name taken off a Trilinos PR comment or CDash
#                      (minus the PR prefix 'PR-<prnum>-test-' and the Jenkins build
#                      ID suffix '-<jobid>' at the end)
#
# This has the following effect:
#
#   * Runs a new bash shell (unless pasing in --ci-mode)
#   * Loads the correct env matching (modules and env vars)
#   * Writes a cmake fragment file GenConfigSettings.cmake
#   * Sets the env var TRILINOS_DIR and any other needed vars
#
# Once in the new shell, just (re)configure and (re)build with:
#
#   rm -rf CMake*
#
#   cmake -C GenConfigSettings.cmake \
#     -D Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES=OFF \
#     -D Trilinos_ENABLE_TESTS=ON \
#     -D Trilinos_ENABLE_<pkg1>=ON \
#     [other options] \
#     $TRILINOS_DIR
#
#   make
#
# When finished building, running tests, etc., exist the new shell by typing
# 'exit[RETURN]' (unless pasing in --ci-mode).
#

# Assert this script is sourced, not run!
if [ "${BASH_SOURCE[0]}" == "${0}" ] ; then
  echo "This script '${0}' is being called.  Instead, it must be sourced!"
  exit 1
fi

# Get the base dir for the sourced script
SCRIPT_BASE_DIR=$(readlink -f \
  $(echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"))

# Get the build name
if [ "$1" == "" ] ; then
  echo "Error, the first argument must be the build name!"
  return
fi
FULL_BUILD_NAME=$1 ; shift

export TRILINOS_DIR=$(realpath $SCRIPT_BASE_DIR/../..)

export WORKSPACE=$(realpath $TRILINOS_DIR/..)  # Needed for some platforms

if [[ -e GenConfigSettings.cmake ]] ; then
  echo "Removing existing file GenConfigSettings.cmake ..."
  rm GenConfigSettings.cmake
fi

# Run new shell, set env, generate fragment file
source $TRILINOS_DIR/packages/framework/GenConfig/gen-config.sh \
--cmake-fragment GenConfigSettings.cmake \
$FULL_BUILD_NAME \
--force \
"$@" \
$TRILINOS_DIR

# NOTE: The above command will not complete until the user types 'exit' after
# doing the configure, build, testing, etc. (unless --ci-mode is passed
# through).

unset FULL_BUILD_NAME
unset WORKSPACE

# NOTE Above, it is an undocumented requirement that the env vars TRILINOS_DIR
# and WORKSPACE be set.
