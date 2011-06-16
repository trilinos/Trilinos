
# For building CASL VRI add-on packages, you must explicitly list them
# here so that the native Trilinos packages will not have their tests
# and examples enabled and built.  That is important because these
# builds are sent to the CDash server on casl-dev which Trilinos
# developers do not have access to.  However, we want the libraries in
# upstream Trilinos packages to be enabled and built independently so
# that if there are errors in the up-stream packages, they will not
# show up in CASL-specific package CDash results.

#
#    source /opt/casldev/env/casl_dev_env.sh
#

SET(Trilinos_EXTRAREPOS_FILE "${TRILINOS_HOME_DIR}/cmake/ctest/drivers/pu241/ExtraExternalRepositories.casl_vri.cmake")
SET(VERA_COUPLED_BOA  OFF CACHE BOOL "")
SET(VERA_COUPLED_RAVE OFF CACHE BOOL "")
