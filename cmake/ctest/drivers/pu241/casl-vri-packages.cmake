
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

SET(Trilinos_EXTRAREPOS_FILE "${CTEST_SCRIPT_DIRECTORY}/ExtraExternalRepositories.casl_vri.cmake")
SET(Trilinos_PACKAGES CASLBOA LIME VRIPSS)
# CASLRAVE 
SET(EXTRA_CONFIGURE_OPTIONS
  "-DTPL_ENABLE_JDK:BOOL=ON"
  "-DJDK_LIBRARY_DIRS='/usr/lib/jvm/java/jre/lib/amd64/server'"
  "-DJDK_INCLUDE_DIRS='/usr/lib/jvm/java/include'"
  "-DTPL_ENABLE_OpenSSL:BOOL=ON"
  "-DTPL_ENABLE_LIBXML2:BOOL=ON"
  "-DTPL_ENABLE_HDF5CPP:BOOL=ON"
  "-DTPL_ENABLE_Zlib:BOOL=ON"
  "-DTPL_ENABLE_TCL:BOOL=ON"
  "-DTCL_LIBRARY_NAMES='tcl8.5'"
  )

# I can't get PVM to configure with these
#  "-DTPL_ENABLE_PVM:BOOL=ON"
#  "-DPVM_LIBRARY_DIRS:PATH='/opt/intel-11.1.064/tpls/pvm3/lib'"
#  "-DPVM_INCLUDE_DIRS:PATH='/opt/intel-11.1.064/tpls/pvm3/include'"

