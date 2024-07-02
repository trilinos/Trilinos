
# NOTE: ADOL-C has a bug in the installtion process.  It fails to
# install the file "config.h" into the install directory.  This is
# required to compile the library, so you must manually copy the file
# over into the installation directory.

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( ADOLC
  REQUIRED_HEADERS adolc/adolc.h adolc/config.h
  REQUIRED_LIBS_NAMES adolc
  )
