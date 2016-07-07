########################################################################
# See associated tribits/Copyright.txt file for copyright and license! #
########################################################################

#
# First, set up the variables for the (backward-compatible) TriBITS way of
# finding HDF5.  These are used in case FIND_PACKAGE(HDF5 ...) is not called
# or does not find HDF5.  Also, these variables need to be non-null in order
# to trigger the right behavior in the function
# TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES().
#

SET(REQUIRED_HEADERS hdf5.h)
SET(REQUIRED_LIBS_NAMES hdf5)

IF (HDF5_REQUIRE_FORTRAN)
  SET(REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} hdf5_fortran)
ENDIF()

IF (TPL_ENABLE_MPI)
  SET(REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} z)
ENDIF()

IF (TPL_ENABLE_Netcdf)
  SET(REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES} hdf5_hl)
ENDIF()

#
# Second, search for HDF5 components (if allowed) using the standard
# FIND_PACKAGE(HDF5 ...).
#
TRIBITS_TPL_ALLOW_PRE_FIND_PACKAGE(HDF5  HDF5_ALLOW_PREFIND)
IF (HDF5_ALLOW_PREFIND)

  MESSAGE("-- Using FIND_PACKAGE(HDF5 ...) ...")

  SET(HDF5_COMPNENTS C)
  IF (HDF5_REQUIRE_FORTRAN)
    LIST(APPEND HDF5_COMPNENTS Fortran)
  ENDIF()

  FIND_PACKAGE(HDF5 COMPONENTS ${HDF5_COMPNENTS})

  # Make sure that HDF5 is parallel.
  IF (TPL_ENABLE_MPI AND NOT HDF5_IS_PARALLEL)
    MESSAGE(FATAL_ERROR "Trilinos is configured for MPI, HDF5 is not.
    Did CMake find the correct libraries?
    Try setting HDF5_INCLUDE_DIRS and/or HDF5_LIBRARY_DIRS explicitly.
    ")
  ENDIF()

  IF (HDF5_FOUND)
    # Tell TriBITS that we found HDF5 and there no need to look any further!
    SET(TPL_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} CACHE PATH
      "HDF5 include dirs")
    SET(TPL_HDF5_LIBRARIES ${HDF5_LIBRARIES} CACHE FILEPATH
      "HDF5 libraries")
    SET(TPL_HDF5_LIBRARY_DIRS ${HDF5_LIBRARY_DIRS} CACHE PATH
      "HDF5 library dirs")
  ENDIF()

ENDIF()

#
# Third, call TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES()
#
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( HDF5
  REQUIRED_HEADERS ${REQUIRED_HEADERS}
  REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )
# NOTE: If FIND_PACKAGE(HDF5 ...) was called and successfully found HDF5, then
# TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES() will use the already-set
# variables TPL_HDF5_INCLUDE_DIRS and TPL_HDF5_LIBRARIES and then print them
# out (and set some other standard variables as well).  This is the final
# "hook" into the TriBITS TPL system.
