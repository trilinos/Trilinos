INCLUDE(CombinedOption)
INCLUDE(TribitsETISupport)

#
# Declare the subpackage
#
TRIBITS_SUBPACKAGE(Core)

# Define ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION,
# which defaults to ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION.
TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

# Add this subpackage to the list of packages that do explicit
# template instantiation (ETI).
TRIBITS_ADD_ETI_SUPPORT()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
# FIXME (mfh 17 Dec 2014) We really just want a macro, not a CMake option.
# This is mainly just for backwards compatibility.
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_EXPLICIT_INSTANTIATION
  HAVE_TPETRA_EXPLICIT_INSTANTIATION
  "Enable explicit template instantiation (ETI) in Tpetra"
  ${${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION}
  )


SET(Tpetra_MACHINE_XML_FILE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/machine_files
    CACHE INTERNAL "")

# Set variables which were implicitly set when the Kokkos packages were optional dependencies
GLOBAL_SET(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT ON)
GLOBAL_SET(HAVE_TPETRACORE_TEUCHOSKOKKOSCOMM ON)
GLOBAL_SET(HAVE_TPETRACORE_KOKKOS ON)

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_SS_TESTING
  HAVE_TPETRA_ENABLE_SS_TESTING
  "Enabling testing of Secondary Stable Code."
  ${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE}
  )

# If this variable isn't defined, then the TpetraTSQR subpackage of
# Tpetra isn't where it should be.
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_TpetraTSQR)
IF (${PROJECT_NAME}_ENABLE_TpetraTSQR)
  SET(${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT ON)
ELSE ()
  SET(${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT OFF)
ENDIF ()

# TpetraCore_ENABLE_TSQR says whether TpetraCore should build its TSQR
# interface.  You may turn this interface off if you like, even if the
# TpetraTSQR subpackage of Tpetra is enabled.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_TSQR
  HAVE_TPETRA_TSQR
  "Whether to build Tpetra with TSQR support.  This is enabled by default if the KokkosTSQR subpackage is enabled.  You may turn off TSQR support if you like, but leaving it enabled is fine.  Please note that Epetra's TSQR adapter lives in the Tpetra package, for various historical reasons.  Therefore, in order to use TSQR with Epetra, you must enable Tpetra as well."
  ${${PACKAGE_NAME}_ENABLE_TSQR_DEFAULT}
  )

# mfh 20 Oct 2017: This option only exists for backwards
# compatibility.  It should always be true.
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_BUGTESTS
  HAVE_TPETRA_BUGTESTS
  "Enable testing of previous Tpetra bugs.  This option only exists for backwards compatibility, and no longer has an effect.  The tests in question always build and run, if Tpetra's tests are enabled."
  ON
  )

# ==================================================================
# Explicit template instantiation (ETI) and test instantiation logic
# ==================================================================

# tpetra/CMakeLists.txt (the package's CMake logic) defines these
# variables.  Despite "ETI" in their names, they exist whether or not
# ETI is defined.  If ETI is defined, these variables govern the set
# of template parameter combinations over which Tpetra instantiates.
# Whether or not ETI is defined, the variables govern the set of
# template parameter combinations over which Tpetra runs tests.

ASSERT_DEFINED (Tpetra_ETI_SCALARS)
ASSERT_DEFINED (Tpetra_ETI_SCALARS_NO_ORDS)
ASSERT_DEFINED (Tpetra_ETI_LORDS)
ASSERT_DEFINED (Tpetra_ETI_GORDS)
ASSERT_DEFINED (Tpetra_ETI_NODES)

SET(${PACKAGE_NAME}_ETI_SCALARS "${Tpetra_ETI_SCALARS}")
# Exclude all ordinal types (GlobalOrdinal and int).
SET(${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS "${Tpetra_ETI_SCALARS_NO_ORDS}")
SET(${PACKAGE_NAME}_ETI_LORDS "${Tpetra_ETI_LORDS}")
SET(${PACKAGE_NAME}_ETI_GORDS "${Tpetra_ETI_GORDS}")
SET(${PACKAGE_NAME}_ETI_NODES "${Tpetra_ETI_NODES}")

# "Export" the names for use in the ETI system.
# If we don't do this, ETI won't see these variables.

GLOBAL_SET(${PACKAGE_NAME}_ETI_SCALARS ${${PACKAGE_NAME}_ETI_SCALARS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS ${${PACKAGE_NAME}_ETI_SCALARS_NO_ORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_GORDS   ${${PACKAGE_NAME}_ETI_GORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_LORDS   ${${PACKAGE_NAME}_ETI_LORDS})
GLOBAL_SET(${PACKAGE_NAME}_ETI_NODES   ${${PACKAGE_NAME}_ETI_NODES})


# ==================================================================
# Other Tpetra options
# ==================================================================

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_Reduced_ETI
  HAVE_TPETRA_REDUCED_ETI
  "Uses reduced ETI of Tpetra to speed compilation."
  OFF
  )


TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Threaded_MKL
  HAVE_TPETRA_THREADED_MKL
  "Indicates that the linked BLAS is a threaded version of the MKL."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_THROW_Warnings
  HAVE_TPETRA_THROW_WARNINGS
  "Enable exception throwing for a number of warnings in Tpetra."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_PRINT_Warnings
  HAVE_TPETRA_PRINT_WARNINGS
  "Enable printing of a number of warnings in Tpetra."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_PRINT_Efficiency_Warnings
  HAVE_TPETRA_PRINT_EFFICIENCY_WARNINGS
  "Enable printing of Tpetra efficiency warnings."
  ${Tpetra_PRINT_Warnings}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_THROW_Abuse_Warnings
  HAVE_TPETRA_THROW_ABUSE_WARNINGS
  "Enable exception throwing for potential Tpetra abuse warnings."
  ${Tpetra_THROW_Warnings}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_PRINT_Abuse_Warnings
  HAVE_TPETRA_PRINT_ABUSE_WARNINGS
  "Enable printing of potential Tpetra abuse warnings."
  ${Tpetra_PRINT_Warnings}
  )

#
# FIXME (mfh 09 Jan 2015) This shouldn't even be an option.  We should
# examine the Map at run time to determine the best hash function to
# use.
#
TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_USE_MURMUR_HASH
  TPETRA_USE_MURMUR_HASH
  "Use the Murmur hash function in Tpetra::Map for global-to-local index lookups, rather than the default hash function.  Murmur hash is more expensive to evaluate, but does a better job of avoiding hash table collisions for Maps that do not look like local permutations of contiguous Maps."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_MMM_Timings
  HAVE_TPETRA_MMM_TIMINGS
  "Enable the timings for the MMM kernels"
  OFF
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_Distributor_Timings
  HAVE_TPETRA_DISTRIBUTOR_TIMINGS
  "Enable the timings for the Distributor"
  OFF
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_ENABLE_MMM_Statistics
  HAVE_TPETRA_MMM_STATISTICS
  "Enable statistics for the MMM kernels.  Warning: This involves extra communication and should only be enabled for diagnostic purposes."
  OFF
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  Tpetra_BCRS_Point_Import
  HAVE_TPETRA_BCRS_DO_POINT_IMPORT
  "Do MV import instead of BMV import."
  OFF
  )


#
# Add libraries, tests, and examples
#

ADD_SUBDIRECTORY(compat)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(inout)
ADD_SUBDIRECTORY(ext)

TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example guide)

#
# Do standard subpackage postprocessing
#
TRIBITS_SUBPACKAGE_POSTPROCESS()
