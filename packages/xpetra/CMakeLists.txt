INCLUDE(CombinedOption)
INCLUDE(DualScopeSet)

#
# A) Define the package
#

TRIBITS_PACKAGE(Xpetra ENABLE_SHADOWING_WARNINGS)

# Add directory xpetra/cmake/ to the cmake module search path
DUAL_SCOPE_SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PACKAGE_SOURCE_DIR}/cmake)

#
# B) Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

# ETI requires that Tpetra is enabled.
IF( DEFINED ${PROJECT_NAME}_ENABLE_Tpetra AND ${${PACKAGE_NAME}_ENABLE_Tpetra} )
  TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
ELSE()
  MESSAGE( STATUS "Xpetra: Skipping ETI check because Tpetra is not enabled." )
  SET( HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION OFF )
ENDIF()

# Deprecation stuff
SET (${PACKAGE_NAME}_ENABLE_DEPRECATED_CODE_DEFAULT ON)
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEPRECATED_CODE
  ${PACKAGE_NAME_UC}_ENABLE_DEPRECATED_CODE
  "Whether Xpetra enables deprecated code (that is, anything marked with the XPETRA_DEPRECATED macro) at compile time.  Default is ON (deprecated code enabled).  If OFF, then deprecated code does not exist."
  ${${PACKAGE_NAME}_ENABLE_DEPRECATED_CODE_DEFAULT}
)

ASSERT_DEFINED (${PACKAGE_NAME}_ENABLE_DEPRECATED_CODE)
IF (${PACKAGE_NAME}_ENABLE_DEPRECATED_CODE)
  MESSAGE(STATUS "Xpetra: Enabling deprecated code")
ELSE ()
  MESSAGE(STATUS "Xpetra: Disabling deprecated code")
ENDIF ()

# Enable experimental code
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_${PACKAGE_NAME_UC}_EXPERIMENTAL
  "Enable experimental code."
  NO)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Thyra)

GLOBAL_SET(HAVE_XPETRA_THYRA ${${PACKAGE_NAME}_ENABLE_Thyra})

# Does Tpetra provide LocalOrdinal = int and GlobalOrdinal = long long?
ASSERT_DEFINED(Tpetra_INST_INT_LONG_LONG)
ASSERT_DEFINED(Tpetra_INST_INT_INT)
GLOBAL_SET (${PACKAGE_NAME}_INT_LONG_LONG ${Tpetra_INST_INT_LONG_LONG})

##################################

TRIBITS_ADD_OPTION_AND_DEFINE(XPETRA_ENABLE_PROFILING
  HAVE_XPETRA_PROFILING
  "Enable detailed timing of Xpetra."
  OFF)


IF (DEFINED ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor)
  MESSAGE(WARNING "Xpetra no longer uses the option \'Xpetra_ENABLE_Kokkos_Refactor\' and will always build with Kokkos code enabled. \'Xpetra_ENABLE_Kokkos_Refactor\' can safely be removed from the CMake invocation.")
ENDIF()

#
#
#

ADD_SUBDIRECTORY(src)

IF (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/CMakeLists.txt)
  ADD_SUBDIRECTORY(doc)
ENDIF()

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(research)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
