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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Epetra)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Thyra)

IF(NOT ${PACKAGE_NAME}_ENABLE_Epetra AND NOT ${PACKAGE_NAME}_ENABLE_Tpetra)
  MESSAGE(FATAL_ERROR "Neither Tpetra nor Epetra support is enabled in Xpetra."
    " Set Xpetra_ENABLE_Tpetra and/or Xpetra_ENABLE_Epetra to ON,"
    " or just disable Xpetra.")
ENDIF()

# If Epetra is enabled, does it provide 64-bit global indices?
IF(${PACKAGE_NAME}_ENABLE_Epetra)

  # Set default value of whether Xpetra should allow for 32 bit Epetra data structures or not.
  IF(Trilinos_NO_32BIT_GLOBAL_INDICES)
    SET(${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES_DEFAULT ON)
  ELSE()
    SET(${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES_DEFAULT OFF)
  ENDIF()

  # Set default value of whether Xpetra should allow for 64 bit Epetra data structures or not.
  IF(Trilinos_NO_64BIT_GLOBAL_INDICES)
    SET(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT ON)
  ELSE()
    IF(Tpetra_INST_INT_LONG_LONG)
      SET(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT OFF)
    ELSE()
      SET(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT ON)
    ENDIF()
  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Tpetra AND NOT Kokkos_ENABLE_OPENMP AND NOT Kokkos_ENABLE_SERIAL)
    MESSAGE(FATAL_ERROR "If you use Xpetra with both the Epetra and Tpetra stack enabled you have to enable either/or the Serial or OpenMP stack for Kokkos (and Tpetra).")
  ENDIF()

ELSE()

  SET(${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES_DEFAULT ON)
  SET(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT ON)

  # Make sure EpetraExt is turned off
  DUAL_SCOPE_SET(${PACKAGE_NAME}_ENABLE_EpetraExt OFF)
ENDIF()

#
# Decide whether it makes sense to enable LocalOrdinal = int and
# GlobalOrdinal = long long.  This requires support for long long in
# Teuchos, and at least one of Epetra or Tpetra must also support
# LocalOrdinal = int and GlobalOrdinal = long long.
#

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES
  ${PACKAGE_NAME_UC}_EPETRA_NO_32BIT_GLOBAL_INDICES
  "If option is enabled, 32-bit Epetra global indices are disabled in Xpetra.  One must use GlobalOrdinal = long long only."
  ${${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES_DEFAULT} )

IF(${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES_DEFAULT AND (NOT ${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES))
  MESSAGE(FATAL_ERROR "There are no 32 bit Epetra classes available but user requested Xpetra use them by setting ${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES OFF.")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES
  ${PACKAGE_NAME_UC}_EPETRA_NO_64BIT_GLOBAL_INDICES
  "If option is enabled, 64-bit Epetra global indices are disabled in Xpetra.  One must use GlobalOrdinal = int only."
  ${${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT} )

IF(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT AND (NOT ${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES))
  MESSAGE(FATAL_ERROR "There are no 64 bit Epetra classes available possibly because they are disabled or long long support is not present in Tpetra but user requested Xpetra use them by setting ${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES OFF.")
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  IF(NOT ${${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES})
    MESSAGE("-- Xpetra support for 32 bit Epetra enabled.")
    MESSAGE("--    ${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES=" ${${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES})

    # Tpetra and Epetra need to match GO types
    IF(${PACKAGE_NAME}_ENABLE_Tpetra)
      IF(NOT Tpetra_INST_INT_INT)
        MESSAGE(FATAL_ERROR "If Xpetra Epetra support is enabled and Epetra 32-bit global indices are enabled, Xpetra requires that you enable Tpetra_INST_INT_INT")
      ENDIF()
    ENDIF()
  ENDIF()
  IF(NOT ${${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES})
    MESSAGE("-- Xpetra support for 64 bit Epetra enabled.")
    MESSAGE("--    ${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES=" ${${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES})

    # Tpetra and Epetra need to match GO types
    IF(${PACKAGE_NAME}_ENABLE_Tpetra)
      IF(NOT Tpetra_INST_INT_LONG_LONG)
        MESSAGE(FATAL_ERROR "If Xpetra Epetra support is enabled and Epetra 64-bit global indices are enabled, Xpetra requires that you enable Tpetra_INST_INT_LONG_LONG")
      ENDIF()
    ENDIF()

  ENDIF()
ENDIF()

GLOBAL_SET(HAVE_XPETRA_EPETRA ${${PACKAGE_NAME}_ENABLE_Epetra})
GLOBAL_SET(HAVE_XPETRA_EPETRAEXT ${${PACKAGE_NAME}_ENABLE_EpetraExt})
GLOBAL_SET(HAVE_XPETRA_THYRA ${${PACKAGE_NAME}_ENABLE_Thyra})

# Check whether Epetra or Tpetra support long long.
IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  # Does Tpetra provide LocalOrdinal = int and GlobalOrdinal = long long?
  ASSERT_DEFINED(Tpetra_INST_INT_LONG_LONG)
  ASSERT_DEFINED(Tpetra_INST_INT_INT)
  GLOBAL_SET (${PACKAGE_NAME}_INT_LONG_LONG ${Tpetra_INST_INT_LONG_LONG})
ELSE () # Tpetra NOT enabled
  IF (${PACKAGE_NAME}_ENABLE_Epetra)
    # Does Epetra provide LocalOrdinal = int and GlobalOrdinal = long long?
    IF (${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES)
      GLOBAL_SET (${PACKAGE_NAME}_INT_LONG_LONG OFF)
    ELSE ()
      GLOBAL_SET (${PACKAGE_NAME}_INT_LONG_LONG ON)
    ENDIF ()
   ELSE () # Neither Epetra nor Tpetra are enabled
    GLOBAL_SET (${PACKAGE_NAME}_INT_LONG_LONG OFF)
  ENDIF ()
ENDIF ()

##################################

TRIBITS_ADD_OPTION_AND_DEFINE(XPETRA_ENABLE_PROFILING
  HAVE_XPETRA_PROFILING
  "Enable detailed timing of Xpetra."
  OFF)

#
#
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(sup)

IF (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/CMakeLists.txt)
  ADD_SUBDIRECTORY(doc)
ENDIF()

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(research)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
