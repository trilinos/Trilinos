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

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

# Enable experimental code
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental
  HAVE_${PACKAGE_NAME_UC}_EXPERIMENTAL
  "Enable experimental code."
  NO)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Tpetra)
GLOBAL_SET(HAVE_XPETRA_TPETRA ${${PACKAGE_NAME}_ENABLE_Tpetra})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Epetra)
GLOBAL_SET(HAVE_XPETRA_EPETRA ${${PACKAGE_NAME}_ENABLE_Epetra})

GLOBAL_SET(HAVE_XPETRA_EPETRAEXT ${${PACKAGE_NAME}_ENABLE_EpetraExt})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Thyra)
GLOBAL_SET(HAVE_XPETRA_THYRA ${${PACKAGE_NAME}_ENABLE_Thyra})

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
    IF(Teuchos_ENABLE_LONG_LONG_INT AND Tpetra_INST_INT_LONG_LONG)
      SET(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT OFF)
    ELSE()
      SET(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT ON)
    ENDIF()
  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Tpetra AND NOT Kokkos_ENABLE_OpenMP AND NOT Kokkos_ENABLE_Serial)
    MESSAGE(FATAL_ERROR "If you use Xpetra with both the Epetra and Tpetra stack enabled you have to enable either/or the Serial or OpenMP stack for Kokkos (and Tpetra).")
  ENDIF()

ELSE()

  SET(${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES_DEFAULT ON)
  SET(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT ON)

ENDIF()

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
  MESSAGE(FATAL_ERROR "There are no 64 bit Epetra classes available possibly because they are disabled or Teuchos_ENABLE_LONG_LONG_INT is OFF or long long support is not present in Tpetra but user requested Xpetra use them by setting ${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES OFF.")
ENDIF()

#
# Decide whether it makes sense to enable LocalOrdinal = int and
# GlobalOrdinal = long long.  This requires support for long long in
# Teuchos, and at least one of Epetra or Tpetra must also support
# LocalOrdinal = int and GlobalOrdinal = long long.
#
ASSERT_DEFINED (Teuchos_ENABLE_LONG_LONG_INT)
IF (NOT Teuchos_ENABLE_LONG_LONG_INT)
  # If Teuchos doesn't have support for long long, we won't get anywhere.
  GLOBAL_SET (${PACKAGE_NAME}_INT_LONG_LONG OFF)
ELSE ()
  # Teuchos supports long long.  How about Epetra or Tpetra?
  IF (${PACKAGE_NAME}_ENABLE_Tpetra)
    # Does Tpetra provide LocalOrdinal = int and GlobalOrdinal = long long?
    ASSERT_DEFINED(Tpetra_INST_INT_LONG_LONG)
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
ENDIF ()

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  IF(NOT ${${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES})
    MESSAGE("-- Xpetra support for 32 bit Epetra enabled.")
    MESSAGE("--    ${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES=" ${${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES})
  ENDIF()
  IF(NOT ${${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES})
    MESSAGE("-- Xpetra support for 64 bit Epetra enabled.")
    MESSAGE("--    ${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES=" ${${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES})
  ENDIF()
ENDIF()

##################################

TRIBITS_ADD_OPTION_AND_DEFINE(XPETRA_ENABLE_PROFILING
  HAVE_XPETRA_PROFILING
  "Enable detailed timing of Xpetra."
  OFF)

# Define option to enable code that is refactored to use Kokkos
# This is ON by default when Kokkos is enabled.
IF (${PACKAGE_NAME}_ENABLE_KokkosCore)
  SET(Kokkos_Refactor_On_By_Default ON)
ELSE()
  SET(Kokkos_Refactor_On_By_Default OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor
  HAVE_${PACKAGE_NAME_UC}_KOKKOS_REFACTOR
  "Enable Xpetra code refactored to use Kokkos."
  ${Kokkos_Refactor_On_By_Default})

#
#
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(sup)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(research)
TRIBITS_ADD_TEST_DIRECTORIES(test)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  test/BlockedCrsOperator
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()

