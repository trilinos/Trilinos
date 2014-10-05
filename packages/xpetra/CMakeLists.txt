INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
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

GLOBAL_SET(HAVE_XPETRA_TPETRA ${${PACKAGE_NAME}_ENABLE_Tpetra})

GLOBAL_SET(HAVE_XPETRA_EPETRA ${${PACKAGE_NAME}_ENABLE_Epetra})

GLOBAL_SET(HAVE_XPETRA_EPETRAEXT ${${PACKAGE_NAME}_ENABLE_EpetraExt})

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

ELSE()

  SET(${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES_DEFAULT ON)
  SET(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT ON)

ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES
  ${PACKAGE_NAME_UC}_EPETRA_NO_32BIT_GLOBAL_INDICES
  "If option is enabled, 32 bit Epetra global indices are disabled in Xpetra.  One must use long long global indices only."
  ${${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES_DEFAULT} )

IF(${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES_DEFAULT AND (NOT ${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES))
  MESSAGE(FATAL_ERROR "There are no 32 bit Epetra classes available but user requested Xpetra use them by setting ${PACKAGE_NAME}_Epetra_NO_32BIT_GLOBAL_INDICES OFF.")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES
  ${PACKAGE_NAME_UC}_EPETRA_NO_64BIT_GLOBAL_INDICES
  "If option is enabled, 64 bit Epetra global indices are disabled in Xpetra.  One must use int global indices only."
  ${${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT} )

IF(${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES_DEFAULT AND (NOT ${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES))
  MESSAGE(FATAL_ERROR "There are no 64 bit Epetra classes available possibly because they are disabled or Teuchos_ENABLE_LONG_LONG_INT is OFF or long long support is not present in Tpetra but user requested Xpetra use them by setting ${PACKAGE_NAME}_Epetra_NO_64BIT_GLOBAL_INDICES OFF.")
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(XPETRA_ENABLE_PROFILING
  HAVE_XPETRA_PROFILING
  "Enable detailed timing of Xpetra."
  OFF)

#
#
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(ext)
ADD_SUBDIRECTORY(sup)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
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

