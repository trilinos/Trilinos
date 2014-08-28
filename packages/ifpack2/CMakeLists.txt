
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE( Ifpack2 )

#
# B) Set up package-specific options
#

ASSERT_DEFINED(Tpetra_ETI_SCALARS)
SET(Ifpack2_ETI_SCALARS "double")
SET(Ifpack2_ETI_GORDS   "int")
SET(Ifpack2_ETI_LORDS   "int")

IF(Tpetra_INST_INT_LONG_LONG) 
  LIST(APPEND Ifpack2_ETI_GORDS   "long long")
ENDIF()


TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_Trilinos
  HAVE_${PACKAGE_NAME}_Trilinos
  "Building as a Trilinos package"
  ON
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_${PACKAGE_NAME}_DEBUG
  "Enable runtime debug checking."
  ${${PROJECT_NAME}_ENABLE_DEBUG}
  )

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

# we use this for testing as well as library eti, so do it regardless
TRIBITS_ADD_ETI_SUPPORT()

#Adding a define for when the TPL QD is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_${PACKAGE_NAME}_QD "${PACKAGE_NAME}_ENABLE_QD")

#Adding a define for when MPI is enabled.
#This replaces an option that was reusing an existing name.
SET(HAVE_${PACKAGE_NAME}_MPI "${${PACKAGE_NAME}_ENABLE_MPI}")

IF(${PACKAGE_NAME}_ENABLE_Lemon)
IF(${PACKAGE_NAME}_ENABLE_Amesos2)
IF(${PACKAGE_NAME}_ENABLE_Cholmod)
  SET(HAVE_IFPACK2_SUPPORTGRAPH ON )
ENDIF()
ENDIF()
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_Experimental	
  HAVE_IFPACK2_EXPERIMENTAL
  "Enable building and installation of experimental Ifpack2 features."
  NO )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_IFPACK2_TIMER_BARRIER
  HAVE_IFPACK2_TIMER_BARRIER
  "Enables a barrier in Ifpack2 Additive Schwarz so that timing information is useful"
  NO )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(adapters)

TRIBITS_ADD_TEST_DIRECTORIES(test)

GLOBAL_SET(Ifpack2_ETI_SCALARS ${Ifpack2_ETI_SCALARS})
GLOBAL_SET(Ifpack2_ETI_GORDS   ${Ifpack2_ETI_GORDS})
GLOBAL_SET(Ifpack2_ETI_LORDS   ${Ifpack2_ETI_LORDS})

#
# Exclude files for source package.
#


#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
