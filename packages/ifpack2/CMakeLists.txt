#
# Define the package
#
TRIBITS_PACKAGE( Ifpack2 )

#
# Set up package-specific options
#

#
# Do explicit template instantiation (ETI) and testing over the same
# Scalar, LocalOrdinal, and GlobalOrdinal types as for Tpetra objects.
# However, by default, exclude all of Tpetra's GlobalOrdinal (GO)
# types from the ETI list.  This is because Tpetra adds all GO ETI
# types to the list of Scalar ETI types.  (This lets it implement
# certain communication patterns involving e.g., Vector<GO, LO, GO>.)
#

# Make sure that Tpetra actually defined these variables, even if they
# are empty.
ASSERT_DEFINED(TpetraCore_ETI_SCALARS)
ASSERT_DEFINED(TpetraCore_ETI_LORDS)
ASSERT_DEFINED(TpetraCore_ETI_GORDS)
ASSERT_DEFINED(TpetraCore_ETI_NODES)

# Exclude GlobalOrdinal types from the list of Scalar types over which
# Ifpack2 does ETI.
#
# Look at each entry in Tpetra's list of ETI Scalar types.  If it is
# also a Tpetra ETI GlobalOrdinal type, don't append it to Ifpack2's
# list of ETI Scalar types.  Otherwise, do.
SET (Ifpack2_ETI_SCALARS_TEMP ) # start with the empty list
FOREACH (scalarName IN LISTS TpetraCore_ETI_SCALARS)
  LIST (FIND TpetraCore_ETI_GORDS ${scalarName} scalarNameIndex)
  IF (-1 EQUAL scalarNameIndex) # not in the list of integer types
    LIST (APPEND Ifpack2_ETI_SCALARS_TEMP ${scalarName})
  ENDIF ()
ENDFOREACH ()

SET(Ifpack2_ETI_SCALARS ${Ifpack2_ETI_SCALARS_TEMP})
SET(Ifpack2_ETI_LORDS   ${TpetraCore_ETI_LORDS})
SET(Ifpack2_ETI_GORDS   ${TpetraCore_ETI_GORDS})
SET(Ifpack2_ETI_NODES   ${TpetraCore_ETI_NODES})

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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Lemon)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos2)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Cholmod)
IF(${PACKAGE_NAME}_ENABLE_Lemon AND ${PACKAGE_NAME}_ENABLE_Amesos2 AND ${PACKAGE_NAME}_ENABLE_Cholmod)
  SET(HAVE_IFPACK2_SUPPORTGRAPH ON )
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

# Promote these variables to be visible outside of the "scope" of this
# directory -- e.g., to other packages.
GLOBAL_SET(Ifpack2_ETI_SCALARS ${Ifpack2_ETI_SCALARS})
GLOBAL_SET(Ifpack2_ETI_LORDS   ${Ifpack2_ETI_LORDS})
GLOBAL_SET(Ifpack2_ETI_GORDS   ${Ifpack2_ETI_GORDS})
GLOBAL_SET(Ifpack2_ETI_NODES   ${Ifpack2_ETI_NODES})

#
# Exclude files for source package.
#


#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
