
TRIBITS_SUBPACKAGE(TpetraAdapters)

ASSERT_DEFINED(
  Thyra_ENABLE_EXPLICIT_INSTANTIATION
  Tpetra_INST_FLOAT
  Tpetra_INST_COMPLEX_DOUBLE
  Tpetra_INST_COMPLEX_FLOAT
  Teuchos_INST_FLOAT
  Teuchos_INST_COMPLEX_DOUBLE
  Teuchos_INST_COMPLEX_FLOAT
)

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

TRIBITS_ADD_ETI_SUPPORT()

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
ASSERT_DEFINED(TpetraCore_ETI_SCALARS_NO_ORDS)
ASSERT_DEFINED(TpetraCore_ETI_LORDS)
ASSERT_DEFINED(TpetraCore_ETI_GORDS)
ASSERT_DEFINED(TpetraCore_ETI_NODES)

# Promote these variables to be visible outside of the "scope" of this
# directory -- e.g., to other packages.
GLOBAL_SET(ThyraTpetraAdapters_ETI_SCALARS ${TpetraCore_ETI_SCALARS_NO_ORDS})
GLOBAL_SET(ThyraTpetraAdapters_ETI_LORDS   ${TpetraCore_ETI_LORDS})
GLOBAL_SET(ThyraTpetraAdapters_ETI_GORDS   ${TpetraCore_ETI_GORDS})
GLOBAL_SET(ThyraTpetraAdapters_ETI_NODES   ${TpetraCore_ETI_NODES})

IF (Thyra_ENABLE_EXPLICIT_INSTANTIATION)

  IF ((Tpetra_INST_FLOAT AND NOT Teuchos_INST_FLOAT) OR (NOT Tpetra_INST_FLOAT AND Teuchos_INST_FLOAT))
    MESSAGE(FATAL_ERROR "Tpetra_INST_FLOAT=${Tpetra_INST_FLOAT} != Teuchos_INST_FLOAT=${Teuchos_INST_FLOAT}. The two need to match. Consider setting Trilinos_ENABLE_FLOAT instead of setting Tpetra_INST_FLOAT and/or Teuchos_INST_FLOAT individually.")
  ENDIF()

  IF ((Tpetra_INST_COMPLEX_FLOAT AND NOT Teuchos_INST_COMPLEX_FLOAT) OR (NOT Tpetra_INST_COMPLEX_FLOAT AND Teuchos_INST_COMPLEX_FLOAT))
    MESSAGE(FATAL_ERROR "Tpetra_INST_COMPLEX_FLOAT=${Tpetra_INST_COMPLEX_FLOAT} != Teuchos_INST_COMPLEX_FLOAT=${Teuchos_INST_COMPLEX_FLOAT}. The two need to match. Consider setting Trilinos_ENABLE_COMPLEX_FLOAT instead of setting Tpetra_INST_COMPLEX_FLOAT and/or Teuchos_INST_COMPLEX_FLOAT individually.")
  ENDIF()

  IF ((Tpetra_INST_COMPLEX_DOUBLE AND NOT Teuchos_INST_COMPLEX_DOUBLE) OR (NOT Tpetra_INST_COMPLEX_DOUBLE AND Teuchos_INST_COMPLEX_DOUBLE))
    MESSAGE(FATAL_ERROR "Tpetra_INST_COMPLEX_DOUBLE=${Tpetra_INST_COMPLEX_DOUBLE} != Teuchos_INST_COMPLEX_DOUBLE=${Teuchos_INST_COMPLEX_DOUBLE}. The two need to match. Consider setting Trilinos_ENABLE_COMPLEX_DOUBLE instead of setting Tpetra_INST_COMPLEX_DOUBLE and/or Teuchos_INST_COMPLEX_DOUBLE individually.")
  ENDIF()

ENDIF()

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_TEST_DIRECTORIES(example)

TRIBITS_SUBPACKAGE_POSTPROCESS()
