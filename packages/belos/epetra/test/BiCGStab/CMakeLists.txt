
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Epetra)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

# Dependency on Triutils fixes Bug 6449.
IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BiCGStab_file
    SOURCES BiCGStabExFile.cpp
    COMM serial mpi
    ARGS "--verbose"
         "--verbose --num-rhs=2"
    STANDARD_PASS_OUTPUT
    )

IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BiCGStab_prec_file
    SOURCES PrecBiCGStabExFile.cpp
    COMM serial mpi
    ARGS "--verbose"
    STANDARD_PASS_OUTPUT
    )

ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBiCGSTABFiles
    SOURCE_DIR ${${PACKAGE_NAME}_SOURCE_DIR}/epetra/test/BiCGStab
    SOURCE_FILES cage4.hb
    EXEDEPS BiCGStab_file
    )

ENDIF()
