INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  ${PACKAGE_NAME}_ENABLE_Amesos
  )


IF (${PACKAGE_NAME}_ENABLE_Epetra AND
    ${PACKAGE_NAME}_ENABLE_EpetraExt AND
    ${PACKAGE_NAME}_ENABLE_Amesos)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Structure2D_Test #FIXME: not the same name as the cpp file. Confusing
  SOURCES Structure2D_epetra.cpp
  COMM serial mpi
  DEPLIBS muelu xpetra
)

IF(Trilinos_ENABLE_AztecOO)
	# build only if AztecOO available, too
	TRIBITS_ADD_EXECUTABLE_AND_TEST(
	  Structure2D_MLmode_Test
	  SOURCES Structure2D_MLmode_epetra.cpp
	  COMM serial mpi
	  DEPLIBS muelu xpetra
	)

	TRIBITS_ADD_EXECUTABLE_AND_TEST(
	  Structure2D_visualization
	  SOURCES Structure2D_visualization.cpp
	  COMM serial mpi
	  DEPLIBS muelu xpetra
	)

	TRIBITS_ADD_EXECUTABLE_AND_TEST(
	  Structure2D_strided
	  SOURCES Structure2D_strided.cpp
	  COMM serial mpi
	  DEPLIBS muelu xpetra
	)
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Structure2d_test_cp
SOURCE_FILES stru2d_A.txt stru2d_b.txt stru2d_ns.txt nodes0.txt MueLu_Agg2VTK.py
)

ENDIF()
