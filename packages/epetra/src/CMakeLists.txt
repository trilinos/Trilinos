# $Header$

# Package-specific configuration options
IF(${PROJECT_NAME}_ENABLE_FORTRAN)
  SET(HAVE_FORTRAN_SUPPORT TRUE)
ELSE(${PROJECT_NAME}_ENABLE_FORTRAN)
  SET(HAVE_FORTRAN_SUPPORT FALSE)
ENDIF(${PROJECT_NAME}_ENABLE_FORTRAN)

# Create <PACKAGE>_config.h
CONFIGURE_FILE(
  ${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}_config.h.in 
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h)
INSTALL(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h
  DESTINATION ${TRILINOS_INSTALL_LIB_INCLUDE_DIR}
  COMPONENT ${PROJECT_NAME}
  )

# Setup directories based on current directory
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Set headers
SET(HEADERS
	Epetra_BLAS.h
	Epetra_BLAS_wrappers.h
	Epetra_BasicDirectory.h
	Epetra_BasicRowMatrix.h
	Epetra_BlockMap.h
	Epetra_BlockMapData.h
	Epetra_C_wrappers.h
	Epetra_CombineMode.h
	Epetra_Comm.h
	Epetra_CompObject.h
	Epetra_ConfigDefs.h
	Epetra_CrsGraph.h
	Epetra_CrsGraphData.h
	Epetra_CrsMatrix.h
	Epetra_CrsSingletonFilter.h
	Epetra_Data.h
	Epetra_DataAccess.h
	Epetra_Directory.h
	Epetra_DistObject.h
	Epetra_Distributor.h
	Epetra_Export.h
	Epetra_FECrsGraph.h
	Epetra_FECrsMatrix.h
	Epetra_FEVbrMatrix.h
	Epetra_FEVector.h
#	Epetra_FastCrsMatrix.h
	Epetra_Flops.h
	Epetra_Fortran_wrappers.h
	Epetra_HashTable.h
	Epetra_Import.h
	Epetra_IntSerialDenseMatrix.h
	Epetra_IntSerialDenseVector.h
	Epetra_IntVector.h
	Epetra_InvOperator.h
	Epetra_JadMatrix.h
	Epetra_LAPACK.h
	Epetra_LAPACK_wrappers.h
	Epetra_LinearProblem.h
#	Epetra_LinearProblemRedistor.h
	Epetra_LocalMap.h
	Epetra_Map.h
	Epetra_MapColoring.h
	Epetra_MultiVector.h
	Epetra_Object.h
	Epetra_OffsetIndex.h
	Epetra_Operator.h
	Epetra_RowMatrix.h
	Epetra_RowMatrixTransposer.h
	Epetra_SerialComm.h
	Epetra_SerialCommData.h
	Epetra_SerialDenseMatrix.h
	Epetra_SerialDenseOperator.h
	Epetra_SerialDenseSVD.h
	Epetra_SerialDenseSolver.h
	Epetra_SerialDenseVector.h
	Epetra_SerialDistributor.h
	Epetra_SerialSpdDenseSolver.h
	Epetra_SerialSymDenseMatrix.h
	Epetra_SrcDistObject.h
	Epetra_Time.h
	Epetra_Util.h
	Epetra_VbrMatrix.h
	Epetra_VbrRowMatrix.h
	Epetra_Vector.h
	Epetra_Version.h
	)

# Set sources
SET(SOURCES
	Epetra_BLAS.cpp
	Epetra_BasicDirectory.cpp
	Epetra_BasicRowMatrix.cpp
	Epetra_BlockMap.cpp
	Epetra_BlockMapData.cpp
	Epetra_C_wrappers.cpp
	Epetra_CompObject.cpp
	Epetra_CrsGraph.cpp
	Epetra_CrsGraphData.cpp
	Epetra_CrsMatrix.cpp
	Epetra_CrsSingletonFilter.cpp
	Epetra_Data.cpp
	Epetra_DistObject.cpp
	Epetra_Export.cpp
	Epetra_FECrsGraph.cpp
	Epetra_FECrsMatrix.cpp
	Epetra_FEVbrMatrix.cpp
	Epetra_FEVector.cpp
#	Epetra_FastCrsMatrix.cpp
	Epetra_Flops.cpp
	Epetra_Fortran_wrappers.cpp
	Epetra_Import.cpp
	Epetra_IntSerialDenseMatrix.cpp
	Epetra_IntSerialDenseVector.cpp
	Epetra_IntVector.cpp
#	Epetra_InvOperator.cpp
	Epetra_JadMatrix.cpp
	Epetra_LAPACK.cpp
	Epetra_LinearProblem.cpp
#	Epetra_LinearProblemRedistor.cpp
	Epetra_LocalMap.cpp
	Epetra_Map.cpp
	Epetra_MapColoring.cpp
	Epetra_MultiVector.cpp
	Epetra_Object.cpp
	Epetra_OffsetIndex.cpp
	Epetra_RowMatrixTransposer.cpp
	Epetra_SerialComm.cpp
	Epetra_SerialCommData.cpp
	Epetra_SerialDenseMatrix.cpp
	Epetra_SerialDenseSVD.cpp
	Epetra_SerialDenseSolver.cpp
	Epetra_SerialDenseVector.cpp
	Epetra_SerialDistributor.cpp
	Epetra_SerialSpdDenseSolver.cpp
	Epetra_SerialSymDenseMatrix.cpp
	Epetra_Time.cpp
	Epetra_Util.cpp
	Epetra_VbrMatrix.cpp
	Epetra_Vector.cpp
	)

SET(FORTRAN_SOURCES
    Epetra_dcrsmm.f
    Epetra_dcrsmv.f
    Epetra_dcrssm.f
    Epetra_dcrssv.f
)

# Set addition headers and source based on package-specific 
# configuration options
IF(TRILINOS_ENABLE_MPI)
	SET(HEADERS
		${HEADERS}
		Epetra_MpiComm.h
		Epetra_MpiCommData.h
		Epetra_MpiDistributor.h
	#	Epetra_MpiSmpComm.h
	#	Epetra_MpiSmpCommData.h

	)

	SET(SOURCES
		${SOURCES}
		Epetra_MpiComm.cpp
		Epetra_MpiCommData.cpp
		Epetra_MpiDistributor.cpp
	#	Epetra_MpiSmpComm.cpp
	#	Epetra_MpiSmpCommData.cpp
	)

ENDIF(TRILINOS_ENABLE_MPI)

IF(${PROJECT_NAME}_ENABLE_TESTS)
  INCLUDE_DIRECTORIES(${${PROJECT_NAME}_SOURCE_DIR}/test/src)

  SET(HEADERS
    ${HEADERS}
    ${${PROJECT_NAME}_SOURCE_DIR}/test/src/Epetra_test_functions.h
    ${${PROJECT_NAME}_SOURCE_DIR}/test/src/Epetra_matrix_data.h
  )

  SET(SOURCES
    ${SOURCES}
    ${${PROJECT_NAME}_SOURCE_DIR}/test/src/Epetra_test_functions.cpp
    ${${PROJECT_NAME}_SOURCE_DIR}/test/src/Epetra_matrix_data.cpp
  )
ENDIF(${PROJECT_NAME}_ENABLE_TESTS)

# Target library setup
ADD_LIBRARY(epetra ${HEADERS} ${SOURCES})
LIST(APPEND TARGET_LIBS epetra)
LIST(APPEND LIB_DEPS ${TRILINOS_LAPACK_LIBRARY} ${TRILINOS_BLAS_LIBRARY})
TARGET_LINK_LIBRARIES(epetra ${LIB_DEPS})

# Auxilary library setup
IF(${PROJECT_NAME}_ENABLE_FORTRAN)  
  ADD_LIBRARY(epetra_fortran ${FORTRAN_SOURCES})
  LIST(APPEND TARGET_LIBS epetra_fortran) 
  TARGET_LINK_LIBRARIES(epetra epetra_fortran)
  INSTALL(TARGETS epetra_fortran
	  RUNTIME DESTINATION bin
	  LIBRARY DESTINATION lib
	  ARCHIVE DESTINATION lib
	  )
ENDIF(${PROJECT_NAME}_ENABLE_FORTRAN)

# MPI setup
IF(TRILINOS_ENABLE_MPI)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  LIST(APPEND LIB_DEPS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
  TARGET_LINK_LIBRARIES(epetra ${LIB_DEPS})
ENDIF(TRILINOS_ENABLE_MPI)

INSTALL(
  TARGETS epetra
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
  COMPONENT ${PROJECT_NAME}
	)

INSTALL(
  FILES ${HEADERS}
  DESTINATION ${TRILINOS_INSTALL_INCLUDE_DIR}
  COMPONENT ${PROJECT_NAME}
  )

# Define current package dependency variables
GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${INCLUDE_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} include directories")
GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT LINK_DIRECTORIES)
SET(${PROJECT_NAME}_LIBRARY_DIRS ${${PROJECT_NAME}_LIBRARY_DIRS}  ${LIBRARY_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} library directories")
SET(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES} ${TARGET_LIBS} ${LIB_DEPS} CACHE INTERNAL "${PROJECT_NAME} libraries")

