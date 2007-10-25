ENABLE_LANGUAGE(Fortran)

OPTION(BUILD_SHARED_LIBS OFF)

SET(CMAKE_MODULE_PATH "${epetra_SOURCE_DIR}/cmake")
#FIND_PACKAGE(LAPACK REQUIRED)

FIND_LIBRARY(BLAS_LIBRARY blas DOC "Path to the BLAS implementation")
FIND_LIBRARY(LAPACK_LIBRARY lapack DOC "Path to the LAPACK implementation")

# Setup configuration options

IF(WIN32)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) name ## __")
ELSE(WIN32)
	SET(F77_FUNC "(name,NAME) name ## _")
	SET(F77_FUNC_ "(name,NAME) name ## _")
ENDIF(WIN32)

SET(HAVE_CASSERT TRUE)
SET(HAVE_CMATH TRUE)
SET(HAVE_CSTDIO TRUE)
SET(HAVE_CSTDLIB TRUE)
SET(HAVE_IOMANIP TRUE)
SET(HAVE_IOSTREAM TRUE)
SET(HAVE_STRING TRUE)

CONFIGURE_FILE(${epetra_SOURCE_DIR}/cmake/Epetra_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/Epetra_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_DEFINITIONS(-DHAVE_CONFIG_H)

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

# The following currently only works on Linux with changes to Epetra_BLAS_wrappers.h
# #ifdef __GNUC__ /* standard Epetra implementation */ 
# #define EPETRA_DCRSMV_F77   F77_FUNC_(epetra_dcrsmv,EPETRA_DCRSMV)
# #define EPETRA_DCRSMM_F77   F77_FUNC_(epetra_dcrsmm,EPETRA_DCRSMM)
# #define EPETRA_DCRSSV_F77   F77_FUNC_(epetra_dcrssv,EPETRA_DCRSSV)
# #define EPETRA_DCRSSM_F77   F77_FUNC_(epetra_dcrssm,EPETRA_DCRSSM)
#
# should be changed to
# #ifdef __GNUC__ /* standard Epetra implementation */ 
# #define EPETRA_DCRSMV_F77   F77_FUNC(epetra_dcrsmv,EPETRA_DCRSMV)
# #define EPETRA_DCRSMM_F77   F77_FUNC(epetra_dcrsmm,EPETRA_DCRSMM)
# #define EPETRA_DCRSSV_F77   F77_FUNC(epetra_dcrssv,EPETRA_DCRSSV)
# #define EPETRA_DCRSSM_F77   F77_FUNC(epetra_dcrssm,EPETRA_DCRSSM)
SET(FORTRAN_SOURCES
    Epetra_dcrsmm.f
    Epetra_dcrsmv.f
    Epetra_dcrssm.f
    Epetra_dcrssv.f
)

IF(ENABLE_MPI)
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

	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(ENABLE_MPI)

ADD_LIBRARY(epetra_fortran ${FORTRAN_SOURCES})

ADD_LIBRARY(epetra ${HEADERS} ${SOURCES})
TARGET_LINK_LIBRARIES(epetra epetra_fortran ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

IF(ENABLE_MPI)
	TARGET_LINK_LIBRARIES(epetra ${MPI_LIBRARY})
ENDIF(ENABLE_MPI)

