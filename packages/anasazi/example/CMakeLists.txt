# $Header$

INCLUDE(Trilinos_Add_Executable)

INCLUDE_DIRECTORIES(${Anasazi_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Anasazi_SOURCE_DIR}/src)
#INCLUDE_DIRECTORIES(${Anasazi_SOURCE_DIR}/util/ModeLaplace)
INCLUDE_DIRECTORIES(${Teuchos_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Teuchos_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${Epetra_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Epetra_SOURCE_DIR}/src)

# 2008/06/09: rabartl: Above: This is exactly the same list of include
# directories in ../test/CMakeLists.txt (see comment there).

LINK_LIBRARIES(anasazi teuchos epetra ${TRILINOS_LAPACK_LIBRARY} ${TRILINOS_BLAS_LIBRARY})

IF(TRILINOS_ENABLE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	LINK_LIBRARIES(${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
ENDIF(TRILINOS_ENABLE_MPI)

TRILINOS_ADD_EXECUTABLE(BlockDavidson NAME BlockDavidsonEpetraEx EXAMPLE SOURCES BlockDavidsonEpetraEx.cpp)
TRILINOS_ADD_EXECUTABLE(BlockKrylovSchur NAME BlockKrylovSchurEpetraEx EXAMPLE SOURCES BlockKrylovSchurEpetraEx.cpp)
TRILINOS_ADD_EXECUTABLE(BlockKrylovSchur NAME BlockKrylovSchurEpetraExSVD EXAMPLE SOURCES BlockKrylovSchurEpetraExSVD.cpp)
TRILINOS_ADD_EXECUTABLE(LOBPCG NAME LOBPCGEpetraEx EXAMPLE SOURCES LOBPCGEpetraEx.cpp)
TRILINOS_ADD_EXECUTABLE(LOBPCG NAME LOBPCGEpetraExSimple EXAMPLE SOURCES LOBPCGEpetraExSimple.cpp)
TRILINOS_ADD_EXECUTABLE(MVOPTester NAME EXAMPLE SOURCES MVOPTesterEx.cpp)

# Requires Trilinos packages not yet setup to use CMake
#TRILINOS_ADD_EXECUTABLE(BlockDavidson NAME BlockDavidsonEpetraExFile EXAMPLE SOURCES BlockDavidsonEpetraExFile.cpp)
#TRILINOS_ADD_EXECUTABLE(BlockKrylovSchur NAME BlockKrylovSchurEpetraExFile EXAMPLE SOURCES BlockKrylovSchurEpetraExFile.cpp)
#TRILINOS_ADD_EXECUTABLE(LOBPCG NAME LOBPCGEpetraExFile EXAMPLE SOURCES LOBPCGEpetraExFile.cpp)
#TRILINOS_ADD_EXECUTABLE(BlockKrylovSchur NAME BlockKrylovSchurEpetraExGenAmesos EXAMPLE SOURCES BlockKrylovSchurEpetraExGenAmesos.cpp)
#TRILINOS_ADD_EXECUTABLE(BlockKrylovSchur NAME BlockKrylovSchurEpetraExGenAztecOO EXAMPLE SOURCES BlockKrylovSchurEpetraExGenAztecOO.cpp)
#TRILINOS_ADD_EXECUTABLE(BlockKrylovSchur NAME BlockKrylovSchurEpetraExGenBelos EXAMPLE SOURCES BlockKrylovSchurEpetraExGenBelos.cpp)
#TRILINOS_ADD_EXECUTABLE(BlockDavidson NAME BlockDavidsonEpetraExGen EXAMPLE SOURCES BlockDavidsonEpetraExGen.cpp)
#TRILINOS_ADD_EXECUTABLE(BlockDavidson NAME BlockDavidsonEpetraExGenPrecIfpack EXAMPLE SOURCES BlockDavidsonEpetraExGenPrecIfpack.cpp)
#TRILINOS_ADD_EXECUTABLE(LOBPCG NAME LOBPCGEpetraExGen EXAMPLE SOURCES LOBPCGEpetraExGen.cpp)
#TRILINOS_ADD_EXECUTABLE(LOBPCG NAME LOBPCGEpetraExGenPrecIfpack EXAMPLE SOURCES LOBPCGEpetraExGenPrecIfpack.cpp)
