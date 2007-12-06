INCLUDE(Trilinos_Add_Executable)

INCLUDE_DIRECTORIES(${Anasazi_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Anasazi_SOURCE_DIR}/src)
#INCLUDE_DIRECTORIES(${Anasazi_SOURCE_DIR}/util/ModeLaplace)
INCLUDE_DIRECTORIES(${Teuchos_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Teuchos_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${Epetra_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Epetra_SOURCE_DIR}/src)

LINK_LIBRARIES(anasazi teuchos epetra ${TRILINOS_LAPACK_LIBRARY} ${TRILINOS_BLAS_LIBRARY})

TRILINOS_ADD_EXECUTABLE(BlockDavidson NAME BlockDavidsonEpetra EXAMPLE SOURCES BlockDavidsonEpetraEx.cpp)
#TRILINOS_ADD_EXECUTABLE(BlockDavidson NAME BlockDavidsonEpetraGen EXAMPLE SOURCES BlockDavidsonEpetraExGen.cpp)
#TRILINOS_ADD_EXECUTABLE(BlockDavidson NAME BlockDavidsonEpetraFile EXAMPLE SOURCES BlockDavidsonEpetraExFile.cpp)
#TRILINOS_ADD_EXECUTABLE(BlockDavidson NAME BlockDavidsonEpetraGenPrecIfpack EXAMPLE SOURCES BlockDavidsonEpetraExGenPrecIfpack.cpp)

