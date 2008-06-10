# $Header$ 

INCLUDE(Trilinos_Add_Executable)

INCLUDE_DIRECTORIES(${Teuchos_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${Teuchos_SOURCE_DIR}/src)

# 2008/06/09: rabartl: Above: You should not be listing include directories
# here!  The list of header files should be defined in one place and one place
# only and reused everywhere!  NO DUPLICATION!

LINK_LIBRARIES(teuchos ${TRILINOS_LAPACK_LIBRARY} ${TRILINOS_BLAS_LIBRARY})

# 2008/06/09: rabartl: Above: The link libraries required by Teuchos should
# not be defined here!  Insteady, they should be defined in one file and
# reused everwhere!  NO DUPLICATION!

IF(TRILINOS_ENABLE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	LINK_LIBRARIES(${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
ENDIF(TRILINOS_ENABLE_MPI)

TRILINOS_ADD_EXECUTABLE(CommandLineProcessor EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(DenseMatrix EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(LAPACK EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(ParameterList EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(RefCountPtr EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(TimeMonitor EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(XML EXAMPLE SOURCES XmlToParameterList.cpp)
TRILINOS_ADD_EXECUTABLE(hilbert EXAMPLE SOURCES cxx_main.cpp)

