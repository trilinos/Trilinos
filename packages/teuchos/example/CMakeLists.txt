
INCLUDE(Trilinos_Add_Executable_And_Test)

TRILINOS_ADD_EXECUTABLE(
  BLAS_example
  DIRECTORY BLAS
  COMM serial mpi  
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  CLP_example
  DIRECTORY CommandLineProcessor
  COMM serial mpi  
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  DenseMatrix_example
  DIRECTORY DenseMatrix
  COMM serial mpi  
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  SymDenseMatrix_example
  DIRECTORY DenseMatrix
  COMM serial mpi  
  SOURCES cxx_main_sym.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  hilbert
  DIRECTORY hilbert
  COMM serial mpi  
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  LAPACK_example
  DIRECTORY LAPACK
  COMM serial mpi  
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  ParameterList_example
  DIRECTORY ParameterList
  COMM serial mpi  
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  RefCountPtr_example
  DIRECTORY RefCountPtr
  COMM serial mpi  
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  TimeMonitor_example
  DIRECTORY TimeMonitor
  COMM serial mpi  
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE(
  XmlToParameterList
  DIRECTORY XML
  COMM serial mpi  
  SOURCES XmlToParameterList.cpp
  )

ADD_SUBDIRECTORY(ExplicitInstantiation/four_files)

# Copy all files over needed for testing

ADD_CUSTOM_TARGET(
  TeuchosCopyExampleFiles ALL
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/XML/_Example1.xml ${CMAKE_CURRENT_BINARY_DIR}/XML/Example1.xml
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/XML/_Example2.xml ${CMAKE_CURRENT_BINARY_DIR}/XML/Example2.xml
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/XML/_Example1.xml ${CMAKE_CURRENT_SOURCE_DIR}/XML/_Example2.xml
  )

# 2008/09/11: rabartl: NOTE: Above, as long as the user types 'make' to build
# the default target in this directory, all of the test input files will be
# copied over correctly.  In addition, if the file in the source directory
# changes, it will be copied over also.
