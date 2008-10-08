
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  BLAS_example
  DIRECTORY BLAS
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  CLP_example
  DIRECTORY CommandLineProcessor
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  DenseMatrix_example
  DIRECTORY DenseMatrix
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  SymDenseMatrix_example
  DIRECTORY DenseMatrix
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main_sym.cpp
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  hilbert
  DIRECTORY hilbert
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  LAPACK_example
  DIRECTORY LAPACK
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  ParameterList_example
  DIRECTORY ParameterList
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  RefCountPtr_example
  DIRECTORY RefCountPtr
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  TimeMonitor_example
  DIRECTORY TimeMonitor
  COMM serial mpi  
  NUM_MPI_PROCS 1
  SOURCES cxx_main.cpp
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  XmlToParameterList
  DIRECTORY XML
  COMM serial mpi
  ARGS
    "--xml-in-file=XML/Example1.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  SOURCES XmlToParameterList.cpp
  )

COPY_FILES_TO_BINARY_DIR(TeuchosCopyExampleXmlFiles
  DEST_FILES Example1.xml Example2.xml
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/XML
  SOURCE_PREFIX "_"
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/XML
  )

ADD_SUBDIRECTORY(ExplicitInstantiation/four_files)
