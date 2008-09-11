
# Required modules
#INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Trilinos_Add_Executable)
INCLUDE(Trilinos_Add_Test)

# Setup include directories and library dependencies
#INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDE_DIRS})
#LINK_LIBRARIES(${${PROJECT_NAME}_LIBRARIES})

# Example definitions
#TRILINOS_ADD_EXECUTABLE_AND_TEST(CommandLineProcessor EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  CommandLineProcessor  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY CommandLineProcessor
  )
#TRILINOS_ADD_EXECUTABLE_AND_TEST(DenseMatrix EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  DenseMatrix  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY DenseMatrix
  )
#TRILINOS_ADD_EXECUTABLE_AND_TEST(LAPACK EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  LAPACK  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY LAPACK
  )
#TRILINOS_ADD_EXECUTABLE_AND_TEST(ParameterList EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  ParameterList  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY ParameterList
  )
#TRILINOS_ADD_EXECUTABLE_AND_TEST(RefCountPtr EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  RefCountPtr  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY RefCountPtr
  )
#TRILINOS_ADD_EXECUTABLE_AND_TEST(TimeMonitor EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  TimeMonitor  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY TimeMonitor
  )
#TRILINOS_ADD_EXECUTABLE_AND_TEST(XML EXAMPLE SOURCES XmlToParameterList.cpp)
TRILINOS_ADD_EXECUTABLE(
  XML  
  COMM serial mpi  
  SOURCES XmlToParameterList.cpp
  DIRECTORY XML
  )
#TRILINOS_ADD_EXECUTABLE_AND_TEST(hilbert EXAMPLE SOURCES cxx_main.cpp)
TRILINOS_ADD_EXECUTABLE(
  hilbert  
  COMM serial mpi  
  SOURCES cxx_main.cpp
  DIRECTORY hilbert
  )

ADD_SUBDIRECTORY(ExplicitInstantiation/four_files)
