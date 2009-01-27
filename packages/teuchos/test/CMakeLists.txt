
# Required modules
INCLUDE(PackageAddExecutableAndTest)

ADD_SUBDIRECTORY(UnitTest)
ADD_SUBDIRECTORY(TabularOutputter)
ADD_SUBDIRECTORY(MemoryManagement)
ADD_SUBDIRECTORY(Comm)
ADD_SUBDIRECTORY(DenseMatrix)

# Test definitions

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  BLAS_test  
  DIRECTORY BLAS
  SOURCES cxx_main.cpp
  ARGS -v
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Containers_test  
  DIRECTORY Containers
  SOURCES cxx_main.cpp
  COMM serial mpi  
  NUM_MPI_PROCS 1
  )

PACKAGE_ADD_EXECUTABLE(
  FancyOutputting_test  
  COMM serial mpi  
  SOURCES FancyOutputting_test.cpp
  DIRECTORY FancyOutputting
  )
PACKAGE_ADD_TEST(
  FancyOutputting_test
  DIRECTORY FancyOutputting
  NAME FancyOutputting_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

#PACKAGE_ADD_EXECUTABLE(
#  Handle_test  
#  COMM serial mpi  
#  SOURCES cxx_main.cpp
#  DIRECTORY Handle
#  )
#PACKAGE_ADD_TEST(
#  Handle_test
#  DIRECTORY Handle
#  NAME Handle_test
#  NUM_MPI_PROCS 1
#  COMM serial mpi
#  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  LAPACK_test  
  DIRECTORY LAPACK
  SOURCES cxx_main.cpp
  NUM_MPI_PROCS 1
  COMM serial mpi
  ARGS -v
  STANDARD_PASS_OUTPUT
  )
  
PACKAGE_ADD_EXECUTABLE(
  ScalarTraits_test  
  COMM serial mpi  
  SOURCES ScalarTraits_test.cpp 
  DIRECTORY ScalarTraits
  )
PACKAGE_ADD_TEST(
  ScalarTraits_test
  DIRECTORY ScalarTraits
  NAME ScalarTraits_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE(
  Time_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY Time
  )
PACKAGE_ADD_TEST(
  Time_test
  DIRECTORY Time
  NAME Time_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE(
  Workspace_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY Workspace
  )
PACKAGE_ADD_TEST(
  Workspace_test
  DIRECTORY Workspace
  NAME Workspace_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  XML_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY XML
  )
PACKAGE_ADD_TEST(
  XML_test
  DIRECTORY XML
  NAME XML_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  #STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE(
  XMLParser_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY XMLParser
  )
PACKAGE_ADD_TEST(
  XMLParser_test
  DIRECTORY XMLParser
  NAME XMLParser_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE(
  dyn_cast_test  
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY dyn_cast
  )
PACKAGE_ADD_TEST(
  dyn_cast_test
  DIRECTORY dyn_cast
  NAME dyn_cast_test
  NUM_MPI_PROCS 1
  COMM serial mpi
  )

PACKAGE_ADD_EXECUTABLE(
  exit_test  
  DIRECTORY exit
  COMM serial mpi  
  SOURCES cxx_main.cpp some_c_func.c 
  )
PACKAGE_ADD_TEST(
  exit_test  
  DIRECTORY exit
  NUM_MPI_PROCS 1
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Polling_test  
  DIRECTORY Polling
  SOURCES cxx_main.cpp 
  COMM mpi  
  NUM_MPI_PROCS 2
)
# 2008/12/20: rabartl: Above is commented out for the same reason that
# Containers_test is commented out as described above

PACKAGE_ADD_EXECUTABLE(
  ParameterList_test 
  COMM serial mpi  
  SOURCES cxx_main.cpp 
  DIRECTORY ParameterList
  )
PACKAGE_ADD_TEST(
  ParameterList_test
  DIRECTORY ParameterList
  NAME ParameterList_test
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_EXECUTABLE(
  Polynomial_test 
  COMM serial
  SOURCES cxx_main.cpp 
  DIRECTORY Polynomial
  )
PACKAGE_ADD_TEST(
  Polynomial_test
  DIRECTORY Polynomial
  NAME Polynomial_test
  COMM serial 
  STANDARD_PASS_OUTPUT
)
