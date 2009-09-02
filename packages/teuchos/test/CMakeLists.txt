
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  BLAS
  Comm
  Containers
  DenseMatrix
  FancyOutputting
  LAPACK
  MemoryManagement
  ObjectBuilder
  ParameterList
  Polling
  Polynomial
  ScalarTraits
  TabularOutputter
  Time
  UnitTest
  Workspace
  XML
  XMLParser
  dyn_cast
  exit
  )


ASSERT_DEFINED(CMAKE_HOST_SYSTEM_NAME)
IF (NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
  ADD_SUBDIRECTORY(CTestScriptsUnitTests)
ENDIF()
