
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  BLAS
  Comm
  Containers
  CTestScriptsUnitTests
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
