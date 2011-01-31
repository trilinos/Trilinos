
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  BLAS
  Comm
  Containers
  Conditions
  DenseMatrix
  Dependencies
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
  Validators
  Workspace
  XML
  XMLParser
  dyn_cast
  )
