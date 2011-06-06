
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  ParameterList # Put this first to define a test-only library
  BLAS
  Comm
  Containers
  Conditions
  DenseMatrix
  Dependencies
  FancyOutputting
  FunctionObjects
  LAPACK
  MemoryManagement
  ObjectBuilder
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
  TwoDArray
  )
