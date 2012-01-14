
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  ParameterList # Put this first to define a test-only library
  BLAS
  Comm
  Containers
  Conditions
  CWrapperSupport
  DenseMatrix
  Dependencies
  FancyOutputting
  FunctionObjects
  LAPACK
  MemoryManagement
  ObjectBuilder
  Polynomial
  ScalarTraits
  TabularOutputter
  Time
  TwoDArray
  UnitTest
  Validators
  Workspace
  XML
  XMLParser
  dyn_cast
  )
