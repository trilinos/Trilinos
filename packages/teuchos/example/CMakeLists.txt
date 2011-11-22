
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  BLAS
  CommandLineProcessor
  DenseMatrix
  ExplicitInstantiation/four_files
  LAPACK
  ParameterList
  RefCountPtr
  TimeMonitor
  XML
  hilbert
  show_stack
  )
