
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  UnitTesting # // Tests parallel reductions across all tests ...
  Comm
  ParameterList
  Time
  TimeMonitor
  )
