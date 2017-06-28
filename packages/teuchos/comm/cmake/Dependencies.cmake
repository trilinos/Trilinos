TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES TeuchosCore TeuchosParameterList
  TEST_OPTIONAL_TPLS  yaml-cpp
  )

# NOTE: The optional direct dependency on yaml-cpp is needed to make warnings
# from the yaml-cpp header files go away.  See Trilinos GitHub #1458.

# ToDo: Ross Bartlett: Make TeuchosComm (this package) only optionally depend
# on TeuchosParameterList.  However, to do so will require some work and I
# don't have the time to do that right now.
