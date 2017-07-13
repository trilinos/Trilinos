TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES  TeuchosCore
  LIB_OPTIONAL_TPLS  yaml-cpp
  TEST_OPTIONAL_TPLS  yaml-cpp
  )
# NOTE: The optional direct dependency on yaml-cpp is needed to make warnings
# from the yaml-cpp header files go away.  See Trilinos GitHub #1458.
