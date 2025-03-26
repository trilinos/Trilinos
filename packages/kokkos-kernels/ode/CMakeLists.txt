LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/ode/src)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/ode/impl)

IF (NOT KokkosKernels_ENABLE_COMPONENT_BATCHED)
  MESSAGE("blas enabled and batched not enabled, we need to include some headers manually!")
  LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched)
  LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/src)
  LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/batched/dense/impl)
ENDIF()


# Adding unit-tests
# Note BMK: Since ODE has no auto-generated ETI files, this directory does not exist in a build without unit tests.
# This causes configure errors when building an app against a Trilinos install, and the unit test build dir doesn't contain any headers that need to be found.
# But uncomment the next line if ETI headers are ever added.
# KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/ode)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/ode)
