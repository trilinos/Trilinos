INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsAddExecutable)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_HOST
  IMPLICIT_THERMAL_GMRES_HOST
  "Enable the host test for the Implicit Thermal GMRES usecase."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_PTHREAD
  IMPLICIT_THERMAL_GMRES_PTHREAD
  "Enable the pthread test for the Implicit Thermal GMRES usecase."
  OFF
  )

IF ( Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_HOST)
  MESSAGE(STATUS "Including Host test in Implicit Thermal GMRES usecase.")
  TRIBITS_ADD_EXECUTABLE(
    host_mini_test
    SOURCES mainhost.cpp testHost.cpp
    DEPLIBS kokkosdevicehost kokkos
  )
ENDIF (Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES_HOST)


