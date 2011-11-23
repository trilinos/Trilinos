INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_IMPLICIT_THERMAL_GMRES_HOST
  IMPLICIT_THERMAL_GMRES_TEST_KOKKOS_HOST
  "Enable the host test for the Implicit Thermal GMRES usecase."
  OFF
  )

SET(SRCS "main.cpp")
SET(deplibs "kokkosnodeapi")

IF ( Kokkos_IMPLICIT_THERMAL_GMRES_HOST)
  MESSAGE(STATUS "Including Host test in Implicit Thermal GMRES usecase."
  APPEND_SET(SRCS testHost.cpp)
  APPEND_SET(deplibs kokkosdevicehost)
ENDIF ( Kokkos_IMPLICIT_THERMAL_GMRES_HOST)

TRIBITS_ADD_EXECUTABLE(
  mini_test
  SOURCES ${SRCS}
  DEPLIBS ${deplibs}
)
