tribits_get_package_enable_status(Kokkos  KokkosEnable "")

macro(disable_warnings_for_deprecated_packages)
    message(STATUS "Disabling all warnings/errors for deprecated packages")
    foreach(package ${DEPRECATED_PACKAGES})
        set(${package}_CXX_FLAGS "-w ${${package}_CXX_FLAGS}")
    endforeach()
endmacro()


macro(enable_warnings warnings)
    message(STATUS "Trilinos warnings enabled: ${warnings}")
    foreach(warning ${warnings})
        set(CMAKE_CXX_FLAGS "-W${warning} -Wno-error=${warning} ${CMAKE_CXX_FLAGS}")
    endforeach()
endmacro()


macro(enable_errors errors)
    message(STATUS "Trilinos warnings-as-errors enabled: ${errors}")
    foreach(error ${errors})
        set(CMAKE_CXX_FLAGS "-Werror=${error} ${CMAKE_CXX_FLAGS}")
    endforeach()
endmacro()


IF (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  MESSAGE("-- " "Adding '-fp-model=precise' to C++ compiler flags because Trilinos needs it when using the Intel OneAPI C++ compiler.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model=precise")
ENDIF()

IF (KokkosEnable)
  MESSAGE("-- " "Skip adding flags for OpenMP because Kokkos flags does that ...")
  SET(OpenMP_CXX_FLAGS_OVERRIDE " ")

  # There is a top-level CMAKE_CXX_FLAGS. Warnings and other flags get added
  # in the sub-scope of each individual package
  # Grab those variables here so we know what was the original top-level flags
  # and what are the CMAKE_CXX_FLAGS added afterwards for an individual package
  SET(TRILINOS_TOPLEVEL_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  # NOTE: Setting TRILINOS_TOPLEVEL_CXX_FLAGS only has any impact if Kokkos is
  # being treated as an internal package.
ENDIF()

set(upcoming_warnings shadow ${Trilinos_ADDITIONAL_WARNINGS})
set(promoted_warnings parentheses sign-compare unused-variable)

if("${Trilinos_WARNINGS_MODE}" STREQUAL "WARN")
    enable_warnings("${upcoming_warnings}")
    enable_errors("${promoted_warnings}")
    disable_warnings_for_deprecated_packages()
elseif("${Trilinos_WARNINGS_MODE}" STREQUAL "ERROR")
    enable_errors("${promoted_warnings};${upcoming_warnings}")
    disable_warnings_for_deprecated_packages()
endif()
