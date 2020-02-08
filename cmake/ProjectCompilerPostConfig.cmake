IF (${PROJECT_NAME}_ENABLE_Kokkos)

  MESSAGE("-- " "Skip adding flags for C++11 because Kokkos flags does that ...")
  # Set this to empty to trick Tribits into passing the C++11 flag check
  SET(${PROJECT_NAME}_CXX11_FLAGS)
  
  MESSAGE("-- " "Skip adding flags for OpenMP because Kokkos flags does that ...")
  SET(OpenMP_CXX_FLAGS_OVERRIDE " ")

  # There is a top-level CMAKE_CXX_FLAGS. Warnings and other flags get added
  # in the sub-scope of each individual package
  # Grab those variables here so we know what was the original top-level flags
  # and what are the CMAKE_CXX_FLAGS added afterwards for an individual package
  SET(TRILINOS_TOPLEVEL_CXX_FLAGS ${CMAKE_CXX_FLAGS})
ENDIF()
