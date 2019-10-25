IF (${Trilinos_ENABLE_Kokkos})

  MESSAGE("-- " "Skip adding flags for C++11 because Kokkos flags does that ...")
  SET(${PROJECT_NAME}_CXX11_FLAGS " ")
  
  MESSAGE("-- " "Skip adding flags for OpenMP because Kokkos flags does that ...")
  SET(OpenMP_CXX_FLAGS_OVERRIDE " ")

ENDIF()
