IF(${Trilinos_ENABLE_Kokkos})

  # This is where to generate the gen_kokkos.cmake and KokkosCore_config.h 
  # that we will use in the configuration
  set(Kokkos_GEN_DIR ${CMAKE_BINARY_DIR})

  # Basic initialization (Used in KOKKOS_SETTINGS)
  set(KOKKOS_SRC_PATH ${Kokkos_SOURCE_DIR})
  set(KOKKOS_PATH ${KOKKOS_SRC_PATH})

  #------------ COMPILER AND FEATURE CHECKS ------------------------------------
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_functions.cmake)
  set_kokkos_cxx_compiler()
  set_kokkos_compiler_standard()
  set_kokkos_cxx_standard()
  
  #------------ GET OPTIONS ----------------------------------------------------
  set(KOKKOS_CMAKE_VERBOSE True)
  set(KOKKOS_HAS_TRILINOS True)
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_options.cmake)

  #------------ COMPUTE KOKKOS_SETTINGS ----------------------------------------
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_settings.cmake)

  #------------ GENERATE HEADER AND SOURCE FILES -------------------------------
  execute_process(
    COMMAND ${KOKKOS_SETTINGS} make -f ${KOKKOS_SRC_PATH}/core/src/Makefile build-makefile-cmake-kokkos
    WORKING_DIRECTORY "${Kokkos_GEN_DIR}"
    OUTPUT_FILE ${Kokkos_GEN_DIR}/core_src_make.out
    RESULT_VARIABLE res
  )
  include(${Kokkos_GEN_DIR}/gen_kokkos.cmake)

  set(CXXFLAGS ${CXXFLAGS} ${KOKKOS_CXXFLAGS})
  # TODO -- need to remove the -lkokkos.  Check on LDFlags
  #set(KOKKOS_LINK_DEPENDS libkokkos.a CACHE STRING "")
  #set(KOKKOS_LIBS -lkokkos -ldl -lpthread CACHE STRING "")
  #set(KOKKOS_LDFLAGS -L/scr_gabrielle/kruger/builds-ptsolve/trilinos/par2 --gcc-toolchain=/usr CACHE STRING "")

ENDIF()
