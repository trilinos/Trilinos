IF (${Trilinos_ENABLE_Kokkos})

  PRINT_VAR(KOKKOS_ARCH)

  # This is where to generate the gen_kokkos.cmake and KokkosCore_config.h 
  # that we will use in the configuration
  set(Kokkos_GEN_DIR ${CMAKE_BINARY_DIR})

  # Enable debug checking in Kokkos by default if
  # ${PROJECT_NAME}_ENABLE_DEBUG=ON
  set(KOKKOS_ENABLE_DEBUG ${${PROJECT_NAME}_ENABLE_DEBUG}
    CACHE BOOL
    "Enable debug checking in Kokkos.")
  set(Kokkos_ENABLE_Debug_Bounds_Check ${KOKKOS_ENABLE_DEBUG}
    CACHE BOOL
    "Enable bounds checking in Kokkos array classes.")

  # Basic initialization (Used in KOKKOS_SETTINGS)
  set(KOKKOS_SRC_PATH ${Kokkos_SOURCE_DIR})
  set(KOKKOS_PATH ${KOKKOS_SRC_PATH})

  #------------ COMPILER AND FEATURE CHECKS ------------------------------------
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_functions.cmake)
  set_kokkos_cxx_compiler()
  set_kokkos_cxx_standard()
  
  #------------ GET OPTIONS ----------------------------------------------------
  set(KOKKOS_CMAKE_VERBOSE True)
  set(KOKKOS_HAS_TRILINOS True)
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_options.cmake)

  #------------ COMPUTE KOKKOS_SETTINGS ----------------------------------------
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_settings.cmake)

  #------------ GENERATE HEADER AND SOURCE FILES -------------------------------
  execute_process(
    COMMAND ${KOKKOS_SETTINGS} make -f ${KOKKOS_SRC_PATH}/cmake/Makefile.generate_cmake_settings CXX=${CMAKE_CXX_COMPILER} generate_build_settings
    WORKING_DIRECTORY "${Kokkos_GEN_DIR}"
    OUTPUT_FILE ${Kokkos_GEN_DIR}/core_src_make.out
    RESULT_VARIABLE GEN_SETTINGS_RESULT
  )
  if (GEN_SETTINGS_RESULT)
    message(FATAL_ERROR "Kokkos settings generation failed:\n"
        "${KOKKOS_SETTINGS} make -f ${KOKKOS_SRC_PATH}/cmake/Makefile.generate_cmake_settings CXX=${CMAKE_CXX_COMPILER} generate_build_settings")
  endif()
  include(${Kokkos_GEN_DIR}/kokkos_generated_settings.cmake)

  IF (NOT KOKKOS_ARCH STREQUAL "None")

    # Convert CMakeList into string for CXX_FLAGS
    set(CMAKE_CXX_FLAGSl "")
    foreach(opt ${KOKKOS_CXX_FLAGS})
      set(CMAKE_CXX_FLAGSl "${CMAKE_CXX_FLAGSl} ${opt}")
    endforeach()
  
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGSl}")
  
    # TODO -- need to remove the -lkokkos.  Check on LDFlags
    #set(KOKKOS_LINK_DEPENDS libkokkos.a CACHE STRING "")
    #set(KOKKOS_LIBS -lkokkos -ldl -lpthread CACHE STRING "")
    #set(KOKKOS_LDFLAGS -L/scr_gabrielle/kruger/builds-ptsolve/trilinos/par2 --gcc-toolchain=/usr CACHE STRING "")
  
    MESSAGE("-- " "Skip adding flags for C++11 because Kokkos flags does that ...")
    SET(${PROJECT_NAME}_CXX11_FLAGS " ")
  
    MESSAGE("-- " "Skip adding flags for OpenMP because Kokkos flags does that ...")
    SET(OpenMP_CXX_FLAGS_OVERRIDE " ")
  
  ENDIF()

  # Above, It is important not to distrube the default configuraiton of
  # Trilinos if KOKKOS_ARCH is not set.  But the implementation of the new
  # Kokkos TriBITS CMake files requires kokkos_generated_settings.cmake be
  # included.

ENDIF()
