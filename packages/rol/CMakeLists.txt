# Rapid Optimization Library

if( COMMAND TRIBITS_PACKAGE )
  set( STANDALONE_ROL FALSE )
else()
  set( STANDALONE_ROL TRUE )
endif()

if( STANDALONE_ROL )

  cmake_minimum_required(VERSION 3.23)

  set(CMAKE_POSITION_INDEPENDENT_CODE ON)

  project( ROL
           VERSION 2.0
           DESCRIPTION "Rapid Optimization Library"
           LANGUAGES CXX )
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  if(NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Build type" FORCE)
  endif()

  # Set F77_BLAS_MANGLE macro based on Fortran-C interface (similar to Kokkos Kernels)
  if("${F77_BLAS_MANGLE}" STREQUAL "")
    enable_language(C)
    enable_language(Fortran)
    include(FortranCInterface)
    FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/FortranCInterface.h MACRO_NAMESPACE "F77_")
    if(FortranCInterface_GLOBAL_SUFFIX STREQUAL "")
      set(F77_BLAS_MANGLE "(name,NAME) ${FortranCInterface_GLOBAL_PREFIX}name")
    else()
      set(F77_BLAS_MANGLE "(name,NAME) ${FortranCInterface_GLOBAL_PREFIX}name ## ${FortranCInterface_GLOBAL_SUFFIX}")
    endif()
    message(STATUS "Detected Fortran name mangling: ${F77_BLAS_MANGLE}")
  endif()

  # Add pugixml for XML parameter file support
  include(FetchContent)
  FetchContent_Declare(
    pugixml
    GIT_REPOSITORY https://github.com/zeux/pugixml.git
    GIT_TAG        v1.14
  )
  FetchContent_MakeAvailable(pugixml)

  include(ROLUtils)

  # Configure the config header
  configure_file("${PROJECT_SOURCE_DIR}/cmake/ROL_config.h.in" ROL_config.h @ONLY)

  add_library(rol SHARED)

  set(SRC ${PROJECT_SOURCE_DIR}/src)

  # Set up compatibility includes (specific directories to avoid conflicts)
  set(ROL_COMPATIBILITY_INCLUDES ${SRC}/compatibility/teuchos/blas
                                 ${SRC}/compatibility/teuchos/la
                                 ${SRC}/compatibility/simple/lapack
                                 ${SRC}/compatibility/simple/mpi
                                 ${SRC}/compatibility/simple/parameterlist
                                 ${SRC}/compatibility/std/shared_ptr
                                 ${SRC}/compatibility/teuchos/stacktrace
                                 ${SRC}/compatibility/teuchos-lite)

  add_library(Teuchos_BLAS OBJECT ${SRC}/compatibility/teuchos-lite/Teuchos_BLAS.cpp)
  add_library(Teuchos_CompObject OBJECT ${SRC}/compatibility/teuchos-lite/Teuchos_CompObject.cpp)
  add_library(Teuchos_TestForException OBJECT ${SRC}/compatibility/teuchos-lite/Teuchos_TestForException.cpp)

  # Use GLOB_RECURSE to find all directories in core components
  set(ROL_CORE_DIRS ${SRC}/algorithm
                    ${SRC}/elementwise
                    ${SRC}/function
                    ${SRC}/oed
                    ${SRC}/sol
                    ${SRC}/status
                    ${SRC}/step
                    ${SRC}/utils
                    ${SRC}/vector
                    ${SRC}/zoo)

  set(ROL_CORE_INCLUDES)
  foreach(core_dir ${ROL_CORE_DIRS})
    file(GLOB_RECURSE subdirs LIST_DIRECTORIES true "${core_dir}/*")
    list(APPEND ROL_CORE_INCLUDES ${core_dir})
    foreach(subdir ${subdirs})
      if(IS_DIRECTORY ${subdir})
        list(APPEND ROL_CORE_INCLUDES ${subdir})
      endif()
    endforeach()
  endforeach()

  target_include_directories(rol PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
                                        ${ROL_COMPATIBILITY_INCLUDES}
                                        ${ROL_CORE_INCLUDES})

  target_link_libraries(rol PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
                                    pugixml
                                    Teuchos_BLAS Teuchos_CompObject Teuchos_TestForException)

  add_library(rol::rol ALIAS rol)

  option(ENABLE_TESTS OFF)
  option(ENABLE_EXAMPLES OFF)

  if(ENABLE_TESTS)
    enable_testing()
    add_subdirectory(test)
  endif()

  if(ENABLE_EXAMPLES)
    add_subdirectory(example)
  endif()

else() # Build with Trilinos
  include(cmake/TrilinosROL.cmake)
endif()
