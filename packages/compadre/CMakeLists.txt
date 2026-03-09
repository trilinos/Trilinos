CMAKE_MINIMUM_REQUIRED(VERSION 3.24.0 FATAL_ERROR)

# get Compadre_VERSION
file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/cmake/Compadre_Version.txt Compadre_VERSION)

# Called at the top of every CMakeLists.txt file
macro(include_tribits_build)
  if (${PROJECT_NAME}_TRIBITS_DIR)
    include("${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.tribits.cmake")
    return()
  endif()
endmacro()

# cmake helper functions
include(${CMAKE_CURRENT_LIST_DIR}/cmake/bob.cmake)

if (${PROJECT_NAME}_TRIBITS_DIR) # TriBITS CMake Project

  #
  # A) Define the package
  #

  tribits_package(Compadre)

  #
  # B) Set up package-specific options
  #
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
      Compadre_DEBUG
      "Enable debug code in compadre"
      ${Trilinos_ENABLE_DEBUG} )

  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_EXTREME_DEBUG
      Compadre_EXTREME_DEBUG
      "Enable extreme debug code in compadre"
      OFF )

  if ("${TPL_ENABLE_CUDA}" STREQUAL "ON")
      set(Compadre_USE_CUDA   ON  CACHE BOOL "Whether to use CUDA")
  else ()
      set(Compadre_USE_CUDA   OFF CACHE BOOL "Whether to use CUDA")
  endif()

  set(Compadre_USE_MATLAB OFF CACHE BOOL "Don't use Matlab when building in Trilinos")
  set(Compadre_USE_MPI ${TPL_ENABLE_MPI} CACHE BOOL "Use MPI for parallelism")
  set(Compadre_USE_PYTHON OFF CACHE BOOL "Don't use Python when building in Trilinos")

  # breakdown Compadre_VERSION from Version.cmake for storage in Compadre_Config.h
  string(REPLACE "." ";" VLIST ${Compadre_VERSION})
  list(GET VLIST 0 Compadre_VERSION_MAJOR)
  list(GET VLIST 1 Compadre_VERSION_MINOR)
  list(GET VLIST 2 Compadre_VERSION_PATCH)

else() # Raw CMake Project

  ##########
  #
  #   Compadre Details
  #
  ########## 

  # Disable in-source builds to prevent source tree corruption.
  if( "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" )
    message( FATAL_ERROR "FATAL: In-source builds are not allowed. You should create a separate directory for build files and delete CMakeCache.txt." )
  endif()

  PROJECT(Compadre VERSION ${Compadre_VERSION} LANGUAGES CXX)
  if (NOT CMAKE_CXX_STANDARD)
      set(CMAKE_CXX_STANDARD 20)
  endif()
  set(CMAKE_DISABLE_SOURCE_CHANGES ON)
  set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
  set(CMAKE_POLICY_DEFAULT_CMP0128 NEW)
  set(CMAKE_CXX_EXTENSIONS OFF)

  get_directory_property(Compadre_HAS_PARENT PARENT_DIRECTORY)

  if (NOT Compadre_HAS_PARENT)
    include(GNUInstallDirs)
  endif()

  bob_begin_package()

  ##########
  #
  #
  #  CMAKE_BUILD_TYPE and Compadre_DEBUG
  #
  #
  ##########

  # Set to OFF for significantly faster performance and ON for error tracking
  bob_option(Compadre_EXTREME_DEBUG "Run Compadre Toolkit in EXTREME DEBUG mode" OFF)
  # do not print variable now
  bob_option(Compadre_DEBUG "Run Compadre Toolkit in DEBUG mode" OFF OFF)
  if (Compadre_EXTREME_DEBUG)
    set(Compadre_DEBUG ON CACHE BOOL "" FORCE)
    message(STATUS "Compadre_DEBUG: ON (because Compadre_EXTREME_DEBUG is ON)")
  endif()

  # set defaults for build type
  if(NOT CMAKE_BUILD_TYPE)
    if(Compadre_DEBUG)
      set(DEFAULT_BUILD_TYPE "Debug")
    else()
      set(DEFAULT_BUILD_TYPE "RelWithDebInfo")
    endif()
    message(STATUS "No build type specified. Defaulting to: ${DEFAULT_BUILD_TYPE}")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING
        "Choose the type of build: {Debug,Release,RelWithDebInfo,MinSizeRel,None}." FORCE)
  endif()

  # ensure conflicting BUILD_TYPE not specified w.r.t. Compadre_DEBUG
  # Compadre_DEBUG takes precedence
  if(Compadre_DEBUG AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Compadre_DEBUG is ON so CMAKE_BUILD_TYPE set to Debug")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type {Release, RelWithDebInfo, Debug, None}" FORCE)
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT Compadre_DEBUG)
    message(STATUS "CMAKE_BUILD_TYPE is 'Debug' so Compadre_DEBUG is set to ON")
    set(Compadre_DEBUG ON CACHE STRING "" FORCE)
  endif()

  message(STATUS "Compadre_DEBUG: ${Compadre_DEBUG}")
  message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

  ##########
  #
  #
  #  Python and RPATH logicCMAKE_BUILD_TYPE and Compadre_DEBUG
  #
  #
  ##########

  # RPATH related settings
  # https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  SET(CMAKE_INSTALL_RPATH "../${CMAKE_INSTALL_LIBDIR}")

  # RPATH should always include the folder it is called from
  bob_option(PYTHON_CALLING_BUILD "Python setuptools calling build" OFF)
  bob_option(Compadre_USE_PYTHON "Use PYTHON" OFF)
  bob_option(Compadre_USE_MATLAB "Use MATLAB interface for PYTHON" OFF)
  bob_option(Compadre_USE_MPI "Use MPI (not needed for Compadre toolkit)" OFF)
  bob_option(Compadre_USE_CUDA "Require CUDA" OFF)
  set(PYTHON_LIBRARY_PREFIX "..") # relative to examples folder
  bob_input(PYTHON_EXECUTABLE "" PATH "Python executable location")
  if (Compadre_USE_PYTHON)
    IF(PYTHON_EXECUTABLE)
      set(Python_EXECUTABLE "${PYTHON_EXECUTABLE}")
    ENDIF()
    find_package(Python REQUIRED COMPONENTS Interpreter Development)
    if (${Python_VERSION} VERSION_LESS "3.10")
        MESSAGE(FATAL_ERROR "
        Compadre_USE_PYTHON:BOOL=ON but Python version <3.10 is not supported.
        ")
    endif()
    MESSAGE(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")
    # change RPATH for a flat directory structure
    # when installing pycompadre as Python package
    if (PYTHON_CALLING_BUILD)
      if (APPLE)
        SET(CMAKE_BUILD_RPATH ${CMAKE_BUILD_RPATH} "@loader_path/")
        SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} "@loader_path/")
      else()
        SET(CMAKE_BUILD_RPATH ${CMAKE_BUILD_RPATH} "$ORIGIN/")
        SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} "$ORIGIN/")
        SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-z,origin")
      endif()
    endif()
  endif()



  ##########
  #
  #
  #  TPLS
  #
  #  Three choices for providing TPLs:
  #
  #      a.) None
  #          - Trilinos_ROOT, Kokkos_ROOT, and KokkosKernels_ROOT 
  #            not specified
  #      b.) Trilinos installation with Kokkos AND KokkosKernels
  #          - Trilinos_ROOT specified
  #      c.) Kokkos AND KokkosKernels standalone installations
  #          - Kokkos_ROOT and KokkosKernels_ROOT specified
  #
  #  For option a, vendored Kokkos and KokkosKernels built for user
  #
  ##########



  #TRILINOS specify with Trilinos_ROOT

  #KOKKOS specify with Kokkos_ROOT

  #KOKKOS-KERNELS specify with KokkosKernels_ROOT
  set(KOKKOS_IN_TRILINOS OFF)
  set(KOKKOS_EXISTING_ELSEWHERE OFF)
  set(KOKKOS_BUILT_FOR_USER OFF)
  set(KOKKOSKERNELS_IN_TRILINOS OFF)
  set(KOKKOSKERNELS_EXISTING_ELSEWHERE OFF)
  set(KOKKOSKERNELS_BUILT_FOR_USER OFF)
  set(KOKKOS_MIN "4.5.01")

  if (Trilinos_ROOT)
    FIND_PACKAGE(Trilinos COMPONENTS Kokkos KokkosKernels OPTIONAL_COMPONENTS Gtest)

    LIST(REVERSE Trilinos_INCLUDE_DIRS)
    LIST(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
    LIST(REVERSE Trilinos_INCLUDE_DIRS)

    LIST(REVERSE Trilinos_LIBRARIES)
    LIST(REMOVE_DUPLICATES Trilinos_LIBRARIES)
    LIST(REVERSE Trilinos_LIBRARIES)
    
    LIST(REVERSE Trilinos_TPL_LIBRARIES)
    LIST(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
    LIST(REVERSE Trilinos_TPL_LIBRARIES)
    
    MESSAGE("\nFound Trilinos!  Here are the details: ")
    MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
    MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
    MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
    MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
    MESSAGE("   Trilinos_BIN_DIRS = ${Trilinos_BIN_DIRS}")
    MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
    MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
    MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
    MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
    MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
    MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
    
    LIST(REVERSE Trilinos_INCLUDE_DIRS)
    LIST(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
    LIST(REVERSE Trilinos_INCLUDE_DIRS)
    MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_CXX_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
    MESSAGE("End of Trilinos details\n")

    IF (Trilinos_VERSION VERSION_GREATER_EQUAL "16.1")
        set(CompadreKokkosPackageName "Kokkos")
    ELSE ()
        MESSAGE(FATAL_ERROR "Found Trilinos but Compadre requires version 16.1 or greater")
    ENDIF()
    LIST(FIND Trilinos_PACKAGE_LIST ${CompadreKokkosPackageName} ${CompadreKokkosPackageName}ID)
    IF (${CompadreKokkosPackageName}ID GREATER -1 )
      MESSAGE(STATUS "Found ${CompadreKokkosPackageName} inside Trilinos!")
      set(KOKKOS_IN_TRILINOS ON)
    ELSE()
      MESSAGE(FATAL_ERROR "Found Trilinos but could not find ${CompadreKokkosPackageName}.")
    ENDIF()

    LIST(FIND Trilinos_PACKAGE_LIST KokkosKernels KokkosKernelsID)
    IF (KokkosKernelsID GREATER -1 )
      MESSAGE(STATUS "Found KokkosKernels inside Trilinos!")
      set(KOKKOSKERNELS_IN_TRILINOS ON)
      set(KokkosKernels_FOUND ON)
    ELSE()
      MESSAGE(FATAL_ERROR "Found Trilinos but could not find KokkosKernels.")
    ENDIF()

    LIST(FIND Trilinos_PACKAGE_LIST Gtest GtestID)

    if(KOKKOS_IN_TRILINOS AND NOT(KOKKOSKERNELS_IN_TRILINOS))
        MESSAGE(FATAL_ERROR "Kokkos found in Trilinos, but not KokkosKernels")
    endif()

    # Check that CUDA_LAMBDA is enabled if using CUDA in Kokkos
    list(FIND Kokkos_DEVICES CUDA cuda_index)
    if (NOT cuda_index EQUAL -1)
      list(FIND Kokkos_OPTIONS CUDA_LAMBDA cuda_lambda_index)
      if (cuda_lambda_index EQUAL -1)
        message(FATAL_ERROR "Please reconfigure Trilinos with -DKokkos_ENABLE_CUDA_LAMBDA:BOOL=ON")
      endif()
    else()
      if (Compadre_USE_CUDA)
        message(FATAL_ERROR "Compadre_USE_CUDA is ON, but CUDA not found in Trilinos")
      endif()
    endif()

    find_package(Kokkos PATHS "${Trilinos_ROOT}/lib/cmake/Kokkos" NO_DEFAULT_PATH REQUIRED)
    message(STATUS "Kokkos Version: ${Kokkos_VERSION}")
    if (Kokkos_VERSION VERSION_LESS "${KOKKOS_MIN}")
        message(FATAL_ERROR "${Kokkos_VERSION}: Requires Trilinos with Kokkos version ${KOKKOS_MIN} or greater")
    endif()

  ENDIF()

  if (NOT KOKKOS_IN_TRILINOS)
    if (NOT Kokkos_ROOT)
      SET(KOKKOS_BUILT_FOR_USER ON)
      if (APPLE)
        bob_option(Kokkos_ENABLE_CUDA "Whether to use CUDA" OFF)
        bob_option(Kokkos_ENABLE_OPENMP "Whether to use OPENMP" OFF)
        bob_option(Kokkos_ENABLE_THREADS "Whether to use THREADS" ON)
      else()
        bob_option(Kokkos_ENABLE_CUDA "Whether to use CUDA" OFF)
        bob_option(Kokkos_ENABLE_OPENMP "Whether to use OPENMP" ON)
        bob_option(Kokkos_ENABLE_THREADS "Whether to use THREADS" OFF)
      endif()
      if (Kokkos_ENABLE_CUDA)
        set(Compadre_USE_CUDA ON CACHE BOOL "" FORCE)
      endif()
      option(Kokkos_ENABLE_DEBUG "" "${Compadre_DEBUG}")
      option(Kokkos_ENABLE_CUDA_UVM "" OFF)
      option(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK "" "${Compadre_EXTREME_DEBUG}")
      option(Kokkos_ENABLE_CUDA_LAMBDA "" "${Kokkos_ENABLE_CUDA}")
      option(Kokkos_ENABLE_TESTS "" OFF)
      option(Kokkos_ENABLE_EXAMPLES "" OFF)
      if (KokkosKernels_ROOT)
          MESSAGE(FATAL_ERROR "KokkosKernels_ROOT specified but Kokkos_ROOT not specified. \
          Either provide Kokkos_ROOT as well, or remove KokkosKernels_ROOT specification.")
      endif()
      # set Compadre_USE_ to use whatever Kokkos_ENABLE_ specified, since user is building Kokkos
      set(Kokkos_ENABLE_CUDA ${Compadre_USE_CUDA})
    else()
      # Existing Kokkos indicated by user by setting Kokkos_ROOT, so we look for it
      # where specified (only)
      find_package(Kokkos PATHS "${Kokkos_ROOT}" NO_DEFAULT_PATH REQUIRED)
      message(STATUS "Kokkos Version: ${Kokkos_VERSION}")
      if (Kokkos_VERSION VERSION_LESS "${KOKKOS_MIN}")
          message(FATAL_ERROR "${Kokkos_VERSION}: Requires Kokkos version ${KOKKOS_MIN} or greater")
      endif()
      set(KOKKOS_EXISTING_ELSEWHERE ON)
      # set Compadre_USE_ to whatever is enabled in Kokkos build
      if(Compadre_USE_CUDA AND NOT(TARGET Kokkos::CUDA))
        MESSAGE(FATAL_ERROR "Compadre_USE_CUDA set to ON, but Kokkos::CUDA not defined")
      endif()
    endif()
  endif()

  if (NOT(KOKKOS_IN_TRILINOS) AND NOT(Kokkos_ROOT))
      add_subdirectory(kokkos)
      get_directory_property(Kokkos_VERSION DIRECTORY kokkos DEFINITION Kokkos_VERSION)
      message(STATUS "Kokkos Version: ${Kokkos_VERSION}")
  endif()
  set(KOKKOS_THREADS_ARG "--kokkos-num-threads")

  if (KOKKOSKERNELS_IN_TRILINOS)
    # register Trilinos so that it shows up in CompadreConfig.cmake
    set(Compadre_EXT_DEPS Trilinos)
  else()
    if (NOT KokkosKernels_ROOT)
      if (Kokkos_ROOT)
          MESSAGE(FATAL_ERROR "Kokkos_ROOT specified but KokkosKernels_ROOT not specified. \
          Either provide KokkosKernels_ROOT as well, or remove Kokkos_ROOT specification.")
      endif()
      SET(KOKKOSKERNELS_BUILT_FOR_USER ON)
      OPTION(KokkosKernels_ENABLE_EXAMPLES "" OFF)
      OPTION(KokkosKernels_ENABLE_TESTS "" OFF)
      OPTION(KokkosKernels_ENABLE_COMPONENT_BATCHED "" ON)
      OPTION(KokkosKernels_ENABLE_COMPONENT_BLAS "" ON)
      OPTION(KokkosKernels_ENABLE_COMPONENT_LAPACK "" OFF)
      OPTION(KokkosKernels_ENABLE_COMPONENT_GRAPH "" OFF)
      OPTION(KokkosKernels_ENABLE_COMPONENT_Sparse "" OFF)
      OPTION(KokkosKernels_ENABLE_COMPONENT_NODE "" OFF)
      OPTION(KokkosKernels_ENABLE_ALL_COMPONENTS "" OFF)
      OPTION(KokkosKernels_ETI_ONLY "" OFF)
      SET(KokkosKernels_ETI_ONLY OFF)
      OPTION(KokkosKernels_INST_OFFSET_SIZE_T "" OFF)
      SET(KokkosKernels_INST_OFFSET_SIZE_T OFF)
      OPTION(KokkosKernels_INST_OFFSET_INT "" OFF)
      SET(KokkosKernels_INST_OFFSET_INT OFF)
      OPTION(KokkosKernels_INST_LAYOUTLEFT "" OFF)
      SET(KokkosKernels_INST_LAYOUTLEFT OFF)
      OPTION(KokkosKernels_INST_LAYOUTRIGHT "" OFF)
      SET(KokkosKernels_INST_LAYOUTRIGHT OFF)
      OPTION(KokkosKernels_INST_DOUBLE "" OFF)
      SET(KokkosKernels_INST_DOUBLE OFF)
      OPTION(KokkosKernels_INST_ORDINAL_INT "" OFF)
      SET(KokkosKernels_INST_ORDINAL_INT OFF)
      OPTION(KokkosKernels_ENABLE_TPL_CUBLAS "" OFF)
      SET(KokkosKernels_ENABLE_TPL_CUBLAS OFF)
      OPTION(KokkosKernels_ENABLE_TPL_CUSPARSE "" OFF)
      SET(KokkosKernels_ENABLE_TPL_CUSPARSE OFF)
      OPTION(KokkosKernels_ENABLE_TPL_CUSOLVER "" OFF)
      SET(KokkosKernels_ENABLE_TPL_CUSOLVER OFF)
      OPTION(KokkosKernels_INST_MEMSPACE_CUDAUVMSPACE "" OFF)
      SET(KokkosKernels_INST_MEMSPACE_CUDAUVMSPACE OFF)
      OPTION(KokkosKernels_INST_EXECSPACE_THREADS "" ${Kokkos_ENABLE_THREADS})
      SET(KokkosKernels_INST_EXECSPACE_THREADS ${Kokkos_ENABLE_THREADS})
      OPTION(KokkosKernels_INST_EXECSPACE_OPENMP "" ${Kokkos_ENABLE_OPENMP})
      SET(KokkosKernels_INST_EXECSPACE_OPENMP ${Kokkos_ENABLE_OPENMP})
      add_subdirectory(kokkos-kernels)
    else()
      # Existing KokkosKernels indicated by user by setting KokkosKernels_ROOT, so we look for it
      # where specified (only)
      find_package(KokkosKernels PATHS "${KokkosKernels_ROOT}" NO_DEFAULT_PATH REQUIRED)
      set(KOKKOSKERNELS_EXISTING_ELSEWHERE ON)
    endif()
    # register Kokkos and KokkosKernels so that they show up in CompadreConfig.cmake
    set(Compadre_EXT_DEPS Kokkos KokkosKernels)
  endif()



  #MPI (Not really needed, only used so that if a kokkos-tool was built with MPI, it won't segfault)
  if (Compadre_USE_MPI)
    FIND_PACKAGE(MPI REQUIRED)
    bob_option(Compadre_USE_MPI "Use MPI for parallelism" ${MPI_CXX_FOUND})
    message(STATUS "MPI Enabled: ${MPI_CXX_FOUND}")
    if (MPI_CXX_FOUND)
        MESSAGE(STATUS "MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
        MESSAGE(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
    endif()
  endif (Compadre_USE_MPI)



  ##########
  #
  #
  #  OPTIONS
  #
  #
  ##########



  bob_option(Compadre_TESTS "Compile Compadre tests" ON)
  bob_option(Compadre_EXAMPLES "Compile Compadre examples" "${Compadre_TESTS}")
  bob_option(Compadre_COVERAGE "Use code coverage" OFF)


endif() 

##########
#
#
#  PREPROCESSOR DEFINES FOR USE WHEN BUILDING CODE
#
#
##########


# if you want to change optimization options, i.e.
# setting CMAKE_CXX_FLAGS=-O3 with CMAKE_BUILD_TYPE=RelWithDebInfo
# will result in flags: '-O3 -O2 -g' which will disregard -O3
# to address this, set: CMAKE_CXX_FLAGS_RELWITHDEBINFO = '-O3 -g' 
# if you want to have the final effect of '-O3 -g'

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(Compadre_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(Compadre_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(Compadre_BUILD_TYPE "${CMAKE_BUILD_TYPE}")

set(Compadre_KEY_BOOLS
    Compadre_DEBUG
    Compadre_EXTREME_DEBUG
    Compadre_USE_CUDA
    Compadre_USE_MPI
    Compadre_USE_PYTHON
    PYTHON_CALLING_BUILD
   )

set(Compadre_KEY_INTS
    Compadre_VERSION_MAJOR
    Compadre_VERSION_MINOR
    Compadre_VERSION_PATCH
   )

set(Compadre_KEY_STRINGS
    Compadre_SEMVER
    Compadre_COMMIT
    Compadre_CXX_FLAGS
    Compadre_BUILD_TYPE
    Compadre_CMAKE_ARGS
    Compadre_INSTALL_PREFIX
)

if (${PROJECT_NAME}_TRIBITS_DIR) # TriBITS CMake Project

  #
  # C) Add the libraries, tests, and examples
  #

  add_subdirectory(src)
  tribits_add_test_directories(examples)

  #
  # D) Do standard postprocessing
  #

  tribits_package_postprocess()

else() # Raw CMake Project

  if (NOT(PYTHON_CALLING_BUILD))
    bob_get_commit()
    bob_get_semver()
  endif()

  ##########
  #
  #
  #  COMPADRE LIBRARIES
  #
  #
  ##########

  if (Compadre_COVERAGE)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
    include(CodeCoverage)
    append_coverage_compiler_flags()
    setup_target_for_coverage_lcov(
      NAME coverage
      EXECUTABLE ctest
      EXCLUDE "/Library/*" "${PROJECT_SOURCE_DIR}/kokkos/*" 
              "${PROJECT_SOURCE_DIR}/kokkos-kernels/*" 
              "${PROJECT_SOURCE_DIR}/pycompadre/*")
  endif()


  # compadre library
  add_subdirectory(src)

  if (Compadre_EXAMPLES AND Compadre_TESTS)
    # gtest from kokkos
    # based on kokkos-kernels/CMakeLists.txt
    SET(GTEST_SOURCE_DIR ${PROJECT_SOURCE_DIR}/kokkos/tpls/gtest)
    message(STATUS "GTEST:: ${GTEST_SOURCE_DIR}")
    add_library(compadre_gtest
        ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
    )
    target_include_directories(compadre_gtest PUBLIC ${GTEST_SOURCE_DIR}/gtest)
    # Disables pthreads, this is a problem for serial builds in Trilinos & Sierra if it's enabled.
    TARGET_COMPILE_DEFINITIONS(compadre_gtest PUBLIC "-DGTEST_HAS_PTHREAD=0")
    TARGET_INCLUDE_DIRECTORIES(compadre_gtest PUBLIC $<BUILD_INTERFACE:${GTEST_SOURCE_DIR}>)
    #Gtest minimally requires C++11
    TARGET_COMPILE_FEATURES(compadre_gtest PUBLIC cxx_std_11)
  endif()

  if(Compadre_EXAMPLES)
    add_subdirectory(examples)
  endif()

  if(Compadre_USE_PYTHON)
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE nanobind_DIR
    )
    find_package(nanobind CONFIG REQUIRED)
    nanobind_add_module(_pycompadre NB_DOMAIN "pycompadre_domain" pycompadre/pycompadre.cpp)
    target_link_libraries(_pycompadre PRIVATE compadre)
    if (PYTHON_CALLING_BUILD)
      configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pycompadre/__init__.py" "${CMAKE_CURRENT_BINARY_DIR}/__init__.py" COPYONLY)
    else()
      set_target_properties(_pycompadre PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pycompadre")
      configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pycompadre/__init__.py" "${CMAKE_CURRENT_BINARY_DIR}/pycompadre/__init__.py" COPYONLY)
    endif()
    set(CMAKE_POLICY_DEFAULT_CMP0127 NEW)
  endif()


  if (Trilinos_ROOT)
    bob_end_package_no_recurse("${Trilinos_LIBRARIES}")
  else()
    bob_end_package()
  endif()
endif()

## print variables beginning with some string
#function (getListOfVarsStartingWith _prefix _varResult)
#  get_cmake_property(_vars VARIABLES)
#  string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
#  set (${_varResult} ${_matchedVars} PARENT_SCOPE)
#endfunction()
#getListOfVarsStartingWith("SomeString" matchedVars)
#foreach (_var IN LISTS matchedVars)
#    message("${_var}=${${_var}}")
#endforeach()
