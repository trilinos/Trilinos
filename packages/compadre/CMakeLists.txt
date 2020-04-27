CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0 FATAL_ERROR)

set(Compadre_VERSION "1.0.3")

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

  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_LAPACK_DECLARED_THREADSAFE
      LAPACK_DECLARED_THREADSAFE
      "Set to OFF if non-threadsafe version of LAPACK being used"
      ON )

  if ("${TPL_ENABLE_CUDA}" STREQUAL "ON")
      set(Compadre_USE_LAPACK OFF CACHE BOOL "Whether to use LAPACK")
      set(Compadre_USE_CUDA   ON  CACHE BOOL "Whether to use CUDA")
  else ()
      set(Compadre_USE_LAPACK ON  CACHE BOOL "Whether to use LAPACK")
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

  PROJECT(Compadre VERSION ${Compadre_VERSION} LANGUAGES CXX)

  get_directory_property(Compadre_HAS_PARENT PARENT_DIRECTORY)

  bob_begin_package()

  # Set to OFF for significantly faster performance and ON for error tracking
  bob_option(Compadre_DEBUG "Run Compadre Toolkit in DEBUG mode" ON)
  bob_option(Compadre_EXTREME_DEBUG "Run Compadre Toolkit in EXTREME DEBUG mode" OFF)

  # RPATH related settings
  # https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
  SET(CMAKE_CXX_STANDARD 11)
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_LIBDIR} ${CMAKE_INSTALL_RPATH})

  # RPATH should always include the folder it is called from
  bob_option(PYTHON_CALLING_BUILD "Python setuptools calling build" OFF)
  bob_option(Compadre_USE_PYTHON "Use PYTHON" OFF)
  if (Compadre_USE_PYTHON)
    if (APPLE)
      SET(CMAKE_BUILD_RPATH ${CMAKE_BUILD_RPATH} "@loader_path/")
      SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} "@loader_path/")
      if (NOT(PYTHON_CALLING_BUILD))
        # If python is calling the build, it will provide these flags
        SET(CMAKE_PYTHON_SHARED_LINKER_FLAGS "-undefined dynamic_lookup -nodefaultlibs")
      endif()
    else()
      SET(CMAKE_BUILD_RPATH ${CMAKE_BUILD_RPATH} "$ORIGIN/")
      SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} "$ORIGIN/")
      SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-z,origin")
    endif()
  endif()



  ##########
  #
  #
  #  TPLS
  #
  #
  ##########



  #TRILINOS
  bob_input(Trilinos_PREFIX "" PATH "Path to Trilinos install")
  if (Trilinos_PREFIX)
    if(NOT IS_ABSOLUTE ${Trilinos_PREFIX})
      set(Trilinos_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/${Trilinos_PREFIX}")
    endif()
    set(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
  endif()

  #KOKKOS
  bob_input(KokkosCore_PREFIX "" PATH "Path to KokkosCore install")
  if (KokkosCore_PREFIX)
    if(NOT IS_ABSOLUTE ${KokkosCore_PREFIX})
      set(KokkosCore_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/${KokkosCore_PREFIX}")
    endif()
    set(CMAKE_PREFIX_PATH ${KokkosCore_PREFIX} ${CMAKE_PREFIX_PATH})
  endif()

  #KOKKOS-KERNELS
  bob_input(KokkosKernels_PREFIX "" PATH "Path to KokkosKernels install")
  if (KokkosKernels_PREFIX)
    if(NOT IS_ABSOLUTE ${KokkosKernels_PREFIX})
      set(KokkosKernels_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/${KokkosKernels_PREFIX}")
    endif()
    set(CMAKE_PREFIX_PATH ${KokkosKernels_PREFIX} ${CMAKE_PREFIX_PATH})
  endif()
  set(KOKKOS_IN_TRILINOS OFF)
  set(KOKKOS_EXISTING_ELSEWHERE OFF)
  set(KOKKOS_BUILT_FOR_USER OFF)
  set(KOKKOSKERNELS_IN_TRILINOS OFF)
  set(KOKKOSKERNELS_EXISTING_ELSEWHERE OFF)
  set(KOKKOSKERNELS_BUILT_FOR_USER OFF)
  set(Compadre_USE_LAPACK ON) # ON until turned OFF if Compadre_USE_CUDA is eventually set
  bob_option(LAPACK_DECLARED_THREADSAFE "Is LAPACK threadsafe? If not, default to serial calls for solves." ON)

  if (Trilinos_PREFIX)
    FIND_PACKAGE(Trilinos)

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
    MESSAGE("End of Trilinos details\n")
    
    LIST(REVERSE Trilinos_INCLUDE_DIRS)
    LIST(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
    LIST(REVERSE Trilinos_INCLUDE_DIRS)
    MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_CXX_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")

    LIST(FIND Trilinos_PACKAGE_LIST KokkosCore KokkosCoreID)
    IF (KokkosCoreID GREATER -1 )
      MESSAGE(STATUS "Found KokkosCore inside Trilinos!")
      set(KOKKOS_IN_TRILINOS ON)
      set(KokkosCore_FOUND ON)
    ELSE()
      MESSAGE(FATAL_ERROR "Found Trilinos but could not find KokkosCore.")
    ENDIF()

    LIST(FIND Trilinos_PACKAGE_LIST KokkosKernels KokkosKernelsID)
    IF (KokkosKernelsID GREATER -1 )
      MESSAGE(STATUS "Found KokkosKernels inside Trilinos!")
      set(KOKKOSKERNELS_IN_TRILINOS ON)
      set(KokkosKernels_FOUND ON)
    ELSE()
      MESSAGE(FATAL_ERROR "Found Trilinos but could not find KokkosKernels.")
    ENDIF()

    if(KOKKOSKERNELS_IN_TRILINOS AND NOT(KOKKOS_IN_TRILINOS))
        MESSAGE(FATAL_ERROR "KokkosKernels found in Trilinos, but not Kokkos")
    elseif(KOKKOS_IN_TRILINOS AND NOT(KOKKOSKERNELS_IN_TRILINOS))
        MESSAGE(FATAL_ERROR "Kokkos found in Trilinos, but not KokkosKernels")
    endif()

    #  Detect KokkosCore functionality from Trilinos
    include(${CMAKE_CURRENT_LIST_DIR}/cmake/detect_trilinos_opts.cmake)
    detect_trilinos_opts()

    # get kokkos settings, change our settings, error to user if we something they requested is not enabled
    if (KokkosCore_HAS_CUDA)
      if (NOT KokkosCore_HAS_CUDA_LAMBDA)
        message(FATAL_ERROR "Please reconfigure Trilinos with -DKokkos_ENABLE_CUDA_LAMBDA:BOOL=ON")
      endif()
      message(STATUS "CUDA enabled in KokkosCore in Trilinos, setting Compadre_USE_CUDA to ON")
      set(Compadre_USE_CUDA ON)
    endif()

    if (KokkosCore_HAS_OPENMP)
      message(STATUS "OPENMP enabled in KokkosCore in Trilinos, setting Compadre_USE_OPENMP to ON")
      set(Compadre_USE_OPENMP ON)
    endif()

    if (KokkosCore_HAS_PTHREAD)
      message(STATUS "PTHREAD enabled in KokkosCore in Trilinos, setting Compadre_USE_PTHREAD to ON")
      set(Compadre_USE_PTHREAD ON)
    endif()

  ENDIF()

  if (KOKKOS_IN_TRILINOS)
    if (NOT(TARGET kokkoscore))
        message(FATAL_ERROR "Target kokkoscore not found.")
    endif()
    if (NOT(TARGET kokkoskernels))
        message(FATAL_ERROR "Target kokkoskernels not found.")
    endif()
  else()
    if (KokkosCore_PREFIX STREQUAL "")
      SET(KOKKOS_BUILT_FOR_USER ON)
      if (APPLE)
        bob_option(Kokkos_ENABLE_CUDA "Whether to use CUDA" OFF)
        bob_option(Kokkos_ENABLE_OPENMP "Whether to use OPENMP" OFF)
        bob_option(Kokkos_ENABLE_PTHREAD "Whether to use PTHREAD" ON)
      else()
        bob_option(Kokkos_ENABLE_CUDA "Whether to use CUDA" OFF)
        bob_option(Kokkos_ENABLE_OPENMP "Whether to use OPENMP" ON)
        bob_option(Kokkos_ENABLE_PTHREAD "Whether to use PTHREAD" OFF)
      endif()
      option(Kokkos_ENABLE_CUDA_LAMBDA "" "${Kokkos_ENABLE_CUDA}")
      option(Kokkos_ENABLE_TESTS "" OFF)
      option(Kokkos_ENABLE_EXAMPLES "" OFF)
      add_subdirectory(kokkos)
      if (NOT(KokkosKernels_PREFIX STREQUAL ""))
          MESSAGE(FATAL_ERROR "KokkosKernels_PREFIX specified but KokkosCore_PREFIX not specified. \
          Either provide KokkosCore_PREFIX as well, or remove KokkosKernels_PREFIX specification.")
      endif()
      # set Compadre_USE_ to use whatever Kokkos_ENABLE_ specified, since user is building Kokkos
      set(Compadre_USE_OPENMP ${Kokkos_ENABLE_OPENMP})
      set(Compadre_USE_PTHREAD ${Kokkos_ENABLE_PTHREAD})
      set(Compadre_USE_CUDA ${Kokkos_ENABLE_CUDA})
      set(Kokkos_PREFIX ${CMAKE_INSTALL_PREFIX})
    else()
      # Existing KokkosCore indicated by user by setting KokkosCore_PREFIX, so we look for it
      # where specified (only)
      find_package(Kokkos PATHS "${KokkosCore_PREFIX}" NO_DEFAULT_PATH REQUIRED)
      set(KOKKOS_EXISTING_ELSEWHERE ON)
      # set Compadre_USE_ to whatever is enabled in Kokkos build
      if (TARGET Kokkos::OPENMP)
        set(Compadre_USE_OPENMP ON)
      elseif(TARGET Kokkos::PTHREAD)
        set(Compadre_USE_PTHREAD ON)
      endif() 
      if (TARGET Kokkos::CUDA)
        set(Compadre_USE_CUDA ON)
      endif()
      set(Kokkos_PREFIX "${KokkosCore_PREFIX}")
    endif()
  endif()
  # set LAPACK to off if CUDA is ON
  if (Compadre_USE_CUDA)
    set(Compadre_USE_LAPACK OFF)
  endif()
  message(STATUS "Compadre_USE_LAPACK: ${Compadre_USE_LAPACK}")
  message(STATUS "Compadre_USE_CUDA: ${Compadre_USE_CUDA}")

  if (KOKKOSKERNELS_IN_TRILINOS)
    if (Compadre_USE_LAPACK)
      list(FIND Trilinos_TPL_LIST "LAPACK" LAPACKID)
      if (NOT(LAPACKID GREATER -1))
        MESSAGE(FATAL_ERROR "Compadre_USE_LAPACK set to ON, but LAPACK not defined in Trilinos")
      endif()
      list(FIND Trilinos_TPL_LIST "BLAS" BLASID)
      if (NOT(BLASID GREATER -1))
          MESSAGE(FATAL_ERROR "Compadre_USE_LAPACK set to ON, but BLAS not defined in Trilinos")
      endif()
    endif()
    if (Compadre_US_CUDA)
      list(FIND Trilinos_TPL_LIST "CUDA" CUDAID)
      if (NOT(CUDAID GREATER -1))
          MESSAGE(FATAL_ERROR "Compadre_USE_CUDA set to ON, but CUDA not defined in Trilinos")
      endif()
    endif()
    # register Trilinos so that it shows up in CompadreConfig.cmake
    set(Compadre_EXT_DEPS Trilinos)
  else()
    if (KokkosKernels_PREFIX STREQUAL "")
      SET(KOKKOSKERNELS_BUILT_FOR_USER ON)
      # both are required to be sure it finds it
      OPTION(KokkosKernels_ENABLE_TPL_LAPACK "" "${Compadre_USE_LAPACK}")
      SET(KokkosKernels_ENABLE_TPL_LAPACK "${Compadre_USE_LAPACK}")
      OPTION(KokkosKernels_ENABLE_TPL_BLAS "" "${Compadre_USE_LAPACK}")
      SET(KokkosKernels_ENABLE_TPL_BLAS "${Compadre_USE_LAPACK}")
      OPTION(KokkosKernels_ENABLE_TPL_CUSPARSE "" "${Compadre_USE_CUDA}")
      SET(KokkosKernels_ENABLE_TPL_CUSPARSE "${Compadre_USE_CUDA}")
      OPTION(KokkosKernels_ENABLE_TPL_CUBLAS "" "${Compadre_USE_CUDA}")
      SET(KokkosKernels_ENABLE_TPL_CUBLAS "${Compadre_USE_CUDA}")
      SET(KokkosKernels_ENABLE_EXAMPLES OFF)
      SET(KokkosKernels_ENABLE_TESTS OFF)
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
      add_subdirectory(kokkos-kernels)
      if (NOT(KokkosCore_PREFIX STREQUAL ""))
          MESSAGE(FATAL_ERROR "KokkosCore_PREFIX specified but KokkosKernels_PREFIX not specified. \
          Either provide KokkosKernels_PREFIX as well, or remove KokkosCore_PREFIX specification.")
      endif()
      set(KokkosKernels_PREFIX ${CMAKE_INSTALL_PREFIX})
    else()
      # Existing KokkosCore indicated by user by setting KokkosCore_PREFIX, so we look for it
      # where specified (only)
      find_package(KokkosKernels PATHS "${KokkosKernels_PREFIX}" NO_DEFAULT_PATH REQUIRED)
      if(Compadre_USE_LAPACK AND NOT(TARGET Kokkos::LAPACK))
        MESSAGE(FATAL_ERROR "Compadre_USE_LAPACK set to ON, but Kokkos::LAPACK not defined")
      endif()
      if(Compadre_USE_BLAS AND NOT(TARGET Kokkos::BLAS))
        MESSAGE(FATAL_ERROR "Compadre_USE_BLAS set to ON, but Kokkos::BLAS not defined")
      endif()
      if(Compadre_USE_CUDA AND NOT(TARGET Kokkos::CUDA))
        MESSAGE(FATAL_ERROR "Compadre_USE_CUDA set to ON, but Kokkos::CUDA not defined")
      endif()
      set(KOKKOSKERNELS_EXISTING_ELSEWHERE ON)
    endif()
    # register Kokkos and KokkosKernels so that they show up in CompadreConfig.cmake
    set(Compadre_EXT_DEPS Kokkos KokkosKernels)
  endif()



  #MPI (Not really needed, only used so that if a kokkos-tool was built with MPI, it won't segfault)
  FIND_PACKAGE(MPI QUIET)
  bob_option(Compadre_USE_MPI "Use MPI for parallelism" ${MPI_CXX_FOUND})
  message(STATUS "MPI Enabled: ${MPI_CXX_FOUND}")
  if (MPI_CXX_FOUND)
      MESSAGE(STATUS "MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
      MESSAGE(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
  endif()



  #PYTHON
  bob_option(Compadre_USE_MATLAB "Use MATLAB interface for PYTHON" OFF)
  bob_input(PYTHON_PREFIX "" PATH "Path to PYTHON install")
  if (PYTHON_PREFIX)
    set(CMAKE_PREFIX_PATH ${PYTHON_PREFIX} ${CMAKE_PREFIX_PATH})
  endif()
  if (Compadre_USE_PYTHON)

    bob_input(PYTHON_EXECUTABLE "" PATH "Python executable location")
    IF(NOT(PYTHON_EXECUTABLE))
      MESSAGE(STATUS "Python executable location PYTHON_EXECUTABLE not given. Search made using 'which python'")
      EXECUTE_PROCESS(
        COMMAND which "python"
        OUTPUT_VARIABLE PYTHON_EXECUTABLE
        OUTPUT_STRIP_TRAILING_WHITESPACE )
    ENDIF()
    MESSAGE(STATUS "PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}")
    
    EXECUTE_PROCESS(
      COMMAND "${PYTHON_EXECUTABLE}" -c "import site; print(site.USER_SITE)"
      OUTPUT_VARIABLE PYTHON_SITEPACKAGES
      OUTPUT_STRIP_TRAILING_WHITESPACE )
    
    IF(NOT EXISTS ${PYTHON_SITEPACKAGES})
      EXECUTE_PROCESS(
        COMMAND "${PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('platlib'))"
        OUTPUT_VARIABLE PYTHON_SITEPACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE )
    ENDIF()
    MESSAGE(STATUS "PYTHON_SITEPACKAGES: ${PYTHON_SITEPACKAGES}")

    IF(NOT EXISTS ${Numpy_PREFIX})
      EXECUTE_PROCESS(
        COMMAND "${PYTHON_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
        OUTPUT_VARIABLE Numpy_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE )
    ENDIF()
    MESSAGE(STATUS "Numpy_PREFIX: ${Numpy_PREFIX}")
    
    FIND_PATH(Numpy_INCLUDE_DIRS numpy/arrayobject.h HINTS ${Numpy_PREFIX})
    IF (Numpy_INCLUDE_DIRS)
      MESSAGE(STATUS "Numpy_INCLUDE_DIRS: ${Numpy_INCLUDE_DIRS}")
    ELSE()
      MESSAGE(FATAL_ERROR "Numpy headers not found, but needed when Compadre_USE_PYTHON:BOOL=ON. Try setting Numpy_PREFIX.")
    ENDIF()

    set(SWIG_PREFIX "../python") # relative to examples folder
    bob_input(GMLS_Module_DEST "${CMAKE_INSTALL_LIBDIR}/python/dist-packages" PATH "Path to install python module")

    IF(NOT EXISTS ${PYTHON_INCLUDE_DIRS})
      EXECUTE_PROCESS(
        COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())"
        OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS
        OUTPUT_STRIP_TRAILING_WHITESPACE )
    ENDIF()

    EXECUTE_PROCESS(
      COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils import sysconfig; print(sysconfig.get_config_var(\"LIBDIR\"))"
      OUTPUT_VARIABLE PYTHON_PATHS
      OUTPUT_STRIP_TRAILING_WHITESPACE )
    MESSAGE(STATUS "PYTHON_PATHS: ${PYTHON_PATHS}")

    file(GLOB PYTHON_LIBNAME
      "${PYTHON_PATHS}/libpython*${CMAKE_SHARED_LIBRARY_SUFFIX}"
    )
    get_filename_component(PYTHON_LIBNAME_BASE "${PYTHON_LIBNAME}" NAME)
    STRING(REGEX REPLACE "lib(.*)${CMAKE_SHARED_LIBRARY_SUFFIX}" "\\1" PYTHON_LIBNAME "${PYTHON_LIBNAME_BASE}")

    MESSAGE(STATUS "PYTHON_LIBNAME: ${PYTHON_LIBNAME}")
    find_library(PYTHON_LIBRARIES NAME "${PYTHON_LIBNAME}" PATHS "${PYTHON_PATHS}" NO_DEFAULT_PATH)

    MESSAGE(STATUS "PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
    MESSAGE(STATUS "PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")

  endif()



  ##########
  #
  #
  #  OPTIONS
  #
  #
  ##########



  bob_option(Compadre_TESTS "Compile Compadre tests" ON)
  bob_option(Compadre_EXAMPLES "Compile Compadre examples" "${Compadre_TESTS}")


endif() 

##########
#
#
#  PREPROCESSOR DEFINES FOR USE WHEN BUILDING CODE
#
#
##########


set(Compadre_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(Compadre_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

set(Compadre_KEY_BOOLS
    Compadre_DEBUG
    Compadre_EXTREME_DEBUG
    Compadre_USE_CUDA
    Compadre_USE_LAPACK
    Compadre_USE_MPI
    Compadre_USE_PYTHON
    LAPACK_DECLARED_THREADSAFE
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
    Compadre_CMAKE_ARGS
    GMLS_Module_DEST
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


  # compadre library
  add_subdirectory(src)

  if(Compadre_EXAMPLES)
    add_subdirectory(examples)
  endif()

  if(Compadre_USE_PYTHON)
    add_subdirectory(python)
  endif()

  ##########
  #
  #
  #  DOXYGEN
  #
  #
  ##########

  find_package(Doxygen QUIET)
  if (DOXYGEN_FOUND)
    if (NOT(Compadre_HAS_PARENT))
      configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile" @ONLY)
      add_custom_target(doc
          COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc 
          VERBATIM 
      )
      set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)
      # still must be called manually with "make doc"
    endif()
  endif()


  if (Trilinos_PREFIX)
    bob_end_package_no_recurse("${Trilinos_LIBRARIES}")
  else()
    bob_end_package()
  endif()
endif()
