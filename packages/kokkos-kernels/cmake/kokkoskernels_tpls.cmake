function(kokkoskernels_append_config_line LINE)
  global_append(KOKKOSKERNELS_TPL_EXPORTS "${LINE}")
endfunction()

macro(KOKKOSKERNELS_ADD_TPL_OPTION NAME DEFAULT_VALUE DOCSTRING)
  cmake_parse_arguments(TPL "" "DEFAULT_DOCSTRING" "" ${ARGN})
  if(NOT TPL_DEFAULT_DOCSTRING)
    set(TPL_DEFAULT_DOCSTRING "${DEFAULT_VALUE}")
  endif()

  set(_NAME_ORIG ${NAME})
  set(_NAME ${NAME})

  # KokkosKernels uses all uppercase for TPLs while TriBits does not...
  # See https://github.com/kokkos/kokkos-kernels/issues/1059
  if(KOKKOSKERNELS_HAS_TRILINOS)
    # Map KK TPL names to Trilinos TPL names
    if(${_NAME} STREQUAL SUPERLU)
      set(_NAME SuperLU)
    elseif(${_NAME} STREQUAL CHOLMOD)
      set(_NAME Cholmod)
    endif()
  endif()

  kokkoskernels_add_option("ENABLE_TPL_${_NAME}" ${DEFAULT_VALUE} BOOL
    "${DOCSTRING} Default: ${TPL_DEFAULT_DOCSTRING}")

  set(ROOT_DEFAULT $ENV{${_NAME_ORIG}_ROOT})
  kokkoskernels_add_option("${_NAME_ORIG}_ROOT" "${ROOT_DEFAULT}" PATH
    "Location of ${_NAME} install root. Default: None or the value of the environment variable ${_NAME}_ROOT if set")

  if(DEFINED TPL_ENABLE_${_NAME})
    if(TPL_ENABLE_${_NAME} AND NOT KOKKOSKERNELS_ENABLE_TPL_${_NAME})
      message("Overriding KOKKOSKERNELS_ENABLE_TPL_${_NAME_ORIG}=OFF with TPL_ENABLE_${_NAME}=ON")
      set(KOKKOSKERNELS_ENABLE_TPL_${_NAME_ORIG} ON)
      set(KOKKOSKERNELS_ENABLE_TPL_${_NAME} ON)
    elseif(NOT TPL_ENABLE_${_NAME} AND KOKKOSKERNELS_ENABLE_TPL_${_NAME})
      message("Overriding KOKKOSKERNELS_ENABLE_TPL_${_NAME_ORIG}=ON with TPL_ENABLE_${_NAME}=OFF")
      set(KOKKOSKERNELS_ENABLE_TPL_${_NAME_ORIG} OFF)
      set(KOKKOSKERNELS_ENABLE_TPL_${_NAME} OFF)
    endif()
  endif()
  if(KOKKOSKERNELS_ENABLE_TPL_${_NAME})
    list(APPEND KOKKOSKERNELS_TPL_LIST ${_NAME})
  endif()
endmacro()

macro(kokkoskernels_create_imported_tpl NAME)
  cmake_parse_arguments(
    TPL "INTERFACE" "LIBRARY;IMPORTED_NAME"
    "LINK_LIBRARIES;INCLUDES;COMPILE_OPTIONS;LINK_OPTIONS" ${ARGN})

  if(NOT TPL_IMPORTED_NAME)
    set(TPL_IMPORTED_NAME KokkosKernels::${NAME})
  endif()

  set(TPL_${NAME}_IMPORTED_NAME ${TPL_IMPORTED_NAME})

  if(KOKKOSKERNELS_HAS_TRILINOS)
    #TODO: we need to set a bunch of cache variables here
  elseif(TPL_INTERFACE)
    add_library(${NAME} INTERFACE)
    #Give this an importy-looking name
    add_library(${TPL_IMPORTED_NAME} ALIAS ${NAME})
    if(TPL_LIBRARY)
      message(SEND_ERROR "TPL Interface library ${NAME} should not have an IMPORTED_LOCATION")
    endif()
    #Things have to go in quoted in case we have multiple list entries
    if(TPL_LINK_LIBRARIES)
      target_link_libraries(${NAME} INTERFACE ${TPL_LINK_LIBRARIES})
      if(NOT DEFINED ${NAME}_FOUND_INFO)
        set(${NAME}_FOUND_INFO ${TPL_LINK_LIBRARIES})
      endif()
    endif()
    if(TPL_INCLUDES)
      target_include_directories(${NAME} INTERFACE ${TPL_INCLUDES})
      if(NOT DEFINED ${NAME}_FOUND_INFO)
        set(${NAME}_FOUND_INFO ${TPL_INCLUDES})
      endif()
    endif()
    if(TPL_COMPILE_OPTIONS)
      target_compile_options(${NAME} INTERFACE ${TPL_COMPILE_OPTIONS})
    endif()
    if(TPL_LINK_OPTIONS)
      target_link_libraries(${NAME} INTERFACE ${TPL_LINK_OPTIONS})
    endif()
  else()
    add_library(${TPL_IMPORTED_NAME} UNKNOWN IMPORTED)
    if(TPL_LIBRARY)
      set_target_properties(${TPL_IMPORTED_NAME} PROPERTIES IMPORTED_LOCATION ${TPL_LIBRARY})
      if(NOT DEFINED ${NAME}_FOUND_INFO)
        set(${NAME}_FOUND_INFO ${TPL_LIBRARY})
      endif()
    endif()
    #Things have to go in quoted in case we have multiple list entries
    if(TPL_LINK_LIBRARIES)
      set_target_properties(${TPL_IMPORTED_NAME} PROPERTIES INTERFACE_LINK_LIBRARIES "${TPL_LINK_LIBRARIES}")
      if(NOT DEFINED ${NAME}_FOUND_INFO)
        set(${NAME}_FOUND_INFO ${TPL_LINK_LIBRARIES})
      endif()
    endif()
    if(TPL_INCLUDES)
      set_target_properties(${TPL_IMPORTED_NAME} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${TPL_INCLUDES}")
      if(NOT DEFINED ${NAME}_FOUND_INFO)
        set(${NAME}_FOUND_INFO ${TPL_INCLUDES})
      endif()
    endif()
    if(TPL_COMPILE_OPTIONS)
      set_target_properties(${TPL_IMPORTED_NAME} PROPERTIES INTERFACE_COMPILE_OPTIONS "${TPL_COMPILE_OPTIONS}")
    endif()
    if(TPL_LINK_OPTIONS)
      set_target_properties(${TPL_IMPORTED_NAME} PROPERTIES INTERFACE_LINK_LIBRARIES "${TPL_LINK_OPTIONS}")
    endif()
  endif()
endmacro()

macro(kokkoskernels_find_header VAR_NAME HEADER TPL_NAME)
  cmake_parse_arguments(TPL "ALLOW_SYSTEM_PATH_FALLBACK" "" "PATHS" ${ARGN})

  set(${VAR_NAME} "${HEADER}-NOTFOUND")
  set(HAVE_CUSTOM_PATHS FALSE)
  if(NOT ${VAR_NAME} AND ${TPL_NAME}_ROOT)
    #ONLY look in the root directory
    find_path(${VAR_NAME} ${HEADER} PATHS ${${TPL_NAME}_ROOT}/include NO_DEFAULT_PATH)
    set(HAVE_CUSTOM_PATHS TRUE)
  endif()

  if(NOT ${VAR_NAME} AND KOKKOSKERNELS_${TPL_NAME}_ROOT)
    #ONLY look in the root directory
    find_path(${VAR_NAME} ${HEADER} PATHS ${KOKKOSKERNELS_${TPL_NAME}_ROOT}/include NO_DEFAULT_PATH)
    set(HAVE_CUSTOM_PATHS TRUE)
  endif()

  if(NOT ${VAR_NAME} AND TPL_PATHS)
    #we got custom paths
    #ONLY look in these paths and nowhere else
    find_path(${VAR_NAME} ${HEADER} PATHS ${TPL_PATHS} NO_DEFAULT_PATH)
    set(HAVE_CUSTOM_PATHS TRUE)
  endif()

  if(NOT HAVE_CUSTOM_PATHS OR TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    #Now go ahead and look in system paths
    if(NOT ${VAR_NAME})
      find_path(${VAR_NAME} ${HEADER})
    endif()
  endif()
endmacro()

macro(kokkoskernels_find_library VAR_NAME TPL_NAME)
  cmake_parse_arguments(TPL "ALLOW_SYSTEM_PATH_FALLBACK" "" "PATHS;LIBRARY_NAMES" ${ARGN})

  set(${VAR_NAME} "${TPL_NAME}-NOTFOUND")
  set(HAVE_CUSTOM_PATHS FALSE)
  if(NOT ${VAR_NAME} AND ${TPL_NAME}_ROOT)
    find_library(${VAR_NAME}
      NAMES ${TPL_LIBRARY_NAMES}
      PATHS ${${TPL_NAME}_ROOT}/lib ${${TPL_NAME}_ROOT}/lib64
      NO_DEFAULT_PATH)
    set(HAVE_CUSTOM_PATHS TRUE)
  endif()

  if(NOT ${VAR_NAME} AND KOKKOSKERNELS_${TPL_NAME}_ROOT)
    #we got root paths, only look in these paths and nowhere else
    find_library(${VAR_NAME}
      NAMES ${TPL_LIBRARY_NAMES}
      PATHS ${KOKKOSKERNELS_${TPL_NAME}_ROOT}/lib
            ${KOKKOSKERNELS_${TPL_NAME}_ROOT}/lib64
      NO_DEFAULT_PATH)
    set(HAVE_CUSTOM_PATHS TRUE)
  endif()

  if(NOT ${VAR_NAME} AND TPL_PATHS)
    #we got custom paths, only look in these paths and nowhere else
    find_library(${VAR_NAME}
      NAMES ${TPL_LIBRARY_NAMES}
      PATHS ${TPL_PATHS}
      NO_DEFAULT_PATH)
    set(HAVE_CUSTOM_PATHS TRUE)
  endif()

  if(NOT HAVE_CUSTOM_PATHS OR TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    if(NOT ${VAR_NAME})
      #Now go ahead and look in system paths
      find_library(${VAR_NAME} NAMES ${TPL_LIBRARY_NAMES})
    endif()
  endif()

endmacro()

macro(kokkoskernels_find_imported NAME)
  cmake_parse_arguments(
    TPL "INTERFACE;ALLOW_SYSTEM_PATH_FALLBACK" "HEADER;IMPORTED_NAME"
    "LIBRARY;HEADERS;LIBRARIES;HEADER_PATHS;LIBRARY_PATHS" ${ARGN})
  #LIBRARY can be a list of possible library names
  #matching the NAMES keyword to CMake find_library

  if(NOT TPL_MODULE_NAME)
    set(TPL_MODULE_NAME TPL${NAME})
  endif()

  if(TPL_ALLOW_SYSTEM_PATH_FALLBACK)
    set(ALLOW_PATH_FALLBACK_OPT ALLOW_SYSTEM_PATH_FALLBACK)
  else()
    set(ALLOW_PATH_FALLBACK_OPT)
  endif()

  if(NOT TPL_IMPORTED_NAME)
    set(TPL_IMPORTED_NAME KokkosKernels::${NAME})
  endif()

  set(${NAME}_INCLUDE_DIRS)
  if(TPL_HEADER)
    kokkoskernels_find_header(
      ${NAME}_INCLUDE_DIRS ${TPL_HEADER} ${NAME} ${ALLOW_PATH_FALLBACK_OPT}
      PATHS ${TPL_HEADER_PATHS})
  endif()

  foreach(HEADER ${TPL_HEADERS})
    kokkoskernels_find_header(
      HEADER_FIND_TEMP ${HEADER} ${NAME} ${ALLOW_PATH_FALLBACK_OPT}
      PATHS ${TPL_HEADER_PATHS})
    if(HEADER_FIND_TEMP)
      list(APPEND ${NAME}_INCLUDE_DIRS ${HEADER_FIND_TEMP})
    endif()
  endforeach()

  set(${NAME}_LIBRARY)
  if(TPL_LIBRARY)
    kokkoskernels_find_library(${NAME}_LIBRARY ${NAME} ${ALLOW_PATH_FALLBACK_OPT}
      LIBRARY_NAMES ${TPL_LIBRARY}
      PATHS ${TPL_LIBRARY_PATHS})
  endif()

  set(${NAME}_FOUND_LIBRARIES)
  #We must find every library in this list
  foreach(LIB ${TPL_LIBRARIES})
    #we want the actual name, not the name -lblas, etc
    set(LIB_CLEAN ${LIB})
    string(FIND "${LIB}" "-l" PREFIX_IDX)
    if(PREFIX_IDX STREQUAL "0")
      string(SUBSTRING ${LIB} 2 -1 LIB_CLEAN)
    endif()

    kokkoskernels_find_library(${LIB}_LOCATION ${NAME} ${ALLOW_PATH_FALLBACK_OPT}
      LIBRARY_NAMES ${LIB_CLEAN}
      PATHS ${TPL_LIBRARY_PATHS})
    if(${LIB}_LOCATION)
      list(APPEND ${NAME}_FOUND_LIBRARIES ${${LIB}_LOCATION})
    else()
      set(${NAME}_FOUND_LIBRARIES ${${LIB}_LOCATION})
      break()
    endif()
  endforeach()

  include(FindPackageHandleStandardArgs)
  set(TPL_VARS_NEEDED)
  if(TPL_LIBRARY)
    list(APPEND TPL_VARS_NEEDED ${NAME}_LIBRARY)
  endif()
  if(TPL_HEADER)
    list(APPEND TPL_VARS_NEEDED ${NAME}_INCLUDE_DIRS)
  endif()
  if(TPL_HEADERS)
    list(APPEND TPL_VARS_NEEDED ${NAME}_INCLUDE_DIRS)
  endif()
  if(TPL_LIBRARIES)
    list(APPEND TPL_VARS_NEEDED ${NAME}_FOUND_LIBRARIES)
  endif()
  find_package_handle_standard_args(${TPL_MODULE_NAME} REQUIRED_VARS ${TPL_VARS_NEEDED})

  if(${TPL_MODULE_NAME}_FOUND)
    set(IMPORT_TYPE)
    if(TPL_INTERFACE)
      set(IMPORT_TYPE "INTERFACE")
    endif()
    kokkoskernels_create_imported_tpl(${NAME} ${IMPORT_TYPE}
      IMPORTED_NAME ${TPL_IMPORTED_NAME}
      INCLUDES "${${NAME}_INCLUDE_DIRS}"
      LIBRARY  "${${NAME}_LIBRARY}"
      LINK_LIBRARIES "${${NAME}_FOUND_LIBRARIES}")
  endif()
  #This is a macro, clear variables we don't to escape
  set(TPL_MODULE_NAME)
endmacro()

macro(kokkoskernels_export_imported_tpl NAME)
  cmake_parse_arguments(TPL "" "IMPORTED_NAME" "" ${ARGN})
  if(NOT TPL_IMPORTED_NAME)
    set(TPL_IMPORTED_NAME KokkosKernels::${NAME})
  endif()
  if(NOT KOKKOSKERNELS_HAS_TRILINOS)
    get_target_property(LIB_TYPE ${TPL_IMPORTED_NAME} TYPE)
    if(${LIB_TYPE} STREQUAL "INTERFACE_LIBRARY")
      # This is not an imported target
      # This an interface library that we created
      install(TARGETS ${NAME}
        EXPORT KokkosKernelsTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
    else()
      #make sure this also gets "exported" in the config file
      kokkoskernels_append_config_line("IF(NOT TARGET KokkosKernels::${NAME})")
      kokkoskernels_append_config_line("ADD_LIBRARY(KokkosKernels::${NAME} UNKNOWN IMPORTED)")
      kokkoskernels_append_config_line("SET_TARGET_PROPERTIES(KokkosKernels::${NAME} PROPERTIES")

      get_target_property(TPL_LIBRARY ${TPL_IMPORTED_NAME} IMPORTED_LOCATION)
      if(TPL_LIBRARY)
        kokkoskernels_append_config_line("IMPORTED_LOCATION ${TPL_LIBRARY}")
      endif()

      get_target_property(TPL_INCLUDES ${TPL_IMPORTED_NAME} INTERFACE_INCLUDE_DIRECTORIES)
      if(TPL_INCLUDES)
        # remove duplicates to prevent incorrect number of arguments to INTERFACE_INCLUDE_DIRECTORIES
        # see issue #2238
        list(REMOVE_DUPLICATES TPL_INCLUDES)
        kokkoskernels_append_config_line("INTERFACE_INCLUDE_DIRECTORIES ${TPL_INCLUDES}")
      endif()

      get_target_property(TPL_COMPILE_OPTIONS ${TPL_IMPORTED_NAME} INTERFACE_COMPILE_OPTIONS)
      if(TPL_COMPILE_OPTIONS)
        kokkoskernels_append_config_line("INTERFACE_COMPILE_OPTIONS ${TPL_COMPILE_OPTIONS}")
      endif()

      set(TPL_LINK_OPTIONS)
      get_target_property(TPL_LINK_OPTIONS ${TPL_IMPORTED_NAME} INTERFACE_LINK_OPTIONS)
      if(TPL_LINK_OPTIONS)
        kokkoskernels_append_config_line("INTERFACE_LINK_OPTIONS ${TPL_LINK_OPTIONS}")
      endif()

      get_target_property(TPL_LINK_LIBRARIES ${TPL_IMPORTED_NAME} INTERFACE_LINK_LIBRARIES)
      if(TPL_LINK_LIBRARIES)
        kokkoskernels_append_config_line("INTERFACE_LINK_LIBRARIES ${TPL_LINK_LIBRARIES}")
      endif()
      kokkoskernels_append_config_line(")")
      kokkoskernels_append_config_line("ENDIF()")
    endif()
  endif()
endmacro()

macro(kokkoskernels_import_tpl NAME)
  set(${NAME}_LIBRARIES "" CACHE STRING
    "Optional override for the libraries that comprise TPL ${NAME}. Default: None. Default common library names will be searched")
  set(${NAME}_LIBRARY_DIRS "" CACHE STRING
    "Optional override for the library directories that comprise TPL ${NAME}. Default: None. Default common library locations will be searched")
  set(${NAME}_INCLUDE_DIRS "" CACHE STRING
    "Optional override for the header directories that comprise TPL ${NAME}. Default: None. Default common header locations will be searched")

  cmake_parse_arguments(TPL "NO_EXPORT" "" "" ${ARGN})

  # Even though this policy gets set in the top-level CMakeLists.txt,
  # I have still been getting errors about ROOT variables being ignored
  # I'm not sure if this is a scope issue - but make sure
  # the policy is set before we do any find_package calls
  cmake_policy(SET CMP0074 NEW)

  if(KOKKOSKERNELS_ENABLE_TPL_${NAME})
    #Tack on a TPL here to make sure we avoid using anyone else's find
    find_package(TPL${NAME} REQUIRED MODULE)
    if(NOT TPL_${NAME}_IMPORTED_NAME)
      message(FATAL_ERROR "Find module did not produce valid IMPORTED_NAME for ${NAME}")
    endif()

    if(NOT TARGET ${TPL_${NAME}_IMPORTED_NAME})
      message(FATAL_ERROR "Find module succeeded for ${NAME}, but did not produce valid target ${TPL_${NAME}_IMPORTED_NAME}")
    endif()
    if(NOT TPL_NO_EXPORT)
      kokkoskernels_export_imported_tpl(${NAME} IMPORTED_NAME ${TPL_${NAME}_IMPORTED_NAME})
    endif()
  endif()
endmacro()

function(kokkoskernels_link_tpl TARGET)
  cmake_parse_arguments(TPL "PUBLIC;PRIVATE;INTERFACE" "IMPORTED_NAME" "" ${ARGN})
  #the name of the TPL
  set(TPL ${TPL_UNPARSED_ARGUMENTS})
  if(KOKKOSKERNELS_HAS_TRILINOS)
    #Do nothing, they will have already been linked
  else()
    if(NOT TPL_IMPORTED_NAME)
      set(TPL_IMPORTED_NAME KokkosKernels::${TPL})
    endif()
    if(KOKKOSKERNELS_ENABLE_TPL_${TPL})
      if(TPL_PUBLIC)
        target_link_libraries(${TARGET} PUBLIC ${TPL_IMPORTED_NAME})
      elseif(TPL_PRIVATE)
        target_link_libraries(${TARGET} PRIVATE ${TPL_IMPORTED_NAME})
      elseif(TPL_INTERFACE)
        target_link_libraries(${TARGET} INTERFACE ${TPL_IMPORTED_NAME})
      else()
        target_link_libraries(${TARGET} ${TPL_IMPORTED_NAME})
      endif()
    endif()
  endif()
endfunction()

kokkoskernels_add_tpl_option(BLAS OFF "Whether to enable BLAS")
#Default on if BLAS is enabled
kokkoskernels_add_tpl_option(
  LAPACK ${KokkosKernels_ENABLE_TPL_BLAS} "Whether to enable LAPACK"
  DEFAULT_DOCSTRING "ON if BLAS is enabled, otherwise OFF")
kokkoskernels_add_tpl_option(MKL OFF "Whether to enable MKL")
kokkoskernels_add_tpl_option(MAGMA OFF "Whether to enable MAGMA")
kokkoskernels_add_tpl_option(CBLAS OFF "Whether to enable CBLAS")
kokkoskernels_add_tpl_option(LAPACKE OFF "Whether to enable LAPACKE")
kokkoskernels_add_tpl_option(ARMPL OFF "Whether to enable ARMPL")
kokkoskernels_add_tpl_option(ACCELERATE   OFF  "Whether to enable ACCELERATE")

# Set F77_BLAS_MANGLE macro based on Fortran-C interface (unless already set
# by Trilinos or user)
if("${F77_BLAS_MANGLE}" STREQUAL "")
  if(KOKKOSKERNELS_ENABLE_TPL_BLAS
     OR KOKKOSKERNELS_ENABLE_TPL_LAPACK
     OR KOKKOSKERNELS_ENABLE_TPL_MKL
     OR KOKKOSKERNELS_ENABLE_TPL_MAGMA
     OR KOKKOSKERNELS_ENABLE_TPL_ARMPL)
    enable_language(C)
    enable_language(Fortran)
    include(FortranCInterface)
    if(FortranCInterface_GLOBAL_SUFFIX STREQUAL "")
      set(F77_BLAS_MANGLE "(name,NAME) ${FortranCInterface_GLOBAL_PREFIX}name")
    else()
      set(F77_BLAS_MANGLE "(name,NAME) ${FortranCInterface_GLOBAL_PREFIX}name ## ${FortranCInterface_GLOBAL_SUFFIX}")
    endif()
  endif()
endif()

kokkoskernels_add_option("NO_DEFAULT_CUDA_TPLS" OFF BOOL
  "Whether CUDA TPLs should be enabled by default. Default: OFF")
set(CUBLAS_DEFAULT ${KOKKOS_ENABLE_CUDA})
set(CUSPARSE_DEFAULT ${KOKKOS_ENABLE_CUDA})
set(CUSOLVER_DEFAULT ${KOKKOS_ENABLE_CUDA})
if(KOKKOSKERNELS_NO_DEFAULT_CUDA_TPLS)
  set(CUBLAS_DEFAULT OFF)
  set(CUSPARSE_DEFAULT OFF)
  set(CUSOLVER_DEFAULT OFF)
endif()
kokkoskernels_add_tpl_option(CUBLAS ${CUBLAS_DEFAULT} "Whether to enable CUBLAS"
  DEFAULT_DOCSTRING "ON if CUDA-enabled Kokkos, otherwise OFF")
kokkoskernels_add_tpl_option(CUSPARSE ${CUSPARSE_DEFAULT} "Whether to enable CUSPARSE"
  DEFAULT_DOCSTRING "ON if CUDA-enabled Kokkos, otherwise OFF")
kokkoskernels_add_tpl_option(CUSOLVER ${CUSOLVER_DEFAULT} "Whether to enable CUSOLVER"
  DEFAULT_DOCSTRING "ON if CUDA-enabled Kokkos, otherwise OFF")

kokkoskernels_add_option("NO_DEFAULT_ROCM_TPLS" OFF BOOL
  "Whether ROCM TPLs should be enabled by default. Default: OFF")
# Unlike CUDA, ROCm does not automatically install these TPLs
set(ROCBLAS_DEFAULT OFF)
set(ROCSPARSE_DEFAULT OFF)
set(ROCSOLVER_DEFAULT OFF)
# Since the default is OFF we do not really need this piece of logic here.
# IF(KOKKOSKERNELS_NO_DEFAULT_ROCM_TPLS)
#   SET(ROCBLAS_DEFAULT   OFF)
#   SET(ROCSPARSE_DEFAULT OFF)
# ENDIF()
kokkoskernels_add_tpl_option(ROCBLAS ${ROCBLAS_DEFAULT} "Whether to enable ROCBLAS"
  DEFAULT_DOCSTRING "OFF even if HIP-enabled Kokkos")
kokkoskernels_add_tpl_option(ROCSPARSE ${ROCSPARSE_DEFAULT} "Whether to enable ROCSPARSE"
  DEFAULT_DOCSTRING "OFF even if HIP-enabled Kokkos")
kokkoskernels_add_tpl_option(ROCSOLVER ${ROCSOLVER_DEFAULT} "Whether to enable ROCSOLVER"
  DEFAULT_DOCSTRING "OFF even if HIP-enabled Kokkos")

if(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
  if(F77_BLAS_MANGLE STREQUAL "(name,NAME) name ## _")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DADD_ -fopenmp -lgfortran")
  elseif(F77_BLAS_MANGLE STREQUAL "(name,NAME) NAME")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUPCASE -fopenmp -lgfortran")
  elseif(F77_BLAS_MANGLE STREQUAL "(name,NAME) name")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOCHANGE -fopenmp -lgfortran")
  else()
    message(FATAL_ERROR
      "F77_BLAS_MANGLE ${F77_BLAS_MANGLE} detected while MAGMA only accepts Fortran mangling that is one of single underscore (-DADD_), uppercase (-DUPCASE), and no change (-DNOCHANGE)")
  endif()
  list(APPEND TPL_LIST "MAGMA")
endif()

kokkoskernels_add_tpl_option(SUPERLU OFF "Whether to enable SUPERLU")
kokkoskernels_add_tpl_option(CHOLMOD OFF "Whether to enable CHOLMOD")
kokkoskernels_add_tpl_option(METIS OFF "Whether to enable METIS")

# We need to do all the import work
if(NOT KOKKOSKERNELS_HAS_TRILINOS)
  if(KokkosKernels_ENABLE_SuperLU5_API)
    set(HAVE_KOKKOSKERNELS_SUPERLU5_API TRUE)
  endif()

  kokkoskernels_import_tpl(BLAS)
  kokkoskernels_import_tpl(LAPACK)
  kokkoskernels_import_tpl(MKL)
  kokkoskernels_import_tpl(CUBLAS)
  kokkoskernels_import_tpl(CUSPARSE)
  kokkoskernels_import_tpl(CUSOLVER)
  kokkoskernels_import_tpl(CBLAS)
  kokkoskernels_import_tpl(LAPACKE)
  kokkoskernels_import_tpl(CHOLMOD)
  kokkoskernels_import_tpl(SUPERLU)
  kokkoskernels_import_tpl(METIS)
  kokkoskernels_import_tpl(ARMPL)
  kokkoskernels_import_tpl(MAGMA)
  kokkoskernels_import_tpl(ROCBLAS)
  kokkoskernels_import_tpl(ROCSPARSE)
  kokkoskernels_import_tpl(ROCSOLVER)
  kokkoskernels_import_tpl(ACCELERATE)
else()
  if(Trilinos_ENABLE_SuperLU5_API)
    set(HAVE_KOKKOSKERNELS_SUPERLU5_API TRUE)
  endif()
endif()

#Convert list to newlines (which CMake doesn't always like in cache variables)
string(REPLACE ";" "\n" KOKKOSKERNELS_TPL_EXPORT_TEMP "${KOKKOSKERNELS_TPL_EXPORTS}")
#Convert to a regular variable
unset(KOKKOSKERNELS_TPL_EXPORTS CACHE)
set(KOKKOSKERNELS_TPL_EXPORTS ${KOKKOSKERNELS_TPL_EXPORT_TEMP})
