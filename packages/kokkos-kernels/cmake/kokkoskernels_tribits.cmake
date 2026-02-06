include(CMakeParseArguments)
include(CTest)

if (KOKKOSKERNELS_HAS_TRILINOS)
  include(TribitsETISupport)
endif()

function(verify_empty CONTEXT)
  if(${ARGN})
    message(FATAL_ERROR "Kokkos does not support all of Tribits. Unhandled arguments in ${CONTEXT}:\n${ARGN}")
  endif()
endfunction()

#MESSAGE(STATUS "The project name is: ${PROJECT_NAME}")

macro(kokkoskernels_package_postprocess)
  if (KOKKOSKERNELS_HAS_TRILINOS)
    tribits_package_postprocess()
  else()
    include(CMakePackageConfigHelpers)
    configure_package_config_file(cmake/KokkosKernelsConfig.cmake.in
      "${KokkosKernels_BINARY_DIR}/KokkosKernelsConfig.cmake"
      INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/KokkosKernels)
    write_basic_package_version_file("${KokkosKernels_BINARY_DIR}/KokkosKernelsConfigVersion.cmake"
      VERSION "${KokkosKernels_VERSION_MAJOR}.${KokkosKernels_VERSION_MINOR}.${KokkosKernels_VERSION_PATCH}"
      COMPATIBILITY AnyNewerVersion)

    install(FILES
      "${KokkosKernels_BINARY_DIR}/KokkosKernelsConfig.cmake"
      "${KokkosKernels_BINARY_DIR}/KokkosKernelsConfigVersion.cmake"
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/KokkosKernels)

    install(EXPORT KokkosKernelsTargets
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/KokkosKernels
      NAMESPACE Kokkos::)
  endif()
endmacro(kokkoskernels_package_postprocess)

macro(kokkoskernels_subpackage NAME)
if (KOKKOSKERNELS_HAS_TRILINOS)
  tribits_subpackage(${NAME})
else()
  set(PACKAGE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  set(PARENT_PACKAGE_NAME ${PACKAGE_NAME})
  set(PACKAGE_NAME ${PACKAGE_NAME}${NAME})
  string(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UC)
  set(${PACKAGE_NAME}_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()
endmacro(kokkoskernels_subpackage)

macro(kokkoskernels_subpackage_postprocess)
if (KOKKOSKERNELS_HAS_TRILINOS)
  tribits_subpackage_postprocess()
endif()
endmacro(kokkoskernels_subpackage_postprocess)

macro(kokkoskernels_process_subpackages)
if (kokkoskernels_has_trilinos)
  tribits_process_subpackages()
endif()
endmacro(kokkoskernels_process_subpackages)

macro(kokkoskernels_package)
if (KOKKOSKERNELS_HAS_TRILINOS)
  tribits_package(KokkosKernels)
else()
  set(PACKAGE_NAME KokkosKernels)
  set(PACKAGE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  string(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UC)
  set(${PACKAGE_NAME}_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()
endmacro(kokkoskernels_package)

function(kokkoskernels_internal_add_library LIBRARY_NAME)
  cmake_parse_arguments(PARSE
    "STATIC;SHARED"
    ""
    "HEADERS;SOURCES"
    ${ARGN})

  if (PARSE_HEADERS)
    list(REMOVE_DUPLICATES PARSE_HEADERS)
  endif()
  if (PARSE_SOURCES)
    list(REMOVE_DUPLICATES PARSE_SOURCES)
  endif()
  if (Kokkos_COMPILE_LANGUAGE)
    foreach(source ${PARSE_SOURCES})
      set_source_files_properties(${source} PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE})
    endforeach()
  endif()

  add_library(
    ${LIBRARY_NAME}
    ${PARSE_HEADERS}
    ${PARSE_SOURCES}
  )
  add_library(Kokkos::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

  install(
    TARGETS ${LIBRARY_NAME}
    EXPORT KokkosKernelsTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

  install(
    FILES  ${PARSE_HEADERS}
    DESTINATION ${KOKKOSKERNELS_HEADER_INSTALL_DIR}
    COMPONENT ${PACKAGE_NAME}
  )

  install(
    FILES  ${PARSE_HEADERS}
    DESTINATION ${KOKKOSKERNELS_HEADER_INSTALL_DIR}
  )

endfunction(kokkoskernels_internal_add_library LIBRARY_NAME)

function(kokkoskernels_add_library LIBRARY_NAME)
  if (KOKKOSKERNELS_HAS_TRILINOS)
    tribits_add_library(${LIBRARY_NAME} ${ARGN})
  else()
    kokkoskernels_internal_add_library(${LIBRARY_NAME} ${ARGN})
  endif()
endfunction()

function(kokkoskernels_add_executable EXE_NAME)
  cmake_parse_arguments(PARSE
    ""
    ""
    "SOURCES;COMPONENTS;TESTONLYLIBS"
    ${ARGN})
  verify_empty(KOKKOSKERNELS_ADD_EXECUTABLE ${PARSE_UNPARSED_ARGUMENTS})

  kokkoskernels_is_enabled(
    COMPONENTS ${PARSE_COMPONENTS}
    OUTPUT_VARIABLE IS_ENABLED
  )

  if (IS_ENABLED)
    if (KOKKOSKERNELS_HAS_TRILINOS)
      tribits_add_executable(${EXE_NAME}
        SOURCES ${PARSE_SOURCES}
        TESTONLYLIBS ${PARSE_TESTONLYLIBS})
    else()
      # Set the correct CMake language on all source files for this exe
      if (Kokkos_COMPILE_LANGUAGE)
        foreach(source ${PARSE_SOURCES})
          set_source_files_properties(${source} PROPERTIES LANGUAGE ${Kokkos_COMPILE_LANGUAGE})
        endforeach()
      endif()
      add_executable(${EXE_NAME} ${PARSE_SOURCES})
      #AJP, BMK altered:
      if (KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE)
        target_link_libraries(${EXE_NAME} PRIVATE common ${PARSE_TESTONLYLIBS})
      endif()

      if (PARSE_TESTONLYLIBS)
        target_link_libraries(${EXE_NAME} PRIVATE Kokkos::kokkoskernels ${PARSE_TESTONLYLIBS})
      else ()
        target_link_libraries(${EXE_NAME} PRIVATE Kokkos::kokkoskernels)
      endif()
    endif()
  else()
    message(STATUS "Skipping executable ${EXE_NAME} because not all necessary components enabled")
  endif()
endfunction()

function(kokkoskernels_add_unit_test ROOT_NAME)
  kokkoskernels_add_executable_and_test(
    ${ROOT_NAME}
    TESTONLYLIBS kokkoskernels_gtest
    ${ARGN}
  )
endfunction()

function(KOKKOSKERNELS_IS_ENABLED)
  cmake_parse_arguments(PARSE
    ""
    "OUTPUT_VARIABLE"
    "COMPONENTS"
    ${ARGN})

  if (KOKKOSKERNELS_ENABLED_COMPONENTS STREQUAL "ALL")
    set(${PARSE_OUTPUT_VARIABLE} TRUE PARENT_SCOPE)
  elseif(PARSE_COMPONENTS)
    set(ENABLED TRUE)
    foreach(comp ${PARSE_COMPONENTS})
      string(TOUPPER ${comp} COMP_UC)
      # make sure this is in the list of enabled components
      if (NOT "${COMP_UC}" IN_LIST KOKKOSKERNELS_ENABLED_COMPONENTS)
        # if not in the list, one or more components is missing
        set(ENABLED FALSE)
      endif()
    endforeach()
    set(${PARSE_OUTPUT_VARIABLE} ${ENABLED} PARENT_SCOPE)
  else()
    # we did not enable all components and no components
    # were given as part of this - we consider this enabled
    set(${PARSE_OUTPUT_VARIABLE} TRUE PARENT_SCOPE)
  endif()
endfunction()

function(kokkoskernels_add_executable_and_test ROOT_NAME)

  cmake_parse_arguments(PARSE
    ""
    ""
    "SOURCES;CATEGORIES;COMPONENTS;TESTONLYLIBS"
    ${ARGN})

  verify_empty(KOKKOSKERNELS_ADD_EXECUTABLE_AND_RUN_VERIFY ${PARSE_UNPARSED_ARGUMENTS})

  kokkoskernels_is_enabled(
    COMPONENTS ${PARSE_COMPONENTS}
    OUTPUT_VARIABLE IS_ENABLED
  )

  if (IS_ENABLED)
    if (KOKKOSKERNELS_HAS_TRILINOS)
      tribits_add_executable_and_test(
        ${ROOT_NAME}
        SOURCES ${PARSE_SOURCES}
        CATEGORIES ${PARSE_CATEGORIES}
        TESTONLYLIBS ${PARSE_TESTONLYLIBS}
        NUM_MPI_PROCS 1
        COMM serial mpi
        )
    else()
      set(EXE_NAME ${PACKAGE_NAME}_${ROOT_NAME})
      kokkoskernels_add_executable(${EXE_NAME} SOURCES ${PARSE_SOURCES})
      if (PARSE_TESTONLYLIBS)
        target_link_libraries(${EXE_NAME} PRIVATE ${PARSE_TESTONLYLIBS})
      endif()
      kokkoskernels_add_test(NAME ${ROOT_NAME} EXE ${EXE_NAME})
    endif()
  else()
    message(STATUS "Skipping executable/test ${ROOT_NAME} because not all necessary components enabled")
  endif()

endfunction()

macro(add_component_subdirectory SUBDIR)
  kokkoskernels_is_enabled(
    COMPONENTS ${SUBDIR}
    OUTPUT_VARIABLE COMP_SUBDIR_ENABLED
  )
  if (COMP_SUBDIR_ENABLED)
    add_subdirectory(${SUBDIR})
  else()
    message(STATUS "Skipping subdirectory ${SUBDIR} because component is not enabled")
  endif()
  unset(COMP_SUBDIR_ENABLED)
endmacro()
