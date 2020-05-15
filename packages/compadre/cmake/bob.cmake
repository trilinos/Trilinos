# MIT License
# 
# Copyright (c) 2016 Dan Ibanez
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

function(bob_always_full_rpath)
  # CMake RPATH "always full" configuration, see:
  # https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH
  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH False PARENT_SCOPE)
  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH False PARENT_SCOPE)
  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
       "${CMAKE_INSTALL_LIBDIR}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_LIBDIR}" PARENT_SCOPE)
  endif()
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True PARENT_SCOPE)
endfunction(bob_always_full_rpath)

function(bob_cmake_arg2 var type default)
if (NOT ${var} STREQUAL "${default}")
  if (${PROJECT_NAME}_CMAKE_ARGS)
    set(sep " ")
  else()
    set(sep "")
  endif()
  set(${PROJECT_NAME}_CMAKE_ARGS "${${PROJECT_NAME}_CMAKE_ARGS}${sep}-D${var}:${type}=\"${${var}}\"" CACHE STRING
      "CMake arguments that would replicate this configuration" FORCE)
endif()
endfunction()

function(bob_cmake_arg var type default)
message(STATUS "${var}: ${${var}}")
bob_cmake_arg2("${var}" "${type}" "${default}")
endfunction()

function(bob_option var desc default)
option(${var} "${desc}" "${default}")
bob_cmake_arg(${var} BOOL "${default}")
endfunction()

function(bob_input var default type desc)
set(${var} "${default}" CACHE ${type} "${desc}")
bob_cmake_arg(${var} ${type} "${default}")
endfunction()

macro(bob_begin_package)
  set(${PROJECT_NAME}_CMAKE_ARGS "" CACHE STRING
      "CMake arguments that would replicate this configuration" FORCE)
  message(STATUS "CMAKE_VERSION: ${CMAKE_VERSION}")
  if (${PROJECT_NAME}_VERSION)
    message(STATUS "${PROJECT_NAME}_VERSION: ${${PROJECT_NAME}_VERSION}")
  endif()
  option(USE_XSDK_DEFAULTS "enable the XDSK v0.3.0 default configuration" OFF)
  bob_cmake_arg(USE_XSDK_DEFAULTS BOOL OFF)
  #try to force BUILD_TESTING to be OFF by default
  set(BUILD_TESTING OFF CACHE BOOL "Build and run tests")
  include(CTest)
  enable_testing()
  option(BUILD_SHARED_LIBS "Build shared libraries" ON)
  #If not building shared libs, then prefer static
  #dependency libs
  if(NOT BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".dylib")
  endif()
  if (USE_XSDK_DEFAULTS)
    if (NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE "Debug")
    endif()
    bob_cmake_arg(CMAKE_BUILD_TYPE STRING "")
  endif()
  bob_always_full_rpath()
  bob_cmake_arg(BUILD_TESTING BOOL OFF)
  bob_cmake_arg(BUILD_SHARED_LIBS BOOL ON)
  bob_cmake_arg(CMAKE_INSTALL_PREFIX BOOL "")
  option(${PROJECT_NAME}_NORMAL_CXX_FLAGS "Allow CMAKE_CXX_FLAGS to follow \"normal\" CMake behavior" ${USE_XSDK_DEFAULTS})
endmacro(bob_begin_package)

function(bob_get_commit)
  execute_process(COMMAND git rev-parse HEAD
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE NO_SHA1
      OUTPUT_VARIABLE SHA1
      ERROR_VARIABLE SHA1_ERROR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  if(NO_SHA1)
    message(WARNING "bob_get_commit: no Git hash!\n" ${SHA1_ERROR})
  else()
    set(${PROJECT_NAME}_COMMIT "${SHA1}" PARENT_SCOPE)
  endif()
endfunction(bob_get_commit)

function(bob_get_semver)
  execute_process(COMMAND git describe --exact-match HEAD
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE NOT_TAG
      OUTPUT_VARIABLE TAG_NAME
      ERROR_VARIABLE TAG_ERROR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  if(NOT_TAG)
    if(${PROJECT_NAME}_VERSION)
      set(SEMVER ${${PROJECT_NAME}_VERSION})
      execute_process(COMMAND git log -1 --format=%h
          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
          RESULT_VARIABLE NO_SHA1
          OUTPUT_VARIABLE SHORT_SHA1
          ERROR_VARIABLE SHA1_ERROR
          OUTPUT_STRIP_TRAILING_WHITESPACE
          )
      if(NO_SHA1)
        message(WARNING "bob_get_semver no Git hash!\n" ${SHA1_ERROR})
      else()
        set(SEMVER "${SEMVER}-sha.${SHORT_SHA1}")
      endif()
    else()
      message(FATAL_ERROR "bob_get_semver needs either ${PROJECT_NAME}_VERSION or a Git tag\n" ${TAG_ERROR})
    endif()
  else()
    if(TAG_NAME MATCHES "^v([0-9]+[.])?([0-9]+[.])?([0-9]+)$")
      string(SUBSTRING "${TAG_NAME}" 1 -1 SEMVER)
      if(${PROJECT_NAME}_VERSION AND (NOT (SEMVER VERSION_EQUAL ${PROJECT_NAME}_VERSION)))
        message(FATAL_ERROR "bob_get_semver: tag is ${TAG_NAME} but ${PROJECT_NAME}_VERSION=${${PROJECT_NAME}_VERSION} !")
      endif()
    else()
      if(${PROJECT_NAME}_VERSION)
        set(SEMVER "${${PROJECT_NAME}_VERSION}-tag.${TAG_NAME}")
      else()
        message(FATAL_ERROR "bob_get_semver needs either ${PROJECT_NAME}_VERSION or a Git tag of the form v1.2.3")
      endif()
    endif()
  endif()
  if(${PROJECT_NAME}_KEY_BOOLS)
    set(SEMVER "${SEMVER}+")
    foreach(KEY_BOOL IN LISTS ${PROJECT_NAME}_KEY_BOOLS)
      if(${KEY_BOOL})
        set(SEMVER "${SEMVER}1")
      else()
        set(SEMVER "${SEMVER}0")
      endif()
    endforeach()
  endif()
  set(${PROJECT_NAME}_SEMVER "${SEMVER}" PARENT_SCOPE)
  message(STATUS "${PROJECT_NAME}_SEMVER = ${SEMVER}")
endfunction(bob_get_semver)

function(bob_begin_cxx_flags)
  if (${PROJECT_NAME}_NORMAL_CXX_FLAGS)
    set(BOB_CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" PARENT_SCOPE)
    bob_cmake_arg2(CMAKE_CXX_FLAGS STRING "")
  else()
    if(CMAKE_BUILD_TYPE)
      message(FATAL_ERROR "can't set CMAKE_BUILD_TYPE and use bob_*_cxx_flags")
    endif()
    option(${PROJECT_NAME}_CXX_OPTIMIZE "Compile C++ with optimization" ON)
    option(${PROJECT_NAME}_CXX_SYMBOLS "Compile C++ with debug symbols" ON)
    option(${PROJECT_NAME}_CXX_WARNINGS "Compile C++ with warnings" ON)
    bob_cmake_arg(${PROJECT_NAME}_CXX_OPTIMIZE BOOL ON)
    bob_cmake_arg(${PROJECT_NAME}_CXX_SYMBOLS BOOL ON)
    #CDash's simple output parser interprets the variable name WARNINGS as a warning...
    message(STATUS "${PROJECT_NAME}_CXX_W**NINGS: ${${PROJECT_NAME}_CXX_WARNINGS}")
    bob_cmake_arg2(${PROJECT_NAME}_CXX_WARNINGS BOOL ON)
    set(FLAGS "")
    if(${PROJECT_NAME}_CXX_OPTIMIZE)
      set(FLAGS "${FLAGS} -O3")
    else()
      set(FLAGS "${FLAGS} -O0")
    endif()
    if(${PROJECT_NAME}_CXX_SYMBOLS)
      set(FLAGS "${FLAGS} -g")
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      if (${PROJECT_NAME}_CXX_WARNINGS)
        #set(FLAGS "${FLAGS} -Werror -Weverything")
        set(FLAGS "${FLAGS} -Wno-padded")
        set(FLAGS "${FLAGS} -Wno-float-equal")
        set(FLAGS "${FLAGS} -Wno-weak-template-vtables")
      endif()
      if(APPLE)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "9.0.0")
          set(FLAGS "${FLAGS} -fcomment-block-commands=file")
        endif()
      else()
        if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "5.0.0")
          set(FLAGS "${FLAGS} -fcomment-block-commands=file")
        endif()
        if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0.0")
          set(FLAGS "${FLAGS} -fcomment-block-commands=file")
        endif()
      endif()
    elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
      if (${PROJECT_NAME}_CXX_WARNINGS)
        #set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
        #set(FLAGS "${FLAGS} -Wdouble-promotion -Wshadow -Wformat=2")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "6.0")
          set(FLAGS "${FLAGS} -Wduplicated-cond -Wnull-dereference")
        endif()
        if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "7.0")
          set(FLAGS "${FLAGS} -Wlogical-op")
          set(FLAGS "${FLAGS} -Wrestrict")
          if(NOT Omega_h_USE_DOLFIN)
            set(FLAGS "${FLAGS} -Wduplicated-branches")
          endif()
        endif()
      endif()
    elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    else()
      message(WARNING "Unexpected compiler type ${CMAKE_CXX_COMPILER_ID}")
    endif()
    set(CMAKE_CXX_FLAGS "${FLAGS}" PARENT_SCOPE)
  endif()
endfunction(bob_begin_cxx_flags)

function(bob_cxx11_flags)
  set(FLAGS "${CMAKE_CXX_FLAGS}")
  set(FLAGS "${FLAGS} --std=c++11")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if (${PROJECT_NAME}_CXX_WARNINGS)
      set(FLAGS "${FLAGS} -Wno-c++98-compat-pedantic -Wno-c++98-compat")
    endif()
  endif()
  set(CMAKE_CXX_FLAGS "${FLAGS}" PARENT_SCOPE)
endfunction(bob_cxx11_flags)

function(bob_end_cxx_flags)
  if (${PROJECT_NAME}_NORMAL_CXX_FLAGS)
    message(STATUS "CMAKE_CXX_FLAGS: ${BOB_CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${BOB_CMAKE_CXX_FLAGS}" PARENT_SCOPE)
  else()
    set(${PROJECT_NAME}_CXX_FLAGS "" CACHE STRING "Override all C++ compiler flags")
    bob_cmake_arg(${PROJECT_NAME}_CXX_FLAGS STRING "")
    set(${PROJECT_NAME}_EXTRA_CXX_FLAGS "" CACHE STRING "Extra C++ compiler flags")
    bob_cmake_arg(${PROJECT_NAME}_EXTRA_CXX_FLAGS STRING "")
    set(FLAGS "${CMAKE_CXX_FLAGS}")
    if(${PROJECT_NAME}_CXX_FLAGS)
      set(FLAGS "${${PROJECT_NAME}_CXX_FLAGS}")
    else()
      set(FLAGS "${FLAGS} ${${PROJECT_NAME}_EXTRA_CXX_FLAGS}")
    endif()
    message(STATUS "CMAKE_CXX_FLAGS: ${FLAGS}")
    set(CMAKE_CXX_FLAGS "${FLAGS}" PARENT_SCOPE)
  endif()
endfunction(bob_end_cxx_flags)

macro(bob_add_dependency)
  set(options PUBLIC PRIVATE)
  set(oneValueArgs NAME)
  set(multiValueArgs COMPONENTS TARGETS INCLUDE_DIR_VARS LIBRARY_VARS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if (NOT ARG_NAME)
    message(FATAL_ERROR "bob_add_dependency: no NAME argument given")
  endif()
  if (ARG_PUBLIC AND ARG_PRIVATE)
    message(FATAL_ERROR "bob_add_dependency: can't specify both PUBLIC and PRIVATE")
  endif()
  if (ARG_COMPONENTS)
    set(ARG_COMPONENTS COMPONENTS ${ARG_COMPONENTS})
  endif()
  if (USE_XSDK_DEFAULTS)
    option(TPL_ENABLE_${ARG_NAME} "Whether to use ${ARG_NAME}"
           "${${PROJECT_NAME}_USE_${ARG_NAME}_DEFAULT}")
    bob_cmake_arg(TPL_ENABLE_${ARG_NAME} BOOL "${${PROJECT_NAME}_USE_${ARG_NAME}_DEFAULT}")
    set(${PROJECT_NAME}_USE_${ARG_NAME} "${TPL_ENABLE_${ARG_NAME}}")
    if(TPL_ENABLE_${ARG_NAME})
      set(TPL_${ARG_NAME}_LIBRARIES "" CACHE STRING "${ARG_NAME} libraries")
      bob_cmake_arg(TPL_${ARG_NAME}_LIBRARIES STRING "")
      set(TPL_${ARG_NAME}_INCLUDE_DIRS "" CACHE STRING "${ARG_NAME} include directories")
      bob_cmake_arg(TPL_${ARG_NAME}_INCLUDE_DIRS STRING "")
      set(tgt "${PROJECT_NAME}-${ARG_NAME}")
      add_library(${tgt} INTERFACE)
      target_include_directories(${tgt} INTERFACE "${TPL_${ARG_NAME}_INCLUDE_DIRS}")
      target_link_libraries(${tgt} INTERFACE "${TPL_${ARG_NAME}_LIBRARIES}")
    endif()
  else()
    option(${PROJECT_NAME}_USE_${ARG_NAME} "Whether to use ${ARG_NAME}"
           ${${PROJECT_NAME}_USE_${ARG_NAME}_DEFAULT})
    bob_cmake_arg(${PROJECT_NAME}_USE_${ARG_NAME} BOOL "${${PROJECT_NAME}_USE_${ARG_NAME}_DEFAULT}")
    if(${PROJECT_NAME}_USE_${ARG_NAME})
      set(${ARG_NAME}_PREFIX "${${ARG_NAME}_PREFIX_DEFAULT}"
          CACHE PATH "${ARG_NAME} install directory")
      bob_cmake_arg(${ARG_NAME}_PREFIX PATH "${${ARG_NAME}_PREFIX_DEFAULT}")
      if (${ARG_NAME}_PREFIX)
        #if ${ARG_NAME}_PREFIX is set, don't find it anywhere else:
        set(ARG_PREFIX PATHS "${${ARG_NAME}_PREFIX}" NO_DEFAULT_PATH)
      else()
        #allow CMake to search other prefixes if ${ARG_NAME}_PREFIX is not set
        set(ARG_PREFIX)
      endif()
      set(${ARG_NAME}_find_package_args
          "${${ARG_NAME}_REQUIRED_VERSION}"
          ${ARG_COMPONENTS}
          ${ARG_PREFIX})
      find_package(${ARG_NAME} ${${ARG_NAME}_find_package_args} REQUIRED)
      if(${ARG_NAME}_CONFIG)
        message(STATUS "${ARG_NAME}_CONFIG: ${${ARG_NAME}_CONFIG}")
      endif()
      if(${ARG_NAME}_VERSION)
        message(STATUS "${ARG_NAME}_VERSION: ${${ARG_NAME}_VERSION}")
      endif()
      set(tgt "${PROJECT_NAME}-${ARG_NAME}")
      add_library(${tgt} INTERFACE)
      if (ARG_TARGETS)
        target_link_libraries(${tgt} INTERFACE ${ARG_TARGETS})
      endif()
      if (ARG_LIBRARY_VARS)
        foreach (library_var IN LISTS ARG_LIBRARY_VARS)
          target_link_libraries(${tgt} INTERFACE ${${library_var}})
        endforeach()
      endif()
      if (ARG_INCLUDE_DIR_VARS)
        foreach (include_dir_var IN LISTS ARG_INCLUDE_DIR_VARS)
          foreach (include_dir IN LISTS ${include_dir_var})
            get_filename_component(abs_include_dir "${include_dir}" ABSOLUTE)
            target_include_directories(${tgt} INTERFACE "${abs_include_dir}")
          endforeach()
        endforeach()
      endif()
      install(TARGETS ${tgt} EXPORT
          ${tgt}-target
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR})
      install(EXPORT ${tgt}-target
              DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
      set(${PROJECT_NAME}_EXPORTED_TARGETS
          ${${PROJECT_NAME}_EXPORTED_TARGETS}
          ${tgt})
      if (ARG_PUBLIC)
        set(${PROJECT_NAME}_DEPS ${${PROJECT_NAME}_DEPS} ${ARG_NAME})
      endif()
    endif()
  endif()
endmacro(bob_add_dependency)

function(bob_link_dependency tgt type dep)
  if (${PROJECT_NAME}_USE_${dep})
    target_link_libraries(${tgt} ${type} ${PROJECT_NAME}-${dep})
  endif()
endfunction(bob_link_dependency)

macro(bob_private_dep pkg_name)
  bob_add_dependency(PRIVATE NAME "${pkg_name}")
endmacro(bob_private_dep)

macro(bob_public_dep pkg_name)
  bob_add_dependency(PUBLIC NAME "${pkg_name}")
endmacro(bob_public_dep)

function(bob_target_includes lib_name)
  #find local headers even with #include <>
  target_include_directories(${lib_name}
      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
  #find generated configuration headers
  target_include_directories(${lib_name}
      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
endfunction(bob_target_includes)

function(bob_library_includes lib_name)
  bob_target_includes("${lib_name}")
  #ensure downstream users include installed headers
  target_include_directories(${lib_name} INTERFACE $<INSTALL_INTERFACE:include>)
endfunction(bob_library_includes)

function(bob_export_target tgt_name)
  get_target_property(tgt_type "${tgt_name}" TYPE)
  if (${tgt_type} STREQUAL "EXECUTABLE")
    install(TARGETS ${tgt_name} DESTINATION bin)
  else()
    if (USE_XSDK_DEFAULTS)
      install(TARGETS ${tgt_name} DESTINATION ${CMAKE_INSTALL_LIBDIR})
    else()
      install(TARGETS ${tgt_name} EXPORT ${tgt_name}-target DESTINATION ${CMAKE_INSTALL_LIBDIR})
      install(EXPORT ${tgt_name}-target NAMESPACE ${PROJECT_NAME}::
              DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
      set(${PROJECT_NAME}_EXPORTED_TARGETS
          ${${PROJECT_NAME}_EXPORTED_TARGETS} ${tgt_name} PARENT_SCOPE)
    endif()
  endif()
endfunction(bob_export_target)

macro(bob_end_subdir)
  set(${PROJECT_NAME}_EXPORTED_TARGETS
      ${${PROJECT_NAME}_EXPORTED_TARGETS} PARENT_SCOPE)
  set(${PROJECT_NAME}_DEPS ${${PROJECT_NAME}_DEPS} PARENT_SCOPE)
  set(${PROJECT_NAME}_DEP_PREFIXES ${${PROJECT_NAME}_DEP_PREFIXES} PARENT_SCOPE)
endmacro(bob_end_subdir)

function(bob_config_header HEADER_PATH)
  if (ARGC GREATER 1)
    set(MODULE_NAME "${ARGV1}")
  else()
    set(MODULE_NAME ${PROJECT_NAME})
  endif()
  get_filename_component(HEADER_NAME "${HEADER_PATH}" NAME)
  string(REPLACE "." "_" INCLUDE_GUARD "${HEADER_NAME}")
  string(TOUPPER "${INCLUDE_GUARD}" INCLUDE_GUARD)
  set(HEADER_CONTENT
"#ifndef ${INCLUDE_GUARD}
#define ${INCLUDE_GUARD}
")
  if (${MODULE_NAME}_KEY_BOOLS)
    foreach(KEY_BOOL IN LISTS ${MODULE_NAME}_KEY_BOOLS)
      if (${KEY_BOOL})
        string(TOUPPER "${KEY_BOOL}" MACRO_NAME)
        set(HEADER_CONTENT
"${HEADER_CONTENT}
#define ${MACRO_NAME}")
      endif()
    endforeach()
  endif()
  if (${MODULE_NAME}_KEY_INTS)
    foreach(KEY_INT IN LISTS ${MODULE_NAME}_KEY_INTS)
      string(TOUPPER "${KEY_INT}" MACRO_NAME)
      set(HEADER_CONTENT
"${HEADER_CONTENT}
#define ${MACRO_NAME} ${${KEY_INT}}")
    endforeach()
  endif()
  if (${MODULE_NAME}_KEY_STRINGS)
    foreach(KEY_STRING IN LISTS ${MODULE_NAME}_KEY_STRINGS)
      string(TOUPPER "${KEY_STRING}" MACRO_NAME)
      set(val "${${KEY_STRING}}")
#escape escapes
      string(REPLACE "\\" "\\\\" val "${val}")
#escape quotes
      string(REPLACE "\"" "\\\"" val "${val}")
      set(HEADER_CONTENT
"${HEADER_CONTENT}
#define ${MACRO_NAME} \"${val}\"")
    endforeach()
  endif()
  set(HEADER_CONTENT
"${HEADER_CONTENT}

#endif
")
  file(WRITE "${HEADER_PATH}" "${HEADER_CONTENT}")
endfunction()

function(bob_get_link_libs tgt var trilinos_components_list)
  get_target_property(tgt_type "${tgt}" TYPE)
  set(sublibs)
  if (NOT tgt_type STREQUAL "INTERFACE_LIBRARY")
    get_target_property(tgt_libs "${tgt}" LINK_LIBRARIES)
    if (tgt_libs)
      set(sublibs ${sublibs} ${tgt_libs})
    endif()
  endif()
  get_target_property(tgt_iface_libs "${tgt}" INTERFACE_LINK_LIBRARIES)
  if (tgt_iface_libs)
    set(sublibs ${sublibs} ${tgt_iface_libs})
  endif()
  set(link_libs)

  # search to see if it is a Trilinos package
  list (FIND trilinos_components_list "${tgt}" _index)
  if (${_index} GREATER -1)
    SET(found_trilinos_component ON)
  else()
    SET(found_trilinos_component OFF)
  endif()
  
  foreach(lib IN LISTS sublibs)
    # do not dig deeper if a trilinos package
    if ((TARGET ${lib}) AND NOT(found_trilinos_component))
      get_target_property(subtgt_type "${lib}" TYPE)
      #if (subtgt_type MATCHES "STATIC_LIBRARY|SHARED_LIBRARY")
      #  get_target_property(sublibtgt_loc "${lib}" LOCATION)
      #  if (sublibtgt_loc)
      #    set(link_libs ${link_libs} ${sublibtgt_loc})
      #  endif()
      #endif()
      #if (subtgt_type MATCHES "UNKNOWN_LIBRARY")
      #  foreach(prop in ITEMS IMPORTED_LOCATION IMPORTED_LOCATION_RELEASE IMPORTED_LOCATION_DEBUG)
      #    get_target_property(sublibtgt_import_loc "${lib}" ${prop})
      #    if (sublibtgt_import_loc)
      #      set(link_libs ${link_libs} ${sublibtgt_import_loc})
      #    endif()
      #  endforeach()
      #endif()
      bob_get_link_libs(${lib} subtgt_link_libs "${trilinos_components_list}")
      set(link_libs ${link_libs} ${subtgt_link_libs})
    else()
      set(link_libs ${link_libs} ${lib})
    endif()
  endforeach()

  if (link_libs)
    list(REVERSE link_libs)
    list(REMOVE_DUPLICATES link_libs)
    list(REVERSE link_libs)
  endif()
  set(${var} ${link_libs} PARENT_SCOPE)
endfunction()

function(bob_install_provenance)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_cmake_args.txt
       "${${PROJECT_NAME}_CMAKE_ARGS}")
  install(FILES
      ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_cmake_args.txt
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_upper)
  foreach(lang IN LISTS languages)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${lang}_compile_line.txt
         "${CMAKE_${lang}_COMPILER} ${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${build_type_upper}}")
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${lang}_compile_line.txt
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  endforeach()
  foreach(tgt IN LISTS ${PROJECT_NAME}_EXPORTED_TARGETS)
    get_target_property(tgt_type "${tgt}" TYPE)
    if (tgt_type MATCHES "STATIC_LIBRARY|SHARED_LIBRARY")
      bob_get_link_libs(${tgt} link_libs "")
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${tgt}_libs.txt
           "${link_libs}")
      install(FILES
          ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${tgt}_libs.txt
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
    endif()
  endforeach()

endfunction(bob_install_provenance)

function(bob_install_provenance_no_recurse no_recurse_components_list)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_cmake_args.txt
       "${${PROJECT_NAME}_CMAKE_ARGS}")
  install(FILES
      ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_cmake_args.txt
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_upper)
  foreach(lang IN LISTS languages)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${lang}_compile_line.txt
         "${CMAKE_${lang}_COMPILER} ${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${build_type_upper}}")
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${lang}_compile_line.txt
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  endforeach()
  foreach(tgt IN LISTS ${PROJECT_NAME}_EXPORTED_TARGETS)
    get_target_property(tgt_type "${tgt}" TYPE)
    if (tgt_type MATCHES "STATIC_LIBRARY|SHARED_LIBRARY")
      bob_get_link_libs(${tgt} link_libs "${no_recurse_components_list}")
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${tgt}_libs.txt
           "${link_libs}")
      install(FILES
          ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${tgt}_libs.txt
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
    endif()
  endforeach()

endfunction(bob_install_provenance_no_recurse)

function(bob_end_package)
  include(CMakePackageConfigHelpers)
  set(INCLUDE_INSTALL_DIR include)
  set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
  set(LATEST_FIND_DEPENDENCY
"#The definition of this macro is really inconvenient prior to CMake
#commit ab358d6a859d8b7e257ed1e06ca000e097a32ef6
#we'll just copy the latest code into our Config.cmake file
macro(latest_find_dependency dep)
  if (NOT \${dep}_FOUND)
    set(cmake_fd_quiet_arg)
    if(\${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      set(cmake_fd_quiet_arg QUIET)
    endif()
    set(cmake_fd_required_arg)
    if(\${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      set(cmake_fd_required_arg REQUIRED)
    endif()

    get_property(cmake_fd_alreadyTransitive GLOBAL PROPERTY
      _CMAKE_\${dep}_TRANSITIVE_DEPENDENCY
    )

    find_package(\${dep} \${ARGN}
      \${cmake_fd_quiet_arg}
      \${cmake_fd_required_arg}
    )

    if(NOT DEFINED cmake_fd_alreadyTransitive OR cmake_fd_alreadyTransitive)
      set_property(GLOBAL PROPERTY _CMAKE_\${dep}_TRANSITIVE_DEPENDENCY TRUE)
    endif()

    if (NOT \${dep}_FOUND)
      set(\${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE \"\${CMAKE_FIND_PACKAGE_NAME} could not be found because dependency \${dep} could not be found.\")
      set(\${CMAKE_FIND_PACKAGE_NAME}_FOUND False)
      return()
    endif()
    set(cmake_fd_required_arg)
    set(cmake_fd_quiet_arg)
    set(cmake_fd_exact_arg)
  endif()
endmacro(latest_find_dependency)"
       )
  # PK added this
  set(FIND_EXT_DEPS_CONTENT)
  foreach(dep IN LISTS ${PROJECT_NAME}_EXT_DEPS)
      set(FIND_EXT_DEPS_CONTENT "${FIND_EXT_DEPS_CONTENT}find_package(${dep} REQUIRED NO_DEFAULT_PATH HINTS ${${dep}_PREFIX})\n")
  endforeach()
  set(FIND_DEPS_CONTENT)
  foreach(dep IN LISTS ${PROJECT_NAME}_DEPS)
    string(REPLACE ";" " " FIND_DEP_ARGS "${${dep}_find_package_args}")
    set(FIND_DEPS_CONTENT
"${FIND_DEPS_CONTENT}
latest_find_dependency(${dep} ${FIND_DEP_ARGS})"
       )
  endforeach()
  set(CONFIG_CONTENT
"set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION})
${LATEST_FIND_DEPENDENCY}
${FIND_EXT_DEPS_CONTENT}
${FIND_DEPS_CONTENT}
set(${PROJECT_NAME}_EXPORTED_TARGETS \"${${PROJECT_NAME}_EXPORTED_TARGETS}\")
foreach(tgt IN LISTS ${PROJECT_NAME}_EXPORTED_TARGETS)
  include(\${CMAKE_CURRENT_LIST_DIR}/\${tgt}-target.cmake)
endforeach()"
  )
  foreach(TYPE IN ITEMS "BOOL" "INT" "STRING")
    if (${PROJECT_NAME}_KEY_${TYPE}S)
      foreach(KEY_${TYPE} IN LISTS ${PROJECT_NAME}_KEY_${TYPE}S)
        set(val "${${KEY_${TYPE}}}")
        #escape escapes
        string(REPLACE "\\" "\\\\" val "${val}")
        #escape quotes
        string(REPLACE "\"" "\\\"" val "${val}")
        set(CONFIG_CONTENT
"${CONFIG_CONTENT}
set(${KEY_${TYPE}} \"${val}\")")
      endforeach()
    endif()
  endforeach()
  set(CONFIG_CONTENT
"${CONFIG_CONTENT}
")
  install(FILES
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  if(PROJECT_VERSION)
    file(WRITE
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        "${CONFIG_CONTENT}")
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)
    install(FILES
      "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  endif()
  bob_install_provenance()
endfunction(bob_end_package)

function(bob_end_package_no_recurse no_recurse_components_list)
  include(CMakePackageConfigHelpers)
  set(INCLUDE_INSTALL_DIR include)
  set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
  set(LATEST_FIND_DEPENDENCY
"#The definition of this macro is really inconvenient prior to CMake
#commit ab358d6a859d8b7e257ed1e06ca000e097a32ef6
#we'll just copy the latest code into our Config.cmake file
macro(latest_find_dependency dep)
  if (NOT \${dep}_FOUND)
    set(cmake_fd_quiet_arg)
    if(\${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      set(cmake_fd_quiet_arg QUIET)
    endif()
    set(cmake_fd_required_arg)
    if(\${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      set(cmake_fd_required_arg REQUIRED)
    endif()

    get_property(cmake_fd_alreadyTransitive GLOBAL PROPERTY
      _CMAKE_\${dep}_TRANSITIVE_DEPENDENCY
    )

    find_package(\${dep} \${ARGN}
      \${cmake_fd_quiet_arg}
      \${cmake_fd_required_arg}
    )

    if(NOT DEFINED cmake_fd_alreadyTransitive OR cmake_fd_alreadyTransitive)
      set_property(GLOBAL PROPERTY _CMAKE_\${dep}_TRANSITIVE_DEPENDENCY TRUE)
    endif()

    if (NOT \${dep}_FOUND)
      set(\${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE \"\${CMAKE_FIND_PACKAGE_NAME} could not be found because dependency \${dep} could not be found.\")
      set(\${CMAKE_FIND_PACKAGE_NAME}_FOUND False)
      return()
    endif()
    set(cmake_fd_required_arg)
    set(cmake_fd_quiet_arg)
    set(cmake_fd_exact_arg)
  endif()
endmacro(latest_find_dependency)"
       )
  set(FIND_DEPS_CONTENT)
  foreach(dep IN LISTS ${PROJECT_NAME}_DEPS)
    string(REPLACE ";" " " FIND_DEP_ARGS "${${dep}_find_package_args}")
    set(FIND_DEPS_CONTENT
"${FIND_DEPS_CONTENT}
latest_find_dependency(${dep} ${FIND_DEP_ARGS})"
       )
  endforeach()
  set(CONFIG_CONTENT
"set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION})
${LATEST_FIND_DEPENDENCY}
${FIND_DEPS_CONTENT}
set(${PROJECT_NAME}_EXPORTED_TARGETS \"${${PROJECT_NAME}_EXPORTED_TARGETS}\")
foreach(tgt IN LISTS ${PROJECT_NAME}_EXPORTED_TARGETS)
  include(\${CMAKE_CURRENT_LIST_DIR}/\${tgt}-target.cmake)
endforeach()"
  )
  foreach(TYPE IN ITEMS "BOOL" "INT" "STRING")
    if (${PROJECT_NAME}_KEY_${TYPE}S)
      foreach(KEY_${TYPE} IN LISTS ${PROJECT_NAME}_KEY_${TYPE}S)
        set(val "${${KEY_${TYPE}}}")
        #escape escapes
        string(REPLACE "\\" "\\\\" val "${val}")
        #escape quotes
        string(REPLACE "\"" "\\\"" val "${val}")
        set(CONFIG_CONTENT
"${CONFIG_CONTENT}
set(${KEY_${TYPE}} \"${val}\")")
      endforeach()
    endif()
  endforeach()
  set(CONFIG_CONTENT
"${CONFIG_CONTENT}
")
  install(FILES
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  if(PROJECT_VERSION)
    file(WRITE
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        "${CONFIG_CONTENT}")
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)
    install(FILES
      "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  endif()
  bob_install_provenance_no_recurse("${no_recurse_components_list}")
endfunction(bob_end_package_no_recurse)
