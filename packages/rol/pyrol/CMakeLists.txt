# -*- cmake -*-

IF(NOT BUILD_SHARED_LIBS)
  MESSAGE(FATAL_ERROR "PyROL can only be built with shared libraries. Building of shared libraries is currently set to OFF. To enable shared libraries please set the cache variable \"BUILD_SHARED_LIBS\" to ON")
ENDIF()

# Set the package version number
SET(PyROL_VERSION ${Trilinos_VERSION})

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_PIP_INSTALL
  PYROL_SCIKIT
  "Logic specific to the scikit-build-core backend."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_ENABLE_BINDER
  PYROL_GENERATE_SRC
  "Enable regeneration of source files with Binder."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_BINDER_SUPPRESS_ERRORS
  PYROL_SUPPRESS_ERRORS
  "Enable the suppress errors option of Binder."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_BINDER_CMAKE_ERROR
  PYROL_CMAKE_ERROR
  "Stop the configuration if binder fails."
  ON )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_BINDER_VERBOSE
  PYROL_B_VERBOSE
  "Increase the verbosity of binder"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PYROL_ENABLE_BINDER_UPDATE
  PYROL_UPDATE_GENERATED_SRC
  "Enable the update of the generated source files with Binder."
  OFF )

MESSAGE("-- PYTHON_EXECUTABLE:")
IF(NOT DEFINED ${PYTHON_EXECUTABLE})
  find_program(PYTHON_EXECUTABLE
      NAMES python3 python
      )
  MESSAGE("  -- CMake has set: PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
ELSE()
  MESSAGE("  -- User has set: PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
ENDIF()


function(get_all_include_dirs  LIBRARY_NAME all_include_dirs all_visited_libs)
  if (TARGET ${LIBRARY_NAME})
    get_property(depend_libs TARGET ${LIBRARY_NAME} PROPERTY INTERFACE_LINK_LIBRARIES)
    foreach(depend_lib IN LISTS depend_libs)
      if (TARGET ${depend_lib} AND (NOT ${depend_lib} IN_LIST all_visited_libs))
        list(APPEND all_visited_libs "${depend_lib}") # Update list in the current scope only
        get_property(current_includes TARGET ${depend_lib} PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
        foreach(include IN LISTS current_includes)
          STRING(REPLACE "$<BUILD_INTERFACE:" "" new_tmp_include ${include})
          STRING(REPLACE ">" "" new_include ${new_tmp_include})
          list(APPEND all_include_dirs "${new_include}") # Update list in the current scope only
        endforeach()
        get_all_include_dirs(${depend_lib} "${all_include_dirs}" "${all_visited_libs}")
      endif()
    endforeach()
    set(all_include_dirs ${all_include_dirs} PARENT_SCOPE)
    set(all_visited_libs ${all_visited_libs} PARENT_SCOPE)
  endif()
endfunction()

# Python files to install
FILE(GLOB PyROLPyFiles ${CMAKE_CURRENT_SOURCE_DIR}/python/*.py)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pyrol)

GET_PROPERTY( PTR_IMPL             GLOBAL PROPERTY PTR_IMPL  )
GET_PROPERTY( PTR_DIR              GLOBAL PROPERTY PTR_DIR  )
GET_PROPERTY( STACKTRACE_DIR       GLOBAL PROPERTY STACKTRACE_DIR  )
GET_PROPERTY( PARAMETERLIST_DIR    GLOBAL PROPERTY PARAMETERLIST_DIR  )
GET_PROPERTY( LA_DIR               GLOBAL PROPERTY LA_DIR  )
GET_PROPERTY( LAPACK_DIR           GLOBAL PROPERTY LAPACK_DIR  )
GET_PROPERTY( BLAS_DIR             GLOBAL PROPERTY BLAS_DIR  )

set(ROL_all_include_dirs "")
list(APPEND ROL_all_include_dirs "${CMAKE_CURRENT_BINARY_DIR}/../src")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeB")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeE")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeG")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeU")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeB/pqn")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeG/augmentedlagrangian/")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeG/fletcher")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeG/interiorpoint")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeG/moreauyosida")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeG/stabilizedlcl")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeU/bundle")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeU/linesearch")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeU/linesearch/descent")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/algorithm/TypeU/trustregion")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/elementwise")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/dynamic")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/boundconstraint")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/constraint")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/nlls")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/polyproj")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/objective")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/operator")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/simopt")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/sketching")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/function/std")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/distribution")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/expectationquad")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/randvarfunctional")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/randvarfunctional/risk")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/randvarfunctional/risk/fdivergence")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/randvarfunctional/risk/spectral")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/randvarfunctional/error")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/randvarfunctional/regret")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/randvarfunctional/deviation")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/randvarfunctional/probability")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/function/progressivehedging")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/algorithm")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/algorithm/PrimalDual")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/vector")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/sampler")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/sampler/SROM")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/sol/status")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/status")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/augmentedlagrangian")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/bundle")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/fletcher")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/interiorpoint")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/krylov")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/linesearch")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/moreauyosidapenalty")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/nonlinearcg")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/secant")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/step/trustregion")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/utils")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/utils/function_bindings")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/vector")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/zoo")
list(APPEND ROL_all_include_dirs "${PTR_DIR}")
list(APPEND ROL_all_include_dirs "${STACKTRACE_DIR}")
list(APPEND ROL_all_include_dirs "${PARAMETERLIST_DIR}")
list(APPEND ROL_all_include_dirs "${LA_DIR}")
list(APPEND ROL_all_include_dirs "${LAPACK_DIR}")
list(APPEND ROL_all_include_dirs "${BLAS_DIR}")

list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/oed")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/oed/constraint")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/oed/factors")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/oed/objective")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/oed/objective/het")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/oed/objective/hom")
list(APPEND ROL_all_include_dirs "${${PACKAGE_NAME}_SOURCE_DIR}/src/oed/utilities")

IF (PYROL_GENERATE_SRC)
  MESSAGE("-- PyROL_BINDER_EXECUTABLE:")
  IF(NOT DEFINED PyROL_BINDER_EXECUTABLE)
    find_program(PyROL_BINDER_EXECUTABLE
        NAMES binder
        )
    MESSAGE("  -- CMake has set: PyROL_BINDER_EXECUTABLE = ${PyROL_BINDER_EXECUTABLE}")
  ELSE()
    MESSAGE("  -- User has set: PyROL_BINDER_EXECUTABLE = ${PyROL_BINDER_EXECUTABLE}")
  ENDIF()
  file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/include_tmp)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include_tmp)
  file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/binder)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/binder)

  MESSAGE("PTR_IMPL = ${PTR_IMPL}")
  IF( PTR_IMPL STREQUAL "std::shared_ptr")
    SET( BINDER_CFG "${CMAKE_CURRENT_SOURCE_DIR}/scripts/PyROL_shared_ptr.cfg" )
  ELSE()
    SET( BINDER_CFG "${CMAKE_CURRENT_SOURCE_DIR}/scripts/PyROL_RCP.cfg" )
  ENDIF()

  set(binder_include_name "${CMAKE_CURRENT_BINARY_DIR}/trilinos_includes.hpp")
  set(all_header_with_dir_list "${CMAKE_CURRENT_BINARY_DIR}/list_with_dir.txt")
  set(all_header_without_dir_list "${CMAKE_CURRENT_BINARY_DIR}/list_without_dir.txt")
  set(all_ETI_classes_list "${CMAKE_CURRENT_BINARY_DIR}/all_ETI_classes_list.txt")
  set(all_include_list "${CMAKE_CURRENT_BINARY_DIR}/list_include.txt")

  set(all_include_dirs "")
  set(all_visited_libs "")
  foreach(depPkg IN LISTS ROL_LIB_ENABLED_DEPENDENCIES)
    get_all_include_dirs(${depPkg}::all_libs "${all_include_dirs}" "${all_visited_libs}")
  endforeach()

  list(REMOVE_DUPLICATES all_include_dirs)
  list(REMOVE_ITEM all_include_dirs "")

  FOREACH(include_dir IN LISTS ROL_all_include_dirs)
    list(APPEND all_include_dirs "${include_dir}")
  ENDFOREACH(include_dir)
  #MESSAGE("all_include_dirs = ${all_include_dirs}")

  set(PyROL_all_include_files_with_dir "")
  set(PyROL_all_include_files_without_dir "")
  foreach(include_dir IN LISTS all_include_dirs)
    file(GLOB include_files
        "${include_dir}/*.hpp"
        "${include_dir}/*.h"
        "${include_dir}/*/*.hpp"
        "${include_dir}/*/*.h"
        "${include_dir}/*/*/*.hpp"
        "${include_dir}/*/*/*.h"
        "${include_dir}/*/*/*/*.hpp"
        "${include_dir}/*/*/*/*.h"
        "${include_dir}/*/*/*/*.inc"
        "${include_dir}/*/*/*/*.inc_predicate"
    )
    foreach(include_file IN LISTS include_files)
      list(APPEND PyROL_all_include_files_with_dir "${include_file}")
      string(REPLACE "${include_dir}/" "" include_file_without_dir "${include_file}")
      list(APPEND PyROL_all_include_files_without_dir "${include_file_without_dir}")
    endforeach()
  endforeach()


  #list(REMOVE_DUPLICATES PyROL_all_include_files_without_dir)
  #list(REMOVE_ITEM PyROL_all_include_files_without_dir "")

  #list(REMOVE_DUPLICATES PyROL_all_include_files_with_dir)
  #list(REMOVE_ITEM PyROL_all_include_files_with_dir "")

  #MESSAGE("PyROL_all_include_files_with_dir = ${PyROL_all_include_files_with_dir}")

  SET(CONTENTS "")
  FOREACH(line IN LISTS all_include_dirs)
    SET(CONTENTS "${CONTENTS}${line}\n")
  ENDFOREACH(line)
  file(WRITE ${all_include_list} ${CONTENTS})

  SET(CONTENTS "")
  FOREACH(line IN LISTS PyROL_all_include_files_with_dir)
    SET(CONTENTS "${CONTENTS}${line}\n")
  ENDFOREACH(line)
  file(WRITE ${all_header_with_dir_list} ${CONTENTS})

  SET(CONTENTS "")
  FOREACH(line IN LISTS PyROL_all_include_files_without_dir)
    SET(CONTENTS "${CONTENTS}${line}\n")
  ENDFOREACH(line)
  file(WRITE ${all_header_without_dir_list} ${CONTENTS})

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)
      
  file (GLOB PyROLPyFiles2 "${CMAKE_CURRENT_BINARY_DIR}/python/*.py")
  list (APPEND PyROLPyFiles ${PyROLPyFiles2})

  EXECUTE_PROCESS(COMMAND
  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/gather_includes.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${all_header_with_dir_list} ${all_header_without_dir_list} ${binder_include_name}
  )

  set(BINDER_OPTIONS "")
  list(APPEND BINDER_OPTIONS ${binder_include_name})
  list(APPEND BINDER_OPTIONS --root-module pyrol)
  list(APPEND BINDER_OPTIONS --prefix ${CMAKE_CURRENT_BINARY_DIR}/binder)
  list(APPEND BINDER_OPTIONS -max-file-size=1000000)
  list(APPEND BINDER_OPTIONS --bind Teuchos)
  list(APPEND BINDER_OPTIONS --bind ROL)
  list(APPEND BINDER_OPTIONS --bind pyrol)
  IF(PYROL_B_VERBOSE)
    list(APPEND BINDER_OPTIONS -v)
  ENDIF()
  list(APPEND BINDER_OPTIONS --config ${BINDER_CFG})
  IF(PYROL_SUPPRESS_ERRORS)
    list(APPEND BINDER_OPTIONS --suppress-errors)
  ENDIF()  
  list(APPEND BINDER_OPTIONS --)
  IF(TPL_ENABLE_CUDA)
    list(APPEND BINDER_OPTIONS -x cuda --cuda-host-only)
  ENDIF()
  list(APPEND BINDER_OPTIONS ${PyROL_BINDER_FLAGS})
  list(APPEND BINDER_OPTIONS -std=c++${CMAKE_CXX_STANDARD})
  list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp)
  list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/src)
  list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR}/src)

  list(APPEND BINDER_OPTIONS -I${PyROL_BINDER_clang_include_dirs})
  list(APPEND BINDER_OPTIONS -iwithsysroot${PyROL_BINDER_LibClang_include_dir})
  IF(PyROL_BINDER_GCC_TOOLCHAIN)
    list(APPEND BINDER_OPTIONS --gcc-toolchain=${PyROL_BINDER_GCC_TOOLCHAIN})
  ENDIF()
  list(APPEND BINDER_OPTIONS -DNDEBUG)

  message("BINDER_OPTIONS='${BINDER_OPTIONS}'")

  EXECUTE_PROCESS(COMMAND
    ${PyROL_BINDER_EXECUTABLE} ${BINDER_OPTIONS}
    RESULT_VARIABLE STATUS
    OUTPUT_VARIABLE OUTPUT_BINDER 
  )

  if(STATUS AND NOT STATUS EQUAL 0)
    message("${OUTPUT_BINDER}")
    if(PYROL_CMAKE_ERROR)
      message(FATAL_ERROR "BINDER FAILED: ${STATUS}")
    else()
      message("BINDER FAILED: ${STATUS}")
    endif()
  else()
    message(STATUS "BINDER SUCCESS:") 
    message("${OUTPUT_BINDER}")
  endif()

  IF(PYROL_UPDATE_GENERATED_SRC)
    FILE(REMOVE_RECURSE ${CMAKE_CURRENT_SOURCE_DIR}/binder)
    FILE(COPY ${CMAKE_CURRENT_BINARY_DIR}/binder DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
  ENDIF()
ELSE()
  # Copy binder files to ${CMAKE_CURRENT_BINARY_DIR}
  FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/binder DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
ENDIF()

# Get the python version
EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c
                        "import sys; print(sys.version_info.major)"
  OUTPUT_VARIABLE PYTHON_MAJOR_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c
                        "import sys; print(sys.version_info.minor)"
  OUTPUT_VARIABLE PYTHON_MINOR_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

SET(PYTHON_VERSION ${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION})

SET(PYBIND11_PYTHON_VERSION ${PYTHON_VERSION})

# Determine the install directory
IF(PYROL_SCIKIT)
  SET(PyROL_INSTALL_DIR
    "."
    CACHE STRING "The path where PyROL will be installed."
  )
ELSE()
  SET(PyROL_INSTALL_DIR
    lib/python${PYTHON_VERSION}/site-packages/pyrol
    CACHE STRING "The path where PyROL will be installed"
    )
ENDIF()

cmake_path(IS_RELATIVE PyROL_INSTALL_DIR is_pyrol_install_relative)

IF(NOT is_pyrol_install_relative)
  MESSAGE(FATAL_ERROR "PyROL install dir is not relative.")
ENDIF()

cmake_path(GET PyROL_INSTALL_DIR RELATIVE_PART PyROL_INSTALL_DIR_RELATIVE)

MESSAGE(STATUS "PyROL installation path: ${PyROL_INSTALL_DIR_RELATIVE}")

INSTALL(FILES
  ${PyROLPyFiles}
  DESTINATION ${PyROL_INSTALL_DIR})

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(PyROL_SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/python)

FOREACH(PyROL_subdir ${PyROL_SUBDIRS})
  INSTALL(DIRECTORY
    ${PyROL_subdir}
    DESTINATION ${PyROL_INSTALL_DIR})
ENDFOREACH()

# Find the pybind11 CMake module
EXECUTE_PROCESS(COMMAND
  ${PYTHON_EXECUTABLE} -c "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE pybind11_DIR
  ERROR_VARIABLE  pybind11_CMAKE_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
MESSAGE(STATUS "pybind11 CMake path: ${pybind11_DIR}")

find_package(pybind11 REQUIRED)

EXECUTE_PROCESS(COMMAND
  ${PYTHON_EXECUTABLE} -c "import mpi4py; print(mpi4py.get_include())"
  OUTPUT_VARIABLE Mpi4Py_INCLUDE_DIR
  ERROR_VARIABLE  Mpi4Py_INCLUDE_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

ADD_SUBDIRECTORY( src )

file(COPY ${PyROLPyFiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pyrol/.)

FOREACH(PyROL_subdir ${PyROL_SUBDIRS})
  FILE(COPY ${PyROL_subdir} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pyrol)
ENDFOREACH()

SET(PyROL_PYTHONPATH "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

TRIBITS_ADD_TEST_DIRECTORIES(test)
