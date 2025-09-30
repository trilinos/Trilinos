# Set the package name
TRIBITS_PACKAGE(PyTrilinos2 DISABLE_STRONG_WARNINGS)

###################################################################################################
# LLVM

find_package(LLVM REQUIRED CONFIG)
message(STATUS "Using LLVM includes in: ${LLVM_INCLUDE_DIRS}")

###################################################################################################
# binder

SET(POTENTIAL_BINDER_BIN "$ENV{BINDER_ROOT}/bin")

FIND_PROGRAM(BINDER_EXECUTABLE REQUIRED
  NAMES binder
  HINTS ${POTENTIAL_BINDER_BIN} ${BINDER_PATH}
)

message(STATUS "Found Binder executable: ${BINDER_EXECUTABLE}")


###################################################################################################
# pybind11

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

###################################################################################################

ASSERT_DEFINED(
  BUILD_SHARED_LIBS
  Python3_EXECUTABLE
  Python3_VERSION_MAJOR
  Python3_VERSION_MINOR
  ${PACKAGE_NAME}_ENABLE_MueLu
  Tpetra_INST_OPENMP
  BINDER_EXECUTABLE
)

###################################################################################################

IF(NOT BUILD_SHARED_LIBS)
  MESSAGE(WARNING "PyTrilinos2 can be built with static linking, but might require building with position-independent code. If you encounter build errors please add \"-fPIC\" to CMAKE_C_FLAGS, CMAKE_C_FLAGS and CMAKE_Fortran_FLAGS.")
ENDIF()

# Set the package version number
SET(PyTrilinos2_VERSION ${Trilinos_VERSION})

TRIBITS_ADD_OPTION_AND_DEFINE(PyTrilinos2_BINDER_SUPPRESS_ERRORS
  PYTRILINOS2_SUPPRESS_ERRORS
  "Enable the suppress errors option of binder."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PyTrilinos2_BINDER_USE_ONE_FILE
  PYTRILINOS2_USE_ONE_FILE
  "Enable the use of one file by binder."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PyTrilinos2_BINDER_CMAKE_ERROR
PYTRILINOS2_CMAKE_ERROR
  "Stop the configuration if binder fails."
  ON )

TRIBITS_ADD_OPTION_AND_DEFINE(PyTrilinos2_BINDER_VERBOSE
  PYTRILINOS2_VERBOSE
  "Increase the verbosity of binder."
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(PyTrilinos2_PYBIND11_VERBOSE
  PYBIND11_DETAILED_ERROR_MESSAGES
  "Increase the verbosity of pybind11."
  OFF )


###################################################################################################

SET(PyTrilinos2_BINDER_NUM_FILES "150" CACHE STRING "Maxinum number of generated files by binder.")

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
FILE(GLOB PyTrilinos2PyFiles ${CMAKE_CURRENT_SOURCE_DIR}/python/*.py)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/PyTrilinos2)

file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/include_tmp)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include_tmp)
file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/binder)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/binder)

set(binder_include_name "${CMAKE_CURRENT_BINARY_DIR}/trilinos_includes.hpp")
set(all_header_with_dir_list "${CMAKE_CURRENT_BINARY_DIR}/list_with_dir.txt")
set(all_header_without_dir_list "${CMAKE_CURRENT_BINARY_DIR}/list_without_dir.txt")
set(all_ETI_files_list "${CMAKE_CURRENT_BINARY_DIR}/all_ETI_files_list.txt")
set(all_ETI_classes_list "${CMAKE_CURRENT_BINARY_DIR}/all_ETI_classes_list.txt")
set(all_include_list "${CMAKE_CURRENT_BINARY_DIR}/list_include.txt")

set(all_include_dirs "")
set(all_visited_libs "")
foreach(depPkg IN LISTS PyTrilinos2_LIB_ENABLED_DEPENDENCIES)
  get_all_include_dirs(${depPkg}::all_libs "${all_include_dirs}" "${all_visited_libs}")
endforeach()

foreach(all_include_dir IN LISTS all_include_dirs)
  STRING(FIND ${all_include_dir} "kokkos/algorithms/src" IS_KOKKOS_ALGORITHMS_SRC)
  IF(${IS_KOKKOS_ALGORITHMS_SRC} GREATER -1)
    list(APPEND all_include_dirs "${all_include_dir}/std_algorithms")
    list(APPEND all_include_dirs "${all_include_dir}/std_algorithms/impl")
    list(APPEND all_include_dirs "${all_include_dir}/sorting")
    list(APPEND all_include_dirs "${all_include_dir}/sorting/impl")
  ENDIF()
endforeach()

list(REMOVE_DUPLICATES all_include_dirs)
list(REMOVE_ITEM all_include_dirs "")
#MESSAGE("all_include_dirs = ${all_include_dirs}")

set(PyTrilinos2_all_include_files_with_dir "")
set(PyTrilinos2_all_include_files_without_dir "")
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
    list(APPEND PyTrilinos2_all_include_files_with_dir "${include_file}")
    string(REPLACE "${include_dir}/" "" include_file_without_dir "${include_file}")
    list(APPEND PyTrilinos2_all_include_files_without_dir "${include_file_without_dir}")
  endforeach()
endforeach()


#list(REMOVE_DUPLICATES PyTrilinos2_all_include_files_without_dir)
#list(REMOVE_ITEM PyTrilinos2_all_include_files_without_dir "")

#list(REMOVE_DUPLICATES PyTrilinos2_all_include_files_with_dir)
#list(REMOVE_ITEM PyTrilinos2_all_include_files_with_dir "")

#MESSAGE("PyTrilinos2_all_include_files_with_dir = ${PyTrilinos2_all_include_files_with_dir}")

SET(CONTENTS "")
FOREACH(line IN LISTS all_include_dirs)
  SET(CONTENTS "${CONTENTS}${line}\n")
ENDFOREACH(line)
file(WRITE ${all_include_list} ${CONTENTS})

set(eti_files_with_dir "")
set(eti_files_without_dir "")

###################################################################################################
# Get ETI files for Tpetra:
file(GLOB tpetra_ETI_files
    "${CMAKE_CURRENT_BINARY_DIR}/../tpetra/core/src/*.cpp"
)
list(APPEND tpetra_ETI_files "${CMAKE_CURRENT_BINARY_DIR}/../tpetra/core/src/TpetraCore_ETIHelperMacros.h")

foreach(tpetra_ETI_file IN LISTS tpetra_ETI_files)
  list(APPEND eti_files_with_dir "${tpetra_ETI_file}")
  get_filename_component(tpetra_ETI_file_without_dir "${tpetra_ETI_file}" NAME)
  list(APPEND eti_files_without_dir "${tpetra_ETI_file_without_dir}")
endforeach()

SET(CONTENTS "")
FOREACH(line IN LISTS PyTrilinos2_all_include_files_with_dir)
  SET(CONTENTS "${CONTENTS}${line}\n")
ENDFOREACH(line)
FOREACH(line IN LISTS eti_files_with_dir)
  SET(CONTENTS "${CONTENTS}${line}\n")
ENDFOREACH(line)
file(WRITE ${all_header_with_dir_list} ${CONTENTS})

SET(CONTENTS "")
FOREACH(line IN LISTS PyTrilinos2_all_include_files_without_dir)
  SET(CONTENTS "${CONTENTS}${line}\n")
ENDFOREACH(line)
FOREACH(line IN LISTS eti_files_without_dir)
  SET(CONTENTS "${CONTENTS}${line}\n")
ENDFOREACH(line)
file(WRITE ${all_header_without_dir_list} ${CONTENTS})

SET(CONTENTS "")
FOREACH(line IN LISTS eti_files_without_dir)
  SET(CONTENTS "${CONTENTS}${line}\n")
ENDFOREACH(line)
file(WRITE ${all_ETI_files_list} ${CONTENTS})

SET(ETI_classes "Tpetra_CrsMatrix;Tpetra_Vector;Tpetra_MultiVector;Tpetra_FEMultiVector;Tpetra_FECrsMatrix")
SET(CONTENTS "")
FOREACH(line IN LISTS ETI_classes)
  SET(CONTENTS "${CONTENTS}${line}\n")
ENDFOREACH(line)
file(WRITE ${all_ETI_classes_list} ${CONTENTS})

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/src/PyTrilinos2_Tpetra_ETI.hpp ${CMAKE_CURRENT_BINARY_DIR}/python/getTpetraTypeName.py
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/gather_ETI.py ${CMAKE_CURRENT_BINARY_DIR} ${all_ETI_files_list} ${all_ETI_classes_list} "src/PyTrilinos2_Tpetra_ETI.hpp"
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/include_tmp
)
add_custom_target(generate_ETI_name DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/src/PyTrilinos2_Tpetra_ETI.hpp ${CMAKE_CURRENT_BINARY_DIR}/include_tmp)

file (GLOB PyTrilinos2PyFiles2 "${CMAKE_CURRENT_BINARY_DIR}/python/*.py")
list (APPEND PyTrilinos2PyFiles ${PyTrilinos2PyFiles2})

add_custom_command(
  OUTPUT  ${binder_include_name}
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/gather_includes.py ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${all_header_with_dir_list} ${all_header_without_dir_list} ${binder_include_name}
  DEPENDS generate_ETI_name ${CMAKE_CURRENT_BINARY_DIR}/include_tmp
)
add_custom_target(generate_include_name DEPENDS ${binder_include_name})
add_dependencies(generate_include_name generate_ETI_name)

set(BINDER_OPTIONS "")
list(APPEND BINDER_OPTIONS --root-module PyTrilinos2)
list(APPEND BINDER_OPTIONS --prefix ${CMAKE_CURRENT_BINARY_DIR}/binder)
IF(PYTRILINOS2_USE_ONE_FILE)
  list(APPEND BINDER_OPTIONS -single-file)
ELSE()
  list(APPEND BINDER_OPTIONS -max-file-size=1000000)
  list(APPEND BINDER_OPTIONS -flat)
ENDIF()
IF(${PACKAGE_NAME}_ENABLE_MueLu)
  list(APPEND BINDER_OPTIONS --bind MueLu)
ENDIF()
IF(PYTRILINOS2_VERBOSE)
  list(APPEND BINDER_OPTIONS -v)
ENDIF()
IF(PYTRILINOS2_SUPPRESS_ERRORS)
  list(APPEND BINDER_OPTIONS --suppress-errors)
ENDIF()
list(APPEND BINDER_OPTIONS --config ${CMAKE_CURRENT_SOURCE_DIR}/scripts/PyTrilinos2_config.cfg)
list(APPEND BINDER_OPTIONS --)
IF(TPL_ENABLE_CUDA)
  list(APPEND BINDER_OPTIONS -x cuda --cuda-host-only)
ENDIF()
list(APPEND BINDER_OPTIONS ${PyTrilinos2_BINDER_FLAGS})
list(APPEND BINDER_OPTIONS -std=c++17)
if (NOT(MPI_BASE_DIR STREQUAL ""))
  list(APPEND BINDER_OPTIONS -I${MPI_BASE_DIR}/include)
ENDIF()
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp/mdspan)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp/View/MDSpan)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp/experimental)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp/experimental/__p0009_bits)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp/experimental/__p1684_bits)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp/experimental/__p2389_bits)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp/experimental/__p2630_bits)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/include_tmp/experimental/__p2642_bits)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_BINARY_DIR}/src)
list(APPEND BINDER_OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR}/src)
IF(NOT DEFINED PyTrilinos2_BINDER_GCC_TOOLCHAIN)
  list(APPEND BINDER_OPTIONS -I${LLVM_INCLUDE_DIRS})
  list(APPEND BINDER_OPTIONS -iwithsysroot${LLVM_INSTALL_PREFIX}/lib/clang/${LLVM_VERSION}/include)
ELSE()
  list(APPEND BINDER_OPTIONS --gcc-toolchain=${PyTrilinos2_BINDER_GCC_TOOLCHAIN})
ENDIF()
IF (Tpetra_INST_OPENMP)
  list(APPEND BINDER_OPTIONS "-fopenmp")
ENDIF()
IF ("${TPL_MPI_INCLUDE_DIRS}" MATCHES "")
  find_package(MPI REQUIRED COMPONENTS CXX)
  EXECUTE_PROCESS(COMMAND ${MPI_CXX_COMPILER} --showme:compile OUTPUT_VARIABLE PYTRILONOS2_MPI_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  SEPARATE_ARGUMENTS(PYTRILONOS2_MPI_CXX_FLAGS)
  list(APPEND BINDER_OPTIONS ${PYTRILONOS2_MPI_CXX_FLAGS})
ELSE()
  list(APPEND BINDER_OPTIONS "-I${TPL_MPI_INCLUDE_DIRS}")
ENDIF()
list(APPEND BINDER_OPTIONS -DNDEBUG)
list(APPEND BINDER_OPTIONS -DPYBIND11_DETAILED_ERROR_MESSAGES)

message("BINDER_OPTIONS='${BINDER_OPTIONS}'")

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/binder/PyTrilinos2.cpp
  COMMAND ${BINDER_EXECUTABLE} ${binder_include_name} ${BINDER_OPTIONS}
  DEPENDS ${binder_include_name} generate_include_name
)
add_custom_target(binder_call DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/binder/PyTrilinos2.cpp)
add_dependencies(binder_call generate_ETI_name generate_include_name)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(PyTrilinos2_DEFAULT_INSTALL_PREFIX ${PYTHON_PREFIX})
ELSE()
  SET(PyTrilinos2_DEFAULT_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
ENDIF()

###################################################################################################
# Set the PyTrilinos2 install prefix
SET(PyTrilinos2_INSTALL_PREFIX ${PyTrilinos2_DEFAULT_INSTALL_PREFIX}
  CACHE PATH "The path prefix for where PyTrilinos2 will be installed, e.g. /usr/local")

# Determine the install directory
SET(PyTrilinos2_INSTALL_DIR
  ${PyTrilinos2_INSTALL_PREFIX}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/PyTrilinos2
  CACHE PATH "The path where PyTrilinos2 will be installed"
  )
MESSAGE(STATUS "PyTrilinos2 installation path: ${PyTrilinos2_INSTALL_DIR}")

INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/python/getTpetraTypeName.py
  DESTINATION ${PyTrilinos2_INSTALL_DIR})

###################################################################################################
# mpi4py

EXECUTE_PROCESS(COMMAND
  ${Python3_EXECUTABLE} -c "import mpi4py; print(mpi4py.get_include())"
  OUTPUT_VARIABLE Mpi4Py_INCLUDE_DIR
  ERROR_VARIABLE  Mpi4Py_INCLUDE_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

IF(NOT Mpi4Py_VERSION_ERROR)
  MESSAGE("  -- Mpi4Py Enabled.")
ELSE()
  MESSAGE(FATAL_ERROR "Mpi4Py_VERSION_ERROR is defined; the python executable cannot access mpi4py.")
ENDIF()

EXECUTE_PROCESS(COMMAND
  ${Python3_EXECUTABLE} -c "import mpi4py; print(mpi4py.get_config()['mpicxx'])"
  OUTPUT_VARIABLE Mpi4Py_MPICXX
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

IF(NOT ${Mpi4Py_MPICXX} STREQUAL ${CMAKE_CXX_COMPILER})
  MESSAGE(WARNING "the cpp compiler used to compile mpi4py ${Mpi4Py_MPICXX} is not consistent with CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
ENDIF()

###################################################################################################

ADD_SUBDIRECTORY( src )

#file (GLOB PyTrilinos2PyFilesSo "${CMAKE_CURRENT_BINARY_DIR}/src/*.so")
file(COPY ${PyTrilinos2PyFiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/PyTrilinos2/.)
#file(COPY ${PyTrilinos2PyFilesSo} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/PyTrilinos2/.)

SET(PyTrilinos2_PYTHONPATH "${CMAKE_CURRENT_BINARY_DIR}")

TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

TRIBITS_ADD_TEST_DIRECTORIES(test)

# Execute the package postprocessing
TRIBITS_PACKAGE_POSTPROCESS()
