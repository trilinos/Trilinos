# @HEADER
# ***********************************************************************
#
#          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
#                 Copyright (2022) Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia
# Corporation, the U.S. Government retains certain rights in this
# software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Kim Liegeois (knliege@sandia.gov)
#
# ***********************************************************************
# @HEADER

FILE(GLOB PYTRILINOS2_SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
MESSAGE("PYTRILINOS2_SRC = ${PYTRILINOS2_SRC}")
FILE(COPY ${PYTRILINOS2_SRC} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

list(APPEND PYTRILINOS2_SRC ${CMAKE_CURRENT_BINARY_DIR}/../binder/PyTrilinos2.cpp)

IF(NOT PYTRILINOS2_USE_ONE_FILE)
  MATH(EXPR NUMBER_FILE "${PyTrilinos2_BINDER_NUM_FILES}")

  foreach(index RANGE 0 ${NUMBER_FILE})
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/../binder/PyTrilinos2_${index}.cpp "")
    list(APPEND PYTRILINOS2_SRC ${CMAKE_CURRENT_BINARY_DIR}/../binder/PyTrilinos2_${index}.cpp)
  endforeach()

  MATH(EXPR NUMBER_FILE "${NUMBER_FILE}+1")
  add_custom_target(
    checkBinderNumFiles
    COMMAND "${CMAKE_COMMAND}"
    -D "NUMBER_FILE=${NUMBER_FILE}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/checkNumberFiles.cmake"
    DEPENDS binder_call
  )
ENDIF()

MESSAGE("PYTRILINOS2_SRC with binder = ${PYTRILINOS2_SRC}")

pybind11_add_module(PyTrilinos2 ${PYTRILINOS2_SRC})

add_dependencies(PyTrilinos2 binder_call generate_ETI_name generate_include_name checkBinderNumFiles)
target_include_directories(PyTrilinos2 PUBLIC ${Mpi4Py_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(PyTrilinos2 PUBLIC cxx_std_14)

foreach(depPkg IN LISTS PyTrilinos2_LIB_ENABLED_DEPENDENCIES) 
    target_link_libraries(PyTrilinos2 PUBLIC ${depPkg}::all_libs)
endforeach()
target_link_libraries(PyTrilinos2 PUBLIC ${Trilinos_EXTRA_LINK_FLAGS})
set_target_properties(PyTrilinos2 PROPERTIES SUFFIX ".so")

INSTALL(TARGETS PyTrilinos2
  DESTINATION ${PyTrilinos2_INSTALL_DIR})

add_custom_command(TARGET PyTrilinos2 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/PyTrilinos2.so ${CMAKE_CURRENT_BINARY_DIR}/../PyTrilinos2/.
    COMMENT "Copy ${PROJECT_BINARY_DIR}/src/PyTrilinos2.so"
)