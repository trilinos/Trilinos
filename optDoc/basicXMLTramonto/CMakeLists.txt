cmake_minimum_required(VERSION 2.8)

set(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos REQUIRED)

if(NOT Trilinos_FOUND)
  Message(FATAL_ERROR "Could not find trilinos. Please make sure to set the TRILINOS_PATH variable correctly at configure time. You gave \"" ${TRILINOS_PATH} "\" as the value for TRILINOS_PATH")
endif(NOT Trilinos_FOUND)


SET(CMAKE_C_COMPILIER ${Trilinos_C_COMPILER})
SET(CMAKE_CXX_COMPILIER ${Trilinos_CXX_COMPILER})
SET(CMAKE_CXX_COMPILIER_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})
SET(CMAKE_C_COMPILIER_FLAGS ${Trilinos_C_COMPILER_FLAGS})

project(basicTramontoXML)
include_directories(${Trilinos_TPL_INCLUDE_DIRS})
include_directories(${Trilinos_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})
link_directories(${Trilinos_TPL_LIBRARY_DIRS})

add_executable(basicTramontoXML main.cpp)
target_link_libraries(basicTramontoXML
  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
