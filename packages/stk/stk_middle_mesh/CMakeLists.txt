cmake_minimum_required(VERSION 3.19.4)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(Triangulator)

#set(BLA_VENDOR OpenBLAS)
#find_package(LAPACK REQUIRED)  # OpenBLAS has both BLAS and LAPACK in a single lib
#find_package(BLAS REQUIRED)
#
#message("BLAS_FOUND = ${BLAS_FOUND}")
#message("BLAS_INCLUDE_DIR  = ${BLAS_INCLUDE_DIR}")
#message("BLAS_LINKER_FLAGS = ${BLAS_LINKER_FLAGS}")
#message("BLAS_LIBRARIES    = ${BLAS_LIBRARIES}")
#
#get_target_property(dirs BLAS::BLAS INTERFACE_LINK_LIBRARIES)
#foreach(dir ${dirs})
#  message("dir='${dir}'")
#endforeach()
#
#get_target_property(dirs BLAS::BLAS INTERFACE_INCLUDE_LIBRARIES)
#foreach(dir ${dirs})
#  message("dir='${dir}'")
#endforeach()

message("CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")

if (NOT CBLAS_INCLUDE_DIRS)
  message(FATAL_ERROR "Please provide CBLAS_INCLUDE_DIRS")
else()
  message("CBLAS_INCLUDE_DIRS = ${CBLAS_INCLUDE_DIRS}")
endif()

if (NOT CBLAS_LIBRARY_DIRS)
  message(FATAL_ERROR "Please provide CBLAS_LIBRARY_DIRS")
else()
  message("CBLAS_LIBRARY_DIRS = ${CBLAS_LIBRARY_DIRS}")
endif()

if (NOT CBLAS_LIBRARIES)
  message(FATAL_ERROR "Please provide CBLAS_LIBRARY_DIRS")
else()
  message("CBLAS_LIBRARIES = ${CBLAS_LIBRARIES}")
endif()

if (ENABLE_STK)
  if (NOT Trilinos_ROOT)
    message(FATAL_ERROR "Please provide Trilinos_ROOT")
  endif()

  message("-- Searching for Trilinos in Trilinos_ROOT: ${Trilinos_ROOT}")
  find_package(Trilinos NO_MODULE NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_PACKAGE_REGISTRY NO_CMAKE_SYSTEM_PATH NO_CMAKE_SYSTEM_PACKAGE_REGISTRY REQUIRED COMPONENTS
    STK
    STKMesh
    STKIO
  )

#set(STK_SOURCES src/create_stk_mesh.cc
#                  src/write_stk_mesh.cc
#                  src/stk_interface.cc)
  file(GLOB STK_SOURCES src/stk/*.cc)
endif()

message("STK_SOURCES = ${STK_SOURCES}")

#include_directories(${CBLAS_INCLUDE_DIRS})
#link_directories(${CBLAS_LIBRARY_DIRS})


add_subdirectory(CDT/CDT)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)


file(GLOB PREDICATE_SRCS src/predicates/*.cc)

file(GLOB SRCS src/*.cc)
add_library(triangle_lib ${SRCS} ${STK_SOURCES} ${PREDICATE_SRCS})
target_include_directories(triangle_lib PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           "${PROJECT_SOURCE_DIR}/include"
)

target_link_libraries(triangle_lib PUBLIC ${CBLAS_LIBRARIES} PRIVATE CDT)
target_include_directories(triangle_lib PUBLIC ${CBLAS_INCLUDE_DIRS})
target_link_directories(triangle_lib PUBLIC ${CBLAS_LIBRARY_DIRS})

if (ENABLE_STK)
  message("Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
  message("Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
  message("Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
  target_include_directories(triangle_lib PUBLIC ${Trilinos_INCLUDE_DIRS})
  #                                            ${Trilinos_TPL_INCLUDE_DIRS})
  target_link_directories(triangle_lib PUBLIC ${Trilinos_LIBRARY_DIRS})
  #                                         ${Trilinos_TPL_LIBRARY_DIRS})
  target_link_libraries(triangle_lib PUBLIC ${Trilinos_LIBRARIES})
    #${Trilinos_TPL_LIBRARIES})
endif()


add_executable(vtest2 src/executables/test2.cc)
target_link_libraries(vtest2 triangle_lib CDT)

add_executable(vtest4 src/executables/test4.cc)
target_link_libraries(vtest4 triangle_lib CDT)

add_executable(vtest8 src/executables/test8.cc)
target_link_libraries(vtest8 triangle_lib CDT)


if (ENABLE_STK)
  add_executable(vtest6 src/executables/test6.cc)
  target_include_directories(vtest6 PRIVATE ${Trilinos_INCLUDE_DIRS})
  #                                            ${Trilinos_TPL_INCLUDE_DIRS})
  target_link_directories(vtest6 PRIVATE ${Trilinos_LIBRARY_DIRS})
  #                                         ${Trilinos_TPL_LIBRARY_DIRS})
  target_link_libraries(vtest6 PRIVATE ${Trilinos_LIBRARIES}
                                       ${Trilinos_TPL_LIBRARIES})

  add_executable(vtest7 src/executables/test7.cc)
  target_include_directories(vtest7 PRIVATE ${Trilinos_INCLUDE_DIRS})
  #                                            ${Trilinos_TPL_INCLUDE_DIRS})
  target_link_directories(vtest7 PRIVATE ${Trilinos_LIBRARY_DIRS})
  #                                         ${Trilinos_TPL_LIBRARY_DIRS})
  target_link_libraries(vtest7 PRIVATE ${Trilinos_LIBRARIES})
  #                                       ${Trilinos_TPL_LIBRARIES})


  add_executable(nc_generator src/executables/nc_generator.cc)
  target_link_libraries(nc_generator triangle_lib)
endif()


add_subdirectory(test)
