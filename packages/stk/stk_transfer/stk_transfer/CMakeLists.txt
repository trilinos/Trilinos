# Currently STKTransfer is a header only library which means without special
# help its headers will not be installed correctly. All this file does is define
# the installation of the headers

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)

if(NOT STK_ENABLE_STKMiddle_mesh)
  message(STATUS "STKTransfer: STKMiddle_mesh not enabled")
  LIST(REMOVE_ITEM HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/ConservativeTransfer.hpp")
  LIST(REMOVE_ITEM HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/ConservativeTransferUser.hpp")
  LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ConservativeTransfer.cpp")
else()
  message(STATUS "STKTransfer: STKMiddle_mesh is enabled, ConservativeTransfer is available")
endif()

FILE(GLOB HEADERS_IMPL ${CMAKE_CURRENT_SOURCE_DIR}/copy_by_id/*.hpp)
FILE(GLOB SOURCES_IMPL  ${CMAKE_CURRENT_SOURCE_DIR}/copy_by_id/*.cpp)

if(HAVE_STK_Trilinos)
  TRIBITS_ADD_LIBRARY(stk_transfer
    NOINSTALLHEADERS ${HEADERS} ${HEADERS_IMPL}
    SOURCES  ${SOURCES} ${SOURCES_IMPL}
  )
else()
  add_library(stk_transfer ${SOURCES} ${SOURCES_IMPL})
  target_link_libraries(stk_transfer PUBLIC stk_search)
  target_link_libraries(stk_transfer PUBLIC stk_util_env)
  target_link_libraries(stk_transfer PUBLIC stk_util_parallel)
  if (STK_ENABLE_STKMiddle_mesh)
    target_link_libraries(stk_transfer PUBLIC stk_middle_mesh)
  endif()
  INSTALL(TARGETS stk_transfer EXPORT stkTargets DESTINATION ${STK_INSTALL_LIBDIR})
endif()

target_include_directories(stk_transfer PUBLIC
  $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_transfer>
  $<INSTALL_INTERFACE:include>
)

INSTALL(FILES ${HEADERS} DESTINATION
${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_transfer/)

INSTALL(FILES ${HEADERS_IMPL} DESTINATION
${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_transfer/copy_by_id)

