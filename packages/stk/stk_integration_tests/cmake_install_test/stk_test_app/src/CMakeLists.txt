
target_sources(test_stk_lib PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/test_stk_coupling.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_stk_search.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_stk_simd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_stk_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_stk_tools.cpp
)

