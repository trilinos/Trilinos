cmake_minimum_required(VERSION 3.23)
project(stk_middle_mesh_util)
add_library(stk_middle_mesh_util)
target_sources(stk_middle_mesh_util PRIVATE ./stk_middle_mesh_util/create_stk_mesh.cpp
    ./stk_middle_mesh_util/exodus_writer.cpp
    ./stk_middle_mesh_util/stk_field_copier.cpp
    ./stk_middle_mesh_util/stk_interface.cpp
    ./stk_middle_mesh_util/write_stk_mesh.cpp)
target_link_libraries(stk_middle_mesh_util PUBLIC stk_io)
target_link_libraries(stk_middle_mesh_util PUBLIC stk_mesh_base)
target_link_libraries(stk_middle_mesh_util PUBLIC stk_middle_mesh)
target_include_directories(stk_middle_mesh_util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(stk_middle_mesh_util PUBLIC STK_BUILT_IN_SIERRA)
add_executable(stk_middle_mesh_util_utest)
target_sources(stk_middle_mesh_util_utest PRIVATE ./../stk_unit_tests/stk_middle_mesh_util/test_create_stk_mesh.cpp
    ./../stk_unit_tests/stk_middle_mesh_util/test_exodus_writer.cpp
    ./../stk_unit_tests/stk_middle_mesh_util/test_stk_field_copier.cpp
    ./../stk_unit_tests/stk_middle_mesh_util/test_stk_interface.cpp)
target_link_libraries(stk_middle_mesh_util_utest PUBLIC stk_mesh_fixtures)
target_link_libraries(stk_middle_mesh_util_utest PUBLIC stk_middle_mesh)
target_link_libraries(stk_middle_mesh_util_utest PUBLIC stk_middle_mesh_util)
target_link_libraries(stk_middle_mesh_util_utest PUBLIC stk_unit_main)
target_include_directories(stk_middle_mesh_util_utest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(stk_middle_mesh_util_utest PUBLIC STK_BUILT_IN_SIERRA)