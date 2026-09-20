# Standalone build of PinT example tanks/example_04 (serial: runs single-process).
# Reached only under STANDALONE_ROL AND ENABLE_MPI; MPI is linked so the real
# ROL::GlobalMPISession resolves. Own subdir => own working dir.
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/../../PinT/tanks)

add_executable(tanks_example_04 ${SRC}/example_04.cpp)
target_include_directories(tanks_example_04 PRIVATE ${SRC})
target_link_libraries(tanks_example_04 PRIVATE rol MPI::MPI_CXX)

configure_file(${SRC}/parameters_ex04.xml ${CMAKE_CURRENT_BINARY_DIR}/parameters_ex04.xml COPYONLY)
configure_file(${SRC}/rol-parameters.xml  ${CMAKE_CURRENT_BINARY_DIR}/rol-parameters.xml  COPYONLY)

if(ENABLE_TESTS)
  add_test(NAME tanks_example_04 COMMAND $<TARGET_FILE:tanks_example_04> PrintItAll)
  set_tests_properties(tanks_example_04 PROPERTIES PASS_REGULAR_EXPRESSION "TEST PASSED")
endif()
