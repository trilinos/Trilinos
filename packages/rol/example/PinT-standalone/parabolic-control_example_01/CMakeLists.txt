# Standalone build of PinT example parabolic-control/example_01 (serial: runs single-process).
# Reached only under STANDALONE_ROL AND ENABLE_MPI; MPI is linked so the real
# ROL::GlobalMPISession resolves. Own subdir => own working dir => no XML-name clash
# (its input_ex01.xml / input_rol.xml differ from nonlinear's).
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/../../PinT/parabolic-control)

add_executable(parabolic-control_example_01 ${SRC}/example_01.cpp)
target_include_directories(parabolic-control_example_01 PRIVATE ${SRC})
target_link_libraries(parabolic-control_example_01 PRIVATE rol MPI::MPI_CXX)

configure_file(${SRC}/input_ex01.xml ${CMAKE_CURRENT_BINARY_DIR}/input_ex01.xml COPYONLY)
configure_file(${SRC}/input_rol.xml  ${CMAKE_CURRENT_BINARY_DIR}/input_rol.xml  COPYONLY)

if(ENABLE_TESTS)
  add_test(NAME parabolic-control_example_01 COMMAND $<TARGET_FILE:parabolic-control_example_01> PrintItAll)
  set_tests_properties(parabolic-control_example_01 PROPERTIES PASS_REGULAR_EXPRESSION "TEST PASSED")
endif()
