# Standalone build of PinT example nonlinear/example_01 (serial: runs single-process).
# Reached only under STANDALONE_ROL AND ENABLE_MPI; MPI is linked so the real
# ROL::GlobalMPISession resolves. Own subdir => own working dir => no XML-name clash
# (this example's input_ex01.xml / input_rol.xml differ from parabolic-control's).
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/../../PinT/nonlinear)

add_executable(nonlinear_example_01 ${SRC}/example_01.cpp)
target_include_directories(nonlinear_example_01 PRIVATE ${SRC})
target_link_libraries(nonlinear_example_01 PRIVATE rol MPI::MPI_CXX)

configure_file(${SRC}/input_ex01.xml ${CMAKE_CURRENT_BINARY_DIR}/input_ex01.xml COPYONLY)
configure_file(${SRC}/input_rol.xml  ${CMAKE_CURRENT_BINARY_DIR}/input_rol.xml  COPYONLY)

if(ENABLE_TESTS)
  add_test(NAME nonlinear_example_01 COMMAND $<TARGET_FILE:nonlinear_example_01> PrintItAll)
  set_tests_properties(nonlinear_example_01 PROPERTIES PASS_REGULAR_EXPRESSION "TEST PASSED")
endif()
