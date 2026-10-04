# Standalone build of PinT example Van_der_Pol/example_01 (serial: runs single-process).
# Reached only under STANDALONE_ROL AND ENABLE_MPI; MPI is linked so the real
# ROL::GlobalMPISession resolves. Own subdir => own working dir => no XML-name clash.
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/../../PinT/Van_der_Pol)

add_executable(Van_der_Pol_example_01 ${SRC}/example_01.cpp)
target_include_directories(Van_der_Pol_example_01 PRIVATE ${SRC})
target_link_libraries(Van_der_Pol_example_01 PRIVATE rol MPI::MPI_CXX)

configure_file(${SRC}/VdP_Parameters.xml ${CMAKE_CURRENT_BINARY_DIR}/VdP_Parameters.xml COPYONLY)

if(ENABLE_TESTS)
  add_test(NAME Van_der_Pol_example_01 COMMAND $<TARGET_FILE:Van_der_Pol_example_01> PrintItAll)
  set_tests_properties(Van_der_Pol_example_01 PROPERTIES PASS_REGULAR_EXPRESSION "TEST PASSED")
endif()
