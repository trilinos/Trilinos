# Standalone build of PinT example tanks/LowerBandedMatrix_test (serial: runs single-process).
# Reached only under STANDALONE_ROL AND ENABLE_MPI; MPI is linked so the real
# ROL::GlobalMPISession resolves. Reads no XML input.
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/../../PinT/tanks)

add_executable(tanks_LowerBandedMatrix_test ${SRC}/LowerBandedMatrix_test.cpp)
target_include_directories(tanks_LowerBandedMatrix_test PRIVATE ${SRC})
target_link_libraries(tanks_LowerBandedMatrix_test PRIVATE rol MPI::MPI_CXX)

if(ENABLE_TESTS)
  add_test(NAME tanks_LowerBandedMatrix_test COMMAND $<TARGET_FILE:tanks_LowerBandedMatrix_test> PrintItAll)
  set_tests_properties(tanks_LowerBandedMatrix_test PROPERTIES PASS_REGULAR_EXPRESSION "TEST PASSED")
endif()
