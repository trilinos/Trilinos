
ROL_ADD_EXECUTABLE(
  example_02
  SOURCES example_02.cpp
  ADD_DIR_TO_NAME
)


ROL_COPY_FILES_TO_BINARY_DIR(
  TopologyOptimizationDataCopy
  SOURCE_FILES
    input.xml
    input_ex02.xml
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
)

string( TOLOWER ${CMAKE_BUILD_TYPE} lower_build_type )

if( ${lower_build_type} STREQUAL "release" )

  ROL_ADD_EXECUTABLE_AND_TEST(
    example_01
    SOURCES example_01.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

else()

  ROL_ADD_EXECUTABLE(
    example_01
    SOURCES example_01.cpp
    ADD_DIR_TO_NAME
  )

endif()
