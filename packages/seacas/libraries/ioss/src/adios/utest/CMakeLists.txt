IF (TPL_ENABLE_gtest)

    TRIBITS_ADD_EXECUTABLE(
    Ut_ioad
    SOURCES Ut_ioad.C
    IMPORTEDLIBS Ionit
    )

    TRIBITS_ADD_TEST(
      Ut_ioad
      NAME Ut_ioad
      NUM_MPI_PROCS 1
    )

    INCLUDE_DIRECTORIES(
      "${CMAKE_CURRENT_SOURCE_DIR}/../"
      "${CMAKE_CURRENT_BINARY_DIR}/../"
      "${CMAKE_CURRENT_SOURCE_DIR}/../../"
    )
ENDIF()
