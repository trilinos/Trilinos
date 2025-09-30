
FILE(COPY test_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
SET(TEST_DATA_DIRECTORY_PATH ${CMAKE_CURRENT_BINARY_DIR}/test_data/)
CONFIGURE_FILE(TestDataDirectoryPath.h.in TestDataDirectoryPath.h @ONLY)

ADD_LIBRARY(catalystiosstest STATIC
            CatalystTestFixture.cxx
            TestDriverMain.cxx)

TARGET_LINK_LIBRARIES(catalystiosstest
                      catalystiossapp
                      VTK::jsoncpp
                      VTK::IOParallelXML
                      MPI::MPI_C)

TARGET_INCLUDE_DIRECTORIES(catalystiosstest PUBLIC
                           ${PROJECT_SOURCE_DIR}/ioss2catalyst
                           ${PROJECT_SOURCE_DIR}/test/Catch2
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${CMAKE_CURRENT_SOURCE_DIR})
