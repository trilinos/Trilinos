
ADD_EXECUTABLE(CameraValidTests
               CameraValidTests.cxx)
TARGET_LINK_LIBRARIES(CameraValidTests
                      catalystiosstest)

addCatalystTest("Parser_Camera_Valid" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/CameraValidTests "exclude:[notWorking]")

ADD_EXECUTABLE(CameraInvalidTests
               CameraInvalidTests.cxx)
TARGET_LINK_LIBRARIES(CameraInvalidTests
                      catalystiosstest)

addCatalystTest("Parser_Camera_Invalid" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/CameraInvalidTests "exclude:[notWorking]")
