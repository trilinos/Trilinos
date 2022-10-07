
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

ADD_EXECUTABLE(ImageSetCameraValidTests
               ImageSetCameraValidTests.cxx)
TARGET_LINK_LIBRARIES(ImageSetCameraValidTests
                      catalystiosstest)

addCatalystTest("Parser_Imageset_Camera_Valid" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/ImageSetCameraValidTests "exclude:[notWorking]")

ADD_EXECUTABLE(ImageSetValidTests
               ImageSetValidTests.cxx)
TARGET_LINK_LIBRARIES(ImageSetValidTests
                      catalystiosstest)

addCatalystTest("Parser_Imageset_Valid" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/ImageSetValidTests "exclude:[notWorking]")

ADD_EXECUTABLE(ImageSetInvalidTests
               ImageSetInvalidTests.cxx)
TARGET_LINK_LIBRARIES(ImageSetInvalidTests
                      catalystiosstest)

addCatalystTest("Parser_Imageset_Invalid" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/ImageSetInvalidTests "exclude:[notWorking]")

ADD_EXECUTABLE(ImageSetOperationValidTests
               ImageSetOperationValidTests.cxx)
TARGET_LINK_LIBRARIES(ImageSetOperationValidTests
                      catalystiosstest)

addCatalystTest("Parser_Imageset_Operation_Valid" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/ImageSetOperationValidTests "exclude:[notWorking]")

ADD_EXECUTABLE(ImageSetRepresentationValidTests
               ImageSetRepresentationValidTests.cxx)
TARGET_LINK_LIBRARIES(ImageSetRepresentationValidTests
                      catalystiosstest)

addCatalystTest("Parser_Imageset_Representation_Valid" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/ImageSetRepresentationValidTests "exclude:[notWorking]")

ADD_EXECUTABLE(OperationsValidTests
               OperationsValidTests.cxx)
TARGET_LINK_LIBRARIES(OperationsValidTests
                      catalystiosstest)

addCatalystTest("Parser_Operations_Valid" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/OperationsValidTests "exclude:[notWorking]")
