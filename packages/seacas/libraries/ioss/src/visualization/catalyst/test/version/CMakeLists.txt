ADD_EXECUTABLE(CatalystVersionTests
               CatalystVersionTests.cxx)
TARGET_LINK_LIBRARIES(CatalystVersionTests
                      catalystiosstest)

addCatalystTest("Catalyst_Version" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/CatalystVersionTests "")
