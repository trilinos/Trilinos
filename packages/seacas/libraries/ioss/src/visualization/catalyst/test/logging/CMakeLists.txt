ADD_EXECUTABLE(CatalystLoggingTests
               CatalystLoggingTests.cxx)
TARGET_LINK_LIBRARIES(CatalystLoggingTests
                      catalystiosstest)

addCatalystTest("Catalyst_Logging" "1"
    ${CMAKE_CURRENT_BINARY_DIR}/CatalystLoggingTests "")
