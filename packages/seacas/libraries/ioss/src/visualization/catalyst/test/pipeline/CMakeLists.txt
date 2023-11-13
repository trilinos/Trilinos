ADD_EXECUTABLE(CatalystLoggingPipelineTests
               CatalystLoggingPipelineTests.cxx)
TARGET_LINK_LIBRARIES(CatalystLoggingPipelineTests
                      catalystiosstest)

ADD_EXECUTABLE(ExodusToCatalystToScriptTests
               ExodusToCatalystToScriptTests.cxx)
TARGET_LINK_LIBRARIES(ExodusToCatalystToScriptTests
                      catalystiosstest)

ADD_EXECUTABLE(CgnsToCatalystToScriptTests
               CgnsToCatalystToScriptTests.cxx)
TARGET_LINK_LIBRARIES(CgnsToCatalystToScriptTests
                      catalystiosstest)

ADD_EXECUTABLE(TwoGridsToCatalystToScriptTests
               TwoGridsToCatalystToScriptTests.cxx)
TARGET_LINK_LIBRARIES(TwoGridsToCatalystToScriptTests
                      catalystiosstest)

ADD_EXECUTABLE(TwoGridsToCatalystToScriptTests2
               TwoGridsToCatalystToScriptTests2.cxx)
TARGET_LINK_LIBRARIES(TwoGridsToCatalystToScriptTests2
                      catalystiosstest)

ADD_EXECUTABLE(TwoGridsToCatalystToScriptTests3
               TwoGridsToCatalystToScriptTests3.cxx)
TARGET_LINK_LIBRARIES(TwoGridsToCatalystToScriptTests3
                      catalystiosstest)

ADD_EXECUTABLE(CatalystV2ScriptFromParaViewGuiTest1
               CatalystV2ScriptFromParaViewGuiTest1.cxx)
TARGET_LINK_LIBRARIES(CatalystV2ScriptFromParaViewGuiTest1
                      catalystiosstest)

ADD_EXECUTABLE(CatalystPartitionedDataSetCollectionTests
               CatalystPartitionedDataSetCollectionTests.cxx)
TARGET_LINK_LIBRARIES(CatalystPartitionedDataSetCollectionTests
                      catalystiosstest)

ADD_EXECUTABLE(TwoGridsTwoPipelinesTests
               TwoGridsTwoPipelinesTests.cxx)
TARGET_LINK_LIBRARIES(TwoGridsTwoPipelinesTests
                      catalystiosstest)

addCatalystTest("Catalyst_Logging_Pipeline" "1;4"
    ${CMAKE_CURRENT_BINARY_DIR}/CatalystLoggingPipelineTests "")

addCatalystTest("Exodus_Catalyst_Script" "1;4"
    ${CMAKE_CURRENT_BINARY_DIR}/ExodusToCatalystToScriptTests "")

addCatalystTest("CGNS_Catalyst_Script" "1;4"
    ${CMAKE_CURRENT_BINARY_DIR}/CgnsToCatalystToScriptTests "")

addCatalystTest("Two_Grids_Catalyst_Script" "1;4"
    ${CMAKE_CURRENT_BINARY_DIR}/TwoGridsToCatalystToScriptTests "")

addCatalystTest("Two_Grids_Catalyst_Script2" "1;4"
    ${CMAKE_CURRENT_BINARY_DIR}/TwoGridsToCatalystToScriptTests2 "")

addCatalystTest("Two_Grids_Catalyst_Script3" "1;4"
    ${CMAKE_CURRENT_BINARY_DIR}/TwoGridsToCatalystToScriptTests3 "")

addCatalystTest("CatalystV2ScriptFromParaViewGuiTest1" "1;4"
    ${CMAKE_CURRENT_BINARY_DIR}/CatalystV2ScriptFromParaViewGuiTest1 "")

addCatalystTest("CatalystPartitionedDataSetCollectionTests" "1;4"
    ${CMAKE_CURRENT_BINARY_DIR}/CatalystPartitionedDataSetCollectionTests "")

addCatalystTest("TwoGridsTwoPipelinesTests" "1;4"
    ${CMAKE_CURRENT_BINARY_DIR}/TwoGridsTwoPipelinesTests "")

