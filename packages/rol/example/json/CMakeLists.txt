
# To build this example, add the following to the cmake call in your trilinos build script
# -D TPL_ENABLE_JSONCPP=ON \
# -D JSONCPP_INCLUDE_DIRS:PATH=/usr/include/jsoncpp \
# -D JSONCPP_LIBRARY_DIRS=/usr/lib/x86_64-linux-gnu \
# -D JSONCPP_LIBRARY_NAMES:STRING="jsoncpp" \
#
# These above paths are default for Ubuntu 64 bit if jsoncpp is installed using 
# sudo apt-get install libjsoncpp-dev  
# 
# To install json-cpp for another OS see https://github.com/open-source-parsers/jsoncpp

IF(TPL_ENABLE_JSONCPP)

    include_directories(${JSONCPP_INCLUDE_DIRS})

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
        example_01
        SOURCES example_01.cpp 
        IMPORTEDLIBS jsoncpp
        ARGS PrintItAll
        NUM_MPI_PROCS 1
        PASS_REGULAR_EXPRESSION "TEST PASSED"
        ADD_DIR_TO_NAME
    )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(
	JsonDataCopy
	SOURCE_FILES
	 parameters.json
	SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
	DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
    )

ENDIF()


