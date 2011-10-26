ADD_SUBDIRECTORY(src)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)
PACKAGE_ADD_TEST_DIRECTORIES(test)

# Add directory xpetra/cmake/ to the cmake module search path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PACKAGE_SOURCE_DIR}/xpetra/cmake)
