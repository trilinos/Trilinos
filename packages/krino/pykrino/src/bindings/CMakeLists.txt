find_package(Python3 REQUIRED Development.Module)
find_package(pybind11 REQUIRED)

TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS "")
SET(SOURCES "")
FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)

pybind11_add_module(_pykrino ${SOURCES})

target_include_directories(_pykrino PUBLIC 
    ${${PACKAGE_NAME}_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR})

target_compile_features(_pykrino PUBLIC cxx_std_17)

target_link_libraries(_pykrino PUBLIC krino_lib)
set_target_properties(_pykrino PROPERTIES SUFFIX ".so")
    
INSTALL(TARGETS _pykrino
  DESTINATION ${PyKrino_INSTALL_DIR})

add_custom_command(TARGET _pykrino POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/_pykrino.so ${CMAKE_CURRENT_BINARY_DIR}/../pykrino/.
    COMMENT "Copy ${PROJECT_BINARY_DIR}/src/_pykrino.so"
)
