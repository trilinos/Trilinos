# ROL Version 2

# ROL2 Requires C++14
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
include( Colors )
include( Dependencies ) 

bold_white("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
set( ROL_HOME  ${PACKAGE_SOURCE_DIR})
set( ROL2_HOME ${CMAKE_CURRENT_SOURCE_DIR})
set( ROL_SRC ${ROL_HOME}/src )

bold_green("\nConfiguring to use ROL2")
#yellow("ROL root directory:  ${ROL_HOME}") 
#yellow("ROL2 root directory: ${ROL2_HOME}")

add_subdirectory(src)
add_subdirectory(test)


