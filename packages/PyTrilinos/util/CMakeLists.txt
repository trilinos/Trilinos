CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/CMakeCacheVariables.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/CMakeCacheVariables.py   )

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/MakefileVariables.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/MakefileVariables.py   )

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/copyWithCMakeSubstitutions.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/copyWithCMakeSubstitutions.py   )

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/pkg_info.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/pkg_info.py   )
