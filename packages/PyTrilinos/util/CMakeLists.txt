# -*- cmake -*-

# Generate the CMakeCacheVariables python script
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/CMakeCacheVariables.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeCacheVariables.py)

# Generate the MakefileVariables python script
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/MakefileVariables.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/MakefileVariables.py)

# Generate the copyWithCMakeSubstitutions python script
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/copyWithCMakeSubstitutions.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/copyWithCMakeSubstitutions.py)
