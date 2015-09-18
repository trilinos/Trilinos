INCLUDE(TrilinosCreateClientTemplateHeaders)

# NTS: The funny thing is that I need to include this in order to
# get the link/include directories correct.  Go figure.
SET(HEADERS_LIB "")
SET(SOURCES_LIB "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

if(TPL_ENABLE_MATLAB)
  INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIRS})

  # Must glob the binary dir last to get all of the auto-generated headers
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})

  APPEND_GLOB(HEADERS_LIB *.hpp)
  APPEND_SET(SOURCES_LIB MueLu_SingleLevelMatlabFactory.cpp MueLu_TwoLevelMatlabFactory.cpp MueLu_MatlabSmoother.cpp MueLu_MatlabUtils.cpp)

  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR} NOSIERRABJAM)
  TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

  # Start with the library
  TRIBITS_ADD_LIBRARY(
    muelu-matlab
    HEADERS ${HEADERS_LIB}
    SOURCES ${SOURCES_LIB}
#    DEPLIBS muelu muelu-adapters
  )
ENDIF()

