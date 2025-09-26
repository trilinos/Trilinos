INCLUDE(TrilinosCreateClientTemplateHeaders)

# NTS: The funny thing is that I need to include this in order to
# get the link/include directories correct.  Go figure.
SET(HEADERS_LIB "")
SET(SOURCES_LIB "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

IF(TPL_ENABLE_MATLAB)
  TRIBITS_INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIRS})

  # Must glob the binary dir last to get all of the auto-generated headers
  TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})

  # Force the code into R2017b compatibility mode
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMX_COMPAT_64 -DMATLAB_MEXCMD_RELEASE=R2017b")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMX_COMPAT_64 -DMATLAB_MEXCMD_RELEASE=R2017b")

  APPEND_GLOB(HEADERS_LIB *.hpp)
  APPEND_SET(
    SOURCES_LIB
    MueLu_SingleLevelMatlabFactory.cpp
    MueLu_TwoLevelMatlabFactory.cpp
    MueLu_MatlabSmoother.cpp
    MueLu_MatlabUtils.cpp
  )

  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR} NOSIERRABJAM)
  TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

  # Start with the library
  TRIBITS_ADD_LIBRARY(
    muelu-matlab
    HEADERS ${HEADERS_LIB}
    SOURCES ${SOURCES_LIB}
    DEPLIBS muelu muelu-adapters
  )
ENDIF()
