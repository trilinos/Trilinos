
SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(SOURCES ${DIR}/*.C)

# TriBITs recommends to find boost on your own,
# see seacas/cmake/tribits/common_tpls/FindTPLBoost.cmake
find_package(Boost REQUIRED
    COMPONENTS
    log
    system
    )

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/../"
  "${CMAKE_CURRENT_BINARY_DIR}/../"
  ${Boost_INCLUDE_DIR}
  )

TRIBITS_ADD_LIBRARY(
  Iodw
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS
  Ioss
  Boost::system
  Boost::log
  )

IF (BUILD_TESTING)
  ENABLE_TESTING()
ENDIF (BUILD_TESTING)
