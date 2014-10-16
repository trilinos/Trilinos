
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)
LIST(APPEND HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h)

#-----------------------------------------------------------------------------

# Can't call the library kokkoslinalg as it conflicts with KokkosClassic
# library of the same name.
TRIBITS_ADD_LIBRARY(
    kokkoscorelinalg
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS
    )

