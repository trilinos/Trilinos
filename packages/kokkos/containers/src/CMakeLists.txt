INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------

SET(HEADERS "")
SET(SOURCES "")

SET(HEADERS_IMPL "")

FILE(GLOB HEADERS *.hpp)
FILE(GLOB HEADERS_IMPL impl/*.hpp)

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR})

INSTALL(FILES ${HEADERS} DESTINATION ${TRILINOS_INCDIR}/)
INSTALL(FILES ${HEADERS_IMPL} DESTINATION ${TRILINOS_INCDIR}/impl/)


TRIBITS_ADD_LIBRARY(
    kokkoscontainers
    HEADERS ${HEADERS}
    NOINSTALLHEADERS ${HEADERS_IMPL}
    SOURCES ${SOURCES}
    DEPLIBS
    )

#-----------------------------------------------------------------------------
