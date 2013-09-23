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

# mfh 20 Sep 2013: This subpackage currently only contains header files.  Thus,
# there is no library to build.  See, for example, this Stack Overflow post:
#
# http://stackoverflow.com/questions/16655705/cmake-cannot-determine-linker-language-for-target
#
# If we add source files to this subpackage, for example, to support explicit
# template instantiation, then we will need to uncomment out the
# TRIBITS_ADD_LIBRARY statement below.

TRIBITS_ADD_LIBRARY(
    kokkoscontainers
    HEADERS ${HEADERS}
    NOINSTALLHEADERS ${HEADERS_IMPL}
    SOURCES ${SOURCES}
    DEPLIBS
    )

#-----------------------------------------------------------------------------
