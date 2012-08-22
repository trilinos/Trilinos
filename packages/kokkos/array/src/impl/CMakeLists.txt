INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

SET(HEADERS "")
SET(SOURCES "")

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)

TRIBITS_ADD_LIBRARY(
    kokkosarray_impl
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS 
    )

