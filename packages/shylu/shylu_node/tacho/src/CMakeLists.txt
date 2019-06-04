IF (Tacho_ENABLE_INT_INT)
  SET(TACHO_USE_INT_INT ON)
  SET(TACHO_USE_INT_SIZE_T OFF)
ELSE()
  SET(TACHO_USE_INT_INT OFF)
  SET(TACHO_USE_INT_SIZE_T ON)
ENDIF()

SET(TACHO_HAVE_MKL ${TPL_ENABLE_MKL})
SET(TACHO_HAVE_LAPACK ${TPL_ENABLE_LAPACK})
SET(TACHO_HAVE_BLAS ${TPL_ENABLE_BLAS})

SET(TACHO_HAVE_METIS ${TPL_ENABLE_METIS})
SET(TACHO_HAVE_METIS_MT ${TPL_ENABLE_MTMETIS})
SET(TACHO_HAVE_SCOTCH ${TPL_ENABLE_Scotch})
SET(TACHO_HAVE_TRILINOS_SS ${Trilinos_ENABLE_TrilinosSS})
SET(TACHO_HAVE_SUITESPARSE ${TPL_ENABLE_Cholmod})
                       
SET(TACHO_HAVE_VTUNE ${TPL_ENABLE_VTune})

#-----------------------------------------------------------------------------

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#-----------------------------------------------------------------------------

SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#-----------------------------------------------------------------------------

SET(HEADERS_PUBLIC "")
SET(HEADERS_PRIVATE "")
SET(SOURCES "")

FILE(GLOB HEADERS_PUBLIC *.hpp)
LIST( APPEND HEADERS_PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h )


APPEND_SET(SOURCES
  Tacho_Util.cpp
  Tacho_Blas_External.cpp
  Tacho_Lapack_External.cpp
)

#-----------------------------------------------------------------------------

TRIBITS_ADD_LIBRARY(
    shylu_nodetacho
    HEADERS ${HEADERS_PUBLIC}
    NOINSTALLHEADERS ${HEADERS_PRIVATE}
    SOURCES ${SOURCES}
    DEPLIBS
    )
