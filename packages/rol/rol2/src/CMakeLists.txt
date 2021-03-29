include(TribitsLibraryMacros)
#
# A) Package-specific configuration options
#

#TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
yellow("ROL source directory:  ${ROL_SRC}")

set(ROL_FUNCTION  ${ROL_SRC}/function       )
set(ROL_OBJECTIVE ${ROL_FUNCTION}/objective )
set(ROL_OPERATOR  ${ROL_FUNCTION}/operator  )


set(DIR ${CMAKE_CURRENT_SOURCE_DIR})
yellow("ROL2 source directory: ${DIR}")


set(HEADERS "")
set(SOURCES "")

include_directories( ${DIR}
                     ${DIR}/algorithm
                     ${DIR}/algorithm/TypeU
                     ${DIR}/algorithm/TypeU/linesearch
                     ${DIR}/algorithm/TypeU/linesearch/descent
                     ${DIR}/algorithm/TypeU/trustregion
                     ${DIR}/elementwise
                     ${DIR}/function
                     ${DIR}/function/objective
                     ${DIR}/function/operator
                     ${DIR}/status
                     ${DIR}/step
                     ${DIR}/step/krylov
                     ${DIR}/step/secant
                     ${DIR}/utils
                     ${DIR}/vector
)

APPEND_GLOBAL_SET( HEADERS ${DIR}/ROL2.hpp )

add_subdirectory( algorithm )
add_subdirectory( elementwise )
add_subdirectory( function )
add_subdirectory( utils )
add_subdirectory( vector )



bold_blue("List of headers")
foreach(H ${HEADERS})
  blue(${H})
endforeach()

bold_green("List of sources")
foreach(S ${SOURCES})
  green(${S})
endforeach()


TRIBITS_ADD_LIBRARY( rol
                     HEADERS ${HEADERS}
                     SOURCES ${SOURCES}
                   )        
