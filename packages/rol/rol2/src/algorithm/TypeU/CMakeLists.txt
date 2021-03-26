set( DIR ${CMAKE_CURRENT_SOURCE_DIR} )

APPEND_GLOBAL_SET( HEADERS ${DIR}/*.hpp )

add_subdirectory( linesearch )
add_subdirectory( trustregion )
