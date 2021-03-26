set( DIR ${CMAKE_CURRENT_SOURCE_DIR} )

add_subdirectory( objective )
add_subdirectory( operator )

APPEND_GLOBAL_SET( HEADERS ${DIR}/*.hpp )
