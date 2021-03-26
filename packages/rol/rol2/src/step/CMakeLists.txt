set( DIR ${CMAKE_CURRENT_SOURCE_DIR} )

add_subdirectory( krylov )
add_subdirectory( secant)

APPEND_GLOBAL_SET( HEADERS ${DIR}/*.hpp )
