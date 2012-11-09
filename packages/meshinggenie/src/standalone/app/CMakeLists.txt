cmake_minimum_required( VERSION 2.4.6 )
mark_as_advanced( CMAKE_BACKWARDS_COMPATIBILITY )
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

# add our generated source to the build 
set( APP_SRCS
	MeshingGenie_main.cpp
   )
   
INCLUDE_DIRECTORIES( 
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${CMAKE_CURRENT_BINARY_DIR} 	   
)


SET(EXECUTABLE_OUTPUT_PATH ../bin)

#specify the sources to compile into our exe
ADD_EXECUTABLE(MeshingGenie ${APP_SRCS})

#link our exe to the rest of the librarires
TARGET_LINK_LIBRARIES(MeshingGenie MeshingGenie_mps)
