INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  graph_color
  SOURCES KokkosGraph_color.cpp       
  )

#Below will probably fail on GPUs.
#TRIBITS_ADD_EXECUTABLE(
#  graph_triangle
#  SOURCES KokkosGraph_triangle.cpp      
#  )

