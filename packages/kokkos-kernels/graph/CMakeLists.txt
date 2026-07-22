list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/graph/src)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/graph/impl)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/graph/eti)
list(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/graph/eti)

# Adding unit-tests
kokkoskernels_include_directories(${CMAKE_CURRENT_BINARY_DIR}/graph)
kokkoskernels_include_directories(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/graph)

##################
#                #
# ETI generation #
#                #
##################

#Build up a list of DECL, AVAIL, and INST macros
#that should be instantiated based on input options
#Generate @X@ variables in the template X.hpp.in and X.cpp.in
#files containing the list of all needed macros
kokkoskernels_generate_eti(Graph_color_d1 color_d1
  COMPONENTS  graph
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS ORDINALS OFFSETS LAYOUTS DEVICES)
