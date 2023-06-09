LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/graph/src)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/graph/impl)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/graph/eti)
LIST(APPEND KK_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/graph/eti)


# Adding unit-tests
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/graph)
KOKKOSKERNELS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/graph)


##################
#                #
# ETI generation #
#                #
##################

#Build up a list of DECL, AVAIL, and INST macros
#that should be instantiated based on input options
#Generate @X@ variables in the template X.hpp.in and X.cpp.in
#files containing the list of all needed macros
KOKKOSKERNELS_GENERATE_ETI(Graph_color_d1 color_d1
  COMPONENTS  graph
  HEADER_LIST ETI_HEADERS
  SOURCE_LIST SOURCES
  TYPE_LISTS  FLOATS ORDINALS OFFSETS LAYOUTS DEVICES
)
