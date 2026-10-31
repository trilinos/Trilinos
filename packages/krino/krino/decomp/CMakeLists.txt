SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

FILE(GLOB HEADERS *.hpp)
FILE(GLOB SOURCES *.cpp)

TRIBITS_ADD_LIBRARY(
  krino_decomp_lib
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS krino_interface_geometry_lib krino_parametric_element_lib krino_eikonal_lib krino_foundation_lib krino_field_calcs_lib krino_support_lib krino_levelset_lib krino_contour_lib krino_mesh_interface_lib krino_phase_lib)

INSTALL(FILES ${HEADERS} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/krino_decomp_lib)
