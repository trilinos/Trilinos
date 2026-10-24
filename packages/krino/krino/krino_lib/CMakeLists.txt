TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# krino_lib is now a source-free umbrella library.  Its former contents have
# been dispersed into decomp, interface_geometry, post_process, mesh_clone and
# levelset_shape_sens; it re-exports them (INTERFACE DEPLIBS) so existing
# consumers (region, parser, rebalance_utils, unit_test_utils) keep resolving
# every symbol transitively through the single krino_lib dependency.  It has no
# sources or headers of its own, so it is declared HEADERONLY (an INTERFACE
# library); TRIBITS_ADD_LIBRARY requires SOURCES to be empty in that case.
TRIBITS_ADD_LIBRARY(
  krino_lib
  HEADERONLY
  DEPLIBS krino_decomp_lib krino_interface_geometry_lib krino_post_process_lib krino_mesh_clone_lib krino_levelset_shape_sens_lib)
