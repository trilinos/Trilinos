add_subdirectory(config)
add_subdirectory(diagwriter)
add_subdirectory(master_element)
add_subdirectory(geometry)
add_subdirectory(math_utils)
add_subdirectory(quality_metric)
add_subdirectory(quality_metric_sens)
add_subdirectory(surface)
add_subdirectory(mesh_surface)
add_subdirectory(mesh_utils)
add_subdirectory(refinement)
add_subdirectory(refinement_rebalance)
add_subdirectory(krino_lib)
add_subdirectory(rol)
add_subdirectory(region)
add_subdirectory(rebalance_utils)
add_subdirectory(parser)
add_subdirectory(unit_main)
TRIBITS_ADD_TEST_DIRECTORIES(unit_tests)

SET(SOURCES_MAIN Apps_krino.cpp)

set(EXE_NAME krino)
TRIBITS_ADD_EXECUTABLE(
    ${EXE_NAME} 
    SOURCES ${SOURCES_MAIN}
    NOEXEPREFIX NOEXESUFFIX INSTALLABLE
    )
