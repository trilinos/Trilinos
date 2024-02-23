add_subdirectory(diagwriter)
add_subdirectory(master_element)
add_subdirectory(geometry)
add_subdirectory(math_utils)
add_subdirectory(surface)
add_subdirectory(mesh_surface)
add_subdirectory(quality_metric)
add_subdirectory(mesh_utils)
add_subdirectory(refinement)
add_subdirectory(krino_lib)
add_subdirectory(region)
add_subdirectory(rebalance_utils)
add_subdirectory(parser)
TRIBITS_ADD_TEST_DIRECTORIES(unit_tests)

SET(SOURCES_MAIN Apps_krino.cpp)

set(EXE_NAME krino)
TRIBITS_ADD_EXECUTABLE(
    ${EXE_NAME} 
    SOURCES ${SOURCES_MAIN}
    NOEXEPREFIX NOEXESUFFIX INSTALLABLE
    )
