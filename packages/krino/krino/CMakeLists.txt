add_subdirectory(interface_geometry_interface)
add_subdirectory(master_element)
add_subdirectory(krino_lib)
add_subdirectory(adaptivity_interface)
add_subdirectory(region)
add_subdirectory(rebalance_utils)
add_subdirectory(parser)
add_subdirectory(unit_tests)

SET(SOURCES_MAIN Apps_krino.cpp)

set(EXE_NAME krino)
TRIBITS_ADD_EXECUTABLE(
    ${EXE_NAME} 
    SOURCES ${SOURCES_MAIN}
    NOEXEPREFIX NOEXESUFFIX INSTALLABLE
    )
